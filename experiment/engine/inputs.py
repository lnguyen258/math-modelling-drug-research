"""Load data/processed/ into arrays indexed [county, year] (Layers 0-1 of REPORT.md)."""
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from .params import XI_NAMES

DATA = Path(__file__).resolve().parents[2] / "data"
PROCESSED = DATA / "processed"


@dataclass
class Inputs:
    fips: np.ndarray          # [C] county FIPS (str)
    state_of: np.ndarray      # [C] state index
    states: list              # state names
    years: np.ndarray         # [Y]
    pop: np.ndarray           # [C, Y] real population P_c(y)
    deaths_obs: np.ndarray    # [C, Y] observed deaths (nan where suppressed)
    suppressed: np.ndarray    # [C, Y] bool, count hidden in 1..9
    xi: np.ndarray            # [C, Y, K] structural indices, order XI_NAMES
    a_rx: np.ndarray          # [C, Y] prescriptions per person per year
    centroid: np.ndarray      # [C, 2] population-weighted centroid (km, EPSG:5070)
    tracts: pd.DataFrame      # tract, county index, pop, x_km, y_km, area_km2
    stringency: np.ndarray    # [states, days] OxCGRT stringency / 100, day 0 = 2020-01-01

    def stringency_at(self, t):
        """Stringency S in [0, 1] per state at time t (years); 0 outside 2020-2022."""
        d = int((t - 2020.0) * 365.25)
        if d < 0 or d >= self.stringency.shape[1]:
            return np.zeros(self.stringency.shape[0])
        return self.stringency[:, d]


def _z(x):
    """Pooled standardization over all counties and years, clipped at +-3 SD (REPORT.md)."""
    return np.clip((x - np.nanmean(x)) / np.nanstd(x), -3, 3)


def _std(x):
    return (x - np.nanmean(x)) / np.nanstd(x)


def build_xi(cy_wide):
    """Composite indices of REPORT.md from wide arrays {column: [C, Y]}."""
    v = cy_wide
    xi = {
        "econ": _std((_z(v["poverty_pct"]) + _z(v["unemployment_pct"]) + _z(v["nohs_pct"])
                      - _z(np.log(v["median_hh_income"])) - _z(v["edu_bachelors_plus_pct"])) / 5),
        "pain": _std((_z(v["disability_pct"]) + _z(v["veteran_pct"])) / 2),
        "unins": _z(v["uninsured_pct"]),
        "age": _z(v["median_age"]),
    }
    shape = v["poverty_pct"].shape
    return np.stack([xi.get(k, np.zeros(shape)) for k in XI_NAMES], axis=-1)


def load_inputs():
    cy = pd.read_csv(PROCESSED / "county_year.csv", dtype={"fips": str, "state_fips": str})
    counties = pd.read_csv(PROCESSED / "counties.csv", dtype={"fips": str, "state_fips": str})
    tracts = pd.read_csv(PROCESSED / "tracts.csv", dtype={"tract": str, "fips": str})

    fips = np.sort(counties["fips"].unique())
    years = np.sort(cy["year"].unique())
    states = sorted(counties["state"].unique())
    st = counties.set_index("fips").loc[fips, "state"]
    state_of = np.array([states.index(s) for s in st])

    def wide(col):
        return cy.pivot(index="fips", columns="year", values=col).loc[fips, years].to_numpy(dtype=float)

    cols = ["poverty_pct", "unemployment_pct", "median_hh_income", "nohs_pct", "edu_bachelors_plus_pct",
            "disability_pct", "veteran_pct", "uninsured_pct", "median_age"]
    xi = build_xi({c: wide(c) for c in cols})

    cidx = {f: k for k, f in enumerate(fips)}
    tracts["c"] = tracts["fips"].map(cidx)
    w = tracts["pop"].to_numpy(float)
    cen = np.zeros((len(fips), 2))
    np.add.at(cen, tracts["c"].to_numpy(), np.c_[tracts["x_km"] * w, tracts["y_km"] * w])
    cen /= np.bincount(tracts["c"], weights=w, minlength=len(fips))[:, None]

    ox = pd.read_csv(DATA / "oxcgrt" / "OxCGRT_US_latest.csv", usecols=["RegionName", "Jurisdiction", "Date",
                                                                       "StringencyIndex"])
    ox = ox[ox["Jurisdiction"] == "STATE_WIDE"]
    ox["Date"] = pd.to_datetime(ox["Date"].astype(str))
    ox = ox.pivot_table(index="RegionName", columns="Date", values="StringencyIndex").reindex(states)
    stringency = ox.ffill(axis=1).fillna(0.0).to_numpy() / 100.0

    return Inputs(fips=fips, state_of=state_of, states=states, years=years,
                  pop=wide("pop"), deaths_obs=wide("deaths_obs"),
                  suppressed=cy.pivot(index="fips", columns="year", values="suppressed").loc[fips, years].to_numpy(bool),
                  xi=xi, a_rx=wide("rx_rate") / 100.0, centroid=cen, tracts=tracts, stringency=stringency)
