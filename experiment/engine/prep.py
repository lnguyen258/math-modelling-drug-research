"""Build the engine's input tables from data/ (Layer 0 of REPORT.md).

Run once:  python -m engine.prep        (from the experiment/ folder)

Writes data/processed/:
  county_year.csv  one row per county-year 2010-2023: population, observed deaths
                   (with suppression flag), ACS covariates, CDC Rx dispensing rate
  counties.csv     one row per county: state, land area
  tracts.csv       2020 census tracts: population, projected centroid (km), land area
"""
import io
import zipfile
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import requests
from bs4 import BeautifulSoup
from pyproj import Transformer

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data"
OUT = DATA / "processed"
YEARS = list(range(2010, 2024))
ACS_YEARS = list(range(2012, 2024))

# CDC national opioid dispensing rate per 100 persons. Used only to extend the
# county series past 2019 (county pages for 2020+ are not downloadable).
# 2021 has no published figure found; it is interpolated between 2020 and 2022.
# The 2019 county page is not used: its county medians halve vs 2018 while the
# national rate fell 9%, which points to a method change, so we extend from 2018.
CDC_NATIONAL_RX = {2018: 51.4, 2019: 46.8, 2020: 43.3, 2021: 41.4, 2022: 39.5, 2023: 37.5}
RX_LAST_COUNTY_YEAR = 2018


def load_deaths():
    d = pd.read_csv(DATA / "SRDDRC CDCWONDER_Opioid_YrStCty.csv", dtype=str, encoding="utf-8-sig")
    d = d.rename(columns={"CountyCode": "fips", "Year": "year", "StateCode": "state_fips", "State": "state"})
    d["year"] = d["year"].astype(int)
    d["pop"] = pd.to_numeric(d["Population"])
    d["suppressed"] = d["Deaths"].eq("Suppressed")
    d["deaths_obs"] = pd.to_numeric(d["Deaths"], errors="coerce")
    return d[["fips", "year", "state", "state_fips", "County", "pop", "deaths_obs", "suppressed"]]


def load_acs():
    frames = []
    for y in ACS_YEARS:
        a = pd.read_csv(DATA / f"acs_{y}_county_yearly.csv", dtype={"CountyCode": str})
        a["year"] = y
        frames.append(a)
    a = pd.concat(frames).rename(columns={"CountyCode": "fips"})
    a["nohs_pct"] = a["edu_lt9th_pct"] + a["edu_9to12_nodiploma_pct"]
    keep = ["poverty_pct", "unemployment_pct", "median_hh_income", "nohs_pct", "edu_bachelors_plus_pct",
            "disability_pct", "uninsured_pct", "veteran_pct", "median_age"]
    a = a[["fips", "year"] + keep]
    # 2010-2011 have no comparable ACS county tables: hold at 2012 values.
    early = [a[a.year == 2012].assign(year=y) for y in (2010, 2011)]
    a = pd.concat([a] + early)
    # A few counties miss single values (e.g. income): fill from the county's other years.
    a = a.sort_values(["fips", "year"])
    a[keep] = a.groupby("fips")[keep].transform(lambda s: s.interpolate().bfill().ffill())
    return a


def parse_rx():
    rows = []
    for y in range(2010, RX_LAST_COUNTY_YEAR + 1):
        soup = BeautifulSoup((DATA / "rx_cdc_html" / f"county{y}.html").read_text("utf-8", "replace"), "html.parser")
        table = soup.find("table", attrs={"title": "US County Prescribing Rates"})
        heads = [th.get_text(strip=True).lower() for th in table.find("thead").find_all("th")]
        i_fips = next(k for k, h in enumerate(heads) if "fips" in h)
        i_rate = next(k for k, h in enumerate(heads) if "rate" in h)
        for tr in table.find("tbody").find_all("tr"):
            td = [c.get_text(strip=True) for c in tr.find_all("td")]
            if len(td) <= max(i_fips, i_rate):
                continue
            rate = pd.to_numeric(td[i_rate], errors="coerce")
            rows.append((td[i_fips].zfill(5), y, rate))
    rx = pd.DataFrame(rows, columns=["fips", "year", "rx_rate"])
    # Extend later years by the national trend relative to the last county year.
    base = rx[rx.year == RX_LAST_COUNTY_YEAR].set_index("fips")["rx_rate"]
    for y in range(RX_LAST_COUNTY_YEAR + 1, 2024):
        ratio = CDC_NATIONAL_RX[y] / CDC_NATIONAL_RX[RX_LAST_COUNTY_YEAR]
        rx = pd.concat([rx, pd.DataFrame({"fips": base.index, "year": y, "rx_rate": base.values * ratio})])
    return rx


def fetch_tracts(state_fips, api_key):
    gaz = DATA / "tracts" / "2020_Gaz_tracts_national.txt"
    if not gaz.exists():
        gaz.parent.mkdir(parents=True, exist_ok=True)
        url = "https://www2.census.gov/geo/docs/maps-data/data/gazetteer/2020_Gazetteer/2020_Gaz_tracts_national.zip"
        zipfile.ZipFile(io.BytesIO(requests.get(url, timeout=120).content)).extractall(gaz.parent)
    g = pd.read_csv(gaz, sep="\t", dtype={"GEOID": str})
    g.columns = g.columns.str.strip()
    pop_file = DATA / "tracts" / "tract_pop_2020.csv"
    if not pop_file.exists():
        parts = []
        for s in state_fips:
            r = requests.get("https://api.census.gov/data/2020/dec/pl",
                             params={"get": "P1_001N", "for": "tract:*", "in": [f"state:{s}", "county:*"],
                                     "key": api_key}, timeout=120)
            r.raise_for_status()
            js = r.json()
            parts.append(pd.DataFrame(js[1:], columns=js[0]))
        p = pd.concat(parts)
        p["GEOID"] = p["state"] + p["county"] + p["tract"]
        p[["GEOID", "P1_001N"]].to_csv(pop_file, index=False)
    p = pd.read_csv(pop_file, dtype={"GEOID": str})
    t = g.merge(p, on="GEOID")
    to_albers = Transformer.from_crs("EPSG:4326", "EPSG:5070", always_xy=True)
    x, y = to_albers.transform(t["INTPTLONG"].values, t["INTPTLAT"].values)
    return pd.DataFrame({"tract": t["GEOID"], "fips": t["GEOID"].str[:5], "pop": t["P1_001N"].astype(int),
                         "x_km": x / 1000, "y_km": y / 1000, "area_km2": t["ALAND"] / 1e6})


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    deaths = load_deaths()
    cy = deaths.merge(load_acs(), on=["fips", "year"], how="left").merge(parse_rx(), on=["fips", "year"], how="left")
    # Counties without a CDC rate (suppressed/unreported): use their state's median that year.
    cy["rx_rate"] = cy["rx_rate"].fillna(cy.groupby(["state", "year"])["rx_rate"].transform("median"))
    cy.to_csv(OUT / "county_year.csv", index=False)

    shp = gpd.read_file(DATA / "shapes" / "counties" / "cb_2022_us_county_20m.shp")
    shp = shp.rename(columns={"GEOID": "fips"})[["fips", "ALAND"]]
    counties = deaths.drop_duplicates("fips")[["fips", "state", "state_fips", "County"]].merge(shp, on="fips", how="left")
    counties["area_km2"] = counties["ALAND"] / 1e6
    counties.drop(columns="ALAND").to_csv(OUT / "counties.csv", index=False)

    key = dict(line.strip().split("=", 1) for line in (ROOT / ".env").read_text().splitlines() if "=" in line)["Census_API"]
    tracts = fetch_tracts(sorted(counties["state_fips"].unique()), key)
    tracts = tracts[tracts.fips.isin(counties.fips) & (tracts["pop"] > 0)]
    tracts.to_csv(OUT / "tracts.csv", index=False)

    print("county_year", cy.shape, "missing:", cy.isna().sum()[cy.isna().sum() > 0].to_dict())
    print("counties", counties.shape, "without area:", counties["area_km2"].isna().sum())
    print("tracts", tracts.shape, "counties covered:", tracts.fips.nunique())


if __name__ == "__main__":
    main()
