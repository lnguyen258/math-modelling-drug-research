"""Run viewer.   streamlit run app.py      (from the experiment/ folder, with the project .venv)"""
import json
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

from engine import scoring
from engine.inputs import load_inputs

HERE = Path(__file__).resolve().parent
RUNS = HERE / "runs"
SHAPES = HERE.parent / "data" / "shapes" / "counties" / "cb_2022_us_county_20m.shp"

# Reference palette (dataviz skill): categorical slots 1-4, chrome, diverging pair
BLUE, ORANGE, AQUA, YELLOW = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
CLASS_COLORS = {"S": BLUE, "E": ORANGE, "I": AQUA, "R": YELLOW}
CLASS_NAMES = {"S": "Susceptible", "E": "Exposed (misuse)", "I": "Infectious (OUD)", "R": "Recovered (treatment)"}
GRID, MUTED = "#e1e0d9", "#898781"
DIVERGING = [[0, "#1c5cab"], [0.25, "#86b6ef"], [0.5, "#f0efec"], [0.75, "#f09a98"], [1, "#b42f2e"]]
SEQUENTIAL = [[0, "#cde2fb"], [0.5, "#3987e5"], [1, "#0d366b"]]


def rgba(hex_, a):
    h = hex_.lstrip("#")
    return f"rgba({int(h[0:2], 16)},{int(h[2:4], 16)},{int(h[4:6], 16)},{a})"


def style(fig, title, ytitle, height=380):
    fig.update_layout(title=dict(text=title, font=dict(size=15)), height=height, margin=dict(l=10, r=10, t=50, b=10),
                      hovermode="x unified", legend=dict(orientation="h", y=-0.15),
                      font=dict(family="system-ui, -apple-system, Segoe UI, sans-serif"))
    fig.update_xaxes(showgrid=False, linecolor="#c3c2b7")
    fig.update_yaxes(title=ytitle, gridcolor=GRID, zeroline=False, rangemode="tozero")
    return fig


@st.cache_resource
def inputs():
    return load_inputs()


@st.cache_data
def load_run(name):
    d = RUNS / name
    res = dict(np.load(d / "results.npz"))
    return res, json.loads((d / "params.json").read_text()), json.loads((d / "summary.json").read_text())


@st.cache_data
def county_geojson(fips):
    g = gpd.read_file(SHAPES)
    g = g[g["GEOID"].isin(fips)][["GEOID", "geometry"]].to_crs(4326)
    return json.loads(g.to_json())


st.set_page_config(page_title="Overdose ABM viewer", layout="wide")
inp = inputs()
runs = sorted([p.name for p in RUNS.iterdir() if (p / "results.npz").exists()]) if RUNS.exists() else []
if not runs:
    st.warning("No runs yet. From experiment/: `python -m engine.run --name baseline`")
    st.stop()

with st.sidebar:
    st.header("Run")
    name = st.selectbox("Run", runs, index=runs.index("v2_baseline") if "v2_baseline" in runs else 0)
    geo = st.selectbox("Geography", ["All 17 states"] + inp.states)
    measure = st.radio("Model deaths", ["Expected (smooth)", "Realized (counted)"],
                       help="Expected = sum of each OUD agent's overdose hazard. Realized = deaths that actually "
                            "fired. Same mean; realized is noisier.")

res, params, summary = load_run(name)
years = inp.years
sel = np.ones(len(inp.fips), bool) if geo == "All 17 states" else inp.state_of == inp.states.index(geo)

deaths_key = "hazard_cy" if measure.startswith("Expected") else "deaths_cy"
Ys = res[deaths_key].shape[-1]
model_years = params["t0"] + np.arange(Ys)                 # may extend past the data years
Yd = len(years)
# population beyond the data is held at the last data year
pop_ext = np.concatenate([inp.pop, np.repeat(inp.pop[:, -1:], max(Ys - Yd, 0), axis=1)], axis=1)[:, :Ys]
exp_real = scoring.expected_real_deaths(res[deaths_key], res["n_c"], pop_ext)        # [M, C, Ys]
has_agents = res["n_c"] > 0
pop_sel = pop_ext[sel & has_agents].sum(0)
model_rate = 1e4 * np.nansum(exp_real[:, sel & has_agents, :], axis=1) / pop_sel     # [M, Ys]
lo, hi = scoring.observed_band(inp, np.where(sel, 0, -1), 1)
obs_lo = 1e4 * lo[0] / inp.pop[sel].sum(0)
obs_hi = 1e4 * hi[0] / inp.pop[sel].sum(0)

# ---- summary row ----
st.title("Overdose ABM: run viewer")
st.caption("Opioid overdose deaths per 10,000 people. Model = spatial stochastic SEIRD agents "
           "(see REPORT.md); data = SRDDRC / CDC WONDER county-years.")
mean_rate = model_rate.mean(0)
inside = np.mean((mean_rate[:Yd] >= obs_lo) & (mean_rate[:Yd] <= obs_hi))
c1, c2, c3, c4 = st.columns(4)
c1.metric("Agents", f"{summary['agents']:,}", f"{summary['members']} runs", delta_color="off")
c2.metric("Model years inside observed range", f"{inside:.0%}", geo, delta_color="off")
c3.metric("Contact reproduction number R_c", f"{summary['R_c']:.2f}")
fit = summary.get("state_loglik")
c4.metric("State fit (log-lik per state-year)", "n/a (old run)" if fit is None else f"{fit:.2f}",
          "higher is better", delta_color="off")

tab_d, tab_c, tab_m, tab_e, tab_p = st.tabs(["Deaths vs data", "Compartments", "Map", "Drug environment", "Parameters"])

with tab_d:
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=np.r_[years, years[::-1]], y=np.r_[obs_hi, obs_lo[::-1]], fill="toself",
                             fillcolor=rgba(ORANGE, 0.18), line=dict(width=0), hoverinfo="skip",
                             name="Observed range (hidden counts = 1 to 9)"))
    fig.add_trace(go.Scatter(x=years, y=obs_lo, mode="lines+markers", line=dict(color=ORANGE, width=2, dash="dot"),
                             marker=dict(size=8), name="Observed, hidden counts as 1"))
    fig.add_trace(go.Scatter(x=years, y=obs_hi, mode="lines+markers", line=dict(color=ORANGE, width=2, dash="dot"),
                             marker=dict(size=8, symbol="diamond"), name="Observed, hidden counts as 9"))
    if model_rate.shape[0] > 1:
        fig.add_trace(go.Scatter(x=np.r_[model_years, model_years[::-1]],
                                 y=np.r_[model_rate.max(0), model_rate.min(0)[::-1]],
                                 fill="toself", fillcolor=rgba(BLUE, 0.18), line=dict(width=0), hoverinfo="skip",
                                 name="Model, range across runs"))
    fig.add_trace(go.Scatter(x=model_years, y=mean_rate, mode="lines+markers", line=dict(color=BLUE, width=2),
                             marker=dict(size=8), name="Model mean"))
    if Ys > Yd:
        fig.add_vline(x=years[-1] + 0.5, line=dict(color=MUTED, width=1, dash="dash"))
        fig.add_annotation(x=years[-1] + 0.5, y=1, yref="paper", text="projection →", showarrow=False,
                           xanchor="left", font=dict(color=MUTED))
    st.plotly_chart(style(fig, f"Opioid deaths per 10,000: model vs data, {geo}", "deaths per 10,000", 460),
                    width="stretch")
    st.caption("CDC hides county counts of 1–9, so the real value lies somewhere in the orange band. "
               "COVID (2020–21) is not in the model on purpose; a gap there is expected.")
    pad = lambda a: np.r_[a, np.full(max(Ys - Yd, 0), np.nan)][:Ys]
    table = pd.DataFrame({"year": model_years, "model mean": mean_rate.round(2),
                          "model min": model_rate.min(0).round(2), "model max": model_rate.max(0).round(2),
                          "observed low": pad(obs_lo).round(2), "observed high": pad(obs_hi).round(2)})
    with st.expander("Table"):
        st.dataframe(table, hide_index=True)

with tab_c:
    so = inp.state_of
    series = res["series"]                                   # [M, T, states, 4]
    if geo == "All 17 states":
        cnt = series.sum(2)
    else:
        cnt = series[:, :, inp.states.index(geo), :]
    frac = 100 * cnt / cnt.sum(-1, keepdims=True)            # [M, T, 4]
    t = res["times"]
    cols = st.columns(2)
    for k, X in enumerate("SEIR"):
        f = go.Figure()
        if frac.shape[0] > 1:
            f.add_trace(go.Scatter(x=np.r_[t, t[::-1]], y=np.r_[frac[:, :, k].max(0), frac[:, :, k].min(0)[::-1]],
                                   fill="toself", fillcolor=rgba(CLASS_COLORS[X], 0.2), line=dict(width=0),
                                   hoverinfo="skip", showlegend=False))
        f.add_trace(go.Scatter(x=t, y=frac[:, :, k].mean(0), line=dict(color=CLASS_COLORS[X], width=2),
                               name=CLASS_NAMES[X], hovertemplate="%{y:.3f}%"))
        style(f, f"{X}: {CLASS_NAMES[X]}", "% of population", 300)
        f.update_layout(showlegend=False)
        if t[-1] > years[-1] + 1:
            f.add_vline(x=years[-1] + 1, line=dict(color=MUTED, width=1, dash="dash"))
        f.update_yaxes(rangemode="normal")
        cols[k % 2].plotly_chart(f, width="stretch")
    st.caption("Weekly class shares, mean of runs with the range across runs shaded. Each panel has its own scale.")

with tab_m:
    yr = st.select_slider("Year", options=list(years), value=int(years[-1]))
    k = int(np.where(years == yr)[0][0])
    m_mean = exp_real.mean(0)[:, k]   # data years start at t0, so index k matches
    with np.errstate(divide="ignore", invalid="ignore"):
        model_c = 1e4 * m_mean / inp.pop[:, k]
        obs_c = 1e4 * inp.deaths_obs[:, k] / inp.pop[:, k]
        ratio = np.log2((model_c + 0.05) / (obs_c + 0.05))
    view = st.radio("Show", ["Model ÷ observed (log2)", "Model rate", "Observed rate"], horizontal=True)
    min_agents = st.slider("Hide counties with fewer agents than", 0, 200, 50,
                           help="Counties with few agents are mostly simulation noise.")
    ok = (res["n_c"] >= min_agents)
    df = pd.DataFrame({"fips": inp.fips, "state": np.array(inp.states)[inp.state_of], "agents": res["n_c"],
                       "model": model_c, "observed": obs_c, "log2 ratio": ratio,
                       "hidden": inp.suppressed[:, k]})
    if geo != "All 17 states":
        ok &= df["state"].eq(geo).to_numpy()
    if view.startswith("Model ÷"):
        df = df[ok & ~df["hidden"] & np.isfinite(df["log2 ratio"])]
        z, cs, zr, cb = df["log2 ratio"], DIVERGING, (-3, 3), "log2(model/obs)"
    else:
        col = "model" if view == "Model rate" else "observed"
        df = df[ok & np.isfinite(df[col])]
        z, cs, zr, cb = df[col], SEQUENTIAL, (0, float(np.nanpercentile(obs_c, 99))), "per 10,000"
    fig = go.Figure(go.Choropleth(geojson=county_geojson(tuple(inp.fips)), featureidkey="properties.GEOID",
                                  locations=df["fips"], z=z, colorscale=cs, zmin=zr[0], zmax=zr[1],
                                  marker_line_width=0.2, marker_line_color="#fcfcfb",
                                  colorbar=dict(title=cb),
                                  customdata=np.c_[df["state"], df["agents"], df["model"].round(2),
                                                   df["observed"].round(2)],
                                  hovertemplate="%{location} %{customdata[0]}<br>agents %{customdata[1]}"
                                                "<br>model %{customdata[2]}<br>observed %{customdata[3]}<extra></extra>"))
    fig.update_geos(fitbounds="locations", visible=False)
    fig.update_layout(height=560, margin=dict(l=0, r=0, t=10, b=0))
    st.plotly_chart(fig, width="stretch")
    st.caption("Ratio view: red = model too high, blue = model too low, gray = matches. "
               "Counties whose observed count is hidden (1–9) are left out of the ratio view.")

with tab_e:
    w = pop_ext[sel]
    a_ext = np.concatenate([inp.a_rx, np.repeat(inp.a_rx[:, -1:], max(Ys - Yd, 0), axis=1)], axis=1)[:, :Ys]
    def popmean(a):  # a: [M, C, Ys] -> [M, Ys]
        return (a[:, sel, :] * w).sum(1) / w.sum(0)
    cols = st.columns(3)
    for col, (key, title, yt) in zip(cols, [("v_cy", "Illicit supply index v (1 = start level)", "index"),
                                            ("phi_cy", "Fentanyl share of supply φ", "share"),
                                            (None, "Prescriptions per person a (data)", "Rx per person")]):
        f = go.Figure()
        if key:
            vals = popmean(res[key])
            if vals.shape[0] > 1:
                f.add_trace(go.Scatter(x=np.r_[model_years, model_years[::-1]],
                                       y=np.r_[vals.max(0), vals.min(0)[::-1]],
                                       fill="toself", fillcolor=rgba(BLUE, 0.2), line=dict(width=0),
                                       hoverinfo="skip", showlegend=False))
            f.add_trace(go.Scatter(x=model_years, y=vals.mean(0), line=dict(color=BLUE, width=2),
                                   mode="lines+markers", marker=dict(size=8), showlegend=False))
        else:
            f.add_trace(go.Scatter(x=model_years, y=(a_ext[sel] * w).sum(0) / w.sum(0), mode="lines+markers",
                                   line=dict(color=BLUE, width=2), marker=dict(size=8), showlegend=False))
        col.plotly_chart(style(f, title, yt, 320), width="stretch")
    st.caption("v and φ are hidden states (no supply data yet), population-weighted over the selected area. "
               "φ follows a Bass diffusion from 2013; the fentanyl death hazard is μ_H·(1−φ) + μ_F·φ.")

with tab_p:
    flat = {k: v for k, v in params.items() if k != "weights"}
    st.dataframe(pd.DataFrame({"parameter": list(flat), "value": [str(v) for v in flat.values()]}),
                 hide_index=True, height=500)
    st.subheader("Structural weights (per SD)")
    st.dataframe(pd.DataFrame(params["weights"]).T.fillna(0.0))
    st.caption("Read-only for now. Editing parameters from this page comes next.")
