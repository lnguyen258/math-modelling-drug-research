# Overdose ABM experiment

Model description: [REPORT.md](REPORT.md) (layers, formulas, flow chart). Full theory and derivations: [GLOSSARY.md](GLOSSARY.md). Papers: [../reference/](../reference/README.md).

All commands are run from this `experiment/` folder with the project venv (`..\.venv\Scripts\python.exe`, written `python` below).

| Step | Command | What it does |
|---|---|---|
| 1 (once) | `python -m engine.prep` | Builds `../data/processed/` from `../data/`, and fetches census tracts once |
| 2 | `python -m engine.run --name baseline --N 300000 --M 3` | Runs an ensemble and saves it to `runs/<name>/` |
| 2b (optional) | `python -m engine.fit_lethality --run <name> --fit_until 2023 [--holdout Texas ...]` | Refits the 4 lethality/COVID parameters on a saved run, in seconds (REPORT section 6) |
| 3 | `python -m streamlit run app.py` | Opens the viewer in the browser. Pick a run in the sidebar |

Options for `engine.run`: `--N` agents, `--M` runs in the ensemble, `--seed`, `--rho` interaction radius in km (default 15). Every other parameter lives in `engine/params.py`.

Runtime: 300k agents takes about 65 s per run (24 simulated years, weekly steps); 20k agents about 25 s.

## Code map

| File | REPORT.md section |
|---|---|
| `engine/prep.py` | Layer 1: raw data → processed tables |
| `engine/inputs.py` | Layers 1–2: arrays, ξ indices, COVID stringency |
| `engine/environment.py` | Layer 2: modifiers F, prescriptions a, supply v, fentanyl φ |
| `engine/model.py` | Layers 3–4: rates λ, agents, neighbourhoods, weekly step |
| `engine/calibrate.py` | Section 5: averaged equations, calibration, elasticity matching, R_c |
| `engine/scoring.py` | Comparison to data (state-year interval likelihood) |
| `engine/fit_lethality.py` | Fit of m_F, p_phi, p_cov, w_iso on a saved run |
| `engine/params.py` | every parameter with its default |
| `app.py` | Streamlit viewer |
