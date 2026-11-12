# Project context (running log)

This is a living summary of the project, updated after each work session. Newest state first, then the history.

## 1. Current state (2026-11-05)

**Goal:** a spatial stochastic agent-based model (SEIRD) of opioid overdose deaths across 17 southern states. Agents live on real population density. Transition rates are a literature baseline × an environment term built from Census and drug-availability data. The model is checked against observed county-year deaths for 2010–2023.

**Folder layout**

| Path | What it is |
|---|---|
| `experiment/REPORT.md` | **Main model description (v2)**: layers Data → Derived environment → Rates λ → Agents, a flow chart, results, limitations |
| `experiment/GLOSSARY.md` | Full theory and derivations (v1–v2 design). Kept on request; does not yet describe the v2 changes |
| `experiment/README.md` | How to run (prep → run → fit → viewer) and the code map |
| `experiment/engine/` | Code: `prep`, `inputs`, `environment`, `model`, `calibrate`, `scoring`, `fit_lethality`, `run`, `params` |
| `experiment/app.py` | Streamlit viewer (`python -m streamlit run app.py` from `experiment/`, http://localhost:8501); default run `v2_baseline` |
| `experiment/runs/` | `v2_baseline`, `v2_projection_2060`, `v3_N1.5M_rho30` (1.5M agents, 30 km) (current); `v2_prefit` (source of the lethality fit); `baseline`, `projection_2060` (v1, for comparison) |
| `data/` | SRDDRC deaths, ACS 2012–2023, shapes, CDC Rx pages, tracts, **`oxcgrt/` COVID stringency**, `processed/` |
| `reference/` | 17 papers plus `README.md` index (tags R01–R19) |
| `outputs/` | Figures, incl. `v2_deaths_by_state.png`, `v2_compartments_2060.png` |

**Model v2, in short** (details in REPORT.md)

- No density factor (η removed).
- Structural weights set by **elasticity matching**: prior signs kept, sizes scaled so +1 SD gives a literature-sized change in long-run deaths (econ +15%, pain +20%, unins +10%; age left at its priors).
- COVID enters as state stringency (OxCGRT) through two channels: a permanent push on fentanyl takeover ($p_{cov}$) and a temporary isolation effect on overdose ($w_{iso}$).
- Four parameters are **fitted to deaths** at the state-year level, without rerunning agents: $m_F$ = 3.88, $p_arphi$ = 0.0128, $p_{cov}$ = 1.33, $w_{iso}$ = 0.14. Holding out 5 states gives the same values and improves the score on those states.
- Score: interval-censored negative binomial per state-year (county-level scores are unusable at 300k agents).

**Results v2**

- Fit score −2.44 per state-year (v1: −70.3).
- No runaway states; none dies out by 2060. Region-wide OUD goes 0.54% (2023) → 0.19% (2059).
- **Finding:** the COVID jump is explained by an accelerated, permanent fentanyl takeover, not by isolation (isolation weight drops from 0.98 to 0.14 once 2022–23 are included).

**Known problems (open)**

1. **States are too alike.** 2021 model values span 1.8–3.8 against an observed 1.0–4.9. Texas is still about 2.5× too high; Alabama and Arkansas are high; Kentucky, Tennessee, South Carolina and New Mexico are too low. Likely cause: the real timing of fentanyl arrival by state, which we have no data for.
2. Region-wide 2018–20 is slightly above the observed range.
3. Four parameters are fitted on deaths, so 2020–23 is not a clean out-of-sample test. The held-out-state test is the check.
4. County-level detail is noise at 300k agents.

**Next steps (proposed)**

1. Fentanyl timing by state from data. **Still asking the user:** does SRDDRC have synthetic-opioid (T40.4) deaths?
2. Parameter editing in the viewer.
3. Optionally add a "v2 changes" section to GLOSSARY.md (asked the user).

**Decisions made by the user**

- Interaction radius ρ = **15 km**.
- Structural weights from sign priors (now sized by elasticity matching); to be made dynamic later.
- Keep R→S and all nine transitions.
- COVID: first left in the error; on 2026-11-05 the user asked to include it, which is done.
- Drug-availability precision is deferred.
- Keep GLOSSARY.md; REPORT.md is added next to it, not a replacement.
- Parameter editing in the viewer is wanted later.

## 2. History

### Session 1: earlier correlation project (before this log)
The first project goal was a correlation heatmap of overdose death rates vs. demographic indicators. It produced Spearman, GWR and smoothed-rate maps in `outputs/`. Its findings, kept as context:

- Disability (about +0.3) and median age (about +0.3 to +0.4) are positively correlated with death rates.
- Uninsured rate is negatively correlated (about −0.25).
- Poverty is weak overall.

The notebook, `src/`, `NOTES.md` and logs from that project were **deleted on request** when the project moved to the ABM.

### Session 2: ABM theory
1. Created `experiment/GLOSSARY.md` v1: simple SEIRD with Poisson clocks and a log-linear environment link.
2. The user asked for a holistic, novel design built on literature. I pulled papers into `reference/`, including Battista 2019, White & Comiskey 2007, Ruhm, Hollingsworth, Monnat, Ciccarone, Sordo (BMJ), Mednick, Deardon, Tiwari, Zibbell. Rewrote the glossary as layers L0–L7:
   - raw data → structural indices ξ → dynamic drug availability (prescriptions, a hidden illicit-supply variable driven by demand, fentanyl share φ via spatial Bass diffusion)
   - → space (agents on tract density, radius kernel)
   - → 9 rates (baseline × availability × exp(w·ξ) + peer terms)
   - → stochastic engine → censored negative-binomial observation model → mean-field checks (R_c, equilibrium)
   - Plus "elasticity matching", which sets the weights from published regression effects.
3. Checked `data/`. Deaths are county-year 2010–2023 for 17 states and 1,538 counties, with about 50% suppressed (counts 1–9). ACS is yearly 2012–2023. Adjusted the indices to the available columns: econ, pain (disability + veterans), uninsured, age; coh and tx are 0 for now.

### Session 3: engine v1 + viewer
1. **Data pulled:**
   - CDC county opioid dispensing rates 2010–2018 from the CDC archive. The 2019 county page was rejected (looks like a method change). 2019–2023 are extended by the national trend, with 2021 interpolated.
   - 2020 census tracts: 33,696 tracts, population and centroids.
   - NFLIS skipped (PDF only).
2. **Engine built:**
   - Agents placed on tracts; 15 km kernel on a 5 km grid; pseudo-count fallback k0 = 20.
   - Weekly competing-risk steps; deaths reborn as S in the same county.
   - β⁰ and χ calibrated to the 0.6% OUD and 0.66%/yr new-misuse targets.
   - σ⁰ lowered from 0.15 to 0.08, because 0.15 made the calibration infeasible.
3. **Found:** supply–demand feedback has a tipping point at θ_d ≈ 0.2 (glossary §12.5). Set θ_d = 0.1.
4. **Baseline run** (300k agents × 3 runs, about 65 s each). The region-wide model mean is inside the observed range in 11 of 14 years and too high in 2017–19. Tennessee is in range 0% of years.
5. **Viewer** (Streamlit) tabs: deaths vs data, compartments, county map, drug environment, parameters.
6. **Projection to 2060.** Region-wide nothing dies out: OUD goes 0.64% → 0.33%, which matches the mean-field theory (no drug-free equilibrium while prescriptions exist). By state, Texas and Florida grow while the others fade, which revealed the η artifact.
7. **COVID question.** Per-state check showed the regional fit is error cancellation (see problems 1–2). Proposed the plan above.

### Session 4 (2026-11-05): Texas fix, COVID, REPORT.md
1. Created this CONTEXT.md at the user's request.
2. **Texas:** η = 0 alone did *not* fix it (Texas 3.7, Mississippi 6.6). The real cause was the weight sizes. The system sits near the threshold (R_c 0.87 on average), and the prior weights implied econ ×4.5 per SD, putting 26% of the population in counties with R_c > 1. Elasticity matching shrank econ about 10×, pain to 0.7× and unins to 0.33×. Now 4% of the population has R_c > 1, and there is no runaway in any state.
3. **COVID:** pulled OxCGRT state stringency for 2020–22 (`data/oxcgrt/`). Added the isolation effect and the fentanyl push. Wrote `engine/fit_lethality.py`, which fits from saved OUD counts in seconds.
   - Fitting on years ≤ 2021 overfit the isolation effect and failed on 2022–23.
   - Fitting all years selected the permanent supply push.
   - Holding out 5 states confirmed the fit.
4. Switched scoring to the state-year interval likelihood and removed the county score.
5. **Simplified** (inactive parts only): data nudging, cohesion/treatment indices, hub exposure, fentanyl spatial term, density factor.
6. Wrote `experiment/REPORT.md`: layered structure, mermaid flow chart, input→rate matrix, parameter status table (L/A/C/F), results, limitations. Briefly moved GLOSSARY.md to an archive; restored it on request.
7. Runs: `v2_baseline` and `v2_projection_2060` (300k × 3). Figures in `outputs/`. Viewer restarted with v2 as the default.

### Session 5 (2026-11-05, continued): glossary v2.1, presentation, scale test
1. **GLOSSARY.md updated to v2.1.** The text now matches the code. Updated: indices (4), illicit supply (no nudging), fentanyl Bass + COVID push (§5.3), new §5.4 COVID, density factor removed (§6.4), I→D with isolation, §9.2 state-interval score + the 4 fitted parameters + held-out tests, §10 prior vs matched weights and how matching is implemented, §11 parameter statuses (F = fitted), §14 open questions, new §15 changelog.
2. **Presentation:** 13-slide deck, https://claude.ai/artifact/RFJMf4i4gXhucX8VcR4RZi (private until shared). Covers: question → SEIRD states → four layers → data → indices/elasticity matching → drug availability → rates table → agents → fitting/COVID finding → results figure → 2060 projection → limitations. Local source: scratchpad `deck/project/`.
3. **Scale and distance test** (user chose **30 km**): run `v3_N1.5M_rho30` = 1.5M agents, ρ = 30 km, M = 3, all other parameters as v2 (not refitted). Old `v2_baseline` kept.
   - Region unchanged (2021: 2.54 vs 2.58).
   - Run-to-run noise 0.172 → 0.075 (the √5 predicted).
   - States pulled together (2021 range 1.77–3.78 → 1.98–3.07; AL −0.85, SC +0.70, both toward the data).
   - Fit −2.44 → −2.23; years in range 35% → 34%.
   - Texas still too high, KY/TN still too low.
   - Caveat: N and ρ changed together; a 1.5M / 15 km run would separate them.
   - Added: REPORT §7.4, deck slide `scale` (14 slides now), `outputs/v3_compare_scale_radius.png`. Viewer lists the run.

### Session 6 (2026-11-12): cleanup, GitHub push, math slides
1. **Cleanup:** `outputs/` held 40+ files from the earlier correlation project (GWR, heatmaps, trend/period maps, smoothed rates). They do not belong to the ABM work, so they were moved out of the project (not deleted) to a scratch backup. `outputs/` now keeps only the 5 ABM figures.
2. **Pushed** the project to https://github.com/lnguyen258/math-modelling-drug-research.git (branch main). Excluded: `.venv/`, `.env` (Census API key), LaTeX build files, the compiled slide PDF.
3. **Math slides:** `experiment/slides/model_slides.tex` (Beamer, 46 pages). Every symbol with a legend; redundancies merged (red boxes) or kept with the reason (green boxes). From the raw-variable slide onward, each computation has a worked numeric example (Barbour Co., AL 2019 for the census chain; calibration, R_c, scoring numbers from the engine). Renamed vs. the code to avoid clashes: radius rho -> h, stringency S -> psi, supply noise sigma_v -> tau_v, spillover l_W -> l. Build: `pdflatex` twice (latexmk needs perl, not installed).
