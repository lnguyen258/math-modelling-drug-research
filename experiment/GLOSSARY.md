# Spatial Stochastic Agent-Based Model of Overdose Mortality: Glossary and Theory

Version **2.1** (2026-11-05). This matches the code in `engine/`. Section 15 lists what changed from version 2. A shorter, presentation-style description is in [REPORT.md](REPORT.md). Citations like [R01] point to [../reference/README.md](../reference/README.md).

The document is built in **layers**. Every symbol is defined only from symbols in earlier layers:

```
L0 Raw data  →  L1 Structural indices  →  L2 Drug availability (dynamic)
     →  L3 Space (agents, kernel, local pressure)  →  L4 Transition rates
     →  L5 Stochastic engine  →  L6 Observation & comparison  →  L7 Analytical checks
```

---

## 0. What is novel here

Most opioid models are well-mixed ODEs at national scale [R01, R02, R08]. The review in [R11] finds that only 16% are agent-based, and almost none place agents on real geography. This design combines four things that are not usually combined:

1. **Agents placed on real population density** across the southern states, interacting only within a fixed radius ρ. The model becomes a spatial individual-level model [R14, R15] applied to drug use instead of infection.
2. **Drug availability is a dynamic, partly hidden state, not a fixed covariate.** It has two channels (prescription supply and illicit supply), plus a potency variable (the fentanyl share) that spreads between counties like a new product in a market (Bass diffusion). Illicit supply is pulled by local demand, so supply and demand form a feedback loop [R04, R09, R10]. COVID policy stringency speeds up the fentanyl switch (§5.4).
3. **Every rate is a literature baseline times an environmental modifier:** θ = θ⁰ · f(environment). The baselines come from cohort and modelling studies. The size of each modifier is set so that the model's *implied elasticity* matches published regression effects ("elasticity matching", section 10).
4. **One observable (deaths) is decomposed into nine flows.** The model makes unobserved flows (initiation, relapse, treatment entry) explicit and testable.

---

## 1. Conventions

| Notation | Meaning |
|---|---|
| $c$ | County index (finest geography). $c \in \mathcal{C}$, the set of counties in the study states |
| $s(c)$ | State that contains county $c$ |
| $t$ | Continuous time in **years**. $y = \lfloor t \rfloor$ is the calendar year |
| $i, j$ | Agent indices, $i = 1..N$ |
| $X_i(t)$ | Compartment of agent $i$ at time $t$, $\in \{S, E, I, R\}$ |
| $\mathbb{1}[\cdot]$ | Indicator function (1 if true, else 0) |
| $\operatorname{clip}(u, a, b)$ | $\min(\max(u, a), b)$ |
| superscript $^0$ | A baseline value: the rate in an "average" county with average drug availability |
| $\tilde{\ }$ | A quantity after environmental modification |

All rates are **per person per year**. "Rate λ" means: a Poisson clock that rings with probability ≈ λ·dt in a short time dt.

---

## 2. Compartments (the Civilian agent)

| State | Name | Precise definition used for calibration | Recruits others? |
|---|---|---|---|
| **S** | Susceptible | No current non-medical opioid use | – |
| **E** | Exposed | Current non-medical use (misuse) without opioid use disorder (OUD) | Weakly (weight $\vartheta_E$) |
| **I** | Infectious | Active OUD, not in treatment | Yes (weight 1) |
| **R** | Recovered | In treatment (medication for OUD or other) or in early remission | No; drives peer recovery |
| **D** | Dead | Opioid overdose death. A **counter**: the agent is immediately reborn as S (§5.3) | – |

The nine allowed transitions:

```
          S→E  initiation by contact            E→S  desistance
          S→I  iatrogenic onset (prescriptions) E→I  progression to OUD
          I→R  treatment entry                   I→S  natural remission
          I→D  overdose death (→ reborn S)       R→S  sustained recovery
          R→I  relapse
```

---

## 3. Layer 0: raw inputs (data)

All county-level and annual unless noted. "Dyn" means it changes every year, and the engine linearly interpolates between years.

### 3.1 Socio-economic (U.S. Census)

**In hand** (`data/acs_{2012..2023}_county_yearly.csv`, ACS 5-year vintage ending in year y, 1,538 counties). Column names are as in the files:

| Symbol | File column | Definition | Unit |
|---|---|---|---|
| $x^{pov}_c(y)$ | `poverty_pct` | Persons below the poverty line | % |
| $x^{unemp}_c(y)$ | `unemployment_pct` | Unemployment rate of the civilian labor force | % |
| $x^{inc}_c(y)$ | `median_hh_income` | Median household income | USD (we use log) |
| $x^{nohs}_c(y)$ | `edu_lt9th_pct + edu_9to12_nodiploma_pct` | Adults 25+ without a high-school diploma | % |
| $x^{bach}_c(y)$ | `edu_bachelors_plus_pct` | Adults 25+ with a bachelor's degree or more | % |
| $x^{dis}_c(y)$ | `disability_pct` | Civilian population with a disability | % |
| $x^{unins}_c(y)$ | `uninsured_pct` | Population without health insurance | % |
| $x^{vet}_c(y)$ | `veteran_pct` | Adult veterans (high chronic pain and prescription exposure) | % |
| $x^{age}_c(y)$ | `median_age` | Median age | years |
| $x^{dens}_c$ | derived | Population per km² of land: SRDDRC population / `ALAND` from `data/shapes/counties` | 1/km² |

Years 2010–2011 have no ACS county tables with these definitions. The engine holds them at the 2012 values.

**Not in hand (optional, I can fetch them with the Census API key in `.env`).** They would feed a cohesion index, which was removed in 2.1 until they are fetched (§4.2).

| Symbol | Definition | Source |
|---|---|---|
| $x^{single}_c(y)$ | Single-parent households among family households (Monnat's "family distress") | ACS DP02 |
| $x^{relig}_c(y)$ | Religious organizations per 10,000 residents (Monnat's protective factor) | County Business Patterns, NAICS 813110 |

Crime stays out. Census has no crime data, and NIBRS on SRDDRC is state-level only.

### 3.2 Drug environment (dynamic, region-specific)

| Symbol | Code name | Definition | Unit | Source | Granularity |
|---|---|---|---|---|---|
| $x^{Rx}_c(y)$ | `rx_rate` | Opioid prescriptions dispensed per 100 residents | Rx/100/yr | CDC opioid dispensing rate maps. **In hand for 2010–2018** (`data/rx_cdc_html/`, from the CDC archive). The 2019 page is not used: its county values halve vs 2018 while the national rate fell 9%, which points to a method change. 2019–2023 = county 2018 value × national rate ratio. Counties without a value use their state median | County |
| $x^{op}_s(y)$ | `x_nflis_op` | Opioid drug reports in forensic labs per 100,000 | /100k | DEA NFLIS-Drug public data | **State** |
| $x^{fent}_s(y)$ | `x_nflis_fent` | Fentanyl and analogs as a share of opioid reports | [0,1] | NFLIS-Drug | **State** |
| $x^{tx}_c(y)$ | `x_tx` | Treatment capacity: opioid treatment programs + buprenorphine prescribers per 100,000 | /100k | SAMHSA treatment locator / N-SUMHSS | County |
| $x^{str}_s(t)$ | `stringency` | **COVID policy stringency index**, daily, 0–100 | index | **In hand:** Oxford COVID-19 Government Response Tracker (OxCGRT), `data/oxcgrt/` (2020–2022; 0 outside) | **State** |

**Status: $x^{Rx}$ and $x^{str}$ are in hand; NFLIS and treatment capacity are not.** NFLIS is published as PDF reports (not machine-readable), and treatment capacity is deferred. Illicit supply and the fentanyl share are therefore **hidden dynamic states** (§5). Their timing and lethality are fitted to deaths (§9.2).

### 3.3 Population and geography

| Symbol | Code name | Definition | Source |
|---|---|---|---|
| $P_c(y)$ | `pop_c` | Resident population of county $c$ | **In hand:** `Population` column of the SRDDRC files |
| — | `ALAND`, geometry | County land area and polygon | **In hand:** `data/shapes/counties/cb_2022_us_county_20m.shp` |
| $P_b,\ \mathbf{b}_b$ | `pop_bg`, `xy_bg` | Population and centroid (projected km, Albers EPSG:5070) of block group $b$ | Census 2020 P1 + TIGER Gazetteer, **I fetch** |
| $A_b$ | `area_bg` | Land area of block group $b$ | TIGER Gazetteer, **I fetch** |

The study region has 17 states (AL, AZ, AR, CO, FL, GA, KS, KY, LA, MS, MO, NM, NC, OK, SC, TN, TX), 1,538 counties, and $P_{tot} \approx$ 136.3 M people (2023).

### 3.4 The observable

| Symbol | Code name | Definition | Source |
|---|---|---|---|
| $Y_c(y)$ | `Deaths` | Opioid overdose death count, county-year, 2010–2023 | **In hand:** `data/SRDDRC CDCWONDER_Opioid_YrStCty.csv` (all-drug file kept for comparison) |
| $\mathcal S$ | `Deaths == "Suppressed"` | Set of county-years where CDC hides counts 1–9 (≈ 50% of rows). Treated as **censored** $Y \in [1, 9]$, never imputed (§9.2) | – |
| $m^{obs}(y)$ | `m_obs` | $10^4 \cdot Y / P$, deaths per 10,000 | Derived |

Scale check from the visible counts across the study region: the opioid death rate was about 0.52 per 10k in 2010, 1.0 in 2017, 2.1 in 2021–22 and 1.95 in 2023. Suppressed rows are excluded from these, so they understate the truth. The 2014 and 2020 jumps line up with fentanyl entry and COVID.

---

## 4. Layer 1: structural indices (static or slow)

### 4.1 Standardization

For each raw variable $k$, pool all counties and all years, then compute

$$
z^{k}_c(y) = \operatorname{clip}\!\left(\frac{x^k_c(y) - \bar{x}^k}{s^k},\ -3,\ 3\right)
$$

where $\bar x^k, s^k$ are the pooled mean and standard deviation. Pooling over years (instead of standardizing year by year) keeps real time trends in the signal. Clipping at ±3 SD prevents one outlier county from producing an explosive rate. For income and density, use $\log x$ before standardizing.

### 4.2 Composite indices ξ

Ten correlated covariates would make the weights unidentifiable (§9.3). We group them into four **structural indices**, each an equal-weight signed mean of z-scores, re-standardized to mean 0 and SD 1 (operator $\operatorname{std}[\cdot]$):

| Index | Code name | Formula | Meaning | Literature |
|---|---|---|---|---|
| $\xi^{econ}_c$ | `xi_econ` | $\operatorname{std}\big[\tfrac15(z^{pov} + z^{unemp} + z^{nohs} - z^{\log inc} - z^{bach})\big]$ | Economic precarity | [R05, R06] |
| $\xi^{pain}_c$ | `xi_pain` | $\operatorname{std}\big[\tfrac12(z^{dis} + z^{vet})\big]$ | Pain and injury burden, which drives prescription exposure | [R01, R06] |
| $\xi^{unins}_c$ | `xi_unins` | $z^{unins}$ | Lack of coverage (barrier to treatment) | [R06] |
| $\xi^{age}_c$ | `xi_age` | $z^{age}$ | Population age structure | Earlier correlation work on this data |
Vector form: $\boldsymbol\xi_c(t) = (\xi^{econ}, \xi^{pain}, \xi^{unins}, \xi^{age})_c(t)$.

*Removed in 2.1:* the cohesion index ($\xi^{coh}$, from religious organizations and single-parent households) and the treatment-access index ($\xi^{tx}$). Both were identically 0 because the data were never fetched. The urbanicity index $\xi^{urb}$ only fed the density factor, which was also removed (§6.4). They can come back once the data exist.

Uninsurance is kept separate from disability because the earlier county correlation work on this data found the two pull in **opposite** directions in this region (disability about +0.3, uninsured about −0.25). Every rate uses a **modifier**

$$
F_\theta(c,t) = \exp\!\big(\mathbf{w}_\theta \cdot \boldsymbol\xi_c(t)\big)
$$

with $\mathbf{w}_\theta$ the weight row for rate θ (§10). Because ξ has mean 0, $F_\theta = 1$ in the average county, so the baseline θ⁰ keeps its literature meaning. A weight $w$ means "one SD of this index multiplies the rate by $e^{w}$" (e.g. $w = 0.2$ gives ×1.22).

*Alternative (later):* replace the equal weights with the first principal component of each group, or estimate them from the heatmap regression.

---

## 5. Layer 2: drug availability, the dynamic environment

Following the triple-wave framework [R09], availability has **three** components. Each enters different transitions.

| Symbol | Code name | Meaning | Type | Enters |
|---|---|---|---|---|
| $a_c(t)$ | `a_rx` | Prescription availability: prescriptions per person per year | Observed (exogenous) | S→I |
| $v_c(t)$ | `v_ill` | Illicit supply volume index (1 = study-wide average) | **Hidden**, dynamic | S→E, E→I, R→I |
| $\varphi_c(t)$ | `phi` | Potency: fentanyl share of the illicit opioid supply, ∈ [0,1] | **Hidden**, dynamic | I→D |

### 5.1 Prescription channel (observed)

$$
a_c(t) = x^{Rx}_c(t) / 100
$$

### 5.2 Illicit volume: demand-pull stochastic dynamics

Illicit supply follows local demand with a lag, spills over between neighboring counties (trafficking), and carries market noise:

$$
dv_c = \Big[\underbrace{\kappa_v\big(v^{*}_c - v_c\big)}_{\text{demand pull}} + \underbrace{D_v \textstyle\sum_{c'} W_{cc'}(v_{c'} - v_c)}_{\text{spatial spillover}}\Big]dt + \sigma_v\, v_c\, dB_c(t)
$$

(Version 2 also had a data-nudging term $\kappa_{obs}(\hat v - v)$ toward NFLIS reports. With no NFLIS data it was always 0, so 2.1 drops it.)

with demand target (county demand is first shrunk toward the regional mean with $k_0$ pseudo-agents, like §6.4, so tiny counties do not swing on one agent)

$$
v^{*}_c(t) = \left(\frac{i_c(t) + \vartheta_E\, e_c(t)}{\overline{i + \vartheta_E e}}\right)^{\theta_d}
$$

| Symbol | Code name | Meaning |
|---|---|---|
| $i_c, e_c$ | `i_frac`, `e_frac` | Fractions of county agents in I and E (from L5, so this is the **feedback loop**) |
| $\theta_d$ | `theta_d` | Demand elasticity of supply: 0 = supply ignores demand, 1 = proportional |
| $\kappa_v$ | `kappa_v` | Market adjustment speed (1/years); 1/κ_v is the lag |
| $W_{cc'}$ | `W` | Row-normalized county proximity weights: $W_{cc'} \propto \exp(-d_{cc'}/\ell_W)$ for $c' \ne c$, with $d$ the centroid distance and $\ell_W$ a trafficking length scale |
| $D_v$ | `D_v` | Spillover strength |
| $\sigma_v$, $B_c$ | `sigma_v` | Noise amplitude; independent Brownian motions (multiplicative noise keeps $v > 0$) |

### 5.3 Potency: Bass diffusion of fentanyl with a COVID supply shock

Fentanyl entered regional markets like a new product: spontaneous introduction plus imitation of markets that have already switched [R09, R10, R19]. In 2.1, the pandemic's supply disruption adds an extra introduction push proportional to each state's policy stringency $S_s(t)$ (§5.4):

$$
\frac{d\varphi_s}{dt} = \big(p_\varphi + p_{cov}\,S_s(t) + q_\varphi\,\varphi_s\big)(1 - \varphi_s), \qquad \varphi_s(t) = 0 \text{ for } t < t_\varphi
$$

| Symbol | Code name | Meaning | Value |
|---|---|---|---|
| $p_\varphi$ | `p_phi` | Innovation rate: spontaneous introduction (1/yr) | 0.0128 (fitted, §9.2) |
| $q_\varphi$ | `q_phi` | Imitation rate: spread from markets that already have fentanyl (1/yr) | 0.6 (assumed) |
| $p_{cov}$ | `p_cov` | Extra introduction per unit stringency: the pandemic supply shock (1/yr) | 1.33 (fitted) |
| $t_\varphi$ | `t_phi` | Fentanyl entry year | 2013 |

The push stops when restrictions end, but the switch it caused is **permanent**, because φ never decreases. The Bass curve is S-shaped, which matches the sudden rise of the fentanyl wave.

Every county in a state shares the same $\varphi_s$, because the only thing that differs between places is $S_s$. Version 2's spatial imitation term ($\omega_\varphi$), hub exposure ($h_c$) and NFLIS nudging ($\kappa^\varphi_{obs}$) were removed: with $h_c = 1$ and no NFLIS data they had no effect.

### 5.4 COVID stringency

$$
S_s(t) = x^{str}_s(t)/100 \in [0, 1]
$$

This is the daily state index, 0 before 2020 and after 2022. It enters twice:

1. the fentanyl push $p_{cov} S_s$ above (permanent effect);
2. the overdose hazard: $\times\, e^{w_{iso} S_s}$ (§7, "using alone" during isolation; temporary effect), with $w_{iso} = 0.14$ (fitted).

---

## 6. Layer 3: space, agents and interaction

### 6.1 Population size and scaling

| Symbol | Code name | Meaning |
|---|---|---|
| $N$ | `N_sim` | Total number of simulated agents (user input; larger means less noise, ∝ $1/\sqrt{N}$) |
| $P_{tot}$ | `N_real` | Real total population of the study region |
| $\Upsilon = P_{tot}/N$ | `scale_factor` | Real people per agent |
| $n_c$ | `n_c` | Agents in county $c$: $n_c = \operatorname{round}(N \cdot P_c / P_{tot})$, kept fixed |

### 6.2 Placement following real density

For each county, place its $n_c$ agents by (1) drawing block group $b$ with probability $P_b / \sum_{b \in c} P_b$, then (2) setting the location

$$
\mathbf{x}_i = \mathbf{b}_b + \sqrt{A_b/\pi}\,\sqrt{U_1}\,(\cos 2\pi U_2,\ \sin 2\pi U_2), \quad U_1, U_2 \sim \text{Unif}(0,1)
$$

This is uniform within a disc the size of the block group. The agent density then reproduces the real density at roughly block-group resolution. Locations are fixed for v1 (no migration). $c(i)$ is the county of agent $i$.

### 6.3 Interaction neighborhood

| Symbol | Code name | Meaning |
|---|---|---|
| $\rho$ | `rho_km` | **Interaction radius** (fixed, user input, km). Default **15 km** |
| $\mathcal{N}_i$ | `nbr[i]` | $\{ j \ne i : \lVert \mathbf{x}_i - \mathbf{x}_j \rVert \le \rho \}$, the agents that can interact with $i$ |
| $n_i$ | `n_nbr` | $\lvert\mathcal N_i\rvert$ |
| $k_0$ | `k0` | Rural fallback (pseudo-count): every neighborhood gets $k_0$ extra virtual agents at the county's current class mix. Sparse neighborhoods lean on the county; dense ones are barely affected (§6.4) |
| $h$ | `cell_km` | Implementation: the radius kernel is evaluated on a 5 km grid. Agents are binned into cells, and per-class counts are convolved with a disc of radius ρ. Error ≈ one cell width, and the cost does not depend on N |

The kernel is a step function $K(d) = \mathbb 1[d \le \rho]$, the one-step special case of the piecewise kernels in [R15]. A smooth version, $K(d) = \exp(-d/\rho)$, can be swapped in later.

### 6.4 Local pressures

The local prevalence of class $X$ around agent $i$ is

$$
\pi^X_i(t) = \frac{\sum_{j \in \mathcal N_i} \mathbb 1[X_j(t) = X] \;+\; k_0\, x_{c(i)}(t)}{n_i + k_0}
$$

where $x_c$ is the county's current fraction in class $X$ (the pseudo-count fallback of §6.3).

The **force of initiation** on a susceptible agent $i$:

$$
\Lambda_i(t) = \pi^I_i(t) + \vartheta_E\,\pi^E_i(t)
$$

| Symbol | Code name | Meaning |
|---|---|---|
| $\vartheta_E$ | `theta_E` | Relative recruiting power of misusers (E) vs OUD (I). Misusers often share pills with friends and family [R18] |

*Removed in 2.1: the density factor.* Version 2 multiplied Λ by $g_i = (n_i/\bar n)^{\eta}$ with η = 0.5, so crowded places produced more contacts. In big metros this pushed the local reproduction number above 1 for about 40% of Texas agents, and those clusters grew without bound. Contacts now depend only on the *share* of users nearby (frequency dependence, η = 0).

**Why this scales:** $\pi^X_i$ is a ratio, so their expected values do not change when $N$ changes. Simulating at 1 agent per 100 people gives the same mean dynamics as 1:1, only noisier. The same $\rho$ works at any $N$.

*Optional long-range term:* add $\epsilon_L \sum_{c'} \text{SCI}_{c(i)c'}\, i_{c'}$ using the Facebook Social Connectedness Index [R18] to model friendships beyond radius ρ.

---

## 7. Layer 4: transition rates (the holistic formulas)

Every rate has the form **baseline × availability term × structural modifier $F$ (§4.2) [+ a peer term]**. Availability enters as a power, $v^{\kappa}$, so κ is an elasticity: +1% of supply gives +κ% to the rate.

| # | Transition | Rate for agent $i$ in county $c$ at time $t$ | Mechanism |
|---|---|---|---|
| 1 | S→E | $\lambda_{SE} = \beta^0\ \Lambda_i(t)\ v_c^{\kappa_\beta}\ F_\beta$ | Social initiation: meeting users × how easy drugs are to get [R02, R13, R18] |
| 2 | S→I | $\lambda_{SI} = \gamma_{Rx}\ \chi\ a_c(t)\ F_\gamma$ | Iatrogenic: share of people on opioid therapy ($\chi \cdot a$, with χ calibrated) × addiction rate per therapy-year [R01] |
| 3 | E→I | $\lambda_{EI} = \sigma^0\ v_c^{\kappa_\sigma}\ F_\sigma$ | Progression from misuse to OUD, faster when supply is plentiful |
| 4 | E→S | $\lambda_{ES} = \epsilon^0\ F_\epsilon$ | Desistance: experimenters who stop |
| 5 | I→R | $\lambda_{IR} = \zeta^0\ F_\zeta + \gamma_{p}\ \pi^R_i(t)$ | Treatment entry (capacity, insurance) + peer recovery ("R meets I") [R01] |
| 6 | I→S | $\lambda_{IS} = \nu^0\ F_\nu$ | Natural remission without treatment (no protection afterwards) |
| 7 | I→D | $\lambda_{ID} = \big[\mu_H(1 - \varphi_s) + \mu_F\,\varphi_s\big]\ e^{w_{iso} S_s(t)}\ F_\mu$ | Overdose. A **mixture of hazards**: heroin/Rx-era hazard and fentanyl-era hazard, weighted by market potency [R09, R12, R19], times the COVID isolation effect (§5.4) |
| 8 | R→S | $\lambda_{RS} = \delta^0\ F_\delta$ | Sustained recovery: leaves treatment and returns to the general at-risk pool [R01] |
| 9 | R→I | $\lambda_{RI} = \big[\rho^0\ v_c^{\kappa_\rho} + \rho_{p}\ \pi^I_i(t)\big]\ F_\rho$ | Relapse: spontaneous (supply-driven) + peer-driven [R01, R02] |

Total exit rate of agent $i$ in state $X$: $\lambda^{out}_i = \sum_{Y \ne X} \lambda_{XY,i}$.

### 7.1 Parameter dictionary for L4

| Symbol | Code name | Meaning | Unit |
|---|---|---|---|
| $\beta^0$ | `beta0` | Initiation rate per unit local pressure | 1/yr |
| $\gamma_{Rx}$ | `gamma_rx` | OUD onset per prescription-person-year | 1/yr |
| $\sigma^0$ | `sigma0` | E→I progression | 1/yr |
| $\epsilon^0$ | `eps0` | E→S desistance | 1/yr |
| $\zeta^0$ | `zeta0` | Treatment entry | 1/yr |
| $\gamma_p$ | `gamma_p` | Peer-recovery strength | 1/yr |
| $\nu^0$ | `nu0` | Natural remission | 1/yr |
| $\mu_H$ | `mu_H` | Overdose hazard of OUD without fentanyl | 1/yr |
| $\mu_F$ | `mu_F` | Overdose hazard of OUD in a fully fentanyl market | 1/yr |
| $\delta^0$ | `delta0` | Sustained recovery R→S | 1/yr |
| $\rho^0$ | `rho0` | Spontaneous relapse | 1/yr |
| $\rho_p$ | `rho_p` | Peer-driven relapse | 1/yr |
| $\kappa_\beta, \kappa_\sigma, \kappa_\rho$ | `k_beta`… | Supply elasticities of initiation, progression and relapse | – |
| $\mathbf w_\theta$ | `W_xi[θ]` | Structural weights (§10) | per SD |

---

## 8. Layer 5: the stochastic engine

### 8.1 Competing Poisson clocks per agent

Within a step of length $\Delta t$ (`dt`, default 1/52 yr), all rates are frozen at their values at the start of the step. The exact law of competing exponential clocks then gives:

$$
\Pr[\text{agent } i \text{ leaves } X \text{ in } \Delta t] = 1 - e^{-\lambda^{out}_i \Delta t},
\qquad
\Pr[\text{destination} = Y \mid \text{leaves}] = \frac{\lambda_{XY,i}}{\lambda^{out}_i}
$$

All agents update **synchronously** using neighbor states from the start of the step. The only approximation is this freezing, with error $O(\Delta t)$. An exact event-driven (Gillespie next-reaction) version is kept as a small-$N$ check.

### 8.2 One step, in order

1. Update the environment: interpolate $\boldsymbol\xi_c(t)$ and $a_c(t)$, and advance $v_c, \varphi_c$ with an Euler–Maruyama step of §5.2–5.3.
2. Compute the county fractions $e_c, i_c$ and every agent's local pressures $\pi^I_i, \pi^E_i, \pi^R_i$ (grid-cell or KD-tree neighbor sums).
3. Compute all nine rates for every agent (vectorized), then draw exits and destinations.
4. Apply the transitions, add deaths to $D_c(y)$, and rebirth the dead (§8.3).
5. $t \leftarrow t + \Delta t$. At each year boundary, record the class fractions and reset $D_c$.

### 8.3 Death and regeneration

An agent who fires I→D increments $D_{c(i)}(y)$ by 1. At once, the same agent slot is **reborn as S at a new location** drawn from §6.2 for the same county. County sizes $n_c$ stay fixed and density stays realistic. Total $N$ is constant.

### 8.4 Randomness and reproducibility

| Symbol | Code name | Meaning |
|---|---|---|
| `seed` | `seed` | Master seed. Each ensemble member $m$ uses an independent child stream (`numpy.random.SeedSequence.spawn`) |
| $M$ | `M` | Ensemble size (number of independent runs) |
| $t_0,\ t_{burn}$ | `t0`, `t_burn` | First data year; burn-in years run under the frozen year-$t_0$ environment so the initial state is in equilibrium with it |

### 8.5 Initial conditions (per county, before burn-in)

| Quantity | Value | Basis |
|---|---|---|
| $i_c(t_0)$ | State OUD prevalence (NSDUH; ≈ 0.6% nationally) | [R01] |
| $r_c(t_0)$ | ≈ $i_c \cdot \zeta^0/\delta^0$, refined by burn-in | Equilibrium balance |
| $e_c(t_0)$ | Misuse prevalence (NSDUH; ≈ 2–4%) | NSDUH |
| $s_c(t_0)$ | $1 - e - i - r$ | – |

Agents are assigned their initial states at random within each county.

---

## 9. Layer 6: observation and comparison

### 9.1 Model observable

$$
m^{mod}_c(y) = 10^4 \cdot \frac{D_c(y)}{n_c}
$$

This is deaths per 10,000 per year. It can be compared directly to $m^{obs}$, because it is a rate and does not depend on $\Upsilon$. State values aggregate with population weights: $m_s = \sum_{c\in s} P_c m_c / \sum_{c \in s} P_c$.

Ensemble statistics: mean $\bar m$ and the 5–95% band across $M$ runs.

### 9.2 Fit score and the four fitted parameters

**Model mean.** Expected real deaths per county-year are
$$\hat Y_c(y) = \frac{P_c(y)}{n_c}\sum_{t \in y}\sum_{i \in I,\ c(i)=c} \lambda_{ID,i}(t)\,\Delta t$$
These are the *expected* deaths (the sum of each OUD agent's hazard), which have the same mean as the counted deaths with much less noise.

**Score (2.1): state-year, interval-censored.** In version 2 the score was a county-level censored negative binomial. That turned out to be unusable: at 0.6% prevalence a county with 50 agents holds about 0.3 OUD agents, so its expected deaths are often exactly 0, and one such county dominates the likelihood. Version 2.1 scores state-years instead. Within a state, the observed total is only known to lie in an interval, because each suppressed county contributes 1 to 9:
$$
Y^{lo}_s(y) = \sum_{c \in s,\ c \notin \mathcal S} Y_c + \#\{c \in s \cap \mathcal S\},\qquad
Y^{hi}_s(y) = \sum_{c \in s,\ c \notin \mathcal S} Y_c + 9\,\#\{c \in s \cap \mathcal S\}
$$
$$
\mathcal L = \frac{1}{|\{s,y\}|}\sum_{s,y} \log \Pr\big(Y^{lo}_s \le Y \le Y^{hi}_s\big), \qquad Y \sim \text{NB}\Big(\textstyle\sum_{c\in s}\hat Y_c,\ r = 50\Big)
$$
No value is guessed for suppressed counts.

**What is fitted to deaths.** Only four parameters, all on the lethality side: $m_F$, $p_\varphi$, $p_{cov}$, $w_{iso}$. Everything else is literature, assumption, or calibrated to non-death targets (§11).

**Why they can be fitted without rerunning agents (separability).** φ and S enter only $\lambda_{ID}$. Overdose is a small exit from I (≈ 4%/yr against ≈ 60%/yr of total exits), so the I trajectory is nearly independent of these four parameters. Expected deaths can then be rebuilt from a saved run's yearly I counts:
$$
\hat Y_c(y) \approx \frac{P_c}{n_c}\, \bar I_c(y)\, F_{\mu,c}(y)\, \overline{\big[\mu_H(1-\varphi_s) + m_F\mu_H\varphi_s\big] e^{w_{iso}S_s}}^{\,(y)}
$$
with φ integrated weekly per state. Fitting takes seconds (`engine/fit_lethality.py`, Powell search). A full agent run with the fitted values confirms the result.

**Held-out tests.**

| Fit on | Fitted | Held-out result |
|---|---|---|
| 2010–2021, all states | $w_{iso}$ = 0.98, $p_{cov}$ = 2.7 | Worse on 2022–23: deaths fall when restrictions end, but the data stay high |
| 2010–2023, 12 states | $m_F$ 3.83, $p_\varphi$ 0.017, $p_{cov}$ 1.60, $w_{iso}$ 0.15 | Better on 5 unseen states (TN, TX, CO, GA, LA): −2.75 → −2.59 |
| 2010–2023, all 17 (final) | $m_F$ 3.88, $p_\varphi$ 0.0128, $p_{cov}$ 1.33, $w_{iso}$ 0.14 | – |

Once the years after restrictions are included, the data choose the permanent supply push over the temporary isolation effect.

The viewer also reports the share of years in which the model mean lies inside the observed range.

### 9.3 Identifiability

There is one observable but many parameters, and several combinations give the same death curve. Countermeasures:

1. Baselines come from the literature (§11), not from fitting.
2. Weights are shared by all counties, and they act on only 4 indices.
3. Any free parameters are pinned by *multiple* targets: deaths, plus NSDUH prevalence ($i$, $e$) and treatment counts ($r$), when available.

### 9.4 Output figures

- **Per state:** time series of $s, e, i, r$ fractions; the model line $\bar m$ with its band; and the observed $m^{obs}$ as markers.
- **Map:** county residual $\log(\bar m / m^{obs})$, showing where the hypothesis fails.
- **Sensitivity heatmap:** $\partial \log \bar m / \partial \xi^k$ for each index and state, placed next to the empirical correlation heatmap from the original project goal.

---

## 10. Weights: sign priors and elasticity matching

Rows are indices, columns are rates. Entries are $w$ (per SD). Blank means 0. The signs come from the cited literature.

**Prior weights** (signs and the mix across rates):

| | β (S→E) | γ (S→I) | σ (E→I) | ε (E→S) | ζ (I→R) | ν (I→S) | μ (I→D) | δ (R→S) | ρ (R→I) |
|---|---|---|---|---|---|---|---|---|---|
| $\xi^{econ}$ | +0.15 | | +0.15 | −0.10 | −0.10 | | | −0.10 | +0.10 |
| $\xi^{pain}$ | | +0.30 | | | | | | | |
| $\xi^{unins}$ | | | | | −0.15 | | +0.10 | | |
| $\xi^{age}$ | −0.10 | +0.20 | | | | | +0.10 | | |

**Weights used in 2.1** (after elasticity matching, below). Each row is the prior row × one scale factor $s_k$:

| | $s_k$ | β | γ | σ | ε | ζ | μ | δ | ρ | target $\partial\ln m^*/\partial\xi^k$ |
|---|---|---|---|---|---|---|---|---|---|---|
| $\xi^{econ}$ | 0.104 | +0.016 | | +0.016 | −0.010 | −0.010 | | −0.010 | +0.010 | +0.15 |
| $\xi^{pain}$ | 0.712 | | +0.213 | | | | | | | +0.20 |
| $\xi^{unins}$ | 0.326 | | | | | −0.049 | +0.033 | | | +0.10 |
| $\xi^{age}$ | 1 | −0.10 | +0.20 | | | | +0.10 | | | none (prior kept; −0.04) |

Why it mattered: with the unscaled priors, $\partial\ln m^*/\partial\xi^{econ} = 1.50$ (+1 SD of precarity multiplies long-run deaths ×4.5), and 26% of the 2010 population lived in counties with local $\mathcal R_c > 1$. Those counties (Mississippi and parts of Texas) ran away. After matching, that share is 4% (max $\mathcal R_c$ 1.18), and no state runs away. The cohesion and treatment rows of version 2 are gone with their indices (§4.2).

Supply elasticities: $\kappa_\beta = 1.0$, $\kappa_\sigma = 0.5$, $\kappa_\rho = 0.5$.

Rationale per row:

- **econ:** precarity raises initiation and progression and lowers exits [R05, R06]. Ruhm [R04] finds this channel smaller than the drug environment, so the magnitudes are kept below the supply elasticities.
- **pain:** disability and veteran status drive prescriptions, hence S→I [R01].
- **unins:** in the literature, uninsured people enter treatment less and die more [R06]. **Known tension:** the earlier correlation work in this region found uninsurance *negatively* correlated with death rates (about −0.25, probably because high-uninsured TX/NM/AZ counties had low death rates). The prior is kept so the model *tests* the mechanism instead of absorbing the correlation. If the residual map (§9.4) shows a systematic error along this index, that is evidence that uninsurance is a proxy for something else, such as a later arrival of fentanyl in the west.
- **age:** older populations have more chronic pain and prescriptions (S→I) and higher lethality, but less peer initiation. The earlier work found median age at +0.3 to +0.4 throughout.
- **coh:** cohesion is protective throughout [R06].
- **tx:** treatment capacity raises entry. Medication for OUD lowers overdose and relapse [R12].

**Elasticity matching (how magnitudes become principled).** For a weight vector $\mathbf w$, compute the model's equilibrium death rate $m^*(\boldsymbol\xi)$ from the mean-field equations (§12). Then scale each row so that $\partial \ln m^* / \partial \xi^k$ equals the published effect converted to per-SD units. Example: [R05] gives +3.6% deaths per +1 pp unemployment. Converted to our index: target $\partial \ln m^*/\partial \xi^{econ} = b_{u}\cdot \ln 1.036$, where $b_u$ is the slope (in pp) of unemployment regressed on $\xi^{econ}$ in our county data. With a county unemployment SD of about 2 pp and a correlation of about 0.8, this gives ≈ 1.6 · 0.035 ≈ 0.06 from the unemployment channel alone. Poverty and education channels add to it. This turns regression literature into mechanistic weights without fitting to our own death data. That keeps the death data free for testing.

**As implemented** (`calibrate.match_elasticities`): at the average county (ξ = 0, v = 1, φ = 0) with the calibrated β⁰ and χ, compute $E_k(s) = \partial\ln m^*/\partial\xi^k$ by central difference (h = 0.25) on the ODE equilibrium, using weights $s\cdot w^{prior}_k$. Solve $E_k(s) = \text{target}_k$ for $s$ with Brent's method. Targets: econ 0.15 (unemployment alone ≈ 0.06 [R05]; poverty and education add to it; the channel is weaker than the drug environment [R04]), pain 0.20, unins 0.10 (kept small, because the regional correlation is negative). Rows are independent, because F = 1 at ξ = 0 whatever the other rows are.

---

## 11. Baseline parameter table

Status: **L** = from literature, **A** = assumption to revisit, **C** = calibrated to a non-death target, **F** = fitted to deaths (§9.2).

| Symbol | Value | Status | Basis |
|---|---|---|---|
| $\gamma_{Rx}$ | 0.0074 /Rx-yr | L | 9.3% of chronic users addicted, scaled to yearly [R01] |
| $\beta^0$ | Set so that total S→E flux ≈ 0.66%/yr at $t_0$ | C | 2.1M new prescription-opioid misusers / 320M in 2015 [R01] |
| $\sigma^0,\ \epsilon^0$ | 0.08, 0.92 | A | Mean time in E ≈ 1 yr, of which ~8% progress to OUD. 15% was infeasible: together with relapse it already exceeds the 0.6% OUD target |
| $\chi$ | ≈ 0.032 | C | Share of the population on opioid therapy per prescription per person. Solved with $\beta^0$ so the mean-field equilibrium has $i^* = 0.6\%$ [R01] |
| $\zeta^0$ | 0.5 (range 0.2–2) | L | [R01] |
| $\gamma_p$ | 0.1 | A | Peer recovery; small, same order as $\delta^0$ |
| $\nu^0$ | 0.05 | A | Untreated remission is slow |
| $\mu_H$ | 0.0127 | L | Overdose deaths out of treatment, 12.7 per 1,000 person-years [R12]; compare $\mu^* = 0.0116$ in [R01] |
| $\mu_F$ | $m_F \cdot \mu_H$, $m_F$ = **3.88** | **F** | Fitted to state-year deaths (§9.2) |
| $\delta^0$ | 0.1 | L | Successful treatment rate [R01] |
| $\rho^0,\ \rho_p$ | 0.7, 0.2 | L | ~90% relapse within a year of treatment; split as in [R01] |
| $\vartheta_E$ | 0.3 | A | Misusers recruit less than people with OUD |
| $\rho$ (radius) | **15 km** | User-set | Roughly a commuting and drug-market scale. In rural counties the $k_0$ fallback still applies. Sensitivity sweep 10–30 km |
| $k_0$ | 20 | A | Rural pseudo-count |
| $\theta_d,\ \kappa_v,\ D_v,\ \sigma_v$ | **0.1**, 1.0, 0.2, 0.1 | A | Supply follows demand within about one year. θ_d is kept below the stability threshold ≈ 0.2 (§12.5) |
| $t_\varphi$ | 2013 | A | Fentanyl entry year |
| $p_\varphi,\ q_\varphi$ | **0.0128** (F), 0.6 (A) | F, A | Bass introduction (fitted) and imitation (assumed) |
| $p_{cov}$ | **1.33** | F | COVID push on fentanyl introduction, per unit stringency |
| $w_{iso}$ | **0.14** | F | COVID isolation effect on overdose (×1.07 at S = 0.5) |
| $\ell_W$ | 100 km | A | Trafficking length scale |
| $\Delta t$ | 1/52 yr | – | Weekly step |

**Note on R→D:** people in treatment do still die, at 1.4–2.6 per 1,000 person-years [R12], about 5–10× less than $\mu_H$. Your transition list excludes R→D, so v1 ignores this. It can be added as a tenth transition.

---

## 12. Layer 7: analytical checks (mean-field limit)

### 12.1 Well-mixed equations

In one well-mixed county ($\rho \to \infty$), with constant environment, the fractions obey (rates with tildes include availability and $F$):

$$
\begin{aligned}
\dot s &= -\tilde\beta (i + \vartheta_E e)\, s - \tilde\gamma\, s + \tilde\epsilon\, e + \tilde\nu\, i + \tilde\delta\, r + \tilde\mu\, i \\
\dot e &= \tilde\beta (i + \vartheta_E e)\, s - (\tilde\sigma + \tilde\epsilon)\, e \\
\dot i &= \tilde\gamma\, s + \tilde\sigma\, e + (\tilde\rho + \rho_p i)\, r - (\tilde\zeta + \gamma_p r + \tilde\nu + \tilde\mu)\, i \\
\dot r &= (\tilde\zeta + \gamma_p r)\, i - (\tilde\rho + \rho_p i + \tilde\delta)\, r
\end{aligned}
$$

with $s + e + i + r = 1$. The $+\tilde\mu i$ term in $\dot s$ is the rebirth of the dead.

The mean of the agent simulation must converge to this ODE as $N \to \infty$ (Kurtz's theorem), with noise of order $1/\sqrt N$. The spatial version replaces $i$ by its kernel-smoothed field $(K_\rho * i)(\mathbf x)$, which gives an integro-differential equation.

### 12.2 No drug-free equilibrium

With $\tilde\gamma > 0$ (prescriptions exist), $(s, e, i, r) = (1, 0, 0, 0)$ is **not** an equilibrium: S→I keeps creating cases. This matches the result in [R01]. Only $a_c = 0$ allows elimination.

### 12.3 Contact reproduction number

Set $\tilde\gamma = 0$, $\gamma_p = \rho_p = 0$, and seed one I in an all-S population. Let

- $L_I = \tilde\zeta + \tilde\nu + \tilde\mu$, the total exit rate from I
- $p_{RI} = \tilde\rho / (\tilde\rho + \tilde\delta)$, the probability that a recovery ends in relapse
- $T_I = 1/(L_I - \tilde\zeta\, p_{RI})$, the expected total time spent in I, counting all relapse cycles

Then

$$
\mathcal R_c = \tilde\beta \left[\frac{\tilde\sigma}{\tilde\sigma + \tilde\epsilon}\,T_I + \frac{\vartheta_E}{\tilde\sigma + \tilde\epsilon}\right]
$$

(When $\vartheta_E = 0$ and there is no relapse, this reduces to version 1's $\beta/(\gamma + \mu)$.) $\mathcal R_c > 1$ means social spread can sustain itself even with no prescriptions. This formula is **engine test 1**.

### 12.4 Equilibrium death rate

Solve $\dot{(\cdot)} = 0$ numerically to get $i^*$. Then $m^* = 10^4\, \tilde\mu\, i^*$. This gives

- **engine test 2:** the large-$N$ simulation must match $m^*$;
- the elasticities $\partial \ln m^* / \partial \xi^k$ used in §10.

### 12.5 Supply–demand feedback creates a tipping point (found while building v1)

Add the supply equation of §5.2 (well-mixed, no noise) to the ODE of §12.1 and linearize at the calibrated equilibrium. The leading non-trivial eigenvalue of the Jacobian, as a function of the demand elasticity θ_d:

| θ_d | 0 | 0.1 | 0.2 | 0.3 | 0.5 |
|---|---|---|---|---|---|
| leading eigenvalue (1/yr) | −0.024 | −0.012 | ≈ 0 | +0.012 | +0.039 |

Two consequences:

1. **There is a bifurcation at θ_d ≈ 0.2.** Above it, supply chasing demand makes the endemic state unstable, and the model generates a self-sustaining rise in drug use on its own. That is a candidate mechanism for an epidemic that is *not* driven by exogenous supply. The agent engine confirmed it: with θ_d = 0.5, OUD prevalence grew from 0.6% to about 2% in 10 burn-in years. v1 uses θ_d = 0.1.
2. **Even when stable, the system is slow.** The relaxation time is 1/0.024 ≈ 40 years, because relapse cycling makes the total time spent in I about 8 years. Changes in rates therefore show up in prevalence only gradually.

The density exponent η = 0.5 alone did not destabilize the system once θ_d ≤ 0.1.

---

## 13. Data you need to download

| # | Data | Variables | Status |
|---|---|---|---|
| 1 | SRDDRC opioid + all-drug deaths, county-year 2010–2023 | $Y$, $P_c$ | ✅ `data/` |
| 2 | ACS 5-year county, 2012–2023 | pov, unemp, inc, nohs, bach, dis, unins, vet, age | ✅ `data/` |
| 3 | County and state shapes (2022, 1:20m) | geometry, `ALAND` | ✅ `data/shapes/` |
| 4 | 2020 **census tract** population + centroids (used instead of block groups: 33,696 tracts, about 4,000 people each) | $P_b$, $\mathbf b_b$, $A_b$ | ✅ `data/tracts/`, fetched |
| 5 | ACS DP02 single-parent + CBP religious orgs | single, relig | ⬜ I fetch (optional) |
| 6 | CDC opioid dispensing rate by county | $x^{Rx}$ → $a_c(t)$ | ✅ 2010–2018 county (CDC archive); 2019–2023 extended by the national trend |
| 7 | **NFLIS-Drug annual reports by state, 2010–2023** | $x^{op}$, $x^{fent}$ → $\hat v$, $\varphi$ | ⬜ **You.** DEA Diversion NFLIS public data / annual reports |
| 8 | Treatment facilities and buprenorphine prescribers | $x^{tx}$ | ⬜ Optional (SAMHSA) |
| 9 | Facebook Social Connectedness Index (county–county) | SCI | ⬜ Optional extension |

**If 6 and 7 are not available**, the fallback is:

- **Prescriptions:** a statewide constant prescription trend from the CDC national series (a few numbers I can type in from the published table).
- **Potency:** $\varphi$ runs in free Bass mode with $(p_\varphi, q_\varphi, m_F)$ calibrated on the 2010–2016 deaths only, and 2017–2023 kept as the test period.

The engine works either way. Rows 6 and 7 are what make drug availability data-driven instead of model-driven.

---

## 14. Open questions

Resolved: geography is **county-year**, years are **2010–2023** ($t_0 = 2010$, burn-in under the 2012 environment), and radius **ρ = 15 km**.

Still open:

1. **Fentanyl timing by state:** synthetic-opioid (T40.4) death shares or NFLIS reports per state would replace the fitted, stringency-driven φ. This is the main suspect for the model's states being too alike (Texas too high; Kentucky and Tennessee too low).
2. ~~Data nudging (κ_obs)~~ Removed in 2.1 (no NFLIS data).
3. ~~COVID~~ Included in 2.1 (§5.4) at the user's request.

---

## 15. Changes from version 2 to 2.1 (2026-11-05)

| Change | Why | Where |
|---|---|---|
| Density factor $g_i$ removed (η = 0) | Metros became locally supercritical; Texas ran away (4.7 vs ≈1.0 observed deaths/10k in 2021) | §6.4 |
| Weights sized by elasticity matching | Priors implied econ ×4.5 per SD; 26% of people in counties with $\mathcal R_c > 1$ | §10 |
| COVID stringency (OxCGRT) added, in two channels | Every state's deaths jumped ×1.4–2.2 in 2019–21 | §3.2, §5.3, §5.4, §7 |
| Four lethality parameters fitted to deaths ($m_F, p_\varphi, p_{cov}, w_{iso}$) | Fentanyl timing and lethality were guesses | §9.2, §11 |
| Scoring moved to state-year interval likelihood | The county score was dominated by zero-mean counties | §9.2 |
| Removed with no effect on results: cohesion/treatment/urbanicity indices, data nudging $\kappa_{obs}$, hub exposure $h_c$, fentanyl spatial imitation $\omega_\varphi$ | Inactive (no data) | §4.2, §5 |

Result (300k agents × 3 runs): state-year log-likelihood −2.44 (version 2: −70.3). The region's shape is right, but the states are too alike (§14.1). Details, figures and the run comparison are in [REPORT.md](REPORT.md).
