"""All model parameters, named as in REPORT.md / REPORT.md."""
from dataclasses import dataclass, field, asdict

# Structural indices xi (REPORT section 2.1)
XI_NAMES = ["econ", "pain", "unins", "age"]
RATES = ["beta", "gamma", "sigma", "eps", "zeta", "nu", "mu", "delta", "rho"]


def _prior_weights():
    # Sign priors per SD of each index (REPORT 3.2). Magnitudes are rescaled per index by
    # elasticity matching (calibrate.match_elasticities); only the signs and the mix across rates are kept.
    return {
        "econ":  {"beta": 0.15, "sigma": 0.15, "eps": -0.10, "zeta": -0.10, "delta": -0.10, "rho": 0.10},
        "pain":  {"gamma": 0.30},
        "unins": {"zeta": -0.15, "mu": 0.10},
        "age":   {"beta": -0.10, "gamma": 0.20, "mu": 0.10},
    }


def _elasticity_targets():
    # Target d ln(equilibrium deaths) / d xi, per SD. None = keep the prior weights unscaled.
    # econ: unemployment alone gives ~0.06/SD [R05]; poverty and education add to it; Ruhm [R04] finds the
    # economic channel smaller than the drug environment. pain: prescription exposure channel [R01].
    # unins: literature sign is +, kept small because the regional correlation is negative (the model tests it).
    return {"econ": 0.15, "pain": 0.20, "unins": 0.10, "age": None}


@dataclass
class Params:
    # --- run control ---
    N_sim: int = 300_000          # number of agents
    M: int = 3                    # ensemble size
    seed: int = 2026
    dt: float = 1 / 52            # years per step
    t0: int = 2010                # first data year
    t_end: int = 2024             # simulate up to (exclusive) this year
    t_burn: float = 10.0          # burn-in years under the frozen t0 environment

    # --- space ---
    rho_km: float = 15.0          # interaction radius
    cell_km: float = 5.0          # grid cell used to evaluate the radius kernel
    k0: float = 20.0              # pseudo-count: sparse neighbourhoods borrow from the county mean
    theta_E: float = 0.3          # recruiting power of E relative to I

    # --- transition baselines ---
    beta0: float = float("nan")   # calibrated to the S->E flux target
    gamma_rx: float = 0.0074      # OUD onset per opioid-therapy person-year [R01]
    chi_rx: float = float("nan")  # share on opioid therapy per (Rx per person); calibrated to the I target
    sigma0: float = 0.08          # E->I
    eps0: float = 0.92            # E->S
    zeta0: float = 0.5            # I->R [R01]
    gamma_p: float = 0.1          # peer recovery
    nu0: float = 0.05             # I->S
    mu_H: float = 0.0127          # overdose hazard without fentanyl [R12]
    m_F: float = 3.88             # fentanyl hazard multiplier, mu_F = m_F * mu_H  [FITTED on deaths, fit_lethality]
    delta0: float = 0.1           # R->S [R01]
    rho0: float = 0.7             # spontaneous relapse [R01]
    rho_p: float = 0.2            # peer relapse [R01]
    k_beta: float = 1.0           # supply elasticities
    k_sigma: float = 0.5
    k_rho: float = 0.5

    # --- calibration targets (non-death data) ---
    target_i: float = 0.006       # OUD prevalence [R01]
    target_flux_SE: float = 0.0066  # new misuse per person-year [R01]

    # --- illicit supply v ---
    theta_d: float = 0.1          # demand elasticity; > ~0.2 makes the endemic state unstable
    kappa_v: float = 1.0
    D_v: float = 0.2
    sigma_v: float = 0.1
    ell_W_km: float = 100.0

    # --- fentanyl share phi (Bass diffusion) ---
    t_phi: float = 2013.0         # fentanyl enters the market
    p_phi: float = 0.0128         # spontaneous introduction  [FITTED]
    q_phi: float = 0.6            # imitation

    # --- COVID (OxCGRT stringency S in [0, 1], by state) ---
    w_iso: float = 0.14           # isolation: mu multiplied by exp(w_iso * S)  [FITTED]
    p_cov: float = 1.33           # supply shock: extra fentanyl introduction rate p_cov * S  [FITTED]
    # The four [FITTED] values come from `python -m engine.fit_lethality --run v2_prefit --fit_until 2023`
    # (state-level interval likelihood, all years). Holding out 5 states gives nearly the same values.

    weights_prior: dict = field(default_factory=_prior_weights)
    elasticity_targets: dict = field(default_factory=_elasticity_targets)
    weights: dict = field(default_factory=dict)   # filled by calibrate.match_elasticities

    @property
    def mu_F(self):
        return self.m_F * self.mu_H

    def to_dict(self):
        return asdict(self)
