"""Fit the lethality parameters (m_F, p_phi, w_iso, p_cov) on a saved run.

Why this works without rerunning agents: phi and COVID act only on the overdose hazard mu, and
overdose is a small exit from I (~4%/yr of ~60%/yr total), so the I trajectory is nearly independent
of them. Expected deaths = sum over county-years of I_c(y) * mu_c(y) * P_c / n_c, with mu rebuilt here.
A full agent run with the fitted values is the check (REPORT section 6).

  python -m engine.fit_lethality --run v2_prefit --fit_until 2021
"""
import argparse
import json

import numpy as np
from scipy.optimize import minimize

from .environment import advance_phi
from .inputs import load_inputs
from .params import Params, RATES, XI_NAMES
from .run import RUNS
from . import scoring

FIT_NAMES = ["m_F", "p_phi", "w_iso", "p_cov"]
BOUNDS = {"m_F": (1.0, 10.0), "p_phi": (1e-4, 0.2), "w_iso": (0.0, 4.0), "p_cov": (0.0, 3.0)}


def yearly_hazard_multiplier(inp, p, years):
    """[states, Y]: yearly mean of (mu_H(1-phi) + mu_F phi) * exp(w_iso S) per state, weekly resolution."""
    n_s = len(inp.states)
    phi = np.zeros(n_s)
    out = np.zeros((n_s, len(years)))
    steps = int(round(1 / p.dt))
    for k, y in enumerate(years):
        acc = np.zeros(n_s)
        for w in range(steps):
            t = y + w * p.dt
            S = inp.stringency_at(t)
            if t >= p.t_phi:
                phi = advance_phi(phi, S, p, p.dt)
            acc += (p.mu_H * (1 - phi) + p.mu_F * phi) * np.exp(p.w_iso * S)
        out[:, k] = acc / steps
    return out


def expected_deaths(inp, p, I_cy, n_c):
    """Expected real deaths [C, Y] given mean I counts per county-year from a run."""
    Y = I_cy.shape[1]
    years = inp.years[:Y]
    mult = yearly_hazard_multiplier(inp, p, years)[inp.state_of]            # [C, Y]
    w_mu = np.array([p.weights.get(k, {}).get("mu", 0.0) for k in XI_NAMES])
    F_mu = np.exp(inp.xi[:, :Y, :] @ w_mu)                                   # [C, Y]
    return scoring.expected_real_deaths(I_cy * F_mu * mult, n_c, inp.pop[:, :Y])


def fit(run_name, fit_until=2021, holdout_states=()):
    """Fit on years <= fit_until and on all states except holdout_states; report both held-out scores."""
    inp = load_inputs()
    res = dict(np.load(RUNS / run_name / "results.npz"))
    pd_ = json.loads((RUNS / run_name / "params.json").read_text())
    p = Params(**{k: v for k, v in pd_.items() if k in Params.__dataclass_fields__})
    Yd = len(inp.years)
    I_cy = res["classes_cy"][:, :, :Yd, 2].mean(0)
    n_c = res["n_c"]
    fit_cols = inp.years <= fit_until
    hold = np.isin(inp.states, list(holdout_states))
    train = ~hold

    def loss(z):
        for name, val in zip(FIT_NAMES, z):
            setattr(p, name, float(val))
        return -scoring.state_interval_loglik(inp, expected_deaths(inp, p, I_cy, n_c), fit_cols, states=train)

    z0 = [getattr(p, n) for n in FIT_NAMES]
    r = minimize(loss, z0, method="Powell", bounds=[BOUNDS[n] for n in FIT_NAMES],
                 options=dict(xtol=1e-3, ftol=1e-5, maxfev=3000))
    best = dict(zip(FIT_NAMES, map(float, r.x)))
    for name, val in best.items():
        setattr(p, name, val)
    test_cols = ~fit_cols
    score = lambda m: dict(
        fit=scoring.state_interval_loglik(inp, m, fit_cols, states=train),
        test_years=scoring.state_interval_loglik(inp, m, test_cols, states=train) if test_cols.any() else None,
        test_states=scoring.state_interval_loglik(inp, m, fit_cols, states=hold) if hold.any() else None)
    report = dict(start=dict(zip(FIT_NAMES, z0)), best=best, after=score(expected_deaths(inp, p, I_cy, n_c)))
    for name, val in zip(FIT_NAMES, z0):
        setattr(p, name, val)
    report["before"] = score(expected_deaths(inp, p, I_cy, n_c))
    return report


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", default="v2_prefit")
    ap.add_argument("--fit_until", type=int, default=2021)
    ap.add_argument("--holdout", nargs="*", default=[], help="state names left out of the fit")
    a = ap.parse_args()
    rep = fit(a.run, a.fit_until, a.holdout)
    print(json.dumps(rep, indent=2))
    (RUNS / a.run / "lethality_fit.json").write_text(json.dumps(rep, indent=2))


if __name__ == "__main__":
    main()
