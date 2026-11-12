"""Observation model and comparison to data (REPORT section 6)."""
import numpy as np
from scipy.special import gammaln


def expected_real_deaths(deaths_cy, n_c, pop):
    """Model deaths scaled to real population: D_c * P_c / n_c. Counties with no agents -> nan. [..., C, Y]"""
    with np.errstate(divide="ignore", invalid="ignore"):
        scale = np.where(n_c[:, None] > 0, pop / n_c[:, None], np.nan)
    return deaths_cy * scale


def aggregate(values_cy, groups, n_groups):
    """Sum [..., C, Y] over counties into [..., G, Y] by group index, ignoring nan."""
    v = np.nan_to_num(values_cy)
    out = np.zeros(v.shape[:-2] + (n_groups, v.shape[-1]))
    for g in range(n_groups):
        out[..., g, :] = v[..., groups == g, :].sum(-2)
    return out


def observed_band(inp, groups, n_groups):
    """Observed deaths per group-year as [low, high]: suppressed county-years count as 1 or 9."""
    vis = np.nan_to_num(inp.deaths_obs)
    sup = inp.suppressed.astype(float)
    lo = aggregate(vis + 1 * sup, groups, n_groups)
    hi = aggregate(vis + 9 * sup, groups, n_groups)
    return lo, hi


def nb_logpmf(k, mean, r):
    mean = np.maximum(mean, 1e-9)
    return (gammaln(k + r) - gammaln(r) - gammaln(k + 1)
            + r * np.log(r / (r + mean)) + k * np.log(mean / (r + mean)))


def state_interval_loglik(inp, mean_cy, cols=None, r=50.0, states=None):
    """Mean log-likelihood per state-year: observed state total lies in [lo, hi] (hidden county counts = 1..9),
    model total ~ NB(mean, r). Counties without agents are left out of both sides."""
    from scipy.stats import nbinom
    cols = np.ones(mean_cy.shape[1], bool) if cols is None else cols
    ok = np.isfinite(mean_cy[:, 0])
    n_s = len(inp.states)
    vis = np.where(ok[:, None], np.nan_to_num(inp.deaths_obs), 0.0)
    sup = np.where(ok[:, None], inp.suppressed.astype(float), 0.0)
    lo = aggregate(vis + sup, inp.state_of, n_s)[:, cols]
    hi = aggregate(vis + 9 * sup, inp.state_of, n_s)[:, cols]
    mu = np.maximum(aggregate(np.where(ok[:, None], mean_cy, 0.0), inp.state_of, n_s)[:, cols], 1e-9)
    pr = r / (r + mu)
    prob = nbinom.cdf(hi, r, pr) - nbinom.cdf(lo - 1, r, pr)
    ll = np.log(np.maximum(prob, 1e-300))
    return float((ll if states is None else ll[states]).mean())
