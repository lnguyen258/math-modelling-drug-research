"""Mean-field equations (REPORT section 5) and calibration of beta0, chi_rx to non-death targets.

Calibration point: the average county (xi = 0, F = 1), v = 1, phi = 0, no COVID, at the t0
prescription level. Unknowns beta0 and chi_rx are solved in closed form so that the
equilibrium has i* = target_i and S->E flux = target_flux_SE.
"""
import numpy as np
from scipy.integrate import solve_ivp


def calibrate(p, a_mean):
    i = p.target_i
    f = p.target_flux_SE
    e = f / (p.sigma0 + p.eps0)
    r = p.zeta0 * i / (p.rho0 + p.rho_p * i + p.delta0 - p.gamma_p * i)
    s = 1 - e - i - r
    out_I = (p.zeta0 + p.gamma_p * r + p.nu0 + p.mu_H) * i
    in_other = p.sigma0 * e + (p.rho0 + p.rho_p * i) * r
    chi = (out_I - in_other) / (p.gamma_rx * a_mean * s)
    if chi < 0:
        raise ValueError(f"Infeasible calibration: progression + relapse already exceed the I target (chi={chi:.4f}). "
                         "Lower sigma0 or rho0.")
    p.chi_rx = float(chi)
    p.beta0 = float(f / ((i + p.theta_E * e) * s))
    return dict(s=s, e=e, i=i, r=r)


def rhs(p, a, v=1.0, phi=0.0, F=None):
    """Well-mixed ODE right-hand side for fractions (s, e, i, r)."""
    F = F or {}
    f = lambda k: F.get(k, 1.0)
    mu = (p.mu_H * (1 - phi) + p.mu_F * phi) * f("mu")

    def g(t, x):
        s, e, i, r = x
        lam_SE = p.beta0 * (i + p.theta_E * e) * v ** p.k_beta * f("beta")
        lam_SI = p.gamma_rx * p.chi_rx * a * f("gamma")
        lam_EI = p.sigma0 * v ** p.k_sigma * f("sigma")
        lam_ES = p.eps0 * f("eps")
        lam_IR = p.zeta0 * f("zeta") + p.gamma_p * r
        lam_IS = p.nu0 * f("nu")
        lam_RS = p.delta0 * f("delta")
        lam_RI = (p.rho0 * v ** p.k_rho + p.rho_p * i) * f("rho")
        ds = -(lam_SE + lam_SI) * s + lam_ES * e + (lam_IS + mu) * i + lam_RS * r
        de = lam_SE * s - (lam_EI + lam_ES) * e
        di = lam_SI * s + lam_EI * e + lam_RI * r - (lam_IR + lam_IS + mu) * i
        dr = lam_IR * i - (lam_RS + lam_RI) * r
        return [ds, de, di, dr]
    return g, mu


def equilibrium(p, a, years=200.0, **kw):
    g, mu = rhs(p, a, **kw)
    x0 = [1 - p.target_i - 0.01, 0.01, p.target_i, 0.0]
    sol = solve_ivp(g, (0, years), x0, rtol=1e-9, atol=1e-12)
    x = sol.y[:, -1]
    return dict(s=x[0], e=x[1], i=x[2], r=x[3], deaths_per_10k=1e4 * mu * x[2])


def reproduction_number(p, a=0.0):
    """Contact reproduction number R_c of REPORT.md (prescription onset and peer terms off)."""
    L_I = p.zeta0 + p.nu0 + p.mu_H
    p_RI = p.rho0 / (p.rho0 + p.delta0)
    T_I = 1 / (L_I - p.zeta0 * p_RI)
    return p.beta0 * (p.sigma0 / (p.sigma0 + p.eps0) * T_I + p.theta_E / (p.sigma0 + p.eps0))


def elasticity(p, a, k, weights, h=0.25):
    """d ln(equilibrium deaths) / d xi_k at the average county, for a given weight table."""
    from .params import RATES

    def m_star(x):
        F = {r: float(np.exp(weights.get(k, {}).get(r, 0.0) * x)) for r in RATES}
        return equilibrium(p, a, years=400.0, F=F)["deaths_per_10k"]
    return (np.log(m_star(h)) - np.log(m_star(-h))) / (2 * h)


def match_elasticities(p, a_mean):
    """Scale each index's prior weight row by one factor so its equilibrium elasticity equals the target.

    Keeps the prior signs and the mix across rates; only the overall size per index is set. Fills p.weights.
    Run after calibrate() (needs beta0, chi_rx).
    """
    from scipy.optimize import brentq
    out, report = {}, {}
    for k, row in p.weights_prior.items():
        target = p.elasticity_targets.get(k)
        if target is None:
            out[k] = dict(row)
            report[k] = (1.0, elasticity(p, a_mean, k, {k: row}))
            continue
        f = lambda s: elasticity(p, a_mean, k, {k: {r: s * w for r, w in row.items()}}) - target
        s = brentq(f, 1e-4, 3.0, xtol=1e-3)
        out[k] = {r: s * w for r, w in row.items()}
        report[k] = (s, target)
    p.weights = out
    return report
