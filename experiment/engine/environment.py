"""Layer 2 environment: modifiers F, prescriptions a, illicit supply v, fentanyl share phi, COVID stringency S."""
import numpy as np

from .params import RATES, XI_NAMES


class Environment:
    def __init__(self, inp, p, rng):
        self.inp, self.p, self.rng = inp, p, rng
        C = len(inp.fips)
        # Weight matrix [K, rates] (elasticity-matched weights, see calibrate.match_elasticities)
        self.Wxi = np.array([[p.weights.get(k, {}).get(r, 0.0) for r in RATES] for k in XI_NAMES])
        # County proximity weights W_cc' ~ exp(-d / ell), row-normalized, no self-loops
        d = np.linalg.norm(inp.centroid[:, None, :] - inp.centroid[None, :, :], axis=-1)
        W = np.exp(-d / p.ell_W_km)
        np.fill_diagonal(W, 0.0)
        self.W = W / W.sum(axis=1, keepdims=True)
        self.v = np.ones(C)
        self.phi = np.zeros(C)
        self.demand0 = None  # region-wide (i + theta_E e) when the run starts; fixes the scale of v*

    def _interp(self, arr, t):
        """Annual values sit at mid-year; linear interpolation between them, flat outside the data range."""
        yrs = self.inp.years + 0.5
        tt = np.clip(t, yrs[0], yrs[-1])
        k = min(np.searchsorted(yrs, tt, side="right") - 1, len(yrs) - 2)
        f = (tt - yrs[k]) / (yrs[k + 1] - yrs[k])
        return (1 - f) * arr[:, k] + f * arr[:, k + 1]

    def static(self, t):
        """Modifiers F [C, rates], prescriptions a [C] and COVID stringency S [C] at time t."""
        F = np.exp(self._interp(self.inp.xi, t) @ self.Wxi)
        a = self._interp(self.inp.a_rx, t)
        S = self.inp.stringency_at(t)[self.inp.state_of]
        return F, a, S

    def step(self, t, dt, demand, S):
        """Advance v and phi by one step. demand = smoothed (i + theta_E e) per county."""
        p = self.p
        if self.demand0 is None:
            self.demand0 = max(np.average(demand, weights=self.inp.pop[:, 0]), 1e-9)
        v_star = np.power(np.maximum(demand, 0) / self.demand0, p.theta_d)
        drift = p.kappa_v * (v_star - self.v) + p.D_v * (self.W @ self.v - self.v)
        noise = p.sigma_v * self.v * np.sqrt(dt) * self.rng.standard_normal(len(self.v))
        self.v = np.maximum(self.v + drift * dt + noise, 1e-3)
        if t >= p.t_phi:
            self.phi = advance_phi(self.phi, S, p, dt)


def advance_phi(phi, S, p, dt):
    """Bass diffusion of the fentanyl share with a COVID supply shock: dphi/dt = (p_phi + p_cov*S + q_phi*phi)(1-phi)."""
    return np.clip(phi + (p.p_phi + p.p_cov * S + p.q_phi * phi) * (1 - phi) * dt, 0.0, 1.0)
