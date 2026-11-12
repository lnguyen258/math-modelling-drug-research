"""Spatial stochastic agent engine (REPORT.md and 8).

States are int8 codes S=0, E=1, I=2, R=3. Each step:
  1. environment (F, a, v, phi) for this time
  2. local pressures pi^X within radius rho, evaluated on a grid (cell_km) by convolving
     per-cell counts with a disc kernel; sparse neighbourhoods borrow k0 pseudo-agents
     at the county mean (no density factor: contacts do not scale with crowding)
  3. nine rates per agent, competing-risk draw of exit and destination
  4. deaths are counted and the agent is reborn as S at a new spot in its county
"""
import numpy as np
from scipy.ndimage import convolve

from .environment import Environment
from .params import RATES

S, E, I, R = 0, 1, 2, 3
RI = {r: k for k, r in enumerate(RATES)}


class Space:
    def __init__(self, inp, p, rng):
        self.rng = rng
        C = len(inp.fips)
        pop0 = inp.pop[:, 0]
        self.n_c = np.floor(p.N_sim * pop0 / pop0.sum()).astype(int)
        # distribute the rounding remainder to the largest fractional parts
        rem = p.N_sim - self.n_c.sum()
        frac = p.N_sim * pop0 / pop0.sum() - self.n_c
        self.n_c[np.argsort(-frac)[:rem]] += 1

        tr = inp.tracts.sort_values("c")
        self.t_c = tr["c"].to_numpy()
        self.t_xy = tr[["x_km", "y_km"]].to_numpy()
        self.t_rad = np.sqrt(tr["area_km2"].to_numpy() / np.pi)
        w = tr["pop"].to_numpy(float)
        self.t_start = np.searchsorted(self.t_c, np.arange(C))
        self.t_end = np.searchsorted(self.t_c, np.arange(C), side="right")
        # cumulative tract weights within each county, for O(log n) sampling
        self.t_cum = np.zeros(len(w))
        for c in range(C):
            a, b = self.t_start[c], self.t_end[c]
            if b > a:
                self.t_cum[a:b] = np.cumsum(w[a:b]) / w[a:b].sum()

        self.county = np.repeat(np.arange(C), self.n_c)
        self.xy = self.place(self.county)

        # grid and disc kernel
        self.cell = p.cell_km
        lo = self.t_xy.min(0) - 50
        hi = self.t_xy.max(0) + 50
        self.lo = lo
        self.shape = tuple(np.ceil((hi - lo) / self.cell).astype(int) + 1)
        r = int(np.ceil(p.rho_km / self.cell))
        dx, dy = np.meshgrid(np.arange(-r, r + 1), np.arange(-r, r + 1), indexing="ij")
        self.kernel = ((dx ** 2 + dy ** 2) * self.cell ** 2 <= p.rho_km ** 2).astype(float)
        self.cell_id = self.cell_of(self.xy)

    def place(self, counties):
        """Tract drawn by population within the county, then uniform in a disc of the tract's area (6.2)."""
        u = self.rng.random(len(counties))
        k = np.empty(len(counties), dtype=int)
        for c in np.unique(counties):
            m = counties == c
            a, b = self.t_start[c], self.t_end[c]
            k[m] = a + np.minimum(np.searchsorted(self.t_cum[a:b], u[m]), b - a - 1)
        rad = self.t_rad[k] * np.sqrt(self.rng.random(len(k)))
        ang = 2 * np.pi * self.rng.random(len(k))
        return self.t_xy[k] + np.c_[rad * np.cos(ang), rad * np.sin(ang)]

    def cell_of(self, xy):
        ij = np.floor((xy - self.lo) / self.cell).astype(int)
        return ij[:, 0] * self.shape[1] + ij[:, 1]

    def neighbourhood_counts(self, state):
        """Agents within rho of each cell, per state: array [4, n_cells]."""
        n_cells = self.shape[0] * self.shape[1]
        cnt = np.bincount(state.astype(np.int64) * n_cells + self.cell_id, minlength=4 * n_cells)
        out = np.empty((4, n_cells))
        for x in range(4):
            out[x] = convolve(cnt[x * n_cells:(x + 1) * n_cells].reshape(self.shape).astype(float),
                              self.kernel, mode="constant").ravel()
        return out


class Simulation:
    def __init__(self, inp, p, seed, eq):
        """eq: initial class fractions {s, e, i, r} (the calibrated mean-field equilibrium)."""
        self.inp, self.p = inp, p
        self.rng = np.random.default_rng(seed)
        self.space = Space(inp, p, self.rng)
        self.env = Environment(inp, p, self.rng)
        self.C, self.Y, self.n_states = len(inp.fips), len(inp.years), len(inp.states)
        N = len(self.space.county)
        self.state = self.rng.choice(4, size=N, p=[eq["s"], eq["e"], eq["i"], eq["r"]]).astype(np.int8)

    def county_fractions(self):
        cnt = np.bincount(self.space.county * 4 + self.state, minlength=self.C * 4).reshape(self.C, 4)
        return cnt / np.maximum(cnt.sum(1, keepdims=True), 1)

    def step(self, t, frozen_t=None):
        p, sp, rng = self.p, self.space, self.rng
        te = t if frozen_t is None else frozen_t
        F, a, S_cov = self.env.static(te)
        frac = self.county_fractions()
        # county demand, shrunk toward the regional mean with k0 pseudo-agents (small counties have few agents)
        dem = frac[:, I] + p.theta_E * frac[:, E]
        dem_bar = np.average(dem, weights=sp.n_c)
        dem = (sp.n_c * dem + p.k0 * dem_bar) / (sp.n_c + p.k0)
        self.env.step(te, p.dt, dem, S_cov)
        v, phi = self.env.v, self.env.phi

        # local pressures with pseudo-count fallback to the county mean (6.4)
        K = sp.neighbourhood_counts(self.state)
        cid, cty = sp.cell_id, sp.county
        n_i = K.sum(0)[cid]
        denom = n_i + p.k0
        piE, piI, piR = ((K[X][cid] + p.k0 * frac[cty, X]) / denom for X in (E, I, R))

        st = self.state
        new = st.copy()
        dead = np.zeros(len(st), bool)

        def draw(idx, rates, dests):
            """Competing risks over dt for agents idx with rate columns `rates` to `dests`."""
            if len(idx) == 0:
                return
            lam = np.column_stack(rates)
            tot = lam.sum(1)
            leave = rng.random(len(idx)) < -np.expm1(-tot * p.dt)
            if not leave.any():
                return
            j = idx[leave]
            cum = np.cumsum(lam[leave], axis=1)
            u = rng.random(len(j)) * cum[:, -1]
            k = (u[:, None] > cum).sum(1)
            d = np.asarray(dests)[k]
            new[j] = np.where(d == 4, S, d)
            dead[j[d == 4]] = True

        iS = np.flatnonzero(st == S)
        c = cty[iS]
        draw(iS, [p.beta0 * (piI[iS] + p.theta_E * piE[iS]) * v[c] ** p.k_beta * F[c, RI["beta"]],
                  p.gamma_rx * p.chi_rx * a[c] * F[c, RI["gamma"]]], [E, I])

        iE = np.flatnonzero(st == E)
        c = cty[iE]
        draw(iE, [p.sigma0 * v[c] ** p.k_sigma * F[c, RI["sigma"]],
                  p.eps0 * F[c, RI["eps"]]], [I, S])

        iI = np.flatnonzero(st == I)
        c = cty[iI]
        mu = (p.mu_H * (1 - phi[c]) + p.mu_F * phi[c]) * F[c, RI["mu"]] * np.exp(p.w_iso * S_cov[c])
        self.hazard = np.bincount(c, weights=mu * p.dt, minlength=self.C)   # expected deaths this step
        draw(iI, [p.zeta0 * F[c, RI["zeta"]] + p.gamma_p * piR[iI],
                  p.nu0 * F[c, RI["nu"]],
                  mu], [R, S, 4])

        iR = np.flatnonzero(st == R)
        c = cty[iR]
        draw(iR, [p.delta0 * F[c, RI["delta"]],
                  (p.rho0 * v[c] ** p.k_rho + p.rho_p * piI[iR]) * F[c, RI["rho"]]], [S, I])

        self.state = new
        d_idx = np.flatnonzero(dead)
        if len(d_idx):
            sp.xy[d_idx] = sp.place(sp.county[d_idx])
            sp.cell_id[d_idx] = sp.cell_of(sp.xy[d_idx])
        return np.bincount(cty[d_idx], minlength=self.C)

    def run(self, progress=None):
        """Burn-in, then t0..t_end. Returns a dict of recorded arrays."""
        p, inp = self.p, self.inp
        for k in range(int(round(p.t_burn / p.dt))):
            self.step(p.t0 - p.t_burn + k * p.dt, frozen_t=p.t0)
        self.env.demand0 = None  # re-anchor v* to the post-burn-in demand

        n_steps = int(round((p.t_end - p.t0) / p.dt))
        times = p.t0 + np.arange(n_steps) * p.dt
        Ys = p.t_end - p.t0          # simulated years; may run past the data (inputs held at the last data year)
        deaths_cy = np.zeros((self.C, Ys))
        hazard_cy = np.zeros((self.C, Ys))           # expected deaths (smooth version of deaths_cy)
        classes_cy = np.zeros((self.C, Ys, 4))       # mean class counts over the year
        steps_in_y = np.zeros(Ys)
        v_cy = np.zeros((self.C, Ys))
        phi_cy = np.zeros((self.C, Ys))
        series = np.zeros((n_steps, self.n_states, 4))   # class counts per state per step
        series_deaths = np.zeros((n_steps, self.n_states))
        so = inp.state_of
        for k, t in enumerate(times):
            d = self.step(t)
            y = min(int(t + 1e-9) - p.t0, Ys - 1)
            deaths_cy[:, y] += d
            hazard_cy[:, y] += self.hazard
            cnt = np.bincount(self.space.county * 4 + self.state, minlength=self.C * 4).reshape(self.C, 4)
            classes_cy[:, y] += cnt
            steps_in_y[y] += 1
            v_cy[:, y] += self.env.v
            phi_cy[:, y] += self.env.phi
            np.add.at(series[k], so, cnt)
            np.add.at(series_deaths[k], so, d)
            if progress and k % 52 == 0:
                progress(k / n_steps)
        classes_cy /= steps_in_y[None, :, None]
        v_cy /= steps_in_y
        phi_cy /= steps_in_y
        return dict(times=times, deaths_cy=deaths_cy, hazard_cy=hazard_cy, classes_cy=classes_cy, v_cy=v_cy, phi_cy=phi_cy,
                    series=series, series_deaths=series_deaths, n_c=self.space.n_c)
