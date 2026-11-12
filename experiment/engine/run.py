"""Run an ensemble and save it for the viewer.

  python -m engine.run --name baseline --N 300000 --M 3

Saves experiment/runs/<name>/ with params.json, results.npz, summary.json.
"""
import argparse
import json
import time
from pathlib import Path

import numpy as np

from .calibrate import calibrate, equilibrium, match_elasticities, reproduction_number
from .inputs import load_inputs
from .model import Simulation
from .params import Params
from . import scoring

RUNS = Path(__file__).resolve().parents[1] / "runs"


def run(p, name, log=print):
    inp = load_inputs()
    a_mean = float(np.average(inp.a_rx[:, 0], weights=inp.pop[:, 0]))
    eq0 = calibrate(p, a_mean)
    match_elasticities(p, a_mean)
    mf = equilibrium(p, a_mean)
    log(f"calibrated beta0={p.beta0:.4f} chi_rx={p.chi_rx:.4f} R_c={reproduction_number(p):.3f} "
        f"| mean-field equilibrium deaths/10k={mf['deaths_per_10k']:.3f}, i*={mf['i']:.4f}")

    seeds = np.random.SeedSequence(p.seed).spawn(p.M)
    outs = []
    for m in range(p.M):
        t = time.time()
        sim = Simulation(inp, p, np.random.default_rng(seeds[m]), eq0)
        outs.append(sim.run())
        log(f"member {m + 1}/{p.M} done in {time.time() - t:.0f}s ({len(sim.state):,} agents)")

    res = {k: np.stack([o[k] for o in outs]) for k in outs[0] if k not in ("times", "n_c")}
    res["times"] = outs[0]["times"]
    res["n_c"] = outs[0]["n_c"]

    Yd = len(inp.years)   # score only the years that have data
    exp_real = scoring.expected_real_deaths(res["hazard_cy"][:, :, :Yd], res["n_c"], inp.pop)   # [M, C, Y]
    ll = scoring.state_interval_loglik(inp, exp_real.mean(0))
    summary = dict(name=name, agents=int(res["n_c"].sum()), members=p.M, beta0=p.beta0, chi_rx=p.chi_rx,
                   R_c=reproduction_number(p), mean_field=mf, state_loglik=float(ll))

    out = RUNS / name
    out.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(out / "results.npz", **res)
    (out / "params.json").write_text(json.dumps(p.to_dict(), indent=2))
    (out / "summary.json").write_text(json.dumps(summary, indent=2, default=float))
    log(f"saved {out}  (state-year interval log-lik {ll:.3f})")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", default="baseline")
    ap.add_argument("--N", type=int, default=Params.N_sim)
    ap.add_argument("--M", type=int, default=Params.M)
    ap.add_argument("--seed", type=int, default=Params.seed)
    ap.add_argument("--rho", type=float, default=Params.rho_km)
    ap.add_argument("--t_end", type=int, default=Params.t_end, help="simulate up to this year (exclusive)")
    a = ap.parse_args()
    run(Params(N_sim=a.N, M=a.M, seed=a.seed, rho_km=a.rho, t_end=a.t_end), a.name)


if __name__ == "__main__":
    main()
