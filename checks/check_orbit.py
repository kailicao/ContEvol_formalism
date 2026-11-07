"""The two-body integrators of the paper (2405.05188, §3.4): a check of Nov 6, 2026.

Runs the paper's own integrators, executed from the code cells of
celestial/2024.04.30 ContEvol eccentric timeit.ipynb (LF, LFCC, RK4, RK4CC, CE1, CE1CC),
against the exact Kepler orbit (e = 63/64, a = 16, from pericenter), and
(1) checks symbolically that CE1's derivatives are the acceleration, jerk, snap and
    crackle of the Kepler motion, so that its step is the Taylor series through
    h^5/120 in position and h^4/24 in velocity (§3.2, Eq. 101);
(2) gives each integrator's largest position and velocity error over the first orbit
    (t <= 128 pi) and its run time (best of 5), for the comparison at equal cost (§3.4, Table 1).

Usage, from the repository's root:
    python3 -I checks/check_orbit.py .
"""
import json
import sys
from pathlib import Path
from time import perf_counter

import numpy as np
import sympy as sp

repo = Path(sys.argv[1])
nb = json.load(open(repo / 'celestial' / '2024.04.30 ContEvol eccentric timeit.ipynb'))

# Execute the notebook's definitions: imports, constants, the correction, the
# equation of motion and the six integrators; skip plots, magics and runs.
KEEP = ('from', 'import','def ', '@njit', 'a_test', 'EM_test')
ns = {}
for cell in nb['cells']:
    if cell['cell_type'] != 'code':
        continue
    src = ''.join(cell['source'])
    if '%' in src or 'plt.' in src or 'matplotlib' in src or 'test_integrator(' in src \
            or 'run_integrator(' in src or 'M_exact' in src or 'exact_plot' in src:
        continue
    if src.strip().startswith(KEEP):
        exec(src, ns)

a, e, b = ns['a_test'], ns['e_test'], ns['b_test']
r0, v0 = ns['r0_test'], ns['v0_test']
n_mean = a ** -1.5  # mean motion, 1/64
T = 2 * np.pi / n_mean  # period, 128 pi


def exact(t):
    """Exact position and velocity at times t (Kepler's equation by bisection, then
    two Newton steps; Newton alone from E = M fails near pericenter at e = 63/64)."""
    M = n_mean * t
    lo, hi = M - 1.0, M + 1.0  # E - M = e sin E lies in [-e, e]
    for _ in range(80):
        mid = (lo + hi) / 2
        above = mid - e * np.sin(mid) - M > 0
        hi = np.where(above, mid, hi)
        lo = np.where(above, lo, mid)
    E = (lo + hi) / 2
    for _ in range(2):
        E -= (E - e * np.sin(E) - M) / (1 - e * np.cos(E))
    assert np.max(np.abs(E - e * np.sin(E) - M)) < 1e-13
    Edot = n_mean / (1 - e * np.cos(E))
    r = np.stack([a * (np.cos(E) - e), b * np.sin(E)], axis=1)
    v = np.stack([-a * np.sin(E) * Edot, b * np.cos(E) * Edot], axis=1)
    return r, v


# (1) CE1's derivatives against symbolic differentiation along the motion.
t = sp.symbols('t')
x, y = sp.Function('x')(t), sp.Function('y')(t)
acc = [-x / (x**2 + y**2) ** sp.Rational(3, 2), -y / (x**2 + y**2) ** sp.Rational(3, 2)]
rng = np.random.default_rng(1)
state = {'x': 0.7, 'y': -0.4, 'vx': 0.3, 'vy': 1.1}
derivs = [acc]
for k in range(3):
    nxt = []
    for comp in derivs[-1]:
        d = sp.diff(comp, t)
        # replace x'' and y'' by the acceleration, x' and y' stay as velocities
        d = d.subs({sp.Derivative(x, (t, 2)): acc[0], sp.Derivative(y, (t, 2)): acc[1]})
        nxt.append(d)
    derivs.append(nxt)
vals = {sp.Derivative(x, t): state['vx'], sp.Derivative(y, t): state['vy']}
sym = np.array([[float(c.subs(vals).subs({x: state['x'], y: state['y']})) for c in d] for d in derivs])

# One CE1 step with tiny h exposes h_powers; instead read its derivatives from two steps:
# the step is r + v h + sum_k D_k h^(k+2)/(k+2)!, v + D_0 h + sum_k D_(k+1) h^(k+2)/(k+2)!,
# so fit the step's increments at several h (exact polynomials in h) to recover D_k.
CE1 = ns['CE1_integrator']
rr = np.array([state['x'], state['y']])
vv = np.array([state['vx'], state['vy']])
hs = np.array([0.02, 0.04, 0.06, 0.08, 0.10, 0.12])
dr = np.array([CE1(1, 1, h, rr, vv)[1, :2] - rr - vv * h for h in hs])
dv = np.array([CE1(1, 1, h, rr, vv)[1, 2:] - vv for h in hs])
Vr = np.stack([hs**k for k in range(2, 6)], axis=1)  # h^2 .. h^5
Vv = np.stack([hs**k for k in range(1, 5)], axis=1)  # h^1 .. h^4
cr = np.linalg.solve(Vr[:4], dr[:4]) * np.array([2, 6, 24, 120])[:, None]
cv = np.linalg.solve(Vv[:4], dv[:4]) * np.array([1, 2, 6, 24])[:, None]
res_r = np.max(np.abs(Vr[4:] @ (cr / np.array([2, 6, 24, 120])[:, None]) - dr[4:]))
res_v = np.max(np.abs(Vv[4:] @ (cv / np.array([1, 2, 6, 24])[:, None]) - dv[4:]))
print('(1) CE1 against the Taylor series, at r = (0.7, -0.4), v = (0.3, 1.1)')
print('    symbolic a, j, s, c:\n', sym)
print('    from CE1 position step (h^2/2 .. h^5/120 coefficients):\n', cr)
print('    from CE1 velocity step (h .. h^4/24 coefficients):\n', cv)
print(f'    largest relative difference: position {np.max(np.abs(cr - sym) / np.abs(sym)):.1e},'
      f' velocity {np.max(np.abs(cv - sym) / np.abs(sym)):.1e};'
      f' polynomial in h of degree 5 (position), 4 (velocity): residuals {res_r:.1e}, {res_v:.1e}')

# (2) The six integrators on the orbit.
print('\n(2) Largest error over the first orbit (t <= 128 pi), and run time (best of 5)')
names = ['LF', 'LFCC', 'RK4', 'RK4CC', 'CE1', 'CE1CC']
print(f'{"method":7s} {"h":>7s} {"max|dr|":>10s} {"max|dv|":>10s} {"max|dr| all":>12s} {"time (s)":>9s}')
results = {}
for name in names:
    integ = ns[name + '_integrator']
    for h_res in [16, 64, 256, 1024]:
        Nstep = 432 * min(h_res, 64)
        Nsub = max(h_res // 64, 1)
        h = 1 / h_res
        integ(Nstep, Nsub, h, r0, v0)  # compile
        best = np.inf
        for _ in range(5):
            t0 = perf_counter()
            hist = integ(Nstep, Nsub, h, r0, v0)
            best = min(best, perf_counter() - t0)
        times = np.arange(Nstep + 1) * h * Nsub
        re, ve = exact(times)
        dr_ = np.linalg.norm(hist[:, :2] - re, axis=1)
        dv_ = np.linalg.norm(hist[:, 2:] - ve, axis=1)
        first = times <= T
        results[name, h_res] = (dr_[first].max(), dv_[first].max(), dr_.max(), best)
        print(f'{name:7s} 1/{h_res:<5d} {dr_[first].max():10.2e} {dv_[first].max():10.2e}'
              f' {dr_.max():12.2e} {best:9.4f}')

print('\nCost ratios to RK4 at each step (this computer):')
for name in names:
    print(name, ' '.join(f'{results[name, hr][3] / results["RK4", hr][3]:.2f}' for hr in [16, 64, 256, 1024]))
