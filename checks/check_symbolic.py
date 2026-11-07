"""Symbolic checks of the RK8 benchmark (paper §2.2), the infinite well (§4.2) and the linearized two-body step (§3.2) (Nov 6, 2026).

Run from the repository's root: python -I checks/check_symbolic.py   (sympy, mpmath; a minute or two)
"""
import sympy as sp
import mpmath as mp

h, t, x, z = sp.symbols('h t x z')


def section(title):
    print('\n== ' + title)


# ---------------------------------------------------------------- §2.2: RK8
section('§2.2: the ten-stage RK8 of the paper (Eqs. 28-39), with -234 and with -231')
R = sp.Rational
c = [0, R(4, 27), R(2, 9), R(1, 3), R(1, 2), R(2, 3), R(1, 6), 1, R(5, 6), 1]


def tableau(a70, a80=-127):
    A = [[0] * 10 for _ in range(10)]
    A[1][0] = R(4, 27)
    A[2][:2] = [R(1, 18), R(3, 18)]
    A[3][:3] = [R(1, 12), 0, R(3, 12)]
    A[4][:4] = [R(1, 8), 0, 0, R(3, 8)]
    A[5][:5] = [R(13, 54), 0, R(-27, 54), R(42, 54), R(8, 54)]
    A[6][:6] = [R(389, 4320), 0, R(-54, 4320), R(966, 4320), R(-824, 4320), R(243, 4320)]
    A[7][:7] = [R(a70, 20), 0, R(81, 20), R(-1164, 20), R(656, 20), R(-122, 20), R(800, 20)]
    A[8][:8] = [R(a80, 288), 0, R(18, 288), R(-678, 288), R(456, 288), R(-9, 288), R(576, 288), R(4, 288)]
    A[9][:9] = [R(1481, 820), 0, R(-81, 820), R(7104, 820), R(-3376, 820), R(72, 820), R(-5040, 820),
                R(-60, 820), R(720, 820)]
    b = [R(v, 840) for v in (41, 0, 0, 27, 272, 27, 216, 0, 216, 41)]
    return sp.Matrix(A), sp.Matrix(b)


def stability(A, b):
    """R(z) = 1 + z b^T (I - zA)^{-1} 1 for y' = lambda y, z = lambda h."""
    n = A.shape[0]
    return sp.expand(1 + z * (b.T * (sp.eye(n) - z * A).inv() * sp.ones(n, 1))[0])


for a70 in (-234, -231):
    A, b = tableau(a70)
    rowsum7 = sum(A[7, :])
    Rz = stability(A, b)
    err = sp.expand(Rz - sum(z**k / sp.factorial(k) for k in range(11)))
    lead = sp.Poly(err, z).terms()[-1]
    print(f'a70 = {a70}: row 7 sums to {rowsum7} (node c7 = 1); R(z) - exp(z) through z^10 '
          f'= {sp.factor(err)}; lowest term {lead[1]} z^{lead[0][0]}')
    # oscillator: x'' = -x is y' = J y, J = [[0,1],[-1,0]], eigenvalues +-i: G = Re/Im of R(ih)
    Rih = sp.expand(Rz.subs(z, sp.I * h))
    G = sp.Matrix([[sp.re(Rih), sp.im(Rih)], [-sp.im(Rih), sp.re(Rih)]]).applyfunc(
        lambda e: sp.expand(sp.simplify(e.subs(h, sp.Symbol('hh', positive=True))).subs(sp.Symbol('hh', positive=True), h)))
    detG = sp.expand(G.det())
    print(f'   1 - det G = {sp.series(1 - detG, h, 0, 14).removeO()}')
    f = sp.lambdify(h, 1 - detG, 'mpmath')
    mp.mp.dps = 50
    hmax = mp.findroot(lambda s: f(s) - mp.mpf(2)**-53, 0.08 if a70 == -231 else 0.0009)
    print(f'   1 - det G = 2^-53 at h = {mp.nstr(hmax, 6)}')
    print(f'   coefficient of z^4 in R(z): {Rz.coeff(z, 4)} (paper, for -234: 1397/33600; exact 1/24)')


section('§2.2: order of the corrected RK8 on a nonlinear problem (numerical, 60 digits)')
mp.mp.dps = 60


def rk_step(A, b, f, y, hh):
    k = []
    for i in range(10):
        yi = y + hh * sum(mp.mpf(A[i, j]) * k[j] for j in range(i))
        k.append(f(yi))
    return y + hh * sum(mp.mpf(b[j]) * k[j] for j in range(10))


# y' = 1 + y^2 (y = tan t): nonlinear and non-autonomous-free; exact y(T) = tan(T)
fnl = lambda y: 1 + y**2
T = mp.mpf('0.5')
for a70 in (-234, -231):
    A, b = tableau(a70)
    errs = []
    for n in (4, 8, 16, 32):
        y = mp.mpf(0)
        for _ in range(n):
            y = rk_step(A, b, fnl, y, T / n)
        errs.append(abs(y - mp.tan(T)))
    orders = [mp.nstr(mp.log(errs[i] / errs[i + 1], 2), 4) for i in range(3)]
    print(f'a70 = {a70}: global errors {[mp.nstr(e, 3) for e in errs]}; observed orders {orders}')

# classical order conditions of order 8 that fail (a sample: the bushy and tall trees)
A, b = tableau(-231)
cv = sp.Matrix(c)
print('order-8 conditions, corrected tableau: sum b c^7 =', sp.nsimplify((b.T * cv.applyfunc(lambda s: s**7))[0]),
      '(1/8 needed); b^T A^6 c =', (b.T * A**6 * cv)[0], f'({sp.Rational(1, sp.factorial(8))} needed)')
# count the failed conditions up to order 8 by B-series (rooted trees), using sympy only
from itertools import combinations_with_replacement


def trees(n, _memo={}):
    """Rooted trees with n nodes as sorted tuples of child trees."""
    if n in _memo:
        return _memo[n]
    if n == 1:
        _memo[1] = [()]
        return _memo[1]
    out = set()

    def parts(m, maxpart):
        if m == 0:
            yield []
            return
        for p in range(min(m, maxpart), 0, -1):
            for rest in parts(m - p, p):
                yield [p] + rest
    for part in parts(n - 1, n - 1):
        pools = [trees(p) for p in part]

        def build(i, acc):
            if i == len(part):
                out.add(tuple(sorted(acc)))
                return
            for tr in pools[i]:
                build(i + 1, acc + [tr])
        build(0, [])
    _memo[n] = sorted(out)
    return _memo[n]


def size(tr):
    return 1 + sum(size(s) for s in tr)


def gamma(tr):
    g = size(tr)
    for s in tr:
        g *= gamma(s)
    return g


def phi_vec(tr, A):
    """Vector Phi_i(tree) for the stage values: product over children of A Phi(child)."""
    n = A.shape[0]
    v = sp.ones(n, 1)
    for s in tr:
        w = A * phi_vec(s, A)
        v = sp.Matrix([v[i] * w[i] for i in range(n)])
    return v


for a70 in (-234, -231):
    A, b = tableau(a70)
    fails = {}
    for order in range(1, 9):
        bad = 0
        for tr in trees(order):
            if (b.T * phi_vec(tr, A))[0] != sp.Rational(1, gamma(tr)):
                bad += 1
        fails[order] = (bad, len(trees(order)))
    print(f'a70 = {a70}: failed order conditions by order (failed, of): {fails}')


# ---------------------------------------------- §4.2: the well as Hermite elements
section('§4.2: the infinite well, P = 2M and Q = -2K (Eq. 160)')
s = sp.symbols('s')  # x = s h on one interval
N1 = 1 - 3 * s**2 + 2 * s**3
N2 = 3 * s**2 - 2 * s**3
N3 = h * (s - 2 * s**2 + s**3)
N4 = h * (-s**2 + s**3)
Nb = [N1, N2, N3, N4]  # order (u_i, u_{i+1}, u'_i, u'_{i+1}), as in Eq. 160
psi = sp.symbols('p0 p1 dp0 dp1')
phi = sp.symbols('f0 f1 df0 df1')
Psi = sum(p * n for p, n in zip(psi, Nb))
Phi = sum(p * n for p, n in zip(phi, Nb))
d2 = lambda e: sp.diff(e, s, 2) / h**2
cost = sp.integrate(sp.expand((d2(Psi) + Phi)**2) * h, (s, 0, 1))
grad = [sp.expand(sp.diff(cost, f)) for f in phi]
P = sp.Matrix(4, 4, lambda i, j: grad[i].coeff(phi[j]))
Q = sp.Matrix(4, 4, lambda i, j: grad[i].coeff(psi[j]))
M = sp.Matrix(4, 4, lambda i, j: sp.integrate(Nb[i] * Nb[j] * h, (s, 0, 1)))
K = sp.Matrix(4, 4, lambda i, j: sp.integrate(sp.diff(Nb[i], s) * sp.diff(Nb[j], s) / h, (s, 0, 1)))
P160 = sp.Matrix([[26 * h / 35, 9 * h / 35, 11 * h**2 / 105, -13 * h**2 / 210],
                  [9 * h / 35, 26 * h / 35, 13 * h**2 / 210, -11 * h**2 / 105],
                  [11 * h**2 / 105, 13 * h**2 / 210, 2 * h**3 / 105, -h**3 / 70],
                  [-13 * h**2 / 210, -11 * h**2 / 105, -h**3 / 70, 2 * h**3 / 105]])
Q160 = sp.Matrix([[-12 / (5 * h), 12 / (5 * h), R(-11, 5), R(-1, 5)],
                  [12 / (5 * h), -12 / (5 * h), R(1, 5), R(11, 5)],
                  [R(-1, 5), R(1, 5), -4 * h / 15, h / 15],
                  [R(-1, 5), R(1, 5), h / 15, -4 * h / 15]])
print('P from the cost equals Eq. 160:', sp.simplify(P - P160) == sp.zeros(4))
print('Q from the cost equals Eq. 160:', sp.simplify(Q - Q160) == sp.zeros(4))
print('P - 2M =', sp.simplify(P - 2 * M))
print('Q + 2K =', sp.simplify(Q + 2 * K), '(boundary terms: -2 dpsi_i in row phi_i, +2 dpsi_{i+1} in row phi_{i+1})')


def assemble(N):
    hh = R(1, N)
    n = 2 * (N + 1)  # unknowns ordered (u_0..u_N, u'_0..u'_N), as in Eq. 161
    PP, QQ = sp.zeros(n), sp.zeros(n)
    Pe, Qe = P.subs(h, hh), Q.subs(h, hh)
    for i in range(N):
        idx = [i, i + 1, N + 1 + i, N + 1 + i + 1]
        for a in range(4):
            for bb in range(4):
                PP[idx[a], idx[bb]] += Pe[a, bb]
                QQ[idx[a], idx[bb]] += Qe[a, bb]
    keep = [k for k in range(n) if k not in (0, N)]  # psi_0 = psi_N = 0: drop rows and columns
    return PP.extract(keep, keep), QQ.extract(keep, keep)


for N in (1, 2, 4, 8):
    PP, QQ = assemble(N)
    H = -PP.inv() * QQ
    MM = PP / 2
    sym = sp.simplify(MM * H - (MM * H).T) == sp.zeros(MM.shape[0])
    mp.mp.dps = 40
    Hm = mp.matrix([[mp.mpf(sp.Rational(H[i, j]).p) / sp.Rational(H[i, j]).q for j in range(H.shape[1])]
                    for i in range(H.shape[0])])
    evc = mp.eig(Hm)[0]
    assert max(abs(mp.im(e)) for e in evc) < mp.mpf(10)**-20
    ev = sorted(float(mp.re(e)) for e in evc)
    exact = [(k * mp.pi)**2 for k in range(1, len(ev) + 1)]
    above = all(e >= float(x) * (1 - 1e-12) for e, x in zip(ev, exact))
    print(f'N = {N}: M H symmetric: {sym}; eigenvalues above exact: {above}; '
          f'ground state {ev[0]:.10f}, fractional error {mp.nstr((ev[0] - exact[0]) / exact[0], 3)}'
          + (f'; eigenvalues {[round(e, 6) for e in ev]}' if N == 1 else ''))


# ------------------------------------------------------- §3.2: Hermite corrector
section('§3.2: least-squares fit of the cubic r(t) to the Hermite interpolant of a(t)')
a0, a1, j0, j1, r0, v0, A3, B2 = sp.symbols('a0 a1 j0 j1 r0 v0 A B')
u = t / h
H1 = 1 - 3 * u**2 + 2 * u**3
H2 = 3 * u**2 - 2 * u**3
H3 = h * (u - 2 * u**2 + u**3)
H4 = h * (-u**2 + u**3)
acc = a0 * H1 + a1 * H2 + j0 * H3 + j1 * H4  # cubic interpolant of the acceleration
r = r0 + v0 * t + B2 * t**2 + A3 * t**3
eps = sp.integrate(sp.expand((sp.diff(r, t, 2) - acc)**2), (t, 0, h))
sol = sp.solve([sp.diff(eps, A3), sp.diff(eps, B2)], [A3, B2], dict=True)[0]
rh = sp.expand(r.subs(sol).subs(t, h))
vh = sp.expand(sp.diff(r, t).subs(sol).subs(t, h))
print('r_h =', sp.collect(rh, [a0, a1, j0, j1]))
print('v_h =', sp.collect(vh, [a0, a1, j0, j1]))
# the interpolant integrated once and twice
v_int = v0 + sp.integrate(acc, (t, 0, h))
r_int = r0 + v0 * h + sp.integrate(sp.integrate(acc, (t, 0, t)).subs(t, sp.Symbol('tt')), (sp.Symbol('tt'), 0, h))
print('v_h equals the interpolant integrated once:', sp.simplify(vh - v_int) == 0)
print('r_h equals the interpolant integrated twice:', sp.simplify(rh - sp.expand(r_int)) == 0)
# Makino & Aarseth (1992) corrector, as written in Makino (1991) / Hut, Makino & McMillan (1995):
# v1 = v0 + h(a1 + a0)/2 - h^2 (j1 - j0)/12 ; r1 = r0 + h(v1 + v0)/2 - h^2 (a1 - a0)/12
vMA = v0 + h * (a1 + a0) / 2 - h**2 * (j1 - j0) / 12
rMA = r0 + h * (vMA + v0) / 2 - h**2 * (a1 - a0) / 12
print('v_h equals the time-symmetric Hermite corrector:', sp.simplify(vh - vMA) == 0)
print('r_h minus the time-symmetric corrector r1 = r0 + h(v1+v0)/2 - h^2(a1-a0)/12:',
      sp.factor(sp.simplify(rh - rMA)))
