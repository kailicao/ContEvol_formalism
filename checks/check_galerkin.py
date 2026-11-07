"""Whether least squares or the truncation makes the step linear (first check), and Sharma, Patil & Woolsey (2022,
arXiv:2201.07327, §4.2): the one-step Galerkin step with shifted-Legendre test functions
(P0 = 1, P1 = (2t - h)/h) on ContEvol's cubic trial space, against ContEvol's least-squares step.

(a) The oscillator x'' = -x: the Galerkin step's matrix G, its determinant and its error against
    cos/sin, beside the least-squares step (Eq. 17) and the matrix Sharma et al. print.
(b) The two-body problem r'' = -r/r^3 under ContEvol's adequate truncation (§3.1: the right-hand
    side expanded to t^3 about the step's start, so the residual is a cubic affine in A and B):
    is the Galerkin step linear, and how do it and the least-squares step differ from the
    Taylor series of the motion (Eq. 101)? Exact rational arithmetic at a generic state; the
    linear systems are solved as power series in h.

Run from the repository's root (about a minute):
    python checks/check_galerkin.py
"""
import sympy as sp

h, t = sp.symbols('h t')
TESTS = (sp.Integer(1), 2 * t - h)  # shifted Legendre P0, h P1


def integ(poly):
    """int_0^h of a polynomial in t."""
    p = sp.Poly(sp.expand(poly), t)
    return sum(c * h**(m[0] + 1) / (m[0] + 1) for m, c in zip(p.monoms(), p.coeffs()))


# ------------------------------------------------------------- (a) the oscillator
x0, v0, A, B = sp.symbols('x0 v0 A B')
x = x0 + v0 * t + B * t**2 + A * t**3
R = sp.diff(x, t, 2) + x
exact = sp.Matrix([[sp.cos(h), sp.sin(h)], [-sp.sin(h), sp.cos(h)]])


def step_matrix(sol):
    xh, vh = x.subs(sol).subs(t, h), sp.diff(x, t).subs(sol).subs(t, h)
    return sp.Matrix([[sp.diff(xh, c) for c in (x0, v0)], [sp.diff(vh, c) for c in (x0, v0)]]).applyfunc(sp.factor)


gal = sp.solve([integ(R * w) for w in TESTS], [A, B], dict=True)[0]
eps = integ(R**2)
ls = sp.solve([sp.diff(eps, A), sp.diff(eps, B)], [A, B], dict=True)[0]
ser = lambda M, n=8: M.applyfunc(lambda e: sp.series(e, h, 0, n).removeO())
for name, sol in (('Galerkin', gal), ('least squares', ls)):
    G = step_matrix(sol)
    print(f'(a) {name}: G (rows x, v; columns x0, v0) =', G)
    print('    1 - det G =', sp.series(sp.simplify(1 - G.det()), h, 0, 10), '; G - exact =', ser(G - exact))
J = sp.Matrix([[0, 1], [-1, 0]])
I2 = sp.eye(2)
pade = ((I2 - J * h / 2 + J**2 * h**2 / 12).inv() * (I2 + J * h / 2 + J**2 * h**2 / 12)).applyfunc(sp.factor)
print('    Pade (2,2) (the 4th-order Hermite-Obreshkov step on this problem) =', pade, '; - exact =', ser(pade - exact))
d = h**4 + 16 * h**2 + 240
S = sp.Matrix([[2 * (h**4 - 52 * h**2 + 120) / d, h * (h**4 - 132 * h**2 + 1440) / (6 * d)],
               [-h * (h**4 - 130 * h**2 + 1200) / (5 * d), -(h**6 - 240 * h**4 + 6240 * h**2 - 14400) / (60 * d)]])
print("    Sharma et al.'s printed Galerkin Jacobian: det =", sp.simplify(S.det()),
      '; its diagonal minus cos h:', [sp.series(S[i, i] - sp.cos(h), h, 0, 6).removeO() for i in range(2)])

# ------------------------------------------------------------- (b) the two-body problem
r0 = sp.Matrix([sp.Rational(9, 10), sp.Rational(-2, 7)])
w0 = sp.Matrix([sp.Rational(1, 5), sp.Rational(11, 10)])
unk = sp.symbols('Ax Ay Bx By')
r = r0 + w0 * t + sp.Matrix(unk[2:]) * t**2 + sp.Matrix(unk[:2]) * t**3
N = 3  # adequate expansion: the right-hand side to t^3


def trunc(e, n):
    p = sp.Poly(sp.expand(e), t)
    return sum(p.coeff_monomial(t**k) * t**k for k in range(n + 1))


c0 = r0.dot(r0)


def inv_r3(z, n):
    """|z|^-3 as a series in t, to t^n, for z = r0 + O(t)."""
    u = sp.expand(trunc(z.dot(z), n) / c0 - 1)
    return trunc(sum(sp.binomial(sp.Rational(-3, 2), k) * u**k for k in range(n + 1)), n) * c0**sp.Rational(-3, 2)


F = (-r * inv_r3(r, N)).applyfunc(lambda e: trunc(e, N))
Rv = (sp.diff(r, t, 2) - F).applyfunc(sp.expand)
print('\n(b) truncated residual affine in A, B:',
      all(sp.Poly(Rv[i].coeff(t, k), *unk).total_degree() <= 1 for i in range(2) for k in range(N + 1)))
gal_eqs = [integ(Rv[i] * w) for i in range(2) for w in TESTS]
eps2 = integ(Rv[0]**2) + integ(Rv[1]**2)
ls_eqs = [sp.diff(eps2, c) for c in unk]

# the exact Taylor series of the motion, to t^8, by undetermined coefficients
K = 8
cs = [r0, w0]
for k in range(2, K + 1):
    z = sum((cs[j] * t**j for j in range(k)), sp.zeros(2, 1))
    cs.append((-z * inv_r3(z, k - 2)).applyfunc(lambda e: sp.expand(e).coeff(t, k - 2)) / (k * (k - 1)))
rex = sum((cs[k] * h**k for k in range(K + 1)), sp.zeros(2, 1))
vex = sp.diff(rex, h)


def series_solve(eqs, order=8):
    M, b = sp.linear_eq_to_matrix([sp.expand(e) for e in eqs], unk)
    M = M.applyfunc(sp.expand)
    # A multiplies t^3 and B t^2: solve for (h A, B), so that the columns are of one order in h
    M[:, :2] = (M[:, :2] / h).applyfunc(sp.expand)
    for i in range(4):  # divide each row by its lowest power of h
        lo = min(sp.Poly(e, h).monoms()[-1][0] for e in list(M.row(i)) + [b[i]] if e != 0)
        M[i, :] = (M[i, :] / h**lo).applyfunc(sp.expand)
        b[i] = sp.expand(b[i] / h**lo)
    Mk = [M.applyfunc(lambda e: e.coeff(h, k)) for k in range(order + 1)]
    bk = [b.applyfunc(lambda e: e.coeff(h, k)) for k in range(order + 1)]
    M0inv, X = Mk[0].inv(), []
    for k in range(order + 1):
        X.append(M0inv * (bk[k] - sum((Mk[k - j] * X[j] for j in range(k)), sp.zeros(4, 1))))
    sol = {unk[i]: sum(X[k][i] * h**k for k in range(order + 1)) for i in range(4)}
    for c in unk[:2]:
        sol[c] = sp.expand(sol[c] / h)
    return sol


print('    Galerkin equations linear in A, B:', all(sp.Poly(sp.expand(e), *unk).total_degree() <= 1 for e in gal_eqs))
steps = {}
for name, eqs in (('least squares', ls_eqs), ('Galerkin', gal_eqs)):
    sol = series_solve(eqs)
    rh = r.subs(sol).subs(t, h).applyfunc(sp.expand)
    vh = sp.diff(r, t).subs(sol).subs(t, h).applyfunc(sp.expand)
    steps[name] = (rh, vh)
    for lab, a, e, n in (('position', rh, rex, 8), ('velocity', vh, vex, 7)):
        dd = [[sp.expand(a[i] - e[i]).coeff(h, k) for k in range(n)] for i in range(2)]
        first = min(k for i in range(2) for k in range(n) if dd[i][k] != 0)
        print(f'    {name}, {lab}: step - Taylor starts at h^{first}; coefficients h^{first}..h^{n - 1}:',
              [dd[i][first:] for i in range(2)])
print('    Summary: each step\'s own coefficient as a fraction of the true Taylor coefficient (x, y):')
for name, (rh, vh) in steps.items():
    fr = lambda a, e, k: [sp.N(sp.expand(a[i]).coeff(h, k) / sp.expand(e[i]).coeff(h, k), 6) for i in range(2)]
    print(f'      {name}: position h^5 {fr(rh, rex, 5)}, h^6 {fr(rh, rex, 6)}, h^7 {fr(rh, rex, 7)};'
          f' velocity h^4 {fr(vh, vex, 4)}, h^5 {fr(vh, vex, 5)}')
