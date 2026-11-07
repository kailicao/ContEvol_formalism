"""Whether least squares or the truncation makes the step linear (second check): what the adequate truncation (§3.1) puts into the linear system,
for the comparison with linearly implicit (Rosenbrock and W-) methods.

The two-body problem r'' = f(r), f = -r/r^3, with ContEvol's cubic r = r0 + v0 t + B t^2 + A t^3.
The residual r'' - f(r), expanded to t^3 about the step's start, is affine in A and B; its
coefficient matrices are compared with the Jacobian J = f_r(r0), the Hessian contracted with the
velocity, H[v0, .] = dJ/dt at t = 0, and the third derivative. Exact symbolic arithmetic at a
generic state.

Run from the repository's root (a few seconds):
    python checks/check_truncation.py
"""
import sympy as sp

t = sp.symbols('t')
rx, ry, vx, vy = sp.Rational(7, 5), sp.Rational(-3, 10), sp.Rational(1, 4), sp.Rational(9, 10)
a1, a2, b1, b2 = sp.symbols('a1 a2 b1 b2')
A, B = sp.Matrix([a1, a2]), sp.Matrix([b1, b2])
r0, v0 = sp.Matrix([rx, ry]), sp.Matrix([vx, vy])

X, Y = sp.symbols('X Y')
f = -sp.Matrix([X, Y]) / (X**2 + Y**2)**sp.Rational(3, 2)
J = f.jacobian([X, Y])
at0 = {X: rx, Y: ry}
J0 = J.subs(at0)
dJdt = sum((J.diff(s) * c for s, c in ((X, vx), (Y, vy))), sp.zeros(2, 2)).subs(at0)  # H[v0, .]

r = r0 + v0 * t + B * t**2 + A * t**3
ftr = f.subs({X: r[0], Y: r[1]}, simultaneous=True)
fser = ftr.applyfunc(lambda e: sp.series(e, t, 0, 4).removeO())
res = (r.diff(t, 2) - fser).applyfunc(sp.expand)
C = [res.applyfunc(lambda e: e.coeff(t, k)) for k in range(4)]

unk = [a1, a2, b1, b2]
print('residual through t^3 affine in A, B:',
      all(sp.Poly(c, *unk).total_degree() <= 1 for k in range(4) for c in C[k]))
dA = [C[k].jacobian([a1, a2]) for k in range(4)]
dB = [C[k].jacobian([b1, b2]) for k in range(4)]
I2 = sp.eye(2)
checks = {
    'dC0/dB = 2I': dB[0] == 2 * I2, 'dC1/dA = 6I': dA[1] == 6 * I2,
    'dC2/dB = -J': sp.simplify(dB[2] + J0) == sp.zeros(2, 2),
    'dC3/dA = -J': sp.simplify(dA[3] + J0) == sp.zeros(2, 2),
    'dC3/dB = -H[v0,.]': sp.simplify(dB[3] + dJdt) == sp.zeros(2, 2),
    'C0, C1 free of A at t^0 and B at t^1 beyond the above':
        dA[0] == sp.zeros(2, 2) and dB[1] == sp.zeros(2, 2) and dA[2] == sp.zeros(2, 2),
}
for k, v in checks.items():
    print(f'  {k}: {v}')
# the known part: Taylor coefficients of f along the straight line r0 + v0 t, which need
# f, J v0, H[v0, v0] and the third derivative T[v0, v0, v0]
line = f.subs({X: rx + vx * t, Y: ry + vy * t}, simultaneous=True)
lser = line.applyfunc(lambda e: sp.series(e, t, 0, 4).removeO())
known = res.subs({a1: 0, a2: 0, b1: 0, b2: 0})
print('known part = minus f along r0 + v0 t, to t^3 (needs f, J, H, T at r0):',
      sp.simplify(known + lser) == sp.zeros(2, 1))
# the next order: at t^4 the products appear
f4 = ftr.applyfunc(lambda e: sp.expand(sp.series(e, t, 0, 5).removeO()).coeff(t, 4))
print('at t^4 the coefficients are quadratic in A, B (why §3.1 stops at t^3):',
      max(sp.Poly(c, *unk).total_degree() for c in f4) == 2)
