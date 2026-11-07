"""The linearized two-body step (paper §3.2, Eqs. 90-96) against Dittmann's Hermite corrector (2020, arXiv:2006.02501, Section 2).

Dittmann writes the Hermite corrector of order 2(p+1), after Nitadori (2016), as
the integral over the step of the Hermite interpolant of f (f = a for the velocity,
f = v for the position), with the even scaled derivatives F_0, F_2, ... at the
midpoint recovered from end-point averages and differences through the matrix A of
his Eq. for [F+_0, F-_1, F+_2, ...]; the position gets one more term, weighted by beta.
This script builds his p = 1 corrector symbolically and compares it with ContEvol's
least-squares step (r_h, v_h in terms of end-of-step values),
and checks his beta = 11/3 coefficients (-7/60, 1/60).

Run: python check_dittmann.py   (a few seconds)
"""
import sympy as sp

dt = sp.symbols("Delta_t", positive=True)
beta = sp.symbols("beta")
h = dt / 2
x0, x1, v0, v1, a0, a1, j0, j1 = sp.symbols("x0 x1 v0 v1 a0 a1 j0 j1")


def midpoint_even(fm, fp, dfm, dfp, d2fm=None, d2fp=None):
    """Even scaled midpoint derivatives F_0, F_2, F_4 from end-point data.

    Unknowns: the midpoint Taylor coefficients F_n = h^n f^(n)/n!, n = 0..5.
    End-point data f(+-h), h f'(+-h), (h^2/2) f''(+-h) are matched exactly.
    """
    F = sp.symbols("F0:6")
    s = sp.symbols("s")
    poly = sum(F[n] * s**n for n in range(6))  # f(t_mid + s h)
    eqs = [
        sp.Eq(poly.subs(s, -1), fm), sp.Eq(poly.subs(s, 1), fp),
        sp.Eq(sp.diff(poly, s).subs(s, -1), h * dfm),
        sp.Eq(sp.diff(poly, s).subs(s, 1), h * dfp),
    ]
    unknowns = list(F[:4])
    if d2fm is not None:
        eqs += [sp.Eq(sp.diff(poly, s, 2).subs(s, -1) / 2, h**2 * d2fm / 2),
                sp.Eq(sp.diff(poly, s, 2).subs(s, 1) / 2, h**2 * d2fp / 2)]
        unknowns = list(F)
    else:
        eqs = [e.subs({F[4]: 0, F[5]: 0}) for e in eqs]
    sol = sp.solve(eqs, unknowns, dict=True)[0]
    return [sol.get(F[n], 0) for n in (0, 2, 4)]


# Velocity, p = 1: f = a, data a, j at both ends (cubic interpolant).
F0, F2, _ = midpoint_even(a0, a1, j0, j1)
dv = sp.expand((F0 + F2 / 3) * dt)
# Position: f = v, data v, a, j at both ends (quintic); Dittmann's Delta x.
G0, G2, G4 = midpoint_even(v0, v1, a0, a1, j0, j1)
dx = sp.expand((G0 + G2 / 3 + beta * G4 / 5) * dt)

print("Dittmann p = 1, velocity:", sp.collect(dv, dt))
print("Dittmann p = 1, position:", sp.collect(dx, dt))

# ContEvol's least-squares step, written with end-of-step values:
dv_ce = dt * (a0 + a1) / 2 - dt**2 * (j1 - j0) / 12
dx_ts = dt * (v0 + v1) / 2 - dt**2 * (a1 - a0) / 12          # time-symmetric form
dx_ce = dx_ts + dt**2 * (2 * (a0 - a1) + dt * (j0 + j1)) / 120  # ContEvol's r_h

print("velocity: Dittmann - ContEvol =", sp.simplify(dv - dv_ce))
print("position, beta = 1: Dittmann - ContEvol =", sp.simplify(dx.subs(beta, 1) - dx_ce))
print("position, beta = 1:", sp.collect(sp.expand(dx.subs(beta, 1)), dt))
print("position, beta = 11/3:", sp.collect(sp.expand(dx.subs(beta, sp.Rational(11, 3))), dt))
b = sp.solve(sp.Eq(dx.coeff(dt, 3), 0), beta)
print("beta that removes the jerk term (the -1/12 time-symmetric form):", b,
      "->", sp.collect(sp.expand(dx.subs(beta, b[0])), dt))
