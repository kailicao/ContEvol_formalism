"""Check §2.1-2.2 of Cao (2025), ContEvol (arXiv:2405.05188): the classic harmonic
oscillator's first-order ContEvol step (Eqs. 13-20) and RK4's (Eqs. 21-27).

Run from the repository's root:
    python checks/check_oscillator.py
"""
import sympy as sp
import mpmath as mp

h, t, x0, v0, A, B = sp.symbols('h t x0 v0 A B', real=True)

# ---------------------------------------------------------------- CE1 (§2.1)
x = x0 + v0 * t + B * t**2 + A * t**3
eps = sp.integrate(sp.expand((sp.diff(x, t, 2) + x)**2), (t, 0, h))
sol = sp.solve([sp.diff(eps, A), sp.diff(eps, B)], [A, B], dict=True)[0]
D = 2 * (75600 + 10080*h**2 + 1080*h**4 + 24*h**6 + 5*h**8)
A14 = 7 * (-3600*v0 + 1800*x0*h + 60*v0*h**2 + 120*x0*h**3 + 10*x0*h**5 + 3*v0*h**6) / D
B14 = -15 * (5040*x0 + 1092*x0*h**2 + 168*v0*h**3 + 72*x0*h**4 + 8*v0*h**5 + 5*x0*h**6 + 2*v0*h**7) / D
print('Eq. 15, A and B:', sp.simplify(sol[A] - A14) == 0, sp.simplify(sol[B] - B14) == 0)

xh = sp.simplify(x.subs(sol).subs(t, h))
vh = sp.simplify(sp.diff(x, t).subs(sol).subs(t, h))
G = sp.Matrix([[sp.diff(xh, x0), sp.diff(xh, v0)], [sp.diff(vh, x0), sp.diff(vh, v0)]])
G17 = sp.Matrix([[151200 - 55440*h**2 - 1620*h**4 - 192*h**6 + 5*h**8,
                  151200*h - 5040*h**3 + 60*h**5 - 72*h**7 + h**9],
                 [60*h*(-2520 + 84*h**2 + 6*h**4 + h**6),
                  151200 - 55440*h**2 - 1620*h**4 - 192*h**6 + 13*h**8]]) / D
print('Eq. 17, G:', sp.simplify(G - G17) == sp.zeros(2))
detG = sp.simplify(G.det())
det18 = 1 - 19*h**8 / (302400 + 40320*h**2 + 4320*h**4 + 96*h**6 + 20*h**8)
print('Eq. 18, det G:', sp.simplify(detG - det18) == 0)
mp.mp.dps = 30
f = sp.lambdify(h, 1 - det18 - sp.Rational(1, 2**53), 'mpmath')
print('  1 - det G = 2^-53 at h =', mp.nstr(mp.findroot(f, 0.03), 6), '(paper 0.03396)')
# Eq. 19: series against the exact cos/sin
exact = sp.Matrix([[sp.cos(h), sp.sin(h)], [-sp.sin(h), sp.cos(h)]])
err = (G - exact).applyfunc(lambda e: sp.series(e, h, 0, 10).removeO())
print('Eq. 19, G - exact through h^9:', err)
# Eq. 20: minimized cost
epsmin = sp.simplify(eps.subs(sol))
ser = sp.series(epsmin, h, 0, 9).removeO()
print('Eq. 20, minimized cost through h^8:', sp.expand(ser))
for (a, b, guess, paper) in [(1, 0, 0.0024, '0.002402'), (0, 1, 0.016, '0.01634')]:
    g = sp.lambdify(h, epsmin.subs({x0: a, v0: b}) - sp.Rational(1, 2**53), 'mpmath')
    print(f'  (x0, v0) = ({a}, {b}): cost = 2^-53 at h =', mp.nstr(mp.findroot(g, guess), 6), f'(paper {paper})')

# ---------------------------------------------------------------- RK4 (§2.2)
J = sp.Matrix([[0, 1], [-1, 0]])
I2 = sp.eye(2)
Z = J * h
GRK4 = I2 + Z + Z**2 / 2 + Z**3 / 6 + Z**4 / 24
print('Eq. 25, G_RK4:', GRK4.applyfunc(sp.expand))
detR = sp.expand(GRK4.det())
print('Eq. 26, det G_RK4 =', detR, '(paper 1 - h^6/72 + h^8/576)')
r = mp.findroot(sp.lambdify(h, 1 - detR - sp.Rational(1, 2**53), 'mpmath'), 0.0045)
print('  1 - det G_RK4 = 2^-53 at h =', mp.nstr(r, 6), '(paper 0.004472); CE1/RK4 ratio',
      mp.nstr(mp.findroot(f, 0.03) / r, 5), '(paper 7.594)')
