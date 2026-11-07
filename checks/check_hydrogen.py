"""ContEvol's hydrogen ground state (paper §4.4): why `2024.05.02 ContEvol H failure.nb` failed.

Independent of the notebooks: every element matrix is built here from the cost function
eps = sum over intervals of int (psi'' + 2 psi/r + phi)^2 dr + the tail's integral,
with psi = r x cubic on each interval (the head has no r^4 term, so psi(0) = 0), the tail
psi_N (r/r_N) exp(1 - r/r_N), P = d2 eps / dphi dphi, Q = d2 eps / dphi dpsi (the notebooks'
convention: no factor 1/2), and the unknowns ordered as the notebooks order them:
psi(h), ..., psi(Nh) at 1..N, then psi'(0), ..., psi'((N-1)h) at N+1..2N; psi'(Nh) = 0.

Checks:
1. The notebooks' Tmat (transcribed from the failure notebook) equals the Hermite map here.
2. For every interval i and step h, the 4x4 element matrix with its last row and column
   dropped equals the element matrix rederived from three parameters (the paper's way,
   `Ptmat`/`Qtmat` of the correction notebook): exact, symbolic in h, for i = 1..7.
3. Assembled with `=` at the last interval (the failure notebook's `makePorQ`) against `+=`
   (the solution notebook's): the lowest-|E| eigenvalue for N = 2, 4, 8, against the
   notebooks' saved outputs.
4. One interval (the failure notebook's first part, Figs. 17-19): E1 and psi'_1 against h,
   the h where E1 = -1 and where psi'_1 = 0, against the paper's 1.0469 and 1.0493.

Run from the repository's root: python3 -I checks/check_hydrogen.py   (about half a minute)
"""
import sympy as sp

r, h, x = sp.symbols('r h x', positive=True)


def L(f):
    return sp.diff(f, r, 2) + 2 * f / r


def element(a, b, nodes, head=False):
    """Element matrices (P, Q) on [a, b] in the nodal unknowns `nodes`
    (head: psi(b), psi'(a), psi'(b); body: psi(a), psi(b), psi'(a), psi'(b))."""
    cs = sp.symbols('c1:5')
    mons = [r, r**2, r**3] if head else [r, r**2, r**3, r**4]
    cs = cs[:len(mons)]
    poly = sum(c * m for c, m in zip(cs, mons))
    if head:
        conds = [poly.subs(r, b), sp.diff(poly, r).subs(r, a), sp.diff(poly, r).subs(r, b)]
    else:
        conds = [poly.subs(r, a), poly.subs(r, b), sp.diff(poly, r).subs(r, a),
                 sp.diff(poly, r).subs(r, b)]
    vals = sp.symbols(f'v1:{len(cs) + 1}')
    sol = sp.solve([cnd - v for cnd, v in zip(conds, vals)], cs, dict=True)[0]
    # basis functions: psi = sum_k vals_k * N_k(r)
    basis = [sp.expand(poly.subs(sol).diff(v)) for v in vals]
    n = len(basis)
    P = sp.zeros(n, n)
    Q = sp.zeros(n, n)
    for j in range(n):
        for k in range(n):
            P[j, k] = 2 * sp.integrate(sp.expand(basis[j] * basis[k]), (r, a, b))
            Q[j, k] = 2 * sp.integrate(sp.expand(basis[j] * L(basis[k])), (r, a, b))
    return sp.simplify(P), sp.simplify(Q), basis


def tail_terms(rN):
    """The tail's contributions to P[N, N] and Q[N, N]: d2/dphiN^2 and d2/dphiN dpsiN."""
    psiN, phiN = sp.symbols('psiN phiN')
    f = (r / rN) * sp.exp(1 - r / rN)
    R = psiN * L(f) + phiN * f
    eps = sp.integrate(sp.expand(R**2), (r, rN, sp.oo))
    return sp.simplify(sp.diff(eps, phiN, 2)), sp.simplify(sp.diff(eps, phiN, psiN))


def tmat_notebook(i, hh):
    """`Tmat[i, h]` as the failure notebook prints it: (psi0, psih, psidot0, psidoth) -> (A, B, C, D)."""
    return sp.Matrix([
        [(-1 + 2 * i) / (hh**4 * i**2), (-3 - 2 * i) / (hh**4 * (1 + i)**2), 1 / (hh**3 * i), 1 / (hh**3 * (1 + i))],
        [(2 - 6 * i**2) / (hh**3 * i**2), 2 * (2 + 6 * i + 3 * i**2) / (hh**3 * (1 + i)**2), -(2 + 3 * i) / (hh**2 * i),
         (-1 - 3 * i) / (hh**2 * (1 + i))],
        [(1 + i) * (-1 - 3 * i + 6 * i**2) / (hh**2 * i**2), -i * (8 + 15 * i + 6 * i**2) / (hh**2 * (1 + i)**2),
         (1 + i) * (1 + 3 * i) / (hh * i), i * (2 + 3 * i) / (hh * (1 + i))],
        [-2 * (-1 + i) * (1 + i)**2 / (hh * i), 2 * i**2 * (2 + i) / (hh * (1 + i)), -(1 + i)**2, -i**2]])


def check_tmat():
    psi0, psih, d0, dh = sp.symbols('psi0 psih d0 dh')
    A, B, C, D = sp.symbols('A B C D')
    ok = True
    for i in map(sp.Integer, range(1, 6)):
        a, b = i * h, (i + 1) * h
        poly = D * r + C * r**2 + B * r**3 + A * r**4
        sol = sp.solve([poly.subs(r, a) - psi0, poly.subs(r, b) - psih,
                        sp.diff(poly, r).subs(r, a) - d0, sp.diff(poly, r).subs(r, b) - dh],
                       [A, B, C, D], dict=True)[0]
        mine = sp.Matrix([[sp.diff(sol[s], v) for v in (psi0, psih, d0, dh)] for s in (A, B, C, D)])
        ok &= sp.simplify(mine - tmat_notebook(i, h)) == sp.zeros(4, 4)
    print(f'1. The notebooks\' Tmat equals the Hermite map (i = 1..5, symbolic h): {ok}')


def check_truncation():
    ok = True
    for i in range(1, 8):
        a, b = i * h, (i + 1) * h
        P4, Q4, _ = element(a, b, None)
        # three parameters: psi'(b) = 0 imposed before the element is built
        cs = sp.symbols('c1:5')
        poly = sum(c * r**k for c, k in zip(cs, (1, 2, 3, 4)))
        v = sp.symbols('v1:4')
        sol = sp.solve([poly.subs(r, a) - v[0], poly.subs(r, b) - v[1],
                        sp.diff(poly, r).subs(r, a) - v[2], sp.diff(poly, r).subs(r, b)], cs, dict=True)[0]
        basis = [sp.expand(poly.subs(sol).diff(vv)) for vv in v]
        P3 = sp.Matrix(3, 3, lambda j, k: 2 * sp.integrate(sp.expand(basis[j] * basis[k]), (r, a, b)))
        Q3 = sp.Matrix(3, 3, lambda j, k: 2 * sp.integrate(sp.expand(basis[j] * L(basis[k])), (r, a, b)))
        keep = [0, 1, 2]
        ok &= sp.simplify(P4.extract(keep, keep) - P3) == sp.zeros(3, 3)
        ok &= sp.simplify(Q4.extract(keep, keep) - Q3) == sp.zeros(3, 3)
    print(f'2. Dropping the last row and column = rederiving from three parameters (P and Q, i = 1..7, symbolic h): {ok}')
    # and against the correction notebook's own figures at i = 1, h = 1/2
    P4, Q4, _ = element(h, 2 * h, None)
    print('   at i = 1, h = 1/2, truncated P:', list(P4.extract([0, 1, 2], [0, 1, 2]).subs(h, sp.Rational(1, 2))))
    print('   the notebook\'s Ptmat[1, 1/2]:  [128/315, 617/5040, 187/5040, 617/5040, 3349/10080, 83/4032, 187/5040, 83/4032, 23/5040]')
    print('   truncated Q:', list(Q4.extract([0, 1, 2], [0, 1, 2]).subs(h, sp.Rational(1, 2))))
    print('   the notebook\'s Qtmat[1, 1/2]:  [-149/42, 877/168, -307/140, 877/168, -14011/3360, 43/112, -27/140, 43/112, -173/840]')


def assemble(N, overwrite):
    hh = sp.Rational(1, N)
    n = 2 * N
    P = sp.zeros(n, n)
    Q = sp.zeros(n, n)
    # head: unknowns psi(h) -> 0, psi'(0) -> N, psi'(h) -> N+1 (0-based)
    Ph, Qh, _ = element(0, h, None, head=True)
    idx = [0, N, N + 1]
    for j in range(3):
        for k in range(3):
            P[idx[j], idx[k]] = Ph[j, k].subs(h, hh)
            Q[idx[j], idx[k]] = Qh[j, k].subs(h, hh)
    for i in range(1, N):
        Pe, Qe, _ = element(i * h, (i + 1) * h, None)
        Pe, Qe = Pe.subs(h, hh), Qe.subs(h, hh)
        if i < N - 1:
            idx = [i - 1, i, N + i, N + i + 1]
            for j in range(4):
                for k in range(4):
                    P[idx[j], idx[k]] += Pe[j, k]
                    Q[idx[j], idx[k]] += Qe[j, k]
        else:  # the interval before the tail: psi'(1) = 0, last row and column dropped
            idx = [i - 1, i, N + i]
            for j in range(3):
                for k in range(3):
                    if overwrite:
                        P[idx[j], idx[k]] = Pe[j, k]
                        Q[idx[j], idx[k]] = Qe[j, k]
                    else:
                        P[idx[j], idx[k]] += Pe[j, k]
                        Q[idx[j], idx[k]] += Qe[j, k]
    tp, tq = tail_terms(1)
    P[N - 1, N - 1] += tp
    Q[N - 1, N - 1] += tq
    return P, Q


def lowest(P, Q):
    import mpmath as mp
    mp.mp.dps = 30
    H = -(P.inv() * Q)
    ev = mp.eig(mp.matrix(H.evalf(40).tolist()))[0]
    ev = sorted(ev, key=lambda z: abs(z))
    return ev[0], ev


def check_assembly():
    saved_fail = {2: '-1.43986146871250', 4: '(not shown in the saved output)', 8: '(not shown)'}
    saved_ok = {2: '-0.999996548604010', 4: '-0.999999959796249', 8: '-0.999999999596798'}
    print('3. The lowest-|E| eigenvalue, assembled with "=" (failure notebook) and "+=" (solution notebook):')
    for N in (2, 4, 8):
        Pf, Qf = assemble(N, True)
        Pc, Qc = assemble(N, False)
        ef, allf = lowest(Pf, Qf)
        ec, _ = lowest(Pc, Qc)
        print(f'   N = {N}: "=" gives {sp.N(complex(ef), 15)} (saved: {saved_fail[N]}); '
              f'"+=" gives {sp.N(complex(ec), 15)} (saved: {saved_ok[N]})')
        if N == 2:
            print('          all four with "=":', [sp.N(complex(z), 8) for z in allf],
                  '(saved: 37.388, -37.242, 5.4290, -1.4399)')


def check_one_interval():
    P, Q, basis = element(0, h, None, head=True)
    tp, tq = tail_terms(h)
    print(f'4. One interval: the tail adds {tp} to P and {sp.simplify(tq)} to Q '
          '(the notebook: 5/2 hh and 3 - 1/(2 hh))')
    P = sp.Matrix(P)
    Q = sp.Matrix(Q)
    P[0, 0] += tp
    Q[0, 0] += tq
    Hm = sp.lambdify(h, -(P.inv() * Q), 'mpmath')
    import mpmath as mp
    mp.mp.dps = 30

    def ground(hv):
        E, ER = mp.eig(mp.matrix(Hm(mp.mpf(hv))))
        k = min(range(len(E)), key=lambda j: abs(E[j]))
        vec = [ER[j, k] / ER[0, k] for j in range(3)]  # psi(h) = 1; psi'_1 then has its sign
        return mp.re(E[k]), [mp.re(v) for v in vec]

    for hv in ('0.5', '0.75', '1', '1.25', '1.5'):
        lam, vec = ground(hv)
        print(f'   h = {hv}: E1 = {mp.nstr(lam, 12)}, psi\'(h)/psi(h) = {mp.nstr(vec[2], 8)}')
    hE = mp.findroot(lambda t: ground(t)[0] + 1, 1.05)
    hd = mp.findroot(lambda t: ground(t)[1][2], 1.05)
    print(f'   E1 = -1 at h = {mp.nstr(hE, 8)} (paper 1.0469); psi\'_1 = 0 at h = {mp.nstr(hd, 8)} (paper 1.0493)')
    lam, vec = ground(1)
    print(f'   at h = 1 the saved output is E1 = -1.03436379096583 and (psi_1, psidot_0, psidot_1) '
          f'proportional to (0.736834, 1.943134, 0.046942): here {mp.nstr(lam, 15)}, '
          f'ratios psidot_0/psi_1 = {mp.nstr(vec[1], 9)} (saved {1.94313384696853 / 0.736833858344606:.9f}), '
          f'psidot_1/psi_1 = {mp.nstr(vec[2], 9)} (saved {0.0469416360789075 / 0.736833858344606:.9f})')


def paper_primed(i, hh):
    """P'(i), Q'(i) as the paper prints them (§4.4, after 'we have to plug the two sets of three parameters')."""
    P = sp.Matrix([
        [(234 * i**4 + 42 * i**3 - 17 * i**2 - 4 * i + 1) / (315 * i**4) * hh,
         (162 * i**5 + 324 * i**4 + 154 * i**3 - 8 * i**2 - 15 * i) / (630 * i**3 * (i + 1)**2) * hh,
         (132 * i**3 + 60 * i**2 - i - 4) / (1260 * i**3) * hh**2],
        [0, (468 * i**7 + 1788 * i**6 + 2522 * i**5 + 1560 * i**4 + 360 * i**3) / (630 * i**3 * (i + 1)**4) * hh,
         (78 * i**4 + 180 * i**3 + 127 * i**2 + 30 * i) / (1260 * i**2 * (i + 1)**2) * hh**2],
        [0, 0, (12 * i**2 + 9 * i + 2) / (630 * i**2) * hh**3]])
    P[1, 0], P[2, 0], P[2, 1] = P[0, 1], P[0, 2], P[1, 2]
    Q = sp.Matrix([
        [(-8 * (63 * i**4 + 4 * i**2 - 4 * i + 1) / hh + (312 * i**3 - 16 * i**2 - 20 * i + 3)) / (210 * i**4),
         (8 * (63 * i**4 + 126 * i**3 + 67 * i**2 + 4 * i - 3) / hh + (108 * i**3 + 162 * i**2 + 20 * i - 17)) / (210 * i**2 * (i + 1)**2),
         (-2 * (231 * i**3 + 14 * i**2 + i - 4) + (44 * i**2 + 6 * i - 3) * hh) / (210 * i**3)],
        [0, (-8 * (63 * i**4 + 252 * i**3 + 382 * i**2 + 264 * i + 72) / hh + (312 * i**3 + 952 * i**2 + 948 * i + 305)) / (210 * (i + 1)**4),
         (2 * (21 * i**3 + 56 * i**2 + 50 * i + 12) + (26 * i**2 + 46 * i + 17) * hh) / (210 * i * (i + 1)**2)],
        [(-2 * (21 * i**3 + 14 * i**2 + i - 4) + (44 * i**2 + 6 * i - 3) * hh) / (210 * i**3), 0,
         (-4 * (14 * i**2 + 7 * i + 2) * hh + (8 * i + 3) * hh**2) / (210 * i**2)]])
    Q[1, 0], Q[2, 1] = Q[0, 1], Q[1, 2]
    return P, Q


def check_paper_and_global():
    ok = True
    for i in map(sp.Integer, range(1, 8)):
        P4, Q4, _ = element(i * h, (i + 1) * h, None)
        Pp, Qp = paper_primed(i, h)
        ok &= sp.simplify(P4.extract([0, 1, 2], [0, 1, 2]) - Pp) == sp.zeros(3, 3)
        ok &= sp.simplify(Q4.extract([0, 1, 2], [0, 1, 2]) - Qp) == sp.zeros(3, 3)
    print(f"5. The paper's printed P'(i) and Q'(i) equal the truncated elements (i = 1..7, symbolic h): {ok}")
    # global: keep psi'(Nh) as an unknown (index 2N), assemble everything with +=, then remove its row and column
    for N in (2, 4):
        hh = sp.Rational(1, N)
        n = 2 * N + 1
        P = sp.zeros(n, n)
        Q = sp.zeros(n, n)
        Ph, Qh, _ = element(0, h, None, head=True)
        idx = [0, N, N + 1]
        for j in range(3):
            for k in range(3):
                P[idx[j], idx[k]] += Ph[j, k].subs(h, hh)
                Q[idx[j], idx[k]] += Qh[j, k].subs(h, hh)
        for i in range(1, N):
            Pe, Qe, _ = element(i * h, (i + 1) * h, None)
            idx = [i - 1, i, N + i, N + i + 1]
            for j in range(4):
                for k in range(4):
                    P[idx[j], idx[k]] += Pe[j, k].subs(h, hh)
                    Q[idx[j], idx[k]] += Qe[j, k].subs(h, hh)
        tp, tq = tail_terms(1)
        P[N - 1, N - 1] += tp
        Q[N - 1, N - 1] += tq
        keep = list(range(2 * N))
        Pc, Qc = assemble(N, False)
        same = (P.extract(keep, keep) - Pc == sp.zeros(2 * N, 2 * N)) and (Q.extract(keep, keep) - Qc == sp.zeros(2 * N, 2 * N))
        print(f"   N = {N}: the full matrices (psi'_N kept, all '+=') with psi'_N's row and column removed "
              f"equal the solution notebook's assembly: {same}")


if __name__ == '__main__':
    check_tmat()
    check_truncation()
    check_assembly()
    check_one_interval()
    check_paper_and_global()
