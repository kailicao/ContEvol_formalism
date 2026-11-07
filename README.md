# ContEvol_formalism

Mathematica and Jupyter notebooks for the paper Cao, K. (2025), "ContEvol Formalism: Numerical Methods Based on Hermite Spline Optimization", Mathematics 13, 3981, DOI: [10.3390/math13243981](https://doi.org/10.3390/math13243981), arXiv: [2405.05188](https://arxiv.org/abs/2405.05188). The files are as the paper accessed them, at commit 4c2d155 (May 8, 2024). The folders follow the paper's Sections 2–4: `prototype/` (the classic harmonic oscillator), `celestial/` (the two-body problem) and `quantum/` (the stationary Schrödinger equation). The notebooks are named by the date they were written and their purpose, and each folder's README lists its files.

## Running the notebooks

The Mathematica notebooks (`.nb`) were created with Wolfram Mathematica 11.0. They hold input and output cells only, with comments in the code. Some notebooks use functions defined in an earlier one, as their comments say.

The Jupyter notebooks (`.ipynb`, in `celestial/`) were run with Python 3.11.9. They import NumPy, SciPy (`scipy.optimize.fsolve`), Numba (`njit`) and Matplotlib.

## The paper and the folders

Equation, figure, and table numbers are those of the published version.

| Folder | Paper | What it holds |
| --- | --- | --- |
| [`prototype/`](prototype/README.md) | §2 (§2.1–2.3) | The classic harmonic oscillator: first- and second-order ContEvol, and the Runge–Kutta comparison |
| [`celestial/`](celestial/README.md) | §3 (§3.1–3.4), and Fig. 22 of §5.4 | The two-body problem: the adequate expansion, its equivalences, the conservation correction, and the eccentric-orbit tests |
| [`quantum/`](quantum/README.md) | §4 (§4.1, §4.2, §4.4) | The stationary Schrödinger equation: the infinite potential well and the Coulomb potential |

Sections 1, 3.5, 4.3, and 5 (except Fig. 22) describe the methods without computations of their own. Each folder's README maps its notebooks to the paper's equations, figures, and table, and describes the notebooks the paper does not use (first trials, a failed attempt, and a mistake kept for the record, preparations).

## Notes since publication

These notes rest on the checks below, which I have not checked in detail. Independent checks of the paper's derivations against these notebooks (November 2026) found one error of substance and a few slips in the text. None changes ContEvol's own results; the scripts that show each point are in `checks/` (below).

- **The eighth-order Runge–Kutta benchmark (§2.2).** One coefficient of the stage $k_7$ (Eq. 35), taken from the public implementation the paper cites, should be $-231$, not $-234$ (confirmed by me). With $-234$ the method is of order 3; with $-231$ it is of order 7 on nonlinear problems (no explicit Runge–Kutta method with ten stages is of order 8: Butcher 1985, BIT 25, 521). So the paper's statements about this benchmark (§2.2, and the abstract's comparison with "two commonly used Runge–Kutta methods") rest on the faulty version; against a corrected one, direct second-order ContEvol still does better in the symplecticity measure, $1 - \det G$ ($h^{12}$ against $h^{10}$). The comparisons with fourth-order Runge–Kutta are unaffected. (`check_symbolic.py`)
- **Slips in the text.** Eq. 36 prints $-217 k_0$ where the notebook has $-127$. Two of Eq. 19's coefficients are misprinted ($-18/7$ for $-12/5$ in $x$'s $v_0 h^7$ term, $-18$ for $-84/5$ in $v$'s $v_0 h^8$ term). Eq. 14 drops an $x_0$ in its $h^4$ term. §2 calls the two Runge–Kutta methods "explicit and multistep"; they are single-step (multistage) methods. The intended reference for geometric solvers (§5) is Hairer, Lubich & Wanner, *Geometric Numerical Integration*, not Krantz & Parks. (`check_oscillator.py`)
- **The interval before the tail (§4.4).** The paper's $P'$ and $Q'$ for this interval are right, but its stated reason for rederiving them is not: removing the row and column of $\dot\psi_N = 0$ from the interval's matrices gives the same $P'$ and $Q'$. The attempt with several intervals in `2024.05.02 ContEvol H failure.nb` failed because that interval's block was written into the global matrices with `=` instead of `+=`, which erased the previous interval's share at their common node; with `+=`, its eigenvalues are those of `2024.05.03 ContEvol H solution.nb`. The one-interval part, behind Figs. 17–19, is right. (`check_hydrogen.py`)
- **Related methods.** The paper's parenthetical note on possible reinvention (§1) was apt: in the cases it tests, first-order ContEvol coincides with methods of other fields. For the infinite well it is the cubic Hermite finite-element method (Eq. 160's element matrices are twice the cubic Hermite mass matrix and minus twice the stiffness matrix). For the two-body problem, the linearized step is the fourth-order Hermite corrector of collisional N-body dynamics (Dittmann 2020, MNRAS 496, 1217, with $\beta = 1$), and the step as coded in the timing notebook is a Taylor-series step through $h^5/120$ in position. In numerical analysis, minimizing the integrated squared residual is the least-squares member of the weighted-residual methods (Finlayson 1972); its Galerkin member on the same cubic Hermite form has $\det G = 1$ for the oscillator. By the largest error over the first orbit, first-order ContEvol leads fourth-order Runge–Kutta by 2.5 times at $h = 1/256$ and 1.5 times at $h = 1/1024$ (both with conservation correction). (`check_symbolic.py`, `check_dittmann.py`, `check_orbit.py`, `check_galerkin.py`)

## Checks

The Python scripts in [`checks/`](checks/README.md) check the paper against the notebooks; they were written with AI assistance (Claude) in November 2026, and I have not checked them in detail.
