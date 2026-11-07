# Checks

These scripts were written with AI assistance (Claude) in November 2026, to check the paper against the notebooks; I have not checked them in detail, so treat what they show as leads to verify. Python scripts written for these checks (November 2026), each run from the repository's root; each script's first lines say what it checks and how long it takes. They need SymPy and mpmath; `check_orbit.py` also needs NumPy, SciPy, and Numba, as the notebooks do.

| Script | What it checks |
| --- | --- |
| [check_symbolic.py](check_symbolic.py) | The RK8 benchmark with $-234$ and $-231$, and its order on a nonlinear problem (§2.2); the infinite well as cubic Hermite finite elements (§4.2); the least-squares fit behind the linearized two-body step (§3.2) |
| [check_oscillator.py](check_oscillator.py) | The oscillator's first-order ContEvol and fourth-order Runge–Kutta steps (§2.1–2.2, Eqs. 13–27) |
| [check_dittmann.py](check_dittmann.py) | The linearized two-body step against Dittmann's (2020) Hermite corrector |
| [check_orbit.py](check_orbit.py) | The timing notebook's six integrators against the exact eccentric orbit, and first-order ContEvol against the Taylor series (§3.4); run as `python3 -I checks/check_orbit.py .` |
| [check_galerkin.py](check_galerkin.py) | The Galerkin step on the same cubic form, for the oscillator and the truncated two-body problem, against the least-squares step |
| [check_truncation.py](check_truncation.py) | What the adequate truncation (§3.1) puts into the linear system |
| [check_hydrogen.py](check_hydrogen.py) | The hydrogen element matrices from the cost function, the interval before the tail, and `=` against `+=` (§4.4) |
| [nb_text.py](nb_text.py) | Prints a Mathematica notebook's cells as plain text (it can drop a minus sign before a fraction or a square root, so check what matters against the notebook) |
