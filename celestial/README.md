# Celestial mechanics: two-body and three-body problems

## The paper and the notebooks

Equation, figure, and table numbers are those of the published version.

| Paper | Notebook(s) | What it computes |
| --- | --- | --- |
| §3.1, Eqs. 67–84 | [2024.04.16 ContEvol 2BP adequate.nb](2024.04.16%20ContEvol%202BP%20adequate.nb), [2024.04.18 ContEvol 2BP adequate.nb](2024.04.18%20ContEvol%202BP%20adequate.nb) | Two-body first-order ContEvol with "adequate" expansion: cost function, linear system, solution, Jacobian, and its determinant, minimized cost function |
| §3.1, Eqs. 85–89 | [2024.04.18 ContEvol 2BP adequate.nb](2024.04.18%20ContEvol%202BP%20adequate.nb), [2024.04.17 ContEvol 2BP parabola.nb](2024.04.17%20ContEvol%202BP%20parabola.nb) | Test cases: uniform circular motion and parabolic motion (the exact parabolic solution, Eqs. 87–88, in the parabola notebook) |
| §3.2, Eqs. 90–96 | [2024.04.19 ContEvol 2BP linearization.nb](2024.04.19%20ContEvol%202BP%20linearization.nb) | Equivalence with linearization |
| §3.2, Eqs. 97–101 | [2024.04.29 ContEvol 2BP expansion.nb](2024.04.29%20ContEvol%202BP%20expansion.nb) | Equivalence with fifth-order Taylor expansion |
| §3.3, Eqs. 102–117 | [2024.04.26 ContEvol 2BP conservative.nb](2024.04.26%20ContEvol%202BP%20conservative.nb) | Velocity from position and position from velocity through energy and angular momentum, with the signs resolved |
| §3.3, Eqs. 118–119 | [2024.04.27 ContEvol 2BP Jacobian.nb](2024.04.27%20ContEvol%202BP%20Jacobian.nb) | Jacobian with conservation correction, and its determinant |
| §3.3, Eqs. 120–123 | [2024.04.26 ContEvol 2BP conservative.nb](2024.04.26%20ContEvol%202BP%20conservative.nb) | Conservation correction applied to the two test cases |
| §3.3, Eqs. 124–125 | [2024.04.27 ContEvol 2BP conservative.nb](2024.04.27%20ContEvol%202BP%20conservative.nb) | Conservation laws enforced in the formalism |
| §3.4, Eqs. 126–135; Fig. 1 | [2024.04.28 ContEvol eccentric exact.ipynb](2024.04.28%20ContEvol%20eccentric%20exact.ipynb) | Exact solution of the eccentric orbit through Kepler's equation |
| §3.4, Table 1 | [2024.04.30 ContEvol eccentric timeit.ipynb](2024.04.30%20ContEvol%20eccentric%20timeit.ipynb) | Timing of the six integrators (leapfrog, fourth-order Runge–Kutta, and first-order ContEvol, without and with conservation correction) |
| §3.4, Figs. 2–8 | [2024.04.30 ContEvol eccentric visualization.ipynb](2024.04.30%20ContEvol%20eccentric%20visualization.ipynb) | Predicted orbits and their deviations from the exact solution |
| §5.4, Fig. 22 | [2024.04.30 ContEvol eccentric visualization.ipynb](2024.04.30%20ContEvol%20eccentric%20visualization.ipynb) | Mean anomaly versus eccentric anomaly |

## The notebooks

Notebooks for Section 3 of the paper, in date order. Only the two-body problem is computed; the three-body problem (§3.5) is described without notebooks.

- `2024.04.16 ContEvol 2BP adequate.nb`: First-order ContEvol with "adequate" expansion, from the cost function to the elimination of the linear system, with checks of each coefficient. §3.1, Eqs. 67–76. Continued in the notebook of April 18, which refers back to this one for the full derivation and its checks.
- `2024.04.16 ContEvol 2BP failure.nb`: A failed attempt kept for the record. Keeping the products of the coefficients made the equations nonlinear, which led to the "adequate" truncation of §3.1.
- `2024.04.17 ContEvol 2BP parabola.nb`: Exact solution of the parabolic test case as a series. §3.1, Eqs. 87–88.
- `2024.04.18 ContEvol 2BP adequate.nb`: The "adequate" expansion carried through: solution, Jacobian, and its determinant, minimized cost function, and the results for the two test cases (uniform circular motion and parabolic motion). §3.1, Eqs. 67–86 and 89.
- `2024.04.19 ContEvol 2BP linearization.nb`: Equivalence with linearization. §3.2, Eqs. 90–96.
- `2024.04.25 ContEvol 2BP mistake.nb`: A mistake kept for the record. In solving the conservation laws for the position, it assumed a positive velocity component along the x axis, so it took one branch of the solution. Replaced by the notebook of April 26.
- `2024.04.26 ContEvol 2BP conservative.nb`: Velocity from position and position from velocity through mechanical energy and angular momentum, keeping both branches and resolving the signs, then the conservation correction applied to the two test cases. §3.3, Eqs. 102–117 and 120–123.
- `2024.04.27 ContEvol 2BP Jacobian.nb`: Jacobian with conservation correction, and its determinant. §3.3, Eqs. 118–119.
- `2024.04.27 ContEvol 2BP conservative.nb`: Conservation laws enforced in the formalism. §3.3, Eqs. 124–125.
- `2024.04.28 ContEvol 2BP preparation.nb`: Preparation of the numerical tests: one symbolic step of leapfrog and fourth-order Runge–Kutta, a truncated version of the linear system of April 18, and the resulting coefficients written out. Not used in the paper as such.
- `2024.04.28 ContEvol eccentric exact.ipynb`: Exact solution of the eccentric orbit, through Kepler's equation for the position and the conservation correction for the velocity. §3.4, Eqs. 126–135, Fig. 1.
- `2024.04.29 ContEvol 2BP expansion.nb`: Equivalence with fifth-order Taylor expansion. §3.2, Eqs. 97–101. It ends with an attempt to correct the position using the corrected velocity, marked in the notebook as not successful and not used in the paper.
- `2024.04.30 ContEvol eccentric timeit.ipynb`: The six integrators (leapfrog, fourth-order Runge–Kutta, and first-order ContEvol, without and with conservation correction), each checked on uniform circular motion, then timed. §3.4, Table 1.
- `2024.04.30 ContEvol eccentric visualization.ipynb`: The same integrators run on the eccentric orbit, with the figures. §3.4, Figs. 2–8, and §5.4, Fig. 22. A study of the apocenter and pericenter, switched off in the notebook, is not in the paper.
