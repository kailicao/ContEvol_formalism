# Prototype case: classic harmonic oscillator

## The paper and the notebooks

Equation, figure, and table numbers are those of the published version.

| Paper | Notebook(s) | What it computes |
| --- | --- | --- |
| §2 and §2.1, Eqs. 5–20 | [2024.04.13 ContEvol CHO first order.nb](2024.04.13%20ContEvol%20CHO%20first%20order.nb) | First-order ContEvol for the oscillator: coefficients, time evolution operator, its determinant, series expansion, minimized cost function |
| §2.2, Eqs. 21–27 | Not in the repository; rederived in [checks/check_oscillator.py](../checks/check_oscillator.py) | Fourth-order Runge–Kutta |
| §2.2, Eqs. 28–39 | [2024.04.14 ContEvol CHO Runge-Kutta.nb](2024.04.14%20ContEvol%20CHO%20Runge-Kutta.nb) | Eighth-order Runge–Kutta stages, result, determinant, and step size (one coefficient is wrong: see Notes since publication) |
| §2.3, Eqs. 40–62 | [2024.04.14 ContEvol CHO second order.nb](2024.04.14%20ContEvol%20CHO%20second%20order.nb) | Second-order ContEvol, with and without the equation of motion enforced at the end of the step |

## The notebooks

Notebooks for Section 2 of the paper, in date order:

- `2024.04.11 ContEvol CHO prototype.nb`: First trial, with one step of length 1 from the initial values 1 and 0, comparing first-order ContEvol with Taylor expansion and fourth-order Runge–Kutta. A first trial, not used in the paper as such.
- `2024.04.13 ContEvol CHO first order.nb`: First-order ContEvol for the oscillator. §2 and §2.1, Eqs. 5–20.
- `2024.04.14 ContEvol CHO Runge-Kutta.nb`: Eighth-order Runge–Kutta, with the step size required by fourth-order Runge–Kutta for comparison. §2.2, Eqs. 28–39. One coefficient of the stage $k_7$ is wrong (see Notes since publication in the main README).
- `2024.04.14 ContEvol CHO second order.nb`: Second-order ContEvol, with and without the equation of motion enforced at the end of the step. §2.3, Eqs. 40–62.
