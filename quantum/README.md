# Quantum mechanics: stationary Schrödinger equation

## The paper and the notebooks

Equation, figure, and table numbers are those of the published version.

| Paper | Notebook(s) | What it computes |
| --- | --- | --- |
| §4.1, Eqs. 145–152; Figs. 9–12 | [2024.04.27 ContEvol IPW simple.nb](2024.04.27%20ContEvol%20IPW%20simple.nb) | Infinite potential well, simple discretization |
| §4.2, Eqs. 153–160, 164, 166–168 | [2024.04.20 ContEvol IPW first order.nb](2024.04.20%20ContEvol%20IPW%20first%20order.nb) | Infinite potential well, first-order ContEvol: element matrices, normalization, toy version |
| §4.2, Eqs. 161–164; Figs. 13–16 | [2024.04.21 ContEvol IPW first order.nb](2024.04.21%20ContEvol%20IPW%20first%20order.nb) | Assembled matrices, eigenpairs, and figures for the toy version and the realistic versions |
| §4.4, Eqs. 174–196, 204–205 | [2024.05.01 ContEvol H derivation.nb](2024.05.01%20ContEvol%20H%20derivation.nb) | Coulomb potential: interval matrices, cost function, tail, special version |
| §4.4, Figs. 17–19 | [2024.05.02 ContEvol H failure.nb](2024.05.02%20ContEvol%20H%20failure.nb) | Toy version with varying step size |
| §4.4, Eqs. 197–202 | [2024.05.03 ContEvol H correction.nb](2024.05.03%20ContEvol%20H%20correction.nb) | Treatment of the last interval before the tail |
| §4.4, Eq. 203; Figs. 20–21 | [2024.05.03 ContEvol H solution.nb](2024.05.03%20ContEvol%20H%20solution.nb) | Realistic versions: matrices, ground state, and its energy |

## The notebooks

Notebooks for Section 4 of the paper, in date order. The harmonic oscillator (§4.3) is described without notebooks.

- `2024.04.20 ContEvol IPW first order.nb`: Infinite potential well, first-order ContEvol: the element matrices, the normalization, and the toy version with one interval. §4.2, Eqs. 153–160, 164, and 166–168.
- `2024.04.21 ContEvol IPW first order.nb`: The assembled matrices, eigenpairs, and normalized wavefunctions for one, two, four, and eight intervals. §4.2, Eqs. 161–164, Figs. 13–16.
- `2024.04.27 ContEvol IPW simple.nb`: Infinite potential well, simple discretization, from 2 to 64 points. §4.1, Eqs. 145–152, Figs. 9–12.
- `2024.05.01 ContEvol H derivation.nb`: Coulomb potential: the matrices of each interval, the cost function, the exponential tail, and the special version without intervals. §4.4, Eqs. 174–196 and 204–205.
- `2024.05.02 ContEvol H failure.nb`: The toy version with one interval and varying step size, used in the paper: §4.4, Figs. 17–19. Its versions with more intervals are a failed attempt kept for the record: the block of the interval before the tail is written into the global matrices with `=` instead of `+=`, which erases the previous interval's share at their common node (see Notes since publication in the main README). The notebooks of May 3 rederive that interval's matrices instead.
- `2024.05.03 ContEvol H correction.nb`: The corrected treatment of the interval before the tail. §4.4, Eqs. 197–202.
- `2024.05.03 ContEvol H solution.nb`: The realistic versions, from two to eight intervals: matrices, ground state, and its energy. §4.4, Eq. 203, Figs. 20–21.
