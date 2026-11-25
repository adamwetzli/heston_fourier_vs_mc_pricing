# Martingale Verification of Option Prices under the Heston Model

*Exact-CIR Monte Carlo vs. semi-analytic Fourier (SINC) pricing of European and digital options*

This project checks numerically that, under the risk-neutral measure of the Heston stochastic volatility model, the underlying and the price processes of six European options are **martingales**. Paths are simulated with exact sampling of the CIR variance process. Each option is repriced along every path with a semi-analytic Fourier method, and the Monte Carlo averages are compared with the price at t = 0.

It was written for the Computational / Quantitative Finance exam (Spring session 2026) at the Department of Statistical Sciences, University of Bologna.

---

## The model

$$
\begin{aligned}
S_t &= S_0\, e^{X_t} \\
dX_t &= r\,dt - \tfrac{\nu_t}{2}\,dt + \sqrt{\nu_t}\,dW_t, & X_0 &= 0 \\
d\nu_t &= \lambda(\bar\nu - \nu_t)\,dt + \eta\sqrt{\nu_t}\,dY_t, & \nu_0 &= \nu_o \\
\mathrm{Cov}(dW_t, dY_t) &= \rho\,dt
\end{aligned}
$$

| Parameter | Symbol | Value |
|---|---|---|
| Mean-reversion speed | λ | 1.03179746 |
| Long-run variance | ν̄ | 0.07819284 |
| Vol-of-vol | η | 1.0 |
| Initial variance | ν₀ | 0.04 |
| Spot–vol correlation | ρ | −0.6450755 |
| Spot / Strike | S₀ = K | 1.2 |
| Maturity | T | 1 year |
| Risk-free rate | r | 0 |

With r = 0 the undiscounted prices must be martingales:

$$
\mathbb{E}[S_T \mid \mathcal{F}_t] = S_t, \qquad
\mathbb{E}\big[\mathcal{O}(S_T, \nu_T, 0) \mid \mathcal{F}_t\big] = \mathcal{O}(S_t, \nu_t, T-t)
$$

## Options

| Name | Payoff |
|---|---|
| Put 𝒫 | (K − S_T)⁺ |
| Call 𝒞 | (S_T − K)⁺ |
| Cash-or-nothing put 𝒫_cn | 1[K > S_T] |
| Asset-or-nothing put 𝒫_an | S_T · 1[K > S_T] |
| Cash-or-nothing call 𝒞_cn | 1[K < S_T] |
| Asset-or-nothing call 𝒞_an | S_T · 1[K < S_T] |

---

## Methodology

1. **Fourier reference prices (SINC series).** The Heston characteristic function of X_T = log(S_T/S₀) is evaluated in closed form through the complex-CIR affine solution. The probabilities P(X_T < w) under the risk-neutral measure and the asset measure come from a SINC (Fourier sine) series, which is truncated adaptively at a tolerance of 10⁻⁸. All six prices follow from these two probabilities. The characteristic function is sanity-checked first: φ(0) = 1, |φ| ≤ 1, and conjugate symmetry holds.
2. **Path simulation.** 100,000 trajectories on a daily grid (252 steps):
   - **ν_t:** exact CIR transition, sampled from a scaled non-central χ² distribution.
   - **S_t:** a log-price update conditional on the variance path, using the integrated variance and the correlation term ρ/η (Δν − λ(ν̄Δt − ΔI)).
3. **Underlying martingale test.** At each month t_j, compare ⟨S(t_j)⟩ ± 2σ_MC with S₀.
4. **Option martingale test.** At each month t_j, reprice every option on every path with the SINC method, using the remaining maturity T − t_j and the path's current variance ν(t_j). Then compare ⟨𝒪(S(t_j), T − t_j)⟩ ± 2σ_MC with 𝒪(S₀, T). This avoids a nested Monte Carlo simulation.
5. **Bonus – pricing curves.** Price all six options against maturity T ∈ (0, 1] at t = 0.
6. **Bonus – QE scheme and variance reduction.** Implement Andersen's Quadratic-Exponential (QE) scheme, using a single uniform draw for both QE branches as in Andersen. Compare it with exact CIR, and apply Glasserman's antithetic variates to both schemes.

---

## Results

### Reference prices at t = 0

| Option | Price |
|---|---|
| Put | 0.082747 |
| Call | 0.082747 |
| 𝒫_cn | 0.392665 |
| 𝒫_an | 0.388450 |
| 𝒞_cn | 0.607335 |
| 𝒞_an | 0.811550 |

Consistency checks: 𝒞 − 𝒫 ≈ −1.1·10⁻¹⁶ (put–call parity at the money), 𝒫_cn + 𝒞_cn = 1, and 𝒫_an + 𝒞_an = S₀.

### Martingale property

- **Underlying:** ⟨S(t_j)⟩ stays within ±2σ_MC of S₀ = 1.2 at all 12 monthly dates. For example, it is 1.2002 ± 0.0015 at T.
- **Options:** all six options stay within ±2σ_MC of their Fourier reference price at every date. For the put, deviations range from +0.2σ to +1.3σ.
- **Correlation:** the empirical ρ is −0.6446 ± 0.0006, against −0.6451 in theory.
- **Discretisation:** a plain Euler scheme for CIR stays clearly biased (E[S_T] ≈ 1.188 even with 2000 steps). The exact scheme has no such bias.

### Variance reduction (put, T = 1)

| Method | Put price | Reduction |
|---|---|---|
| Exact CIR, standard | 0.082870 ± 0.001074 | — |
| Exact CIR, antithetic | 0.082761 ± 0.000860 | 1.249× |
| QE, standard | 0.082967 ± 0.001076 | — |
| QE, antithetic | 0.082944 ± 0.000862 | 1.248× |
| Fourier reference | 0.082747 | |

---

## Repository contents

| File | Description |
|---|---|
| `project_code.ipynb` | Full implementation: characteristic function, SINC pricer, exact-CIR and QE simulators, martingale tests, plots |
| `project_report.pdf` | Write-up with algorithms, tables and figures |
| `task.pdf` | Original exam assignment |

## Running it

Requires Python 3.10+.

```bash
pip install numpy pandas scipy matplotlib seaborn jupyter
jupyter notebook project_code.ipynb
```

With 100,000 paths × 252 steps the path arrays take roughly 1 GB of RAM, and the per-path SINC repricing is the slowest step. Lower `N_traj` in the parameters cell for a quicker run. The notebook also writes `fourier_prices.csv`, `call_parity_check.csv` and `prices_from_maturity.csv` to the working directory.

## References

- S. L. Heston (1993). *A Closed-Form Solution for Options with Stochastic Volatility with Applications to Bond and Currency Options.* Review of Financial Studies 6(2).
- L. Andersen (2008). *Simple and Efficient Simulation of the Heston Stochastic Volatility Model.* Journal of Computational Finance 11(3).
- P. Glasserman (2003). *Monte Carlo Methods in Financial Engineering.* Springer.
- P. Rossi, Computational Finance lecture notes and course code (`FT_opt.py`, `cir_evol.py`), University of Bologna.
