# A Fibonacci–Lucas Spectral Law for the Trinomials x^n = x + 1

**Stephanie Alexander, 2026**

The polynomial family $x^n - x - 1$ governs the boundary between Pisot and non-Pisot algebraic integers. We prove that the golden ratio $\varphi = r_2$ acts as a **spectral filter** on this family, with an exact Fibonacci–Lucas structure that distinguishes the Pisot boundary pair $(r_3, r_4)$ from all others.

## Main Results

### Theorem 2 — The Fibonacci–Lucas Spectral Law

For all $n \geq 3$, the minimal polynomial of $\varphi \cdot r_n$, reduced modulo the golden ratio polynomial $t^2 - t - 1$, has remainder coefficients given exactly by **Fibonacci numbers**, and the field norm equals $(-1)^{n+1} L_{2n-1}$, where $L_k$ is the $k$-th **Lucas number**.

| n | Remainder $a\,t + b$ | Norm $N_\varphi$ |
|---|----------------------|-------------------|
| 3 | $-7t - 5$           | $11 = L_5$        |
| 4 | $-22t - 13$         | $-29 = -L_7$      |
| 5 | $-54t - 34$         | $76 = L_9$        |
| 6 | $-145t - 89$        | $-199 = -L_{11}$  |
| 7 | $-376t - 233$       | $521 = L_{13}$    |
| 8 | $-988t - 610$       | $-1364 = -L_{15}$ |

The remainder coefficients $a = -F_{2n} + (-1)^{n+1}$ and $b = -F_{2n-1}$ are exact for all $n$.

### Theorem 4 — The Unit-Norm Identity for (3, 4)

The minimal polynomial of $\rho \cdot r_4$ satisfies the exact division:

$$P_{3,4}(t) = (t^2 - t - 1)\,q(t) + (t - 1)$$

The remainder $t - 1 = \varphi^{-1}$ is the **fundamental unit** of $\mathbb{Z}[\varphi]$, with algebraic norm $-1$.

### Theorem 6 — Uniqueness of the Unit Norm

For **all** products $r_m \cdot r_n$ with $2 \leq m < n \leq 8$, the (3, 4) pair is the **unique** product with unit norm. Every other pair has $|N_\varphi| \geq 11$, and the minimum is achieved at $(2,3)$ with norm $11 = L_5$ — the first term of the spectral law.

| $(m, n)$ | $r_m \cdot r_n$ | $N_\varphi$ |
|-----------|-----------------|-------------|
| **(3, 4)** | **1.61714...** | **−1** |
| (2, 3) | 2.14344... | 11 |
| (2, 4) | 1.97521... | −29 |
| (3, 5) | 1.54635... | −479 |
| (4, 5) | 1.42498... | −7,129 |
| (5, 6) | 1.32457... | −2,002,189 |

### Corollary 8 — Entropy Stacking

The topological entropies satisfy $\ln\varphi = \ln\rho + \ln r_4 + O(10^{-4})$, a surplus of just 0.115%. No other consecutive triple comes within 25% of additivity.

## Quick Verification

The core identity (Theorem 4) in three lines:

```python
from sympy import Symbol, Poly, resultant
t, x = Symbol('t'), Symbol('x')
P = resultant(Poly(x**3-x-1, x), Poly(t**4-t*x**3-x**4, x), x)
print(divmod(Poly(P, t), Poly(t**2-t-1, t)))
```

Output: `(Poly(t**10 + t**9 + 2*t**8 + 2*t**4 + t**3, t), Poly(t - 1, t))`

## Full Verification

```bash
pip install sympy mpmath
python verify_spectral_law.py
```

This verifies all four theorems:
- **Theorem 2**: Spectral law for $n = 3, \ldots, 13$
- **Theorem 4**: Unit-norm remainder identity
- **Proposition 3**: $P_{3,4}(t)$ coefficients and irreducibility checks
- **Theorem 6**: Uniqueness survey for all 21 pairs with $2 \leq m < n \leq 8$
- **Corollary 8**: Entropy stacking values

## Files

- `verify_spectral_law.py` — Full verification of all theorems (symbolic + numeric)
- `verify_unit_norm.py` — Original verification of the (3,4) identity
- `notebooks/` — Colab notebook
- `golden_ratio_unit_norm.tex` — Paper source (LaTeX)
- `golden_ratio_unit_norm.pdf` — Paper (compiled)

## Citation

```bibtex
@article{alexander2026spectral,
  title={A Fibonacci--Lucas Spectral Law for the Trinomials $x^n = x + 1$},
  author={Alexander, Stephanie},
  year={2026},
  doi={10.5281/zenodo.18435677}
}
```

Paper available on Zenodo: [10.5281/zenodo.18435677](https://doi.org/10.5281/zenodo.18435677)

## License

MIT# golden-ratio-unit-norm
Companion code and paper: A Fibonacci–Lucas Spectral Law for the Trinomials x^n = x + 1
