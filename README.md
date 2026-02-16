# Golden Ratio Unit-Norm Identity

**Verification code for:**  
*A Unit-Norm Identity Relating the Product of Pisot Boundary Roots to the Golden Ratio*  
Stephanie Alexander, 2026

## The Result

Let ρ ≈ 1.32472 be the real root of x³ = x + 1 (the plastic constant, smallest Pisot number)  
and r₄ ≈ 1.22074 be the real root of x⁴ = x + 1 (the first non-Pisot root in the family).

Their product ρ·r₄ ≈ 1.61714 is suspiciously close to the golden ratio φ ≈ 1.61803.

**This is not a coincidence.** The minimal polynomial of ρ·r₄ over ℚ, reduced modulo the golden ratio polynomial, gives:

```
P(t) mod (t² - t - 1) = t - 1
```

exactly. The remainder t − 1 = φ⁻¹ is the fundamental unit of ℤ[φ] with algebraic norm −1. No other product r_m · r_n from the family x^n = x + 1 has this property — all others have |norm| ≥ 11.

## Quick Verification

```bash
pip install sympy mpmath
python verify_unit_norm.py
```

The core identity can also be checked in three lines:

```python
from sympy import Symbol, Poly, resultant
t, x = Symbol('t'), Symbol('x')
P = resultant(Poly(x**3-x-1, x), Poly(t**4-t*x**3-x**4, x), x)
print(divmod(Poly(P, t), Poly(t**2-t-1, t)))
```

Output: `(Poly(t**10 + t**9 + 2*t**8 + 2*t**4 + t**3, t), Poly(t - 1, t))`

## Files

- `verify_unit_norm.py` — Full verification: symbolic identity, 50-digit numerics, uniqueness survey
- `notebooks/` — Colab notebook ([![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/stalex444/golden-ratio-unit-norm/blob/main/notebooks/verify_unit_norm.ipynb))
- `golden_ratio_unit_norm.tex` — Paper source (LaTeX)
- `golden_ratio_unit_norm.pdf` — Paper (compiled)

## Citation

Paper available on Zenodo: [10.5281/zenodo.18435677](https://doi.org/10.5281/zenodo.18435677)

## License

MIT
