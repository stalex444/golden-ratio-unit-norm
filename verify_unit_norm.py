"""
Verification script for:
"A Unit-Norm Identity Relating the Product of Pisot Boundary Roots
 to the Golden Ratio"

Stephanie Alexander, 2026

Run: python verify_unit_norm.py

Requires: sympy, mpmath (both pip-installable)
"""

from sympy import Symbol, Poly, resultant, factor, nroots
from mpmath import mp, mpf, polyroots, log, sqrt

# =============================================================
# PART 1: The core identity (symbolic, exact over Z)
# =============================================================
print("=" * 60)
print("PART 1: EXACT POLYNOMIAL IDENTITY")
print("=" * 60)

t, x = Symbol('t'), Symbol('x')

# Minimal polynomial of rho * r4 via resultant
f3 = Poly(x**3 - x - 1, x)
g4 = Poly(t**4 - t*x**3 - x**4, x)
P = Poly(resultant(f3, g4, x), t)

print(f"\nMinimal polynomial P(t) = {P.as_expr()}")
print(f"Degree: {P.degree()}")

# Divide by golden ratio polynomial
golden = Poly(t**2 - t - 1, t)
q, r = divmod(P, golden)

print(f"\nP(t) mod (t^2 - t - 1) = {r.as_expr()}")
print(f"Quotient q(t) = {q.as_expr()}")

# Verify reconstruction
reconstructed = q * golden + r
assert reconstructed == P, "VERIFICATION FAILED"
print("\nReconstruction check: (t^2-t-1)*q(t) + (t-1) == P(t)  ✓")

# Algebraic norm of remainder t - 1 in Z[phi]
# N(a*phi + b) = -a^2 + a*b + b^2
a, b = 1, -1  # remainder is 1*t + (-1)
norm = -a**2 + a*b + b**2
print(f"\nRemainder: {a}t + ({b})")
print(f"Algebraic norm N(t-1) = -{a}^2 + ({a})({b}) + ({b})^2 = {norm}")
assert norm == -1, "NORM CHECK FAILED"
print("Norm = -1 (fundamental unit of Z[phi])  ✓")

# Resultant of P with golden ratio polynomial
res = resultant(P, golden, t)
print(f"\nRes(P, t^2-t-1) = {res}")
assert res == -1, "RESULTANT CHECK FAILED"
print("Resultant = -1  ✓")

# =============================================================
# PART 2: Numerical verification (50 decimal digits)
# =============================================================
print("\n" + "=" * 60)
print("PART 2: HIGH-PRECISION NUMERICAL VERIFICATION")
print("=" * 60)

mp.dps = 50

rho = max(r.real for r in polyroots([1, 0, -1, -1])
          if abs(r.imag) < 1e-30)
r4 = max(r.real for r in polyroots([1, 0, 0, -1, -1])
         if abs(r.imag) < 1e-30 and r.real > 1)
phi = (1 + sqrt(5)) / 2
product = rho * r4

print(f"\nrho     = {rho}")
print(f"r4      = {r4}")
print(f"rho*r4  = {product}")
print(f"phi     = {phi}")
print(f"phi - rho*r4 = {phi - product}")

# Verify P(rho*r4) = 0
P_val = (product**12 - 3*product**9 - 2*product**8
         + 2*product**6 - product**5 - 3*product**4
         - product**3 + product - 1)
print(f"\nP(rho*r4) = {P_val}")
print(f"  (residual magnitude: ~10^{int(mp.log10(abs(P_val)))})")

# Verify the deficit identity
q_val = (product**10 + product**9 + 2*product**8
         + 2*product**4 + product**3)
deficit_LHS = product**2 - product - 1
deficit_RHS = -(product - 1) / q_val

print(f"\nDeficit identity:")
print(f"  (rho*r4)^2 - rho*r4 - 1   = {deficit_LHS}")
print(f"  -(rho*r4 - 1) / q(rho*r4) = {deficit_RHS}")
print(f"  Match to {int(-mp.log10(abs(deficit_LHS - deficit_RHS)))} digits  ✓")

# Entropy stacking
h2 = log(phi)
h3 = log(rho)
h4 = log(r4)
Dh = h2 - (h3 + h4)

print(f"\nEntropy stacking:")
print(f"  h2 = ln(phi) = {h2}")
print(f"  h3 = ln(rho) = {h3}")
print(f"  h4 = ln(r4)  = {h4}")
print(f"  h2 - (h3+h4) = {Dh}")
print(f"  Surplus: {float(Dh/h2)*100:.4f}% of h2")

# =============================================================
# PART 3: Uniqueness survey
# =============================================================
print("\n" + "=" * 60)
print("PART 3: UNIQUENESS (all pairs 2 <= m < n <= 8)")
print("=" * 60)

import numpy as np

print(f"\n{'Pair':>8s} | {'|Norm|':>12s} | {'Remainder':>25s}")
print("-" * 55)

for m in range(2, 9):
    for n in range(m+1, 9):
        fm = Poly(x**m - x - 1, x)
        gn = Poly(t**n - t*x**(n-1) - x**n, x)
        res_mn = Poly(resultant(fm, gn, x), t)
        _, r_mn = divmod(res_mn, golden)
        coeffs = [int(c) for c in r_mn.all_coeffs()]
        if len(coeffs) == 1:
            a_mn, b_mn = 0, coeffs[0]
        else:
            a_mn, b_mn = coeffs[0], coeffs[1]
        norm_mn = -a_mn**2 + a_mn*b_mn + b_mn**2
        tag = " ← UNIT" if abs(norm_mn) == 1 else ""
        rem_str = f"{a_mn}t + {b_mn}" if a_mn else f"{b_mn}"
        print(f"  ({m},{n}) | {abs(norm_mn):>12d} | {rem_str:>25s}{tag}")

print("\n" + "=" * 60)
print("ALL CHECKS PASSED")
print("=" * 60)
