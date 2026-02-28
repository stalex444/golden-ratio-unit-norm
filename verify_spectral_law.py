"""
verify_spectral_law.py

Complete verification of all results in:
"A Fibonacci-Lucas Spectral Law for the Trinomials x^n = x + 1"
by Stephanie Alexander (2026)

Verifies:
  Theorem 2  - Fibonacci-Lucas Spectral Law (n = 3..13)
  Theorem 4  - Unit-Norm Remainder Identity for (3,4)
  Proposition 3 - P_{3,4}(t) coefficients and rational root check
  Theorem 6  - Uniqueness of unit norm for all pairs 2 <= m < n <= 8
  Corollary 8 - Entropy stacking

Requirements: pip install sympy mpmath
"""

from sympy import (Symbol, Poly, resultant, fibonacci, lucas,
                   div, factor_list, ZZ)
from mpmath import mp, mpf, log, fabs

mp.dps = 50  # 50-digit arithmetic

t, x = Symbol('t'), Symbol('x')
golden = Poly(t**2 - t - 1, t)

PASS = 0
FAIL = 0

def check(name, condition):
    global PASS, FAIL
    if condition:
        PASS += 1
        print(f"  ✓ {name}")
    else:
        FAIL += 1
        print(f"  ✗ {name}  *** FAILED ***")

def field_norm(a, b):
    """N_{Q(sqrt5)/Q}(a*phi + b) = -a^2 + a*b + b^2"""
    return -a**2 + a*b + b**2

def get_real_root(n):
    """Largest real root of x^n - x - 1 via mpmath."""
    f = lambda x: x**n - x - 1
    return mp.findroot(f, mpf('1.3'))


# ══════════════════════════════════════════════════════════════
print("=" * 65)
print("THEOREM 4: Unit-Norm Remainder Identity for (3,4)")
print("=" * 65)

P34 = resultant(Poly(x**3 - x - 1, x),
                Poly(t**4 - t*x**3 - x**4, x), x)
P34_poly = Poly(P34, t)

# Check P_{3,4}(t) coefficients match equation (11)
expected_coeffs = [1, 0, 0, -3, -2, 0, 2, -1, -3, -1, 0, 1, -1]
check("P_{3,4}(t) coefficients match eq. (11)",
      P34_poly.all_coeffs() == expected_coeffs)

# Polynomial division
q, r = div(P34_poly, golden)
check("Remainder = t - 1",
      r == Poly(t - 1, t))
check("Quotient = t^10 + t^9 + 2t^8 + 2t^4 + t^3",
      q == Poly(t**10 + t**9 + 2*t**8 + 2*t**4 + t**3, t))

# Norm
check("N_phi(t - 1) = -1",
      field_norm(1, -1) == -1)

# Reconstruct: (t^2-t-1)*q(t) + (t-1) = P_{3,4}(t)
reconstructed = golden * q + Poly(t - 1, t)
check("Reconstruction: (t²-t-1)*q(t) + (t-1) = P_{3,4}(t)",
      reconstructed == P34_poly)


# ══════════════════════════════════════════════════════════════
print()
print("=" * 65)
print("PROPOSITION 3: Rational root check for P_{3,4}")
print("=" * 65)

# P(1) and P(-1) — neither is zero, so no rational roots
P_at_1 = sum(expected_coeffs)
P_at_neg1 = sum(c * ((-1)**(12 - i)) for i, c in enumerate(expected_coeffs))
check("P_{3,4}(1) = -7 (not zero)",
      P_at_1 == -7)
check("P_{3,4}(-1) = 1 (not zero)",
      P_at_neg1 == 1)

# Numeric verification: P_{3,4}(rho * r4) ≈ 0
rho = get_real_root(3)
r4 = get_real_root(4)
phi = get_real_root(2)
prod34 = rho * r4

val = prod34**12 - 3*prod34**9 - 2*prod34**8 + 2*prod34**6 \
      - prod34**5 - 3*prod34**4 - prod34**3 + prod34 - 1
check(f"|P_{{3,4}}(rho*r4)| < 10^-45 (got {float(fabs(val)):.2e})",
      fabs(val) < mpf('1e-45'))


# ══════════════════════════════════════════════════════════════
print()
print("=" * 65)
print("THEOREM 2: Fibonacci-Lucas Spectral Law (n = 3..13)")
print("=" * 65)

for n in range(3, 14):
    P2n = resultant(Poly(x**2 - x - 1, x),
                    Poly(t**n - t*x**(n-1) - x**n, x), x)
    P2n_poly = Poly(P2n, t)
    _, rem = div(P2n_poly, golden)
    coeffs = [int(c) for c in rem.all_coeffs()]
    a_actual = coeffs[0] if len(coeffs) == 2 else 0
    b_actual = coeffs[1] if len(coeffs) == 2 else coeffs[0]

    a_pred = int(-fibonacci(2*n) + (-1)**(n+1))
    b_pred = int(-fibonacci(2*n - 1))
    norm_actual = field_norm(a_actual, b_actual)
    norm_pred = int((-1)**(n+1) * lucas(2*n - 1))

    ok = (a_actual == a_pred and b_actual == b_pred
          and norm_actual == norm_pred)
    check(f"n={n:2d}: a={a_actual:>6d} b={b_actual:>6d}  "
          f"N={norm_actual:>8d} = (-1)^{n+1}*L_{2*n-1}={norm_pred:>8d}",
          ok)


# ══════════════════════════════════════════════════════════════
print()
print("=" * 65)
print("THEOREM 6: Uniqueness — all pairs 2 <= m < n <= 8")
print("=" * 65)

norms = {}
for m in range(2, 9):
    for n in range(m + 1, 9):
        Pmn = resultant(Poly(x**m - x - 1, x),
                        Poly(t**n - t*x**(n-1) - x**n, x), x)
        Pmn_poly = Poly(Pmn, t)
        _, rem = div(Pmn_poly, golden)
        coeffs = [int(c) for c in rem.all_coeffs()]
        if len(coeffs) == 2:
            a, b = coeffs
        else:
            a, b = 0, coeffs[0]
        norm = field_norm(a, b)
        norms[(m, n)] = norm

        rm = get_real_root(m)
        rn = get_real_root(n)
        prod = float(rm * rn)

        tag = " *** UNIT ***" if abs(norm) == 1 else ""
        print(f"  ({m},{n}): r_m*r_n = {prod:.5f}  "
              f"rem = {a}t + ({b})  N = {norm}{tag}")

# Check uniqueness claims
check("(3,4) has norm -1",
      norms[(3, 4)] == -1)
check("(3,4) is the ONLY unit norm",
      sum(1 for v in norms.values() if abs(v) == 1) == 1)
check("Minimum non-unit |norm| = 11 at (2,3)",
      min(abs(v) for v in norms.values() if abs(v) > 1) == 11)
check("|N(P_{2,3})| = 11 = L_5",
      abs(norms[(2, 3)]) == 11)

# Check all non-(3,4) have |norm| >= 11
all_above_11 = all(abs(v) >= 11
                   for k, v in norms.items() if k != (3, 4))
check("All (m,n) != (3,4) have |norm| >= 11",
      all_above_11)


# ══════════════════════════════════════════════════════════════
print()
print("=" * 65)
print("COROLLARY 8: Entropy Stacking")
print("=" * 65)

roots = {n: get_real_root(n) for n in range(2, 9)}

h2 = log(roots[2])
h3 = log(roots[3])
h4 = log(roots[4])
surplus = float(h2 - (h3 + h4))
pct = float(surplus / h2) * 100

check(f"h2 - (h3 + h4) = {surplus:.6f} = {pct:.3f}% of h2 (expect ~0.115%)",
      abs(pct - 0.115) < 0.001)

# Check all other consecutive triples have deficit > 25%
for n in range(3, 7):
    hn = log(roots[n])
    hn1 = log(roots[n + 1])
    hn2 = log(roots[n + 2])
    diff = float(hn - (hn1 + hn2))
    deficit_pct = float(abs(diff / hn)) * 100
    check(f"n={n}: h{n}-(h{n+1}+h{n+2}) = {diff:.4f} "
          f"({deficit_pct:.1f}% of h{n}, expect > 25%)",
          deficit_pct > 25)


# ══════════════════════════════════════════════════════════════
print()
print("=" * 65)
print("NUMERICAL SUMMARY")
print("=" * 65)
print(f"  rho       = {float(rho):.15f}")
print(f"  r4        = {float(r4):.15f}")
print(f"  rho * r4  = {float(prod34):.15f}")
print(f"  phi       = {float(phi):.15f}")
print(f"  phi - rho*r4 = {float(phi - prod34):.15f}")
print(f"  Relative gap = {float((phi - prod34)/phi)*100:.4f}%")
q_val = prod34**10 + prod34**9 + 2*prod34**8 + 2*prod34**4 + prod34**3
print(f"  q(rho*r4) = {float(q_val):.4f}")


# ══════════════════════════════════════════════════════════════
print()
print("=" * 65)
if FAIL == 0:
    print(f"ALL {PASS} CHECKS PASSED")
else:
    print(f"{FAIL} CHECKS FAILED out of {PASS + FAIL}")
print("=" * 65)

if FAIL > 0:
    raise AssertionError(f"{FAIL} checks failed")
