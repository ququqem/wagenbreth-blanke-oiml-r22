"""
wagenbreth_blanke.py
====================
Reference implementation of the Wagenbreth–Blanke equation for the density
of ethanol–water mixtures, as published in:

    OIML International Recommendation R 22 (1975)
    "International Alcoholometric Tables"
    https://www.oiml.org/en/files/pdf_r/r022-e75.pdf

Original polynomial developed by:
    Wagenbreth, H. and Blanke, W. (1971)
    PTB-Mitteilungen 81, pp. 412–415

Usage
-----
    from wagenbreth_blanke import density_kg_m3

    rho = density_kg_m3(p=0.40, t=20.0)  # 40% ethanol by mass at 20°C
    print(f"{rho:.4f} kg/m³")

No dependencies beyond the Python standard library.
Released under CC0 1.0 (public domain).
"""

# ---------------------------------------------------------------------------
# A coefficients — concentration terms at 20°C
# A[k] is the coefficient of p^k (0-indexed, so A[0] = A_1 in the standard)
# Units: kg/m³
# ---------------------------------------------------------------------------
_A = [
     998.20123,      # A_1  — pure water density at 20°C; constant term (p^0)
    -192.9769495,    # A_2  — p^1
     389.1238958,    # A_3  — p^2
   -1668.103923,     # A_4  — p^3
   13522.15441,      # A_5  — p^4
  -88292.78388,      # A_6  — p^5
  306287.4042,       # A_7  — p^6
 -613838.1234,       # A_8  — p^7
  747017.2998,       # A_9  — p^8
 -547846.1354,       # A_10 — p^9
  223446.0334,       # A_11 — p^10
  -39032.85426,      # A_12 — p^11
]

# ---------------------------------------------------------------------------
# B coefficients — thermal expansion along the pure-water axis
# B[k] is the coefficient of (t-20)^(k+1) (0-indexed)
# Units: kg/(m³·°C^k)
# ---------------------------------------------------------------------------
_B = [
  -0.20618513,       # B_1 — (t-20)^1
  -0.0052682542,     # B_2 — (t-20)^2
   3.6130013e-05,    # B_3 — (t-20)^3
  -3.8957702e-07,    # B_4 — (t-20)^4
   7.169354e-09,     # B_5 — (t-20)^5
  -9.9739231e-11,    # B_6 — (t-20)^6
]

# ---------------------------------------------------------------------------
# C coefficients — interaction (cross) terms
# C[i][k] is the coefficient of p^(k+1) * (t-20)^(i+1) (both 0-indexed)
# Units: kg/(m³·°C^i)
#
# The matrix is truncated: row i has m_i terms where
#   m = [11, 10, 9, 4, 2]  (for i = 1..5)
# ---------------------------------------------------------------------------
_C = [
    # C1: (t-20)^1 terms, k = 1..11
    [
       0.1693443461530087,   # C_{1,1}  — p^1  * (t-20)^1
     -10.46914743455169,     # C_{1,2}  — p^2  * (t-20)^1
      71.96353469546523,     # C_{1,3}  — p^3  * (t-20)^1
    -704.7478054272792,      # C_{1,4}  — p^4  * (t-20)^1
    3924.090430035045,       # C_{1,5}  — p^5  * (t-20)^1
   -12101.64659068747,       # C_{1,6}  — p^6  * (t-20)^1
    22486.46550400788,       # C_{1,7}  — p^7  * (t-20)^1
   -26055.62982188164,       # C_{1,8}  — p^8  * (t-20)^1
    18523.73922069467,       # C_{1,9}  — p^9  * (t-20)^1
    -7420.201433430137,      # C_{1,10} — p^10 * (t-20)^1
     1285.617841998974,      # C_{1,11} — p^11 * (t-20)^1
    ],
    # C2: (t-20)^2 terms, k = 1..10
    [
      -0.01193013005057010,  # C_{2,1}  — p^1  * (t-20)^2
       0.2517399633803461,   # C_{2,2}  — p^2  * (t-20)^2
      -2.170575700536993,    # C_{2,3}  — p^3  * (t-20)^2
      13.53034988843029,     # C_{2,4}  — p^4  * (t-20)^2
     -50.29988758547014,     # C_{2,5}  — p^5  * (t-20)^2
     109.6355666577570,      # C_{2,6}  — p^6  * (t-20)^2
    -142.2753946421155,      # C_{2,7}  — p^7  * (t-20)^2
     108.0435942856230,      # C_{2,8}  — p^8  * (t-20)^2
     -44.14153236817392,     # C_{2,9}  — p^9  * (t-20)^2
       7.442971530188783,    # C_{2,10} — p^10 * (t-20)^2
    ],
    # C3: (t-20)^3 terms, k = 1..9
    [
      -0.0006802995733503803, # C_{3,1} — p^1 * (t-20)^3
       0.01876837790289664,   # C_{3,2} — p^2 * (t-20)^3
      -0.2002561813734156,    # C_{3,3} — p^3 * (t-20)^3
       1.022992966719220,     # C_{3,4} — p^4 * (t-20)^3
      -2.895696483903638,     # C_{3,5} — p^5 * (t-20)^3
       4.810060584300675,     # C_{3,6} — p^6 * (t-20)^3
      -4.672147440794683,     # C_{3,7} — p^7 * (t-20)^3
       2.458043105903461,     # C_{3,8} — p^8 * (t-20)^3
      -0.5411227621436812,    # C_{3,9} — p^9 * (t-20)^3
    ],
    # C4: (t-20)^4 terms, k = 1..4
    [
       4.075376675622027e-06,  # C_{4,1} — p^1 * (t-20)^4
      -8.763058573471110e-06,  # C_{4,2} — p^2 * (t-20)^4
       6.515031360099368e-06,  # C_{4,3} — p^3 * (t-20)^4
      -1.515784836987210e-06,  # C_{4,4} — p^4 * (t-20)^4
    ],
    # C5: (t-20)^5 terms, k = 1..2
    [
      -2.788074354782409e-08,  # C_{5,1} — p^1 * (t-20)^5
       1.345612883493354e-08,  # C_{5,2} — p^2 * (t-20)^5
    ],
]


def density_kg_m3(p: float, t: float) -> float:
    """
    Compute the density of an ethanol-water mixture.

    Parameters
    ----------
    p : float
        Ethanol mass fraction, dimensionless, in the range [0.0, 1.0].
        Example: p=0.40 means 40% ethanol by mass.
    t : float
        Temperature in degrees Celsius, in the range [-20.0, 40.0].

    Returns
    -------
    float
        Density in kg/m³. Divide by 1000 for g/mL or g/cm³.

    Notes
    -----
    Implements the Wagenbreth–Blanke polynomial from OIML R 22 (1975) using
    Horner's method for efficient evaluation.

    The polynomial is organised as:
        rho = sum_{k=0}^{11} f_k(t) * p^k
    where f_k(t) is a polynomial in (t-20). Horner's method is applied to
    both the inner f_k polynomials and the outer sum over p.
    """
    tau = t - 20.0  # temperature offset from reference

    # Build f_k(tau) for k = 0..11.
    # f_k is the coefficient of p^k in the full polynomial.
    #
    # f_0(tau) = A_1 + B_1*tau + B_2*tau^2 + ... + B_6*tau^6
    # f_k(tau) = A_{k+1} + C_{1,k}*tau + C_{2,k}*tau^2 + ...  for k >= 1
    #            (C_{i,k} is zero when k > m_i, i.e. beyond truncation limit)

    f = [0.0] * 12

    # f_0: A_1 plus all B terms (Horner from highest power down)
    f[0] = _B[5]
    for b in reversed(_B[:5]):
        f[0] = f[0] * tau + b
    f[0] = f[0] * tau + _A[0]

    # f_k for k = 1..11: A_{k+1} plus C terms where they exist
    for k in range(1, 12):
        # Collect C_{i,k} values for i = 5 down to 1 (where k <= m_i)
        # m = [11, 10, 9, 4, 2] for i = 1..5
        c_col = []
        for i in range(5):           # i index: 0=C1 row, 4=C5 row
            row = _C[i]
            if k <= len(row):        # k is 1-based; row is 0-indexed at k-1
                c_col.append(row[k - 1])
            else:
                c_col.append(0.0)
        # c_col[0] = C_{1,k}, c_col[4] = C_{5,k}
        # Horner from highest temperature power (i=5) down to i=1
        val = c_col[4]
        for i in range(3, -1, -1):
            val = val * tau + c_col[i]
        val = val * tau + _A[k]
        f[k] = val

    # Outer Horner: rho = f_0 + f_1*p + f_2*p^2 + ... + f_11*p^11
    rho = f[11]
    for k in range(10, -1, -1):
        rho = rho * p + f[k]

    return rho


# ---------------------------------------------------------------------------
# Validation — runs automatically when the script is executed directly
# ---------------------------------------------------------------------------
_VALIDATION_CASES = [
    # (p,    t,    expected_rho,  description)
    (0.0,  20.0,  998.20123,    "Pure water at 20°C — equals A_1 exactly"),
    (1.0,  20.0,  789.2391233,  "Pure ethanol at 20°C"),
    (0.0,   0.0,  999.8369332,  "Pure water at 0°C — tests B terms only"),
    (1.0,   0.0,  806.2151206,  "Pure ethanol at 0°C — tests all terms"),
    (0.5,  20.0,  913.7705950,  "50% by mass at 20°C — mid-range check"),
]

if __name__ == "__main__":
    print("Wagenbreth–Blanke equation — validation\n")
    print(f"{'p':>5}  {'t':>6}  {'computed':>14}  {'expected':>14}  {'diff':>12}  result")
    print("-" * 72)
    all_pass = True
    for p, t, expected, desc in _VALIDATION_CASES:
        computed = density_kg_m3(p, t)
        diff = computed - expected
        ok = abs(diff) < 1e-4   # tolerance: 0.1 g/m³
        status = "PASS" if ok else "FAIL"
        if not ok:
            all_pass = False
        print(f"{p:>5.1f}  {t:>6.1f}  {computed:>14.7f}  {expected:>14.7f}  {diff:>+12.2e}  {status}")
        print(f"       {desc}")
    print()
    print("All tests passed." if all_pass else "FAILURES detected — check implementation.")
