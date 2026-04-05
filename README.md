# Wagenbreth–Blanke Equation — Ethanol-Water Mixture Density

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19425655.svg)](https://doi.org/10.5281/zenodo.19425655)

**Repository:** https://github.com/ququqem/wagenbreth-blanke-oiml-r22

Machine-readable coefficients and reference documentation for the
Wagenbreth–Blanke polynomial, as published in
**OIML International Recommendation R 22 (1975)**,
*International Alcoholometric Tables*.

This repository exists because the equation and its 54 coefficients are
widely cited in metrology, laboratory, and alcoholometry software, yet no
clean, openly-licensed, machine-readable version was available online.

---

## What is the Wagenbreth–Blanke equation?

The Wagenbreth–Blanke equation computes the **density of an ethanol–water
mixture** as a function of ethanol mass fraction and temperature. It is the
basis of the official OIML international alcoholometric tables, and is used
in legal metrology, pharmaceutical manufacturing, distillery operations, and
laboratory instrument calibration worldwide.

The equation was originally developed by H. Wagenbreth and W. Blanke at the
*Physikalisch-Technische Bundesanstalt* (PTB), Germany, for pure water density
(PTB-Mitteilungen 81, pp. 412–415, 1971). It was subsequently extended to
ethanol–water mixtures and adopted by the International Organisation of Legal
Metrology (OIML) in Recommendation R 22 (1975).

---

## The Equation

```
ρ(p, t) = A₁
         + Σ_{k=2}^{12}  A_k · p^(k−1)
         + Σ_{k=1}^{6}   B_k · (t − 20)^k
         + Σ_{i=1}^{5}   Σ_{k=1}^{mᵢ}  C_{i,k} · p^k · (t − 20)^i
```

**Variables:**

| Symbol | Description | Range |
|--------|-------------|-------|
| `p` | Ethanol mass fraction (dimensionless) | 0.0 – 1.0 |
| `t` | Temperature (°C) | −20 to +40 |
| `ρ` | Mixture density (kg/m³) | — |

Divide by 1000 to convert to g/mL.

**Term groups:**

- **A terms** (12 coefficients): concentration dependence at the reference
  temperature of 20°C. A₁ = 998.20123 kg/m³ is the density of pure water at
  20°C and is the constant term (p⁰).
- **B terms** (6 coefficients): thermal expansion along the pure-water axis
  (p = 0). Together with A₁ these form a standalone pure-water density
  polynomial.
- **C terms** (36 non-zero coefficients, in a truncated 5×12 matrix):
  interaction (cross) terms coupling concentration and temperature. The matrix
  is truncated at row-dependent limits m₁=11, m₂=10, m₃=9, m₄=4, m₅=2.

**Total non-zero coefficients: 54** (A: 12, B: 6, C: 36).

---

## Coefficients

The full coefficient set is in
[`wagenbreth_blanke_oiml_r22.toml`](wagenbreth_blanke_oiml_r22.toml),
a self-documenting TOML file with inline commentary explaining every
coefficient, its units, and its physical interpretation.

### A coefficients (concentration terms at 20°C)

| Index | Coefficient | Value |
|-------|-------------|-------|
| A₁  | p⁰  (constant) |  998.20123 |
| A₂  | p¹             | −192.9769495 |
| A₃  | p²             |  389.1238958 |
| A₄  | p³             | −1668.103923 |
| A₅  | p⁴             |  13522.15441 |
| A₆  | p⁵             | −88292.78388 |
| A₇  | p⁶             |  306287.4042 |
| A₈  | p⁷             | −613838.1234 |
| A₉  | p⁸             |  747017.2998 |
| A₁₀ | p⁹             | −547846.1354 |
| A₁₁ | p¹⁰            |  223446.0334 |
| A₁₂ | p¹¹            | −39032.85426 |

Units: kg/m³

### B coefficients (thermal expansion, pure water axis)

| Index | Coefficient | Value |
|-------|-------------|-------|
| B₁ | (t−20)¹ | −0.20618513 |
| B₂ | (t−20)² | −0.0052682542 |
| B₃ | (t−20)³ |  3.6130013 × 10⁻⁵ |
| B₄ | (t−20)⁴ | −3.8957702 × 10⁻⁷ |
| B₅ | (t−20)⁵ |  7.169354 × 10⁻⁹ |
| B₆ | (t−20)⁶ | −9.9739231 × 10⁻¹¹ |

Units: kg/(m³·°C^k)

### C coefficients (interaction terms)

See the TOML file for the full double-precision C matrix. Summary:

| Group | (t−20) power | p powers | Count |
|-------|-------------|----------|-------|
| C1 | 1 | 1–11 | 11 |
| C2 | 2 | 1–10 | 10 |
| C3 | 3 | 1–9  |  9 |
| C4 | 4 | 1–4  |  4 |
| C5 | 5 | 1–2  |  2 |

---

## Validation

The following spot-check values can be used to verify a correct
implementation:

| p (mass fraction) | t (°C) | ρ (kg/m³) | Description |
|-------------------|--------|-----------|-------------|
| 0.0 | 20 | 998.20123   | Pure water at 20°C; equals A₁ exactly |
| 1.0 | 20 | 789.2391233 | Pure ethanol at 20°C |
| 0.0 |  0 | 999.8369332 | Pure water at 0°C; tests B terms |
| 1.0 |  0 | 806.2151206 | Pure ethanol at 0°C; tests all terms |
| 0.5 | 20 | 913.7705950 | 50% by mass at 20°C; mid-range check |

---

## Reference Implementation

A minimal Python implementation is provided in
[`wagenbreth_blanke.py`](wagenbreth_blanke.py).

```python
from wagenbreth_blanke import density_kg_m3

rho = density_kg_m3(p=0.40, t=20.0)   # 40% ethanol by mass at 20°C
print(f"{rho:.4f} kg/m³")             # → 935.2151 kg/m³
```

---

## Implementation Notes

### Efficient evaluation (Horner's method)

The triple sum can be reorganised as:

```
ρ = Σ_{k=0}^{11}  f_k(t) · p^k
```

where each `f_k(t)` is a polynomial in `(t − 20)`. Applying Horner's method
to both the inner `f_k` polynomials and the outer sum in `p` minimises
multiplications and avoids computing large integer powers explicitly.

### Converting between mass fraction and % v/v

The equation uses **mass fraction** (`p`), not volume fraction (% v/v).

To convert % v/v to mass fraction at temperature `t`:

```
p = (vv/100) · ρ_mixture(p, t) / ρ_pure_ethanol(t)
```

This equation is implicit in `p` — Newton–Raphson iteration is required.

### Units

- Input `p`: dimensionless mass fraction, range 0.0 to 1.0
- Input `t`: degrees Celsius
- Output `ρ`: kg/m³ (divide by 1000 for g/mL or g/cm³)

---

## Sources

- **OIML R 22 (1975)**: *International Alcoholometric Tables*.
  International Organisation of Legal Metrology.
  [https://www.oiml.org/en/files/pdf_r/r022-e75.pdf](https://www.oiml.org/en/files/pdf_r/r022-e75.pdf)

- **Wagenbreth, H. and Blanke, W. (1971)**: "Die Dichte des Wassers im
  Internationalen Einheitensystem und in der Internationalen Praktischen
  Temperaturskala von 1968". *PTB-Mitteilungen* 81, pp. 412–415.
  (Original paper establishing the polynomial form for pure water density,
  providing the B terms and the A₁ baseline.)

---

## License

The coefficients and equation structure are transcribed from OIML R 22 (1975),
an international standard document. The TOML file, README, Python
implementation, and all original content in this repository are released under
the **Creative Commons Zero (CC0 1.0)** public domain dedication — you may
use, copy, modify, and distribute this material for any purpose without
restriction.

[https://creativecommons.org/publicdomain/zero/1.0/](https://creativecommons.org/publicdomain/zero/1.0/)

---

## Contributing

If you find an error in any coefficient, please open an issue with a reference
to the specific value in OIML R 22 (1975). Precision matters — please supply
the full value from the source document, not a rounded approximation.
