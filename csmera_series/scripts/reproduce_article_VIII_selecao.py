#!/usr/bin/env python3
from __future__ import annotations

import math

phi = (1 + math.sqrt(5)) / 2


def compute() -> dict:
    rows = []
    zero_points = []
    for k in range(1, 8):
        T_k = k * (k + 1) / 2
        C_k = 2 * k
        F_k = T_k - C_k
        rows.append({"k": k, "T_k": T_k, "C_k": C_k, "F_k": F_k})
        if abs(F_k) < 1e-15:
            zero_points.append(k)

    phi_from_k0 = 2 * math.cos(math.pi / 5)
    N_UV = 3**4 - 1
    D_Fib_sq = 1 + phi**2
    e_DE = 2 * N_UV - 3 * 5 + (3 - 1)
    exp_cosmo = 4 * e_DE

    # Assertions used by the audit harness.
    assert zero_points == [3], f"zero points in k=1..7 should be [3], got {zero_points}"
    assert abs(phi**2 - phi - 1) < 1e-10, "phi não satisfaz equação áurea"
    assert abs(phi_from_k0 - phi) < 1e-10, "phi != 2cos(pi/5)"
    assert N_UV == 80, "N_UV != 80"
    assert abs(D_Fib_sq - ((5 + math.sqrt(5)) / 2)) < 1e-10, "D_Fib^2 mismatch"
    assert e_DE == 147, "e_DE != 147"
    assert exp_cosmo == 588, "expoente cosmológico != 588"

    return {
        "frustration_rows": rows,
        "zero_points": zero_points,
        "phi": phi,
        "phi_from_k0": phi_from_k0,
        "N_UV": N_UV,
        "D_Fib_sq": D_Fib_sq,
        "e_DE": e_DE,
        "exp_cosmo": exp_cosmo,
        "exp_ew": -80,
        "all_checks_passed": True,
    }


def main() -> None:
    out = compute()

    print("VERIFICAÇÃO 1 - Frustração nula:")
    for row in out["frustration_rows"]:
        marker = " ← k0" if abs(row["F_k"]) < 1e-15 else ""
        print(
            f"  k={row['k']}: T_k={row['T_k']:.1f}, "
            f"C_k={row['C_k']:.1f}, F(k)={row['F_k']:.1f}{marker}"
        )
    print("PASSOU - k0=3 é o único ponto de frustração nula")

    print(f"PASSOU - phi = 2cos(pi/5) = {out['phi_from_k0']:.10f}")
    print(f"PASSOU - N_UV = 3^4 - 1 = {out['N_UV']}")
    print(f"PASSOU - D_Fib^2 = 1 + phi^2 = {out['D_Fib_sq']:.10f}")
    print(f"PASSOU - e_DE={out['e_DE']}, expoente cosmológico={out['exp_cosmo']}")
    print("\nTodas as verificações do Artigo VIII passaram.")


if __name__ == "__main__":
    main()
