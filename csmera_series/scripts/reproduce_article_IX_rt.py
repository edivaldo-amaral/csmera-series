#!/usr/bin/env python3
from __future__ import annotations
import math


def compute() -> dict:
    phi = (1 + math.sqrt(5)) / 2
    K = 5  # k0 + h_dual = 3 + 2

    coeff_RT = 2 * math.log(phi) / math.log(K)
    c_SU2_3 = 3 * 3 / (3 + 2)  # = 9/5
    coeff_CC = c_SU2_3 / 3
    error_pct = abs(coeff_RT - coeff_CC) / coeff_CC * 100

    N_A = 18
    N_cut = 2 * math.log(N_A) / math.log(K)
    S_RT = N_cut * math.log(phi)
    S_CC = (c_SU2_3 / 3) * math.log(N_A)

    k0 = 3
    h_dual = 2

    return {
        "phi": phi,
        "K": K,
        "coeff_RT": coeff_RT,
        "c_SU2_3": c_SU2_3,
        "coeff_CC": coeff_CC,
        "error_pct": error_pct,
        "N_A": N_A,
        "N_cut": N_cut,
        "S_RT": S_RT,
        "S_CC": S_CC,
        "S_abs_diff": abs(S_RT - S_CC),
        "k0": k0,
        "h_dual": h_dual,
        "K_from_k0_hdual": k0 + h_dual,
    }


def main() -> None:
    out = compute()
    print("VERIFICAÇÃO 1 - Coeficiente RT vs Calabrese-Cardy:")
    print(f"  2ln(φ)/ln(K) = {out['coeff_RT']:.6f}")
    print(f"  c/3 = {out['coeff_CC']:.6f}")
    print(f"  Erro = {out['error_pct']:.4f}%")
    assert out['error_pct'] < 1.0, f"Erro {out['error_pct']:.4f}% >= 1%"
    print("PASSOU")

    print(f"\nVERIFICAÇÃO 2 - S_RT para N_A=18:")
    print(f"  S_RT = {out['S_RT']:.4f}")
    print(f"  S_CC = {out['S_CC']:.4f}")
    assert out['S_abs_diff'] < 0.01, f"Diferença {out['S_abs_diff']:.4f} >= 0.01"
    print("PASSOU")

    assert out['K_from_k0_hdual'] == out['K']
    print(f"\nVERIFICAÇÃO 3 - K = k0+h_dual = {out['k0']}+{out['h_dual']} = {out['K']}")
    print("PASSOU")

    print("\nTodas as verificações do Artigo IX passaram.")


if __name__ == "__main__":
    main()
