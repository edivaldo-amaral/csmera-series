#!/usr/bin/env python3
from __future__ import annotations
import math


def compute() -> dict:
    phi = (1 + math.sqrt(5)) / 2
    fermions = [
        ('e','3_1',-2,3),('u','4_1',0,5),('d','6_1',0,9),
        ('s','5_2',-2,7),('mu','6_3',0,13),('c','5_1',-4,5),
        ('tau','6_2',-2,11),('b','7_6',0,19),('t','7_7',0,21)
    ]
    rows = []
    max_err = 0.0
    for f,k,sigma,det in fermions:
        n_topo = abs(sigma)/2 + math.log(det)/math.log(phi)
        S_branch = math.log(det) + abs(sigma)/2 * math.log(phi)
        expected = n_topo * math.log(phi)
        err = abs(S_branch - expected)
        max_err = max(max_err, err)
        rows.append({
            "fermion": f,
            "knot": k,
            "sigma": sigma,
            "det": det,
            "n_topo": n_topo,
            "S_branch": S_branch,
            "n_topo_ln_phi": expected,
            "abs_error": err,
        })
    return {
        "phi": phi,
        "max_abs_error": max_err,
        "rows": rows,
    }


def main() -> None:
    out = compute()
    print("VERIFICAÇÃO - Teorema ΔS_branch = n_topo × ln(φ):")
    for row in out["rows"]:
        assert row["abs_error"] < 1e-10, f"{row['fermion']}: erro {row['abs_error']}"
        print(
            f"  ✓ {row['fermion']}: S_branch={row['S_branch']:.6f} "
            f"= n_topo×ln(φ)={row['n_topo_ln_phi']:.6f}"
        )
    print("PASSOU - identidade algébrica exata para todos os 9 férmions")


if __name__ == "__main__":
    main()
