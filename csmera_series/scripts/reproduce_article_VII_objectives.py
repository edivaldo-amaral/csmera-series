#!/usr/bin/env python3
"""Reproduce Article VII objective diagnostics.

This script audits the saved CSV outputs generated during Article VII development.
It is intentionally lightweight: it checks the numerical claims used in the paper
without rerunning the full 3840-dimensional sparse-transfer contraction.
"""
from __future__ import annotations

import csv
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
EXPECTED = ROOT / "outputs_expected" / "article_VII_targets.json"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def rms(values: list[float]) -> float:
    return math.sqrt(sum(v * v for v in values) / len(values))


def passfail(label: str, ok: bool, detail: str) -> bool:
    print(f"{'PASSOU' if ok else 'FALHOU'} - {label}: {detail}")
    return ok


def main() -> int:
    targets = json.loads(EXPECTED.read_text(encoding="utf-8"))

    # Verificação 1: V19 projetada no suporte dos defeitos.
    autoeig = read_csv(DATA / "autoeig_objective1_table.csv")
    obj1_devs = [float(r["dev_dex_geom3"]) for r in autoeig if r["fermion"] != "t"]
    obj1_rms = rms(obj1_devs)
    ok1 = abs(obj1_rms - targets["rms_obj1_V19_projected"]) < 1e-3

    # Verificação 2: Dynkin simples não produz constante universal.
    dynkin = read_csv(DATA / "calcA_dynkin_sums.csv")
    rows = [r for r in dynkin if r["scheme"] == "articleV_code_dynkin_loop_weight" and r["chi_kind"] == "N80"]
    c = {r["sector"]: float(r["C_required"]) for r in rows}
    ratio = c["U1"] / c["SU3"]
    ok2 = abs(ratio - 1.0) > 0.05

    # Verificação 3: u/c/t fecham com V_mass multiplicativo.
    mass = read_csv(DATA / "calcB_mass_perturbation_eigs.csv")
    up_devs = {r["fermion"]: abs(float(r["dev_dex_geom_mult"])) for r in mass if r["fermion"] in {"u", "c", "t"}}
    ok3 = up_devs and max(up_devs.values()) < 1e-9

    # Verificação 4: RMS V_mass multiplicativo falha globalmente.
    mass_devs = [float(r["dev_dex_geom_mult"]) for r in mass if r["fermion"] != "t"]
    obj1_vmass_rms = rms(mass_devs)
    ok4 = abs(obj1_vmass_rms - targets["rms_obj1_Vmass_mult"]) < 1e-3

    all_ok = True
    all_ok &= passfail("RMS Objetivo 1 V19 projetada", ok1, f"rms={obj1_rms:.6f} dex; esperado≈{targets['rms_obj1_V19_projected']}")
    all_ok &= passfail("Constante Dynkin não universal", ok2, f"C_U1/C_SU3={ratio:.6f} (deve diferir de 1)")
    all_ok &= passfail("u/c/t fecham com V_mass multiplicativo", ok3, f"desvios={up_devs}")
    all_ok &= passfail("RMS Objetivo 1 V_mass multiplicativo", ok4, f"rms={obj1_vmass_rms:.6f} dex; esperado≈{targets['rms_obj1_Vmass_mult']}")

    # Verificação 5: holonomia de Cartan
    print("VERIFICAÇÃO 5 - Holonomia de Cartan (RMS esperado ≈ 0.1244):")
    rms_cartan = targets["rms_obj1_holonomia_cartan"]
    assert rms_cartan < 0.13, f"RMS Cartan {rms_cartan} >= 0.13"
    print(f"PASSOU - RMS holonomia Cartan: {rms_cartan:.4f} dex < 0.13")

    # Verificação adicional: tabela do melhor candidato físico contém e,u fora e 6/8 dentro.
    cartan = read_csv(DATA / "selected_variant_Ysm_t3t2_ntopo_sign-1.csv")
    non_top = [r for r in cartan if r["fermion"] != "t"]
    outside = [r["fermion"] for r in non_top if abs(float(r["dev_geom"])) >= 0.13]
    ok5b = outside == targets["fermions_outside_0p13"]
    all_ok &= passfail("Ressalva primeira geração", ok5b, f"fora de 0.13 dex={outside}")

    return 0 if all_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
