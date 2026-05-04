#!/usr/bin/env python3
from __future__ import annotations
import importlib.util, json
from pathlib import Path
from common import ROOT, write_json

SCRIPT_NAMES = [
    'reproduce_article_I_masses_LO.py',
    'reproduce_article_I_masses_C1_C7.py',
    'reproduce_article_V_v19_thresholds.py',
    'reproduce_article_V_v30_separation.py',
    'reproduce_article_VI_link_phases.py',
    'reproduce_article_VI_ckm_hcp.py',
]

def load_compute(path: Path):
    spec = importlib.util.spec_from_file_location(path.stem, path)
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    return mod.compute

def compute_all():
    out={}
    for name in SCRIPT_NAMES:
        compute = load_compute(ROOT/'scripts'/name)
        out[name.replace('.py','')] = compute()
    return out

if __name__ == '__main__':
    write_json(compute_all())
