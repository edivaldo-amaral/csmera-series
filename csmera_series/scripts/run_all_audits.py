#!/usr/bin/env python3
from __future__ import annotations
import importlib.util
import subprocess
import sys
from pathlib import Path
from common import ROOT, write_json

SCRIPT_NAMES = [
    'reproduce_article_I_masses_LO.py',
    'reproduce_article_I_masses_C1_C7.py',
    'reproduce_article_V_v19_thresholds.py',
    'reproduce_article_V_v30_separation.py',
    'reproduce_article_VI_link_phases.py',
    'reproduce_article_VI_ckm_hcp.py',
    'reproduce_article_VII_objectives.py',
    'reproduce_article_VIII_selecao.py',
    'reproduce_article_IX_rt.py',
    'reproduce_article_X_branch.py',
]

def load_module(path: Path):
    spec = importlib.util.spec_from_file_location(path.stem, path)
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    return mod

def run_script(path: Path) -> dict:
    result = subprocess.run(
        [sys.executable, str(path)],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode != 0:
        raise RuntimeError(result.stdout + result.stderr)
    return {
        'script_passed': True,
        'stdout': result.stdout,
    }

def compute_all():
    out={}
    for name in SCRIPT_NAMES:
        path = ROOT/'scripts'/name
        mod = load_module(path)
        if hasattr(mod, 'compute'):
            out[name.replace('.py','')] = mod.compute()
        else:
            out[name.replace('.py','')] = run_script(path)
    return out

if __name__ == '__main__':
    write_json(compute_all())
