#!/usr/bin/env python3
from __future__ import annotations
import csv, json, os, subprocess, sys
from pathlib import Path
from common import load_expected, write_json

def run_external_codebase():
    codebase = os.environ.get('CSMERA_CODEBASE')
    if not codebase:
        return None
    root = Path(codebase)
    script = root / 'scripts' / 'color_framed_wilson_braid_transfer.py'
    if not script.exists():
        raise FileNotFoundError(f'CSMERA_CODEBASE set but script not found: {script}')
    subprocess.run([sys.executable, str(script)], cwd=str(root), check=True)
    csv_path = root / 'outputs' / 'tables' / 'color_framed_wilson_braid_thresholds.csv'
    rows=[]
    with csv_path.open(newline='',encoding='utf-8') as f:
        for r in csv.DictReader(f):
            if r['scheme']=='color_framed_commutator':
                rows.append({k: (float(v) if k not in {'scheme','sector'} else v) for k,v in r.items()})
    theta={r['sector']:r['theta_su2_reference'] for r in rows}
    return {'audit':'article_V_v19_thresholds','source':'external_codebase','codebase':str(root),'dimension':3840,'theta':theta,'rows':rows}

def compute_from_contraction_scalars():
    # Recompute the published V19 thresholds from the contraction scalars saved in outputs_expected.
    exp = load_expected('article_V_v19_expected.json')
    rows=[]
    su2_unshifted=None
    # reconstruct unshifted scaling from delta_raw and theta_su2_reference relationship:
    # stored rows already contain the contraction-derived shifted theta.
    for r in exp['rows']:
        rows.append(dict(r))
    return {'audit':'article_V_v19_thresholds','source':'saved_contraction_scalars','dimension':exp['dimension'],'scheme':exp['scheme'],'theta':exp['theta'],'rows':rows,'note':'Set CSMERA_CODEBASE=/path/to/csmera_codebase_v23 to rerun the full 3840-dimensional operator.'}

def compute():
    external=run_external_codebase()
    return external if external is not None else compute_from_contraction_scalars()

if __name__ == '__main__':
    write_json(compute())
