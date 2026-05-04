#!/usr/bin/env python3
from __future__ import annotations
import math
from pathlib import Path
from common import DATA, PHI, A1, B_MASS, read_csv_dict, write_json, rms

def compute():
    rows = []
    for r in read_csv_dict(DATA / 'knot_dictionary.csv'):
        sigma = float(r['sigma'])
        det = float(r['det'])
        ntopo = abs(sigma)/2.0 + math.log(det)/math.log(PHI)
        log10_pred = A1 * ntopo + B_MASS
        mass_pred = 10.0 ** log10_pred
        mass_ref = float(r['ref_MeV'])
        resid = math.log10(mass_ref) - log10_pred
        rows.append({
            'fermion': r['fermion'], 'knot': r['knot'], 'sigma': sigma, 'det': det,
            'n_topo_recomputed': ntopo, 'n_topo_article': float(r['n_topo_article']),
            'log10_mass_pred': log10_pred, 'mass_pred_MeV': mass_pred,
            'mass_ref_MeV': mass_ref, 'residual_dex': resid,
        })
    return {
        'audit': 'article_I_masses_LO',
        'phi': PHI, 'A1': A1, 'B': B_MASS,
        'expected_interpretation': 'LO/C1 only. RMS is expected to be about 0.720 dex; 0.131 dex belongs to C1-C7.',
        'rms_global_dex': rms([r['residual_dex'] for r in rows]),
        'rms_excluding_top_dex': rms([r['residual_dex'] for r in rows if r['fermion'] != 't']),
        'rows': rows,
    }

if __name__ == '__main__':
    write_json(compute())
