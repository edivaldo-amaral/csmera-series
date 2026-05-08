#!/usr/bin/env python3
from __future__ import annotations
import math
from common import DATA, B_MASS, read_csv_dict, write_json, rms

def compute_mode(mode: str):
    if mode not in {'effective','global'}:
        raise ValueError('mode must be effective or global')
    key = 'log10Z_effective_tabulated' if mode == 'effective' else 'log10Z_global_normalized'
    rows = []
    for r in read_csv_dict(DATA / 'knot_dictionary.csv'):
        log10Z = float(r[key])
        log10_mass = B_MASS + log10Z
        mass_pred = 10.0 ** log10_mass
        mass_ref = float(r['ref_MeV'])
        resid = math.log10(mass_ref) - log10_mass
        rows.append({
            'fermion': r['fermion'], 'knot': r['knot'], 'log10Z_C1_C7': log10Z,
            'log10_mass_pred': log10_mass, 'mass_pred_MeV': mass_pred,
            'mass_ref_MeV': mass_ref, 'residual_dex': resid,
        })
    return {
        'mode': mode,
        'B': B_MASS,
        'rms_global_dex': rms([r['residual_dex'] for r in rows]),
        'rms_excluding_top_dex': rms([r['residual_dex'] for r in rows if r['fermion'] != 't']),
        'rows': rows,
    }

def compute():
    return {'audit': 'article_I_masses_C1_C7', 'effective': compute_mode('effective'), 'global': compute_mode('global')}

if __name__ == '__main__':
    write_json(compute())
