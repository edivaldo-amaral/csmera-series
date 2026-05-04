#!/usr/bin/env python3
from __future__ import annotations
import math
from common import load_expected, write_json

def compute():
    exp = load_expected('article_V_v30_expected.json')
    theta_nlo = {k: exp['theta_v11'][k] - exp['theta_v19'][k] for k in exp['theta_v11']}
    k_norm = exp['chi_UV'] / exp['chi_SU2_local']
    return {
        'audit':'article_V_v30_separation',
        'L': exp['L'], 'chi_UV': exp['chi_UV'], 'chi_SU2_local': exp['chi_SU2_local'],
        'k_norm': k_norm,
        'theta_v11': exp['theta_v11'], 'theta_v19': exp['theta_v19'], 'theta_nlo': theta_nlo,
        'expected_note':'theta_nlo = theta_v11 - theta_v19; k_norm = chi_UV / chi_SU2_local.',
    }

if __name__ == '__main__':
    write_json(compute())
