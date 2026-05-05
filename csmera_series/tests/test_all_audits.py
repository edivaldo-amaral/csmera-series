from __future__ import annotations
import importlib.util, math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / 'scripts'
import sys
sys.path.insert(0, str(SCRIPTS))

def load(name: str):
    path = SCRIPTS / name
    spec = importlib.util.spec_from_file_location(path.stem, path)
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    return mod

def approx(a, b, tol):
    assert abs(a-b) <= tol, f'{a} != {b} within {tol}'

def test_article_I_lo_rms_is_not_c1c7():
    mod = load('reproduce_article_I_masses_LO.py')
    out = mod.compute()
    approx(out['rms_global_dex'], 0.7200275861852485, 1e-12)
    assert out['rms_global_dex'] > 0.70

def test_article_I_c1c7_effective():
    mod = load('reproduce_article_I_masses_C1_C7.py')
    out = mod.compute()
    approx(out['effective']['rms_global_dex'], 0.13130631509616475, 1e-12)
    approx(out['effective']['rms_excluding_top_dex'], 0.13927137852963678, 1e-12)

def test_article_V_v19_thresholds():
    mod = load('reproduce_article_V_v19_thresholds.py')
    out = mod.compute()
    approx(out['theta']['SU3'], 0.162322323364, 1e-12)
    approx(out['theta']['SU2'], 0.0, 1e-15)
    approx(out['theta']['U1'], -1.26618486429, 1e-12)

def test_article_V_v30_separation():
    mod = load('reproduce_article_V_v30_separation.py')
    out = mod.compute()
    approx(out['k_norm'], 16.294694713397355, 1e-9)
    approx(out['theta_nlo']['SU3'], 0.024845277053, 1e-12)
    approx(out['theta_nlo']['U1'], 0.05534771811, 1e-12)

def test_article_VI_link_phase_jones_evaluator_no_cache():
    mod = load('reproduce_article_VI_link_phases.py')
    out = mod.compute()
    assert out['cache_used'] is False
    rows = {r['transition']: r for r in out['rows']}
    assert rows['uc']['gamma_recomputed'] == [1,1,1,1,1,2,-1,2,-1]
    assert rows['ds']['gamma_recomputed'] == [1,1,1,2,-1,2,3,-2,3,1,-2,-1,-1]
    approx(rows['uc']['chi_over_pi'], -0.2, 1e-15)
    approx(rows['ds']['chi_over_pi'], 0.2, 1e-15)
    assert rows['uc']['matches_expected_phase']
    assert rows['ds']['matches_expected_phase']
    assert rows['uc']['jones_polynomial_A'] == '(-A**24 + A**20 - 2*A**16 + 2*A**12 - 2*A**8 + A**4 - 1)/A**30'
    assert rows['ds']['jones_polynomial_A'] == '(A**24 + A**16 + A**8 + 1)/A**24'

def test_article_VI_ckm_hcp():
    mod = load('reproduce_article_VI_ckm_hcp.py')
    out = mod.compute()
    approx(out['rho_bar'], 0.1559369711081561, 1e-15)
    approx(out['eta_bar'], 0.3486856676032578, 1e-15)
    approx(out['delta_deg'], 65.94137650808422, 1e-12)
    approx(out['J'], 3.136257074694502e-05, 1e-17)


def test_article_vii_objectives_script_runs():
    import subprocess, sys
    result = subprocess.run(
        [sys.executable, str(ROOT / "scripts" / "reproduce_article_VII_objectives.py")],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert "PASSOU - RMS Objetivo 1 V19 projetada" in result.stdout
    assert "PASSOU - Constante Dynkin não universal" in result.stdout
