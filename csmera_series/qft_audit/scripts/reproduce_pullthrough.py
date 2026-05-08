#!/usr/bin/env python3
"""Standalone SU(2)_3 pull-through audit using NumPy only.

This script reconstructs canonical SU(2) intertwiners for the multiplicity-free
SU(2)_3 fusion channels and verifies the global pull-through equation

    M_out(g) A = A M_in(g)

both infinitesimally and for finite rotations. It also applies a non-constant
local twist and checks that the residual is nonzero in the four-valent channels.

No CS-MERA codebase, SciPy, SymPy, TenPy, quimb, or ITensor is required.
"""
from __future__ import annotations
import json
import math
from itertools import product
from pathlib import Path
import numpy as np

K_LEVEL = 3
LABELS = tuple(range(K_LEVEL + 1))  # label a = 2j, so j = a/2
TOL = 1.0e-9


def fusion_allowed(a: int, b: int, c: int, k: int = K_LEVEL) -> bool:
    """SU(2)_k fusion admissibility for doubled spins a,b,c."""
    if not (0 <= a <= k and 0 <= b <= k and 0 <= c <= k):
        return False
    if (a + b + c) % 2 != 0:
        return False
    if c < abs(a - b):
        return False
    if c > min(a + b, 2 * k - a - b):
        return False
    return True


def spin_matrices(label: int):
    """Hermitian spin matrices Jx,Jy,Jz for spin j=label/2."""
    j = label / 2.0
    dim = label + 1
    mvals = np.asarray([-j + i for i in range(dim)], dtype=float)
    jp = np.zeros((dim, dim), dtype=complex)
    for col, m in enumerate(mvals):
        mp = m + 1.0
        row = int(round(mp + j))
        if 0 <= row < dim:
            jp[row, col] = math.sqrt(max(0.0, (j - m) * (j + m + 1.0)))
    jm = jp.conj().T
    jx = 0.5 * (jp + jm)
    jy = (jp - jm) / (2.0j)
    jz = np.diag(mvals.astype(complex))
    return jx, jy, jz


def kron_all(mats):
    out = mats[0]
    for m in mats[1:]:
        out = np.kron(out, m)
    return out


def kron_sum(gens):
    """Generator on a tensor product, sum over legs."""
    dims = [g.shape[0] for g in gens]
    total = int(np.prod(dims, dtype=int))
    out = np.zeros((total, total), dtype=complex)
    for i, g in enumerate(gens):
        factors = [g if i == j else np.eye(d, dtype=complex) for j, d in enumerate(dims)]
        out += kron_all(factors)
    return out


def single_generators(a):
    return spin_matrices(a)


def pair_generators(a, b):
    ja = spin_matrices(a)
    jb = spin_matrices(b)
    return tuple(kron_sum([ja[i], jb[i]]) for i in range(3))


def nullspace_numpy(mat, rtol=1.0e-12):
    """Numerical nullspace via NumPy SVD."""
    u, s, vh = np.linalg.svd(mat, full_matrices=True)
    if s.size == 0:
        return vh.conj().T
    cutoff = rtol * max(mat.shape) * float(s[0])
    rank = int(np.sum(s > cutoff))
    return vh[rank:].conj().T


def canonical_three_leg_intertwiner(a, b, c):
    """Isometry C: V_a tensor V_b -> V_c satisfying J_c C = C J_ab."""
    if not fusion_allowed(a, b, c):
        raise ValueError(f"Forbidden SU(2)_3 fusion channel: {a},{b}->{c}")
    in_gens = pair_generators(a, b)
    out_gens = single_generators(c)
    din = (a + 1) * (b + 1)
    dout = c + 1
    equations = []
    for go, gi in zip(out_gens, in_gens):
        # vec(go C - C gi) = (I_in kron go - gi.T kron I_out) vec(C)
        equations.append(np.kron(np.eye(din), go) - np.kron(gi.T, np.eye(dout)))
    mat = np.vstack(equations)
    ns = nullspace_numpy(mat)
    if ns.shape[1] < 1:
        raise RuntimeError(f"No intertwiner nullspace for {a},{b}->{c}")
    C = ns[:, 0].reshape((dout, din), order='F')
    norm = math.sqrt(max(0.0, (np.trace(C @ C.conj().T).real / dout)))
    if norm == 0:
        raise RuntimeError("Degenerate intertwiner norm")
    C = C / norm
    # deterministic phase convention
    idx = np.unravel_index(np.argmax(np.abs(C)), C.shape)
    if abs(C[idx]) > 0:
        C = C / (C[idx] / abs(C[idx]))
    return C


def unitary_from_axis(label_or_pair, axis, theta):
    if isinstance(label_or_pair, tuple):
        gens = pair_generators(label_or_pair[0], label_or_pair[1])
    else:
        gens = single_generators(label_or_pair)
    axis = np.asarray(axis, dtype=float)
    axis = axis / np.linalg.norm(axis)
    H = axis[0] * gens[0] + axis[1] * gens[1] + axis[2] * gens[2]
    # H is Hermitian; diagonalize and exponentiate.
    evals, evecs = np.linalg.eigh(H)
    return evecs @ np.diag(np.exp(1j * theta * evals)) @ evecs.conj().T


def three_leg_audit(a, b, c):
    C = canonical_three_leg_intertwiner(a, b, c)
    in_gens = pair_generators(a, b)
    out_gens = single_generators(c)
    iso = float(np.linalg.norm(C @ C.conj().T - np.eye(c + 1), ord='fro'))
    inf_errors = [float(np.linalg.norm(go @ C - C @ gi, ord='fro')) for go, gi in zip(out_gens, in_gens)]
    Uout = unitary_from_axis(c, (0.31, -0.57, 0.76), 0.731)
    Uin = unitary_from_axis((a, b), (0.31, -0.57, 0.76), 0.731)
    finite = float(np.linalg.norm(Uout @ C - C @ Uin, ord='fro'))
    ja = spin_matrices(a)
    jb = spin_matrices(b)
    Jin_var = kron_sum([0.25 * ja[2], -0.70 * jb[2]])
    Jout_var = 1.10 * out_gens[2]
    variable = float(np.linalg.norm(Jout_var @ C - C @ Jin_var, ord='fro'))
    return {
        'a': a, 'b': b, 'c': c,
        'isometry_error': iso,
        'infinitesimal_error_x': inf_errors[0],
        'infinitesimal_error_y': inf_errors[1],
        'infinitesimal_error_z': inf_errors[2],
        'finite_error': finite,
        'constant_residual_error': max(inf_errors),
        'nonconstant_residual_norm': variable,
        'status': 'PASS' if max([iso, finite, *inf_errors]) < TOL else 'FAIL',
    }


def four_valent_intertwiner(a, b, c, d, x):
    Cin = canonical_three_leg_intertwiner(a, b, x)   # V_ab -> V_x
    Cout = canonical_three_leg_intertwiner(c, d, x)  # V_cd -> V_x
    A = Cout.conj().T @ Cin                          # V_ab -> V_cd
    Pin = Cin.conj().T @ Cin
    Pout = Cout.conj().T @ Cout
    return A, Pin, Pout


def four_valent_audit(a, b, c, d, x):
    A, Pin, Pout = four_valent_intertwiner(a, b, c, d, x)
    in_gens = pair_generators(a, b)
    out_gens = pair_generators(c, d)
    inf_errors = [float(np.linalg.norm(go @ A - A @ gi, ord='fro')) for go, gi in zip(out_gens, in_gens)]
    Uout = unitary_from_axis((c, d), (0.21, 0.44, -0.87), -0.913)
    Uin = unitary_from_axis((a, b), (0.21, 0.44, -0.87), -0.913)
    finite = float(np.linalg.norm(Uout @ A - A @ Uin, ord='fro'))
    pin_err = float(np.linalg.norm(Pin @ Pin - Pin, ord='fro'))
    pout_err = float(np.linalg.norm(Pout @ Pout - Pout, ord='fro'))
    ja = spin_matrices(a)
    jb = spin_matrices(b)
    jc = spin_matrices(c)
    jd = spin_matrices(d)
    Jin_var = kron_sum([0.15 * ja[0], -0.50 * jb[0]])
    Jout_var = kron_sum([0.80 * jc[0], 1.35 * jd[0]])
    variable = float(np.linalg.norm(Jout_var @ A - A @ Jin_var, ord='fro'))
    return {
        'a': a, 'b': b, 'c': c, 'd': d, 'x': x,
        'projector_in_error': pin_err,
        'projector_out_error': pout_err,
        'infinitesimal_error_x': inf_errors[0],
        'infinitesimal_error_y': inf_errors[1],
        'infinitesimal_error_z': inf_errors[2],
        'finite_error': finite,
        'constant_residual_error': max(inf_errors),
        'nonconstant_residual_norm': variable,
        'status': 'PASS' if max([pin_err, pout_err, finite, *inf_errors]) < TOL else 'FAIL',
    }


def run_audit():
    three = []
    for a, b, c in product(LABELS, repeat=3):
        if fusion_allowed(a, b, c):
            three.append(three_leg_audit(a, b, c))

    four = []
    for a, b, c, d in product(LABELS, repeat=4):
        for x in LABELS:
            if fusion_allowed(a, b, x) and fusion_allowed(c, d, x):
                four.append(four_valent_audit(a, b, c, d, x))

    max_three_inf = max(max(r['infinitesimal_error_x'], r['infinitesimal_error_y'], r['infinitesimal_error_z']) for r in three)
    max_three_fin = max(r['finite_error'] for r in three)
    max_three_iso = max(r['isometry_error'] for r in three)
    max_four_inf = max(max(r['infinitesimal_error_x'], r['infinitesimal_error_y'], r['infinitesimal_error_z']) for r in four)
    max_four_fin = max(r['finite_error'] for r in four)
    max_four_proj = max(max(r['projector_in_error'], r['projector_out_error']) for r in four)
    max_const = max(max(r['constant_residual_error'] for r in three), max(r['constant_residual_error'] for r in four))
    max_var = max(max(r['nonconstant_residual_norm'] for r in three), max(r['nonconstant_residual_norm'] for r in four))
    nonzero_four = sum(1 for r in four if r['nonconstant_residual_norm'] > 1.0e-12)
    pass_all = all(r['status'] == 'PASS' for r in three) and all(r['status'] == 'PASS' for r in four)

    return {
        'test': 'Standalone SU(2)_3 pull-through verification',
        'status': 'PASS' if pass_all else 'FAIL',
        'dependencies': 'numpy only',
        'k_level': K_LEVEL,
        'num_three_leg_intertwiners': len(three),
        'num_four_valent_intertwiner_channels': len(four),
        'max_three_leg_isometry_error': max_three_iso,
        'max_three_leg_infinitesimal_pullthrough_error': max_three_inf,
        'max_three_leg_finite_pullthrough_error': max_three_fin,
        'max_four_valent_projector_error': max_four_proj,
        'max_four_valent_infinitesimal_pullthrough_error': max_four_inf,
        'max_four_valent_finite_pullthrough_error': max_four_fin,
        'max_constant_twist_residual': max_const,
        'max_nonconstant_twist_residual': max_var,
        'nonzero_nonconstant_four_valent_residuals': nonzero_four,
        'interpretation': (
            'The canonical SU(2)_3 fusion intertwiners satisfy global pull-through '
            'to numerical precision. Constant twists vanish; position-dependent twists '
            'produce nonzero local residuals, identifying the lattice current candidate. '
            'This standalone audit does not require the CS-MERA codebase and does not prove '
            'the global uniform energy bound for an optimized MERA.'
        ),
    }


def main():
    result = run_audit()
    print(json.dumps(result, indent=2, ensure_ascii=False))
    here = Path(__file__).resolve().parent
    out = here.parent / 'results' / 'reproduced_pullthrough.json'
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2, ensure_ascii=False), encoding='utf-8')
    print(f"Wrote {out}")
    return 0 if result['status'] == 'PASS' else 1


if __name__ == '__main__':
    raise SystemExit(main())
