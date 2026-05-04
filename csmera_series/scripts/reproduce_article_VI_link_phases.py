#!/usr/bin/env python3
"""
Independent reproduction of the Article VI link phases chi_uc=-pi/5 and chi_ds=+pi/5.

No phase cache is used.

Method
------
1. Read the braid representatives beta_i for u,c,d,s from data/braid_words_article_VI.csv.
2. Build gamma_ij = beta_j beta_i^{-1}.
3. Evaluate the Jones polynomial of the Markov closure of gamma_ij through the
   Temperley--Lieb/Kauffman-bracket representation.
4. Use the convention
       sigma_i      -> A I + A^{-1} e_i,
       sigma_i^{-1} -> A^{-1} I + A e_i,
       delta        = -A^2 - A^{-2},
       V_L(A)       = (-A^3)^(-w(beta)) <closure(beta)>,
   where the closed-bracket normalization uses delta^(number_of_loops - 1), so
   the unknot has Jones value 1.
5. Evaluate at q = exp(2*pi*i/5) with A = q^{-1/4}=exp(-pi*i/10).
"""
from __future__ import annotations

from collections import defaultdict
import math
import sympy as sp

from common import DATA, read_csv_dict, parse_word, inverse_braid_word, write_json

A = sp.symbols("A")
I = sp.I


def identity_diagram(n: int):
    """Temperley--Lieb identity diagram on n strands."""
    return tuple(sorted(tuple(sorted((i, n + i))) for i in range(n)))


def e_diagram(n: int, i: int):
    """Temperley--Lieb generator e_i, with i zero-based."""
    pairs = [tuple(sorted((i, i + 1))), tuple(sorted((n + i, n + i + 1)))]
    for j in range(n):
        if j not in (i, i + 1):
            pairs.append(tuple(sorted((j, n + j))))
    return tuple(sorted(pairs))


class DSU:
    def __init__(self, n: int):
        self.parent = list(range(n))

    def find(self, x: int) -> int:
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a: int, b: int) -> None:
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[rb] = ra


def compose_diagrams(upper, lower, n: int):
    """Compose two TL diagrams, upper above lower. Return resulting diagram and closed loops."""
    dsu = DSU(3 * n)

    # upper endpoints: top -> 0..n-1; bottom -> n..2n-1
    for a, b in upper:
        aa = a if a < n else n + (a - n)
        bb = b if b < n else n + (b - n)
        dsu.union(aa, bb)

    # lower endpoints: top -> n..2n-1; bottom -> 2n..3n-1
    for a, b in lower:
        aa = n + a if a < n else 2 * n + (a - n)
        bb = n + b if b < n else 2 * n + (b - n)
        dsu.union(aa, bb)

    comps = defaultdict(list)
    for x in range(3 * n):
        comps[dsu.find(x)].append(x)

    boundary = set(range(n)) | set(range(2 * n, 3 * n))
    loops = sum(1 for nodes in comps.values() if not any(x in boundary for x in nodes))

    comp_boundary = defaultdict(list)
    for x in boundary:
        comp_boundary[dsu.find(x)].append(x)

    result_pairs = []
    for nodes in comp_boundary.values():
        if len(nodes) != 2:
            raise RuntimeError(f"Invalid TL composition boundary component: {nodes}")

        def convert(x: int) -> int:
            if x < n:
                return x
            return n + (x - 2 * n)

        result_pairs.append(tuple(sorted((convert(nodes[0]), convert(nodes[1])))))

    return tuple(sorted(result_pairs)), loops


def closure_loop_count(diagram, n: int) -> int:
    """Number of loops after Markov closure of a TL diagram."""
    dsu = DSU(2 * n)
    for a, b in diagram:
        dsu.union(a, b)
    for i in range(n):
        dsu.union(i, n + i)
    return len({dsu.find(x) for x in range(2 * n)})


def add_linear(*terms):
    out = defaultdict(lambda: sp.Integer(0))
    for term in terms:
        for diagram, coeff in term.items():
            out[diagram] += coeff
    return {diagram: sp.simplify(coeff) for diagram, coeff in out.items() if coeff != 0}


def scale_linear(term, scalar):
    return {diagram: sp.simplify(scalar * coeff) for diagram, coeff in term.items()}


def multiply_linear(left, right, n: int, delta):
    out = defaultdict(lambda: sp.Integer(0))
    for d1, c1 in left.items():
        for d2, c2 in right.items():
            d3, loops = compose_diagrams(d1, d2, n)
            out[d3] += c1 * c2 * delta**loops
    return {diagram: sp.simplify(coeff) for diagram, coeff in out.items() if coeff != 0}


def braid_generator_element(n: int, generator: int):
    """Kauffman-bracket representation of sigma_i^{+/-1}."""
    idx = abs(generator) - 1
    ident = {identity_diagram(n): sp.Integer(1)}
    e_i = {e_diagram(n, idx): sp.Integer(1)}
    if generator > 0:
        return add_linear(scale_linear(ident, A), scale_linear(e_i, A**-1))
    return add_linear(scale_linear(ident, A**-1), scale_linear(e_i, A))


def kauffman_bracket_closed_braid(word: list[int], n: int):
    """Kauffman bracket of the Markov closure of a braid word."""
    delta = -A**2 - A**-2
    element = {identity_diagram(n): sp.Integer(1)}
    for g in word:
        element = multiply_linear(element, braid_generator_element(n, g), n, delta)

    bracket = sp.Integer(0)
    for diagram, coeff in element.items():
        loops = closure_loop_count(diagram, n)
        bracket += coeff * delta ** (loops - 1)
    return sp.simplify(sp.expand(bracket))


def jones_polynomial_A(word: list[int], n: int):
    """Jones polynomial in the A-variable for the chosen convention."""
    writhe = sum(1 if g > 0 else -1 for g in word)
    bracket = kauffman_bracket_closed_braid(word, n)
    return sp.simplify(sp.expand((-A**3) ** (-writhe) * bracket))


def exact_A_at_q_5():
    """A = q^{-1/4} = exp(-i*pi/10), q=exp(2*pi*i/5)."""
    return sp.sqrt(10 + 2 * sp.sqrt(5)) / 4 - I * (sp.sqrt(5) - 1) / 4


def phase_over_pi_from_value(z) -> float:
    """Return numeric phase/pi in (-1,1]."""
    return float(sp.N(sp.arg(z) / sp.pi, 40))


def compute():
    beta = {r['label']: parse_word(r['beta_word']) for r in read_csv_dict(DATA / 'braid_words_article_VI.csv')}
    transitions = [
        {'label': 'uc', 'source': 'u', 'target': 'c', 'strands': 3, 'expected_chi_exact': '-pi/5'},
        {'label': 'ds', 'source': 'd', 'target': 's', 'strands': 4, 'expected_chi_exact': 'pi/5'},
    ]
    A0 = exact_A_at_q_5()
    rows = []
    for tr in transitions:
        label, source, target, strands = tr['label'], tr['source'], tr['target'], tr['strands']
        gamma = beta[target] + inverse_braid_word(beta[source])
        V_A = jones_polynomial_A(gamma, strands)
        V_q = sp.simplify(sp.expand(V_A.subs(A, A0)))
        re = sp.simplify(sp.re(V_q))
        im = sp.simplify(sp.im(V_q))
        tan2 = sp.simplify((im / re) ** 2)
        chi_over_pi = phase_over_pi_from_value(V_q)
        expected = -0.2 if tr['expected_chi_exact'] == '-pi/5' else 0.2
        rows.append({
            'transition': label,
            'source': source,
            'target': target,
            'strands': strands,
            'beta_source': beta[source],
            'beta_target': beta[target],
            'gamma_recomputed': gamma,
            'jones_polynomial_A': str(V_A),
            'V_q_re_exact': str(re),
            'V_q_im_exact': str(im),
            'im_over_re_squared_exact': str(tan2),
            'V_q_numeric': str(sp.N(V_q, 18)),
            'chi_over_pi': chi_over_pi,
            'chi_rad': chi_over_pi * math.pi,
            'chi_exact': tr['expected_chi_exact'],
            'matches_expected_phase': abs(chi_over_pi - expected) < 1e-14,
        })
    return {
        'audit': 'article_VI_link_phases',
        'method': 'Temperley-Lieb/Kauffman-bracket Jones evaluator; no phase cache is used.',
        'convention': {
            'sigma_i': 'A I + A^{-1} e_i',
            'sigma_i_inverse': 'A^{-1} I + A e_i',
            'delta': '-A^2 - A^{-2}',
            'jones_normalization': 'V_L(A)=(-A^3)^(-w(beta)) <closure(beta)>',
            'closed_bracket_normalization': 'delta^(number_of_loops - 1); unknot Jones value is 1',
            'q': 'exp(2*pi*i/5)',
            'A_at_q': 'q^{-1/4}=exp(-pi*i/10)',
        },
        'rows': rows,
        'cache_used': False,
    }


if __name__ == '__main__':
    write_json(compute())
