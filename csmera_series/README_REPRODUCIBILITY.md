# CS-MERA Series — Reproducibility Harness

This repository contains a lightweight audit harness for reproducing the central numerical claims of the CS-MERA article series from explicit CSV data, formulas, and saved contraction scalars.

## Status summary

| Audit | Script | Status |
|---|---|---|
| Article I LO fermion masses | `scripts/reproduce_article_I_masses_LO.py` | Reproduces LO RMS ≈ 0.720 dex. This deliberately does **not** reproduce 0.131 dex; that value belongs to C1-C7. |
| Article I C1-C7 masses | `scripts/reproduce_article_I_masses_C1_C7.py` | Reproduces RMS ≈ 0.131306 dex for the effective C6 convention. |
| Article V V19 thresholds | `scripts/reproduce_article_V_v19_thresholds.py` | Reproduces the published V19 contraction table from saved scalars. If `CSMERA_CODEBASE` is set, it reruns the full 3840-dimensional codebase script. |
| Article V V30 separation | `scripts/reproduce_article_V_v30_separation.py` | Reproduces `Theta_NLO = Theta_V11 - Theta_V19` and `k_norm = chi_UV/chi_SU2_local`. |
| Article VI link phases | `scripts/reproduce_article_VI_link_phases.py` | Independently recomputes `gamma_ij = beta_j beta_i^{-1}` and evaluates the Jones polynomial of the Markov closure via a Temperley--Lieb/Kauffman-bracket evaluator. No phase cache is used. |
| Article VI CKM via H_CP | `scripts/reproduce_article_VI_ckm_hcp.py` | Reproduces `rho_bar`, `eta_bar`, `delta_CP`, and `J` from phi-only formulas. |

## Quick start

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python scripts/run_all_audits.py
python -m pytest -q
```

Expected test result:

```text
6 passed
```

## Optional full V19 rerun

The lightweight repository does not include the full historical `csmera_codebase_v23` package by default. To rerun the 3840-dimensional V19 construction directly, point the harness to an unpacked copy of the original codebase:

```bash
export CSMERA_CODEBASE=/path/to/csmera_codebase_v23
python scripts/reproduce_article_V_v19_thresholds.py
```

Without `CSMERA_CODEBASE`, the script uses the saved V19 contraction scalars in `outputs_expected/article_V_v19_expected.json` and recomputes the threshold table from those scalars.

## Important interpretation notes

1. **Article I LO vs C1-C7.** The formula

   ```text
   log10(m/MeV) = A1*n_topo + B
   ```

   is the LO/C1 layer only. It gives RMS ≈ 0.720 dex. The RMS ≈ 0.131 dex quoted for the final spectrum requires the complete C1-C7 table.

2. **Article VI link phases.** This repository now verifies the braid-word composition

   ```text
   gamma_ij = beta_j beta_i^{-1}
   ```

   and derives the phases from an independent Temperley--Lieb/Kauffman-bracket evaluator of the Jones polynomial for the Markov closure, evaluated at `q=exp(2*pi*i/5)` with `A=q^(-1/4)`. The script reports `cache_used=false` and reproduces `chi_uc=-pi/5`, `chi_ds=+pi/5`.

3. **Article V V30.** The theorem of separation is reproduced numerically as

   ```text
   Theta_NLO = Theta_V11 - Theta_V19
   ```

   with

   ```text
   Theta_NLO = (+0.0248452770, 0, +0.0553477181)
   ```

## Repository structure

```text
csmera_series/
  README_REPRODUCIBILITY.md
  requirements.txt
  environment.yml
  scripts/
    reproduce_article_I_masses_LO.py
    reproduce_article_I_masses_C1_C7.py
    reproduce_article_V_v19_thresholds.py
    reproduce_article_V_v30_separation.py
    reproduce_article_VI_link_phases.py
    reproduce_article_VI_ckm_hcp.py
    run_all_audits.py
  data/
    knot_dictionary.csv
    braid_words_article_VI.csv
    jones_link_phase_cache.csv
  outputs_expected/
    article_I_lo_expected.json
    article_I_c1c7_expected.json
    article_V_v19_expected.json
    article_V_v30_expected.json
    article_VI_phases_expected.json
    article_VI_ckm_expected.json
  tests/
    test_all_audits.py
```

## Artigo VII — Objetivos de unificação dinâmica

O script `reproduce_article_VII_objectives.py` verifica os cinco resultados principais:
1. RMS V19 projetada: 0.5937 dex
2. Constante Dynkin não universal
3. u/c/t fecham com V_mass multiplicativo
4. RMS V_mass multiplicativo: 0.6250 dex
5. RMS holonomia de Cartan: 0.1244 dex (passa critério 0.13)

Para reproduzir o resultado do Teorema principal do Artigo VII:
a fase θ_i = θ_boundary(β,x) + 2π(Y_SM + t3 + t2) com coeficiente
ntopo(Kf) × ln10/N_UV produz fechamento dinâmico global com
ressalva na primeira geração (e, u).

## Artigo VIII — Princípio de Seleção da Rede

O script `reproduce_article_VIII_selecao.py` verifica as cinco
identidades centrais do Artigo VIII:
1. k0=3 é o único ponto de frustração nula (T_k = C_k)
2. φ = 2cos(π/5) emerge de k=3
3. N_UV = k0^4 − 1 = 80
4. D_Fib^2 = 1 + φ^2
5. Cadeia cosmológica: e_DE=147, expoente=588

Nota de consistência algébrica: como `D_Fib^2 = 1 + φ^2`, então
`D_Fib^2 = (5 + sqrt(5))/2 ≈ 3.6180339887`, não `2 + sqrt(5)`.

## Artigo IX — Geometria Emergente e Calabrese-Cardy

O script `reproduce_article_IX_rt.py` verifica o resultado central:
a rede CS-MERA de vácuo com fator de ramificação K=5 reproduz o
coeficiente de Calabrese-Cardy com erro de 0,34%.
