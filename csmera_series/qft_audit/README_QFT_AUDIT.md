# CS-MERA QFT Consistency Audit — Summary

## Status: categorical/modular consistency strong; axiomatic QFT local reconstruction conditional

## Key results

| Test | What it checks | Result |
|------|---------------|--------|
| 1  | Interval entropy, c_eff | 1.794 (0.34% from 9/5) |
| 3  | Modular S,T matrices | exact to 10⁻¹⁵ |
| 6  | Entanglement first law | error 10⁻⁹ |
| 10B | WZW entropy scaling | c=9/5 via Z₃×U(1) decomposition |
| 17 | Pull-through global | error < 10⁻¹⁵ (all 104 channels) |
| 17 | Local twist residual | nonzero in 103/104 channels |
| 18 | Uniform energy bound | local bounds finite; global bound conditional on optimized MERA |

## The precise open problem

The SU(2)₃ fusion intertwiners satisfy M_out(g)A = AM_in(g) to machine
precision. A position-dependent twist generates a nonzero local residual
(the candidate Kac-Moody current). The remaining question is whether the
renormalized sum of these residuals satisfies a uniform energy bound in
an optimized SU(2)₃-symmetric MERA.

## Conditional theorem

If the bound holds, the currents converge strongly to the SU(2)₃ WZW
conformal net (see Documento_Tecnico_QFT_CSMERA_v24.tex for full proof).
