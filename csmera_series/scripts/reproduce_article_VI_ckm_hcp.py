#!/usr/bin/env python3
from __future__ import annotations
import math
from common import PHI, write_json

def compute():
    rho = PHI**(-2) / math.sqrt(6)
    eta = math.sqrt(5) * PHI**(-2) / math.sqrt(6)
    lam = PHI**(-3) * (1 - PHI**(-6))
    s23 = PHI**(-6) * (1 - PHI**(-3))
    A = s23 / lam**2
    z = complex(rho, eta)
    s13_complex = A * lam**3 * z * math.sqrt(1 - A*A*lam**4) / (math.sqrt(1-lam*lam) * (1 - A*A*lam**4*z))
    s13 = abs(s13_complex)
    delta = math.atan2(s13_complex.imag, s13_complex.real)
    s12=lam; c12=math.sqrt(1-s12*s12); c23=math.sqrt(1-s23*s23); c13=math.sqrt(1-s13*s13)
    J = s12*s23*s13*c12*c23*c13**2*math.sin(delta)
    return {
        'audit':'article_VI_ckm_hcp','phi':PHI,
        'rho_bar':rho,'eta_bar':eta,'lambda_topo':lam,'s23':s23,'A_topo':A,
        's13':s13,'delta_rad':delta,'delta_deg':math.degrees(delta),'J':J,
    }

if __name__ == '__main__':
    write_json(compute())
