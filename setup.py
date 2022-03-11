import numpy as np
import scipy as sp
import scipy.constants as sc
import warnings
import sys

from numpy              import log, exp
from scipy.optimize     import differential_evolution
from scipy.integrate    import cumtrapz, odeint
from scipy.special      import erfc

from output_redirect    import *

np.set_printoptions(precision=4)
warnings.filterwarnings('error','.*Excess.*',)
warnings.filterwarnings('error','.*invalid value encountered in log.*',)

pi      = sc.pi
NA      = sc.N_A
e       = sc.e                  # elementary charge
eps0    = sc.epsilon_0          # vacuum permittivity
Fc      = NA*e                  # Faraday constant
kB      = sc.k
Accuracy= 1e-6

nc      = @nc                   # carbon number in tail

tem     = 86
T       = 273+tem               # K, temperature 
D       = 87.74*exp(-0.0046*(T-273))
kT      = kB*T
RT      = NA*kT

MH2O    = 15.9994+1.0079*2
MCH3    = 12.0107+1.0079*3
MCH2    = 12.0107+1.0079*2
MCH     = 12.0107+1.0079
MOH     = 15.9994+1.0079
MCOO    = 12.0107+15.9994*2
Mna     = 22.9898
Mk      = 39.0983
Mrb     = 85.4678
Mcs     = 132.9055

MSO3    = 80.063
MSO4    = 96.06                 # g/mol

rho_w   = 0.14395/0.0112**(1+(1-T/649.727)**0.05107) # [g/l], 273K-648K

sv_a    = @sv                   # ml/g, specific volume of surfactant critical liquid
sv_w    = 1000/rho_w            # ml/g, specific volume of water

MAMP    = MCH3+(nc-1)*MCH2+MCOO+M@name
VAMP    = sv_a*(10/6.022)*MAMP
VH2O    = sv_w*(10/6.022)*MH2O

nx      = 1001                  # number of points for PB equaiton integration
B0      = 20                    # [A], height for hexagon, or surface radius for lamellar phase

VCH2    = 26.9+0.0146*(T-298)   # [A]^3
VCH3    = 54.6+0.1240*(T-298)   # [A]^3

aso          = @size            # ion hard sphere size including outer shell of metal cations
ai           = @gauss           # ion gauss size
L_segm       = 4.6              # [A] segment length R.Nagarajan 2002 Langmuir 18, 31-38. Equ(20)
n_CH2        = nc-1             # number of CH2 group in tail

M_amph       = MAMP
V_amph       = VAMP
V_tail       = VCH3+n_CH2*VCH2
l_tail_exd   = 1.5+1.265*(n_CH2+1)

v3           = V_amph           # amphiphile molecule volume  v3=v1+v2
v0           = V_tail           # tail    volume
v2           = 4*pi/3*aso**3    # ion     volume
v1           = v3 - v2          # monomer volume v1 = v0+v_head

mu0_CH2      = (5.85*np.log(T)+ 896/T-36.15-0.00560*T)
mu0_CH3      = (3.38*np.log(T)+4064/T-44.13+0.02595*T)
mu0_transfer = mu0_CH3+n_CH2*mu0_CH2

Q3           = (27/8)*v0*L_segm
Q2           = (20/8)*v0*L_segm
Q1           = (10/8)*v0*L_segm
