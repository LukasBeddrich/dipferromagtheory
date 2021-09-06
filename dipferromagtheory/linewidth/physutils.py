"""
Some equations and parameters for calculation and analysis.

At the same time, I try to figure out the units of these f***ing quantities J, g
"""

from math import pi, sqrt, atan
from scipy.constants import physical_constants

a1 = dict(
    sc=4 * pi,
    bcc=3**(3/2) * pi,
    fcc=2**(5/2) * pi       # Nickl is fcc type structure
)

gLandeNi = 2.18
aNi = 3.524             # angstroem
DNi = 390.0             # meV angstroem^2
a1Ni = a1["fcc"]
Msat = 4.91e5           # A/m --> 55.1 emu/g
muNi = 0.579            # mu_B / atom
SNi = muNi / gLandeNi   # assuming that muNi = gLandeNi * muB * SNi
gNi = g_func(a1Ni, gLandeNi, SNi, aNi)

def g_func(a1key, gL, J, a):
    """

    """
    muB_meV = physical_constants["Bohr magneton in eV/T"] * 1000.0
    return a1[a1key] * (gL * muB_meV)**2 / (2 * J * a**3)

def D_func(J, S, a):
    """
    Note
    ----
    [S]     : ?
    [a]     : angstroem
    [D]     : meV * angstroem**2
    [J]     : ?
    """
    return 2 * J * S * a**2

def qD_func(g, a):
    """
    Parameter
    ---------
    g   : float
        ratio between exchange and dipolar interaction (supposed g << 1)
    a   : float
        lattice constant of the cubic unit cell

    Return
    ------
    qD  : float
        dipolar wavevector of the system
    """
    return g**0.5 / a

def phi(g, xi):
    """
    Parameter
    ---------
    g   : float
        ratio between exchange and dipolar interaction (supposed g << 1)
    xi  : float
        magnetic correlation length (in m ?)

    Return
    ------
    phi : float
        phi angle for parameterization of the dyn. scaling function in radial coordinates
    """

    return atan(g**(3/2) * xi)

def xi_func_theory(xi0, T, Tc, nu):
    """
    Parameter
    ---------
    xi0 : float
        correlation length 0 (in m)
    T   : float
        lattice constant of the cubic unit cell (in K)
    Tc  : float
        critical temperature of the system (in K)
    nu  : float
        critical exponent of the correlation length (nu = gamme_eff/2)

    Return
    ------
    xi  : float
        temperature dependent correlation length (in m)
    """

    return xi0 * (T/Tc - 1.0)**(-nu)

def xi_func_experimental(T, Tc, xi0=0.81, exp=0.701):
    """
    Parameter
    ---------
    T   : float
        lattice constant of the cubic unit cell (in K)
    Tc  : float
        critical temperature of the system (in K)
    xi0 : float
        empirical correlation length 0 (in m)
    exp : float
        empirical critical exponent of the correlation length

    Return
    ------
    xi  : float
        temperature dependent correlation length (in m)

    Note
    ----
    Values for xi0 and exp taken from
    Böni et al., Longitudinal spin fluctuations in nickel, PRB 43, 1 (1991)
    """
    return xi0 * (T/Tc - 1.0)**exp
