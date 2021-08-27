"""
Module for calculating dispersion of a ferromagnet (including dipolar interaction)
with Holstein-Primakoff theory.

For details see Diss Säubert (https://mediatum.ub.tum.de/?id=1448818) 
Chapter 5 (spin waves in Fe)
and
Appendix B.3.2. (dipolar energy)
"""

### IMPORTS
import numpy as np

# Physical constants
from scipy.constants import mu_0, physical_constants
mu_B = physical_constants["Bohr magneton in eV/T"][0]

def D_theo(D0, T, Tc, mu):
    """
    Calculates the theoretical value for spinwave stiffness D as a function of temperature
    
    Parameters
    ----------
    D0   :   float
        spin-wave stiffness at 0 K meV/angstroem^2
    T    :   float
        temperature of the system in K (needs to be <= Tc)
    Tc   :   float
        Curie temperature of the material in K
    beta :   float, optional
        critical exponent from dynamical scaling theory beta = 1/3 is exact within the theory

    Return
    ------

    Note
    ----
    * D0_Fe = 280 meV / angstroem^2
    * D0_Ni = 430 meV / angstroem^2
    """
    
    return D0*(1-T/Tc)**mu

#------------------------------------------------------------------------------

def mag_func_theo(Msat, T, Tc, beta):
    """
    Calculates the theoretical magnetisation as a function of temperature
    from the saturation magnetisation Msat and the Curie temperature Tc
    
    Parameters
    ----------
    Msat :   float
        saturation magnetisation of the material in A/m
    T    :   float
        temperature of the system in K (needs to be <= Tc)
    Tc   :   float
        Curie temperature of the material in K
    mu   :   float, optional
        critical exponent from dynamical scaling theory mu = 1/3 is exact within the theory
    
    Return
    ------
    M    :   float
        value of the temperature dependent magnetisation
    
    Note
    ----
    * Tc_Fe = 1043. K
    * Msat_Fe = 1.75*(10**6) A/m
    """
    
    return Msat*(1-T/Tc)**beta

#------------------------------------------------------------------------------

def e_dipolar(Msat, T, Tc, beta, g_j, sin2mean):
    """
    Calculates the dipolar energy as function of the temperature in meV
    
    Parameters
    ----------
    T         :   float
        temperature of the system in K (needs to be <= Tc)
    Tc        :   float
        Curie temperature of the material in K
    Msat      :   float
        saturation magnetisation of the material in A/m
    beta :   float, optional
        critical component from dynamical scaling theory beta = 1/3 is exact within the theory
    g_j       :   float, optional
        Landé factor of the material
    sin2mean  :   float, optional
        average of the angle between q and the magnetization => 2/3. for randomly distributed domains
    
    Return
    ------
    e_dipolar :   float
        dipolar exchange energy contribution to the spin-wave dispersion in the ferromagnetic system 
    
    Note
    ----
    * g_j_Fe = 2.
    * Tc_Fe = 1043. K
    * Msat_Fe = 1.75*(10**6) A/m
    """
    
    return g_j * mu_B * mu_0 * mag_func_theo(Msat, T, Tc, beta) * sin2mean * 1000.0 # mu_B is in eV/T --> this is meV

#------------------------------------------------------------------------------

def holstein_primakoff_dispersion(q, D0, Msat, T, Tc, beta, g_j, sin2mean):
    """
    Calculates the theoretical value of the spin wave energy
    as a function of q and T
    
    Parameters
    ----------
    q         :   float, ndarray
        wavevector transfer in angstroem
    D0        :   float
        spin-wave stiffness at 0 K meV/angstroem^2
    Msat      :   float
        saturation magnetisation of the material in A/m
    T         :   float
        absolute temperature of the system in K. Needs to be smaller than mat_props["Tc"]
    Tc        :   float
        Curie temperature of the material in K
    beta      :   float
        critical component from dynamical scaling theory beta = 1/3 is exact within the theory
    g_j       :   float
        Landé factor of the material
    sin2mean  :   float
        average of the angle between q and the magnetization => 2/3. for randomly distributed domains

    Return
    ------
    spin_wave_energy :   float, ndarray
    """

    E_dip = e_dipolar(Msat, T, Tc, beta, g_j, sin2mean)
    D = D_theo(D0, T, Tc, beta)
    return np.sqrt(D*q**2*(D*q**2+E_dip))

#------------------------------------------------------------------------------

def generalized_holstein_primakoff_dispersion(q, D, E_dip):
    """
    Calculates a Holstein Primakoff type dispersion

    Parameters
    ----------
    q     :   float, ndarray
        wavevector transfer in 1/angstroem
    D     :   float
        spin-wave stiffness meV/angstroem^
    E_dip :   float
        dipolar exchange energy contribution to the spin-wave dispersion
        in a ferromagnetic system
    
    Return
    ------
          :   ndarray, float
        Spin wave energy
    """

    return np.sqrt(D * q**2 * (D * q**2 + E_dip))

#------------------------------------------------------------------------------

def generalized_holstein_primakoff_dispersion_resoconvo(q, D, E_dip, reso):
    """
    Calculates a Holstein Primakoff type dispersion
    convolved with a q-resolution

    Parameters
    ----------
    q     :   float, ndarray
        wavevector transfer in 1/angstroem
    D     :   float
        spin-wave stiffness meV/angstroem^
    E_dip :   float
        dipolar exchange energy contribution to the spin-wave dispersion
        in a ferromagnetic system
    reso  :   hp-theory.Resolution
        Resolution object for convolution
    
    Return
    ------
          :   ndarray, float
        Spin wave energy after accounting q-resolution

    """
    pass