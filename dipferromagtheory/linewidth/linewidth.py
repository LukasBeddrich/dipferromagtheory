"""
Implementation of the integral equations determining the longitudinal
and transverse linewdith.
"""
import numpy as np
#------------------------------------------------------------------------------
from .vertices import *
from .susceptibility import Chi

class DynScalingFunc():
    """

    """
    def __init__(self, j, g, xi, qs):
        """
        Parameters
        ----------
        j       : float
            (strength of) exchange interaction
        g       : float
            ration of dipolar vs exchange interaction strength
        xi      : float
            magnetic correlation length (no particular unit decided on)
        qs      : ndarray
            q array for which to calculate the dynamical scaling function
        """
        self.settings = {
            "nrho" : 101,
            "neta" : 21
        }

        self.j = j
        self.g = g
        self.xi = xi
        self.qs = qs
        self.Chi = Chi(j, g, xi)

        self.x = np.power(qs * xi, -1)
        self.y = np.sqrt(g) * np.power(qs, -1)

        self.gammaL = 5.1326 * np.ones(qs.shape) # According to paper gamma_0  = 5.1326 ... "\_(°.°)_/"
        self.gammaT = 5.1326 * np.ones(qs.shape)

#------------------------------------------------------------------------------

    def __call__(self, x, y):
        """

        """
        pass

#------------------------------------------------------------------------------

    def _calcL(self):
        """
        Calculates one iteration of the gammaL given the current state.

        Note
        ----
        put the 'temp' calculation in _calcL_integrand?!?
        """
        newgamma = np.zeros(self.gammaL.shape)
        rho = np.logspace(-1, 3.1, self.settings["nrho"])
        eta = np.linspace(-1,1, self.settings["neta"])
        rr, ee = np.meshgrid(rho, eta)
        rrm = rhom(rr, ee)

        for gidx, (ax, ay) in enumerate(zip(self.x, self.y)):
            temp = vLTTscaled(ay, g, rr, ee) * np.power(rrm, -2) * self.Chi.chiTscaled(ax/rr, ay/rr) * self.Chi.chiTscaled(ax/rrm, ay/rrm)
            temp /= (rr**2.5 * self.gammaT[self._select_rho_idx(ax, rho)] + rhom(rr, ee)**2.5 * self.gammaT(self._select_rhom_idx(ax, rrm)))

#------------------------------------------------------------------------------

    def _calcL_integrand(self, rho, eta):
        """
        Calculates one state of the integration along rho and eta.
        """
        pass

#------------------------------------------------------------------------------

    def _select_rho_idx(self, ax, rhos):
        """

        """
        # Avoiding pyhton list iteration -> generate from 1D
        # 2D versions
        x2d = np.tile(np.atleast_2d(self.x).T, (1, len(rhos)))
        rhos2d = np.tile(rhos, (len(self.x), 1))
        # Now divide of both arrays.
        # It means that every x / rho is calculated and afterwards
        # subtract from all x values again.
        x_over_rhos = ax / rhos2d
        # Take the absolute value of the subtracted array and find
        # the minimum value
        # has length 
        return np.tile(np.argmin(np.abs(x2d - x_over_rhos), axis=1), (self.settings["neta"], 1))

#------------------------------------------------------------------------------

    def _select_rhom_idx(self, ax, rhos):
        """

        """
        pass

#------------------------------------------------------------------------------

    def calc(self, x, y):
        """

        """
        pass

###############################################################################
###############################################################################
###############################################################################

def kappaNi(T, Tc):
    """

    """
    return 0.81 * np.abs((T - Tc)/Tc)**0.701

j=0.1
g = 1.0
xi = 1 / kappaNi(T=629, Tc=628)

q = np.linspace(0.0, 0.03, 101)

chi = Chi(j=j, g=g, xi=xi)