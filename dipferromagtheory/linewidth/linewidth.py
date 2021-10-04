"""
Implementation of the integral equations determining the longitudinal
and transverse linewdith.
"""
try:
    import autograd.numpy as np
except ImportError:
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
            "rho_upper_limit" : 3.0,
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

    def update(self):
        """

        """
        self.Chi = Chi(self.j, self.g, self.xi)
        self.x = np.power(self.qs * self.xi, -1)
        self.y = np.sqrt(self.g) * np.power(self.qs, -1)

#------------------------------------------------------------------------------

    def set(self, key, value):
        """

        """
        if key.lower() == "j":
            self.j = value
        elif key.lower() == "g":
            self.g = value
        elif key.lower() == "xi":
            self.xi = value
        elif key.lower() == "nrho":
            self.settings["nrho"] = value
        elif key.lower() == "neta":
            self.settings["neta"] = value
        elif key.lower() == "rho_upper_limit":
            self.settings["rho_upper_limit"] = value
        else:
            print(f"key {key} is not known! Acceptable keys are 'j', 'g', 'xi', 'nrho', 'neta', 'rho_upper_limit'")
        
        self.update()

#------------------------------------------------------------------------------

    def calc(self):
        """

        """
        self.gammaL = self._calcL()
        self.gammaT = self._calcT()

        return self.gammaL, self.gammaT

#------------------------------------------------------------------------------

    def _calcL(self):
        """
        Calculates one iteration of the gammaL given the current state.

        Note
        ----
        put the 'temp' calculation in _calcL_integrand?!?
        """
        newgamma = np.zeros(self.gammaL.shape)
        rho = np.logspace(-1.1, self.settings["rho_upper_limit"], self.settings["nrho"])
        eta = np.linspace(-1,1, self.settings["neta"])
        rr, ee = np.meshgrid(rho, eta)
        rrm = rhom(rr, ee)

        for gidx, (ax, ay) in enumerate(zip(self.x, self.y)):
            # first term alpha = L, beta = sigma = T
            LTT = vLTTscaled(ay, self.g, rr, ee) * np.power(rrm, -2) * self.Chi.chiTscaled(ax/rr, ay/rr) * self.Chi.chiTscaled(ax/rrm, ay/rrm)
            LTT /= (rr**2.5 * self.gammaT[self._select_rho_idx(ax, rho)] + np.power(rrm, 2.5) * self.gammaT[self._select_rhom_idx(ax, rrm)])

            # second term alpha = beta = L, sigma = T
            LLT = vLLTscaled(ay, self.g, rr, ee) * np.power(rrm, -2) * self.Chi.chiLscaled(ax/rr, ay/rr) * self.Chi.chiTscaled(ax/rrm, ay/rrm)
            LLT /= (rr**2.5 * self.gammaL[self._select_rho_idx(ax, rho)] + np.power(rrm, 2.5) * self.gammaT[self._select_rhom_idx(ax, rrm)])

            newgamma[gidx] = np.trapz(np.trapz(LTT + LLT, rho, axis=1), eta)
        return newgamma

#------------------------------------------------------------------------------

    def _calcT(self):
        """
        Calculates one iteration of the gammaT given the current state

        Note
        ----
        """
        newgamma = np.zeros(self.gammaL.shape)
        rho = np.logspace(-1.1, self.settings["rho_upper_limit"], self.settings["nrho"])
        eta = np.linspace(-1,1, self.settings["neta"])
        rr, ee = np.meshgrid(rho, eta)
        rrm = rhom(rr, ee)

        for gidx, (ax, ay) in enumerate(zip(self.x, self.y)):
            # first term alpha = beta = sigma = T
            TTT = vTTTscaled(ay, self.g, rr, ee) * np.power(rrm, -2) * self.Chi.chiTscaled(ax/rr, ay/rr) * self.Chi.chiTscaled(ax/rrm, ay/rrm)
            TTT /= (np.power(rr, 2.5) * self.gammaT[self._select_rho_idx(ax, rho)] + np.power(rrm, 2.5) * self.gammaT[self._select_rhom_idx(ax, rrm)])

            # second term alpha = beta = T, sigma = L
            TLT = vTLTscaled(ay, self.g, rr, ee) * np.power(rrm, -2) * self.Chi.chiLscaled(ax/rr, ay/rr) * self.Chi.chiTscaled(ax/rrm, ay/rrm)
            TLT /= (np.power(rr, 2.5) * self.gammaL[self._select_rho_idx(ax, rho)] + np.power(rrm, 2.5) * self.gammaT[self._select_rhom_idx(ax, rrm)])

            # third term alpha = T, beta = sigma = L
            TLL = vTLLscaled(ay, self.g, rr, ee) * np.power(rrm, -2) * self.Chi.chiLscaled(ax/rr, ay/rr) * self.Chi.chiLscaled(ax/rrm, ay/rrm)
            TLL /= (np.power(rr, 2.5) * self.gammaL[self._select_rho_idx(ax, rho)] + np.power(rrm, 2.5) * self.gammaL[self._select_rhom_idx(ax, rrm)])

            newgamma[gidx] = np.trapz(np.trapz(TTT + TLT + TLL, rho, axis=1), eta)

        return newgamma

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
        return np.tile(np.argmin(np.abs(x2d - x_over_rhos), axis=0), (self.settings["neta"], 1))

#------------------------------------------------------------------------------

    def _select_rhom_idx(self, ax, rhosm):
        """

        """
        # Build 3D versions of the arrays to find all indices without
        # python native iteration
        x3d = np.transpose(np.tile(self.x, list(rhosm.shape) + [1]), (2,0,1))
        rrm3d = np.tile(rhosm, (len(self.x), 1, 1))
        # calculate 
        x_over_rhosm = ax / rrm3d
        return np.argmin(np.abs(x3d - x_over_rhosm), axis=0)

###############################################################################
###############################################################################
###############################################################################