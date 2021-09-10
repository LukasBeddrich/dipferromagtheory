"""
Implements all the susceptibility functions.
"""
import numpy as np

class Chi():
    """

    """
    def __init__(self, j, g, xi):
        """
        Parameters
        ----------
        j   : float
            (strength of) exchange interaction
        g   : float
            ratio of dipolar to exchange interaction
        xi  : float
            mag. correlation length
        """
        self.j = j
        self.g = g
        self.xi = xi

#------------------------------------------------------------------------------

    def set(self, key, value):
        """
        Parameters
        ----------
        key : str
            (strength of) exchange interaction
        g   : float
            ratio of dipolar to exchange interaction
        xi  : float
            mag. correlation length
        """
        if key == "j":
            self.j = value
        elif key == "g":
            self.g = value
        elif key == "xi":
            self.xi = value
        else:
            print(f"key {key} is not known! Acceptable keys are 'j', 'g', 'xi'")

#------------------------------------------------------------------------------
    
    def chiL(self, q):
        """
        Parameters
        ----------
        q       : float, numpy.ndarray
            wave vector / momentum transfer

        Return
        ------
        chiL    : float, numpy.ndarray
            longitudinal mag. susceptibility
        """
        return np.power(self.j * (q**2 + self.xi**(-2) + self.g), -1)

#------------------------------------------------------------------------------

    def chiT(self, q):
        """
        Parameters
        ----------
        q       : float, numpy.ndarray
            wave vector / momentum transfer

        Return
        ------
        chiT    : float, numpy.ndarray
            transversal mag. susceptibility
        """
        return np.power(self.j * (q**2 + self.xi**(-2)), -1)

#------------------------------------------------------------------------------

    def chiLscaled(self, x, y):
        """
        Parameters
        ----------
        x       : float, numpy.ndarray
            1/(xi * q)
        y       : float, numpy.ndarray
            g**(1/2)/ q

        Return
        ------
        chiL    : float, numpy.ndarray
            scaled longitudinal mag. susceptibility
        """
        return np.power(1 + x**2 + y**2, -1)

#------------------------------------------------------------------------------

    def chiTscaled(self, x, y):
        """
        Parameters
        ----------
        x       : float, numpy.ndarray
            1/(xi * q)
        y       : float, numpy.ndarray
            g**(1/2)/ q

        Return
        ------
        chiT    : float, numpy.ndarray
            scaled transversal mag. susceptibility
        """
        return np.power(1 + x**2, -1)

###############################################################################
