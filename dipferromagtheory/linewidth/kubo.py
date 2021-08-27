"""
Implements all the Kubo relaxation functions.
"""

import numpy as np
#------------------------------------------------------------------------------
from .susceptibility import Chi
#------------------------------------------------------------------------------

class Kubo():
    """

    """
    def __init__(self, g, chi=None, gamma=None):
        """

        """
        self.g = g
        self.chi = chi
        self.gamma = gamma

    def __call__(self, dir, q):
        """

        """
        pass

    def phiL(self, q):
        """

        """
        assert isinstance(self.chi, Chi)
#        assert isinstance(self.gamma, ) 
        return self.chi.chiL(q) / self.gamma()

    def phiT(self, q):
        """

        """
        assert isinstance(self.chi, Chi)
#        assert isinstance(self.gamma, ) 
        return self.chi.chiT(q) / self.gamma()