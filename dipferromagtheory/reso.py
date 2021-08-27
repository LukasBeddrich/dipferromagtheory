"""

"""

### IMPORTS
import numpy as np
import json

from dipferromagtheory import resdir



def sigma_lambda(q, q_mean, dlam):
    """
    Standard deviation of the wavelength distribution in q space

    Parameters
    ----------
    q      :   float
        momentum/q transfer value
    q_mean :   float
        average momentum/q transfer value
    dlam   :   float
        wavelength spread in percent

    Return
    ------
           :   float
        standard deviation of the wavelength spread in q space
    
    Note
    ----
    triangular spead to sigma conversion:
    (sig_lam / lam)**2 = 1/6 (delta lam / lam)**2
    """

    return (q_mean - q) * dlam / np.sqrt(6)

#-------------------------------------------------------------------------------------

def sigma_Q_all(q, q_mean, r1, l1, r2, l2, dx, dy, lam, dlam):
    """
    combined q resolution sigma from goemetric and wavelength spread
    
    Parameters
    ----------
    q      :   float
        momentum/q transfer value
    q_mean :   float
        average momentum/q transfer value
    r1     :   float
        circular radius of source aperture
    r2     :   float
        circular radius of sample aperture
    l1     :   float
        source-sample distance
    l2     :   float
        sample-detector distance
    dx     :   float
        horizontal pixel width
    dy     :   float
        vertical pixel width
    lam    :   float
        wavelength
    dlam   :   float
        wavelength spread in percent (delta lambda/lambda)
        for a triangular distr.
        
    Return
    ------
           :   float
        standard deviation of the q resolution is gaussian approximation
    """

    geomcontr = (l2 / l1 * r1 * 0.5)**2 + ((l1 + l2) / l1 * r2 * 0.5)**2
    
    sig_sqr_Q_wave = sigma_lambda(q, q_mean, dlam)**2
    
    sig_sqr_Q_geom = (2 * np.pi / lam / l2)**2 * (2 * geomcontr + 1/12. * (dx**2 + dy**2))
    
    return np.sqrt(sig_sqr_Q_geom + sig_sqr_Q_wave)

#-------------------------------------------------------------------------------------
 
def resolution_Q(q, q_mean, sigma_q):
    """
    Calculates the q resolution approximated in gaussian form
    
    Parameters
    ----------
    q       :   float
        momentum/q transfer value (intefration parameter)
    q_mean  :   float
        average momentum/q transfer value
    sigma_q :   float
        standard deviation of the q resolution
    """
    return np.exp(-0.5 * ((q_mean - q) / sigma_q)**2) / np.sqrt(2 * np.pi * sigma_q**2)

#------------------------------------------------------------------------------

class Maskqdists:
    """

    """
    path = f"{resdir}/qdists_Ni2019.json"
    dists = {}
    # with open(path, "r") as jsonfile:
    #     dists = json.load(jsonfile)

    @classmethod
    def keys(cls):
        return cls.dists.keys()

    @classmethod
    def get(cls, key):
        return cls.dists[key]
    
    @classmethod
    def init(cls):
        with open(cls.path, "r") as jsonfile:
            cls.dists = json.load(jsonfile)
    
    @classmethod
    def update(cls, path):
        cls.path = path

Maskqdists.init()