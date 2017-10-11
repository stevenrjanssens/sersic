"""sersic.py -- functions to work with Sersic profiles

Notes:
    Most equations from Graham & Driver (2005PASA...22..118G).

    Approximation for b from Appendix A in MacArthur et al.
    (2003ApJ...582..689M).
"""

import numpy
from scipy import special

def sersic_b(n):
    """Approx Sersic b param."""
    if numpy.isscalar(n):
        narr = numpy.array([n])
    else:
        narr = numpy.array(n)
    
    big_n = (narr > 0.36)
    small_n = (narr <= 0.36)
    
    b = numpy.zeros(len(narr))
    
    # use expansion for n > 0.36
    b[big_n] = 2*narr[big_n] - 1./3 + 4./405/narr[big_n] + 46./25515/narr[big_n]**2 + 131./1148175/narr[big_n]**3 - 2194697./30690717750/narr[big_n]**4
    
    # use polynomial approx for small n
    b[small_n] = 0.01945 - 0.8902*narr[small_n] + 10.95*narr[small_n]**2 - 19.67*narr[small_n]**3 + 13.43*narr[small_n]**4
    
    if numpy.isscalar(n):
        return b[0]
    else:
        return b
    
def meanmueff_from_mag(mag, R_e, q=1):
    """
    mag: total mag
    R_e: R_e in arcsec
    q:   axis ratio, q=1 for circularized R_e
    """
    return mag + 2.5*numpy.log10(2*numpy.pi*(R_e**2)*q)

def meanmueff_from_mueff(mu_e, sersic_index):
    """
    mu_e: surface brightness at R_e
    sersic_index: sersic index n
    """
    b = sersic_b(sersic_index)
    sersic_f = sersic_index * numpy.exp(b) / numpy.power(b, 2*sersic_index) * special.gamma(2*sersic_index)
    return mu_e - 2.5*numpy.log10(sersic_f)

def meanmueff_from_mu0(mu_0, sersic_index):
    """
    mu_0: central surface brightness
    sersic_index: sersic index n
    """
    mu_e = mueff_from_mu0(mu_0, sersic_index)
    return meanmueff_from_mueff(mu_e, sersic_index)

def mu0_from_mueff(mu_e, sersic_index):
    """
    mu_e: surface brightness at R_e
    sersic_index: sersic index n
    """
    return mu_e - 2.5*sersic_b(sersic_index)/numpy.log(10)

def mueff_from_mu0(mu_0, sersic_index):
    """
    mu_0: central surface brightness
    sersic_index: sersic index n
    """
    return mu_0 + 2.5*sersic_b(sersic_index)/numpy.log(10)

def mag_from_meanmueff(meanmueff, R_e, q=1):
    """
    meanmueff: mean surface brightness within R_e
    R_e:       R_e in arcsec
    q:         axis ratio, q=1 for circularized R_e
    """
    return meanmueff - 2.5*numpy.log10(2*numpy.pi*(R_e**2)*q)

def mag_from_mueff(mu_e, R_e, sersic_index, q=1):
    """
    mueff:        mean surface brightness within R_e
    R_e:          R_e in arcsec
    sersic_index: sersic index n
    q:            axis ratio, q=1 for circularized R_e
    """
    meanmueff = meanmueff_from_mueff(mu_e, sersic_index)
    return mag_from_meanmueff(meanmueff, R_e, q)
