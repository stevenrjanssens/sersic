import numpy

import sersic

def test_b():
    # 2.5*b/ln(10) from Graham & Driver
    assert (2.5*sersic.sersic_b(n=1)/numpy.log(10) - 1.822) < 1e-3
    assert (2.5*sersic.sersic_b(n=4)/numpy.log(10) - 8.327) < 1e-3

def test_meanmueff_from_mueff():
    assert (sersic.meanmueff_from_mueff(mueff=0, n=1) + 0.699) < 1e-2
    assert (sersic.meanmueff_from_mueff(mueff=0, n=4) + 1.393) < 1e-2

def test_mueff_from_meanmueff():
    assert (sersic.mueff_from_meanmueff(meanmueff=0, n=1) - 0.699) < 1e-2
    assert (sersic.mueff_from_meanmueff(meanmueff=0, n=4) - 1.393) < 1e-2

def test_meanmueff_from_mu0():
    # only uses mueff_from_mu0 and meanmueff_from_mueff
    pass

def test_mu0_from_mueff():
    assert (sersic.mu0_from_mueff(mueff=0, n=1) + 1.822) < 1e-3
    assert (sersic.mu0_from_mueff(mueff=0, n=4) + 8.327) < 1e-3

def test_mueff_from_mu0():
    assert (sersic.mueff_from_mu0(mu0=0, n=1) - 1.822) < 1e-3
    assert (sersic.mueff_from_mu0(mu0=0, n=4) - 8.327) < 1e-3
