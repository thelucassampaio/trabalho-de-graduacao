from __future__ import division
from .Potential import Potential, _APY_LOADED
from galpy.util import bovy_conversion
import numpy as np


if _APY_LOADED:
    from astropy import units


class JunqueiraPotential(Potential):
    #""" initialize a spiral arms potential """
    normalize= property() # turn off normalize
    def __init__(self, amp=0.1, ro=None, vo=None, amp_units='velocity2',
                 m=2, alpha=np.radians(14), sigma=4.7, e_s = 0.4, Ri = 8,
                 r_ref=1, phi_ref=0, Rs=0.3, H=0.125, omega=23, Cs=[1]):


        """
        NAME:
            __init__
        PURPOSE:
            initialize a spiral arms potential
        INPUT:
            :amp: amplitude to be applied to the potential (default: 1);
                        can be a Quantity with units of density. (:math:`amp = 4 \\pi G \\rho_0`)
            :ro: distance scales for translation into internal units (default from configuration file)
            :vo: velocity scales for translation into internal units (default from configuration file)
            :N: number of spiral arms
            :alpha: pitch angle of the logarithmic spiral arms in radians (can be Quantity)
            :r_ref: fiducial radius where :math:`\\rho = \\rho_0` (:math:`r_0` in the paper by Cox and Gomez) (can be Quantity)
            :phi_ref: reference angle (:math:`\\phi_p(r_0)` in the paper by Cox and Gomez) (can be Quantity)
            :Rs: radial scale length of the drop-off in density amplitude of the arms (can be Quantity)
            :H: scale height of the stellar arm perturbation (can be Quantity)
            :Cs: list of constants multiplying the :math:`\cos(n \\gamma)` terms
            :omega: rotational pattern speed of the spiral arms (can be Quantity)
        OUTPUT:
            (none)
        HISTORY:
            Started - 2017-05-12  Jack Hong (UBC)

            Completed - 2017-07-04 Jack Hong (UBC)
        """

        Potential.__init__(self, amp=amp, ro=ro, vo=vo, amp_units=amp_units)

        if _APY_LOADED:
            if isinstance(alpha, units.Quantity):
                alpha = alpha.to(units.rad).value
            if isinstance(r_ref, units.Quantity):
                r_ref = r_ref.to(units.kpc).value / self._ro
            if isinstance(phi_ref, units.Quantity):
                phi_ref = phi_ref.to(units.rad).value
            if isinstance(Rs, units.Quantity):
                Rs = Rs.to(units.kpc).value / self._ro
            if isinstance(H, units.Quantity):
                H = H.to(units.kpc).value / self._ro
            if isinstance(omega, units.Quantity):
                omega = omega.to(units.km / units.s / units.kpc).value \
                        / bovy_conversion.freq_in_kmskpc(self._vo, self._ro)
            if isinstance(sigma, units.Quantity):
                sigma = sigma.to(units.kpc).value
            if isinstance(e_s, units.Quantity):
                e_s = e_s.to(1/units.kpc).value
            if isinstance(Ri, units.Quantity):
                Ri = Ri.to(units.kpc).value


        self._m = -m  # trick to flip to left handed coordinate system; flips sign for phi and phi_ref, but also alpha.
        self._alpha = -alpha  # we don't want sign for alpha to change, so flip alpha. (see eqn. 3 in the paper)
        self._sin_alpha = np.sin(-alpha)
        self._tan_alpha = np.tan(-alpha)
        self._sigma = sigma
        self._e_s = e_s        # Epsilon S
        self._Ri = Ri
        self._r_ref = r_ref
        self._phi_ref = phi_ref
        self._omega = omega

        self._Rs = Rs
        self._H = H
        self._Cs = self._Cs0 = np.array(Cs)
        self._ns = self._ns0 = np.arange(1, len(Cs) + 1)
        self._rho0 = 1 / (4 * np.pi)
        self._HNn = self._HNn0 = self._H * self._m * self._ns0

        self.isNonAxi = True   # Potential is not axisymmetric
        self.hasC = False       # Potential has C implementation to speed up orbit integrations
        self.hasC_dxdv = False  # Potential has C implementation of second derivatives

    #""" Evaluate the potential at the given coordinates. (without the amp factor; handled by super class) """
    def _evaluate(self, R, z, phi=0, t=0):
        """
        NAME:
            _evaluate
        PURPOSE:
            Evaluate the potential at the given coordinates. (without the amp factor; handled by super class)
        INPUT:
            :param R: galactocentric cylindrical radius
            :param z: vertical height
            :param phi: azimuth
            :param t: time
        OUTPUT:
            :return: Phi(R, z, phi, t)
        HISTORY:
            2017-05-12  Jack Hong (UBC)
        """


        fm = self._fm(R)
        dfm_dR = self._dfm_dR(R)

        R2_sigma2 = R**2/self._sigma**2
        cos = np.cos(self._m * phi - fm)



        return - R * np.exp(-R2_sigma2 * (1 - cos) - self._e_s * R )

    def _Rforce(self, R, z, phi=0, t=0):
        """
        NAME:
            _Rforce
        PURPOSE:
            Evaluate the radial force for this potential at the given coordinates. (-dPhi/dR)
        INPUT:
            :param R: galactocentric cylindrical radius
            :param z: vertical height
            :param phi: azimuth
            :param t: time
        OUTPUT:
            :return: the radial force
        HISTORY:
            2017-05-12  Jack Hong (UBC)
        """
        fm = self._fm(R)
        dfm_dR = self._dfm_dR(R)

        R2_sigma2 = R**2/self._sigma**2
        sen = np.sin(self._m * phi - fm)
        cos = np.cos(self._m * phi - fm)


        return np.exp(-R2_sigma2 * (1 - cos) - self._e_s * R ) \
                * (1 - R*(2*R/self._sigma**2 * (1 - cos)\
                - R2_sigma2 * sen * dfm_dR + self._e_s ))


    def _zforce(self, R, z, phi=0, t=0):
        """
        NAME:
            _zforce
        PURPOSE:
            Evaluate the vertical force for this potential at the given coordinates. (-dPhi/dz)
        INPUT:
            :param R: galactocentric cylindrical radius
            :param z: vertical height
            :param phi: azimuth
            :param t: time
        OUTPUT:
            :return: the vertical force
        HISTORY:
            2017-05-25  Jack Hong (UBC)
        """
        return 0

    def _phiforce(self, R, z, phi=0, t=0):
        """
        NAME:
            _phiforce
        PURPOSE:
            Evaluate the azimuthal force in cylindrical coordinates. (-dPhi/dphi)
        INPUT:
            :param R: galactocentric cylindrical radius
            :param z: vertical height
            :param phi: azimuth
            :param t: time
        OUTPUT:
            :return: the azimuthal force
        HISTORY:
            2017-05-25  Jack Hong (UBC)
        """
        fm = self._fm(R)

        R2_sigma2 = R**2/self._sigma**2
        sen = np.sin(self._m * phi - fm)
        cos = np.cos(self._m * phi - fm)

        return - R**3/self._sigma**2 * self._m \
                * np.exp(-R2_sigma2 * (1 - cos) - self._e_s * R ) * sen

    def _R2deriv(self, R, z, phi=0, t=0):
        """
        NAME:
            _R2deriv
        PURPOSE:
            Evaluate the second (cylindrical) radial derivative of the potential.
             (d^2 potential / d R^2)
        INPUT:
            :param R: galactocentric cylindrical radius
            :param z: vertical height
            :param phi: azimuth
            :param t: time
        OUTPUT:
            :return: the second radial derivative
        HISTORY:
            2017-05-31  Jack Hong (UBC)
        """
        R2_sigma2 = R**2/self._sigma**2
        R_sigma2 = R/self._sigma**2
        fm = self._fm(R)

        dfm_dR = self._dfm_dR(R)

        sen = np.sin(self._m * phi - fm)
        cos = np.cos(self._m * phi - fm)

        # Potencial/R
        pot_R = - np.exp(-R2_sigma2 * (1 - cos) - self._e_s * R )  # Potencial/R

        dpot_dR = pot_R * (1 - R*(2*R/self._sigma**2 * (1 - cos) \
                - R2_sigma2 * sen * dfm_dR + self._e_s ))


        return dpot_dR * (1 - R * self._e_s + R2_sigma2 * (2 * cos + sen - 2))+\
                pot_R * (R2_sigma2 * (2 * sen * dfm_dR - cos * dfm_dR) + \
                R_sigma2 * (8 * cos + 2 * sen - 4) - self._e_s)

    def _z2deriv(self, R, z, phi=0, t=0):
        """
        NAME:
            _z2deriv
        PURPOSE:
            Evaluate the second (cylindrical) vertical derivative of the potential.
             (d^2 potential / d z^2)
        INPUT:
            :param R: galactocentric cylindrical radius
            :param z: vertical height
            :param phi: azimuth
            :param t: time
        OUTPUT:
            :return: the second vertical derivative
        HISTORY:
            2017-05-26  Jack Hong (UBC)
        """


        return 0

    def _phi2deriv(self, R, z, phi=0, t=0):
        """
        NAME:
            _phi2deriv
        PURPOSE:
            Evaluate the second azimuthal derivative of the potential in cylindrical coordinates.
            (d^2 potential / d phi^2)
        INPUT:
            :param R: galactocentric cylindrical radius
            :param z: vertical height
            :param phi: azimuth
            :param t: time
        OUTPUT:
            :return: d^2 potential / d phi^2
        HISTORY:
            2017-05-29 Jack Hong (UBC)
        """

        R2_sigma2 = R**2/self._sigma**2
        fm = self._fm(R)

        dfm_dR = self._dfm_dR(R)

        sen = np.sin(self._m * phi - fm)
        cos = np.cos(self._m * phi - fm)

        pot = - R * np.exp(-R2_sigma2 * (1 - cos) - self._e_s * R )

        dpot_dphi = - R2_sigma2 * self._m * pot * sen

        return - R2_sigma2 * self._m * (dpot_dphi * sen + pot * cos * self._m)

    def _Rzderiv(self, R, z, phi=0., t=0.):
        """
        NAME:
            _Rzderiv
        PURPOSE:
            Evaluate the mixed (cylindrical) radial and vertical derivative of the potential
            (d^2 potential / dR dz).
        INPUT:
            :param R: galactocentric cylindrical radius
            :param z: vertical height
            :param phi: azimuth
            :param t: time
        OUTPUT:
            :return: d^2 potential / dR dz
        HISTORY:
            2017-05-12  Jack Hong (UBC)
        """

        return 0

    def _Rphideriv(self, R, z, phi=0,t=0):
        """
        NAME:
            _Rphideriv
        PURPOSE:
            Return the mixed radial and azimuthal derivative of the potential in cylindrical coordinates
             (d^2 potential / dR dphi)
        INPUT:
            :param R: galactocentric cylindrical radius
            :param z: vertical height
            :param phi: azimuth
            :param t: time
        OUTPUT:
            :return: the mixed radial and azimuthal derivative
        HISTORY:
            2017-06-09  Jack Hong (UBC)
        """

        R2_sigma2 = R**2/self._sigma**2
        R_sigma2 = R/self._sigma**2
        fm = self._fm(R)

        dfm_dR = self._dfm_dR(R)

        sen = np.sin(self._m * phi - fm)
        cos = np.cos(self._m * phi - fm)

        pot = - R * np.exp(-R2_sigma2 * (1 - cos) - self._e_s * R )

        pot_R = - np.exp(-R2_sigma2 * (1 - cos) - self._e_s * R ) #Potencial/R
        dpot_dR = pot_R * (1 - R*(2*R/self._sigma**2 * (1 - cos) \
                - R2_sigma2 * sen * dfm_dR + self._e_s ))

        return 2 * R_sigma2 * self._m * pot * sen \
                + R2_sigma2 * self._m * dpot_dR * sen \
                + R2_sigma2 * self._m * pot * cos * (-dfm_dR)






    def OmegaP(self):
        """
        NAME:
            OmegaP
        PURPOSE:
            Return the pattern speed. (used to compute the Jacobi integral for orbits).
        INPUT:
            :param self
        OUTPUT:
            :return: the pattern speed
        HISTORY:
            2017-06-09  Jack Hong (UBC)
        """
        return self._omega

    def _fm(self, R):
        """ Return the shape function """
        return self._m / self._tan_alpha * np.log(R/self._Ri) + self._phi_ref

    def _dfm_dR(self, R):
        """ Return the first derivative of the shape function wrt R. """
        return self._m / self._tan_alpha / R





    def _gamma(self, R, phi):
        """Return gamma. (eqn 3 in the paper)"""
        return self._m * (phi - self._phi_ref - np.log(R / self._r_ref) / self._tan_alpha)

    def _dgamma_dR(self, R):
        """Return the first derivative of gamma wrt R."""
        return -self._m / R / self._tan_alpha

    def _K(self, R):
        """Return numpy array from K1 up to and including Kn. (eqn. 5)"""
        return self._ns * self._m / R / self._sin_alpha


    def _dK_dR(self, R):
        """Return numpy array of dK/dR from K1 up to and including Kn."""
        return -self._ns * self._m / R**2 / self._sin_alpha

    def _B(self, R):
        """Return numpy array from B1 up to and including Bn. (eqn. 6)"""
        HNn_R = self._HNn / R

        return HNn_R / self._sin_alpha * (0.4 * HNn_R / self._sin_alpha + 1)

    def _dB_dR(self, R):
        """Return numpy array of dB/dR from B1 up to and including Bn."""
        return -self._HNn / R**3 / self._sin_alpha**2 * (0.8 * self._HNn + R * self._sin_alpha)

    def _D(self, R):
        """Return numpy array from D1 up to and including Dn. (eqn. 7)"""
        return (0.3 * self._HNn**2 / self._sin_alpha / R
                + self._HNn + R * self._sin_alpha) / (0.3 * self._HNn + R * self._sin_alpha)

    def _dD_dR(self, R):
        """Return numpy array of dD/dR from D1 up to and including Dn."""
        HNn_R_sina = self._HNn / R / self._sin_alpha

        return HNn_R_sina * (0.3 * (HNn_R_sina + 0.3 * HNn_R_sina**2. + 1) / R / (0.3 * HNn_R_sina + 1)**2
                             - (1/R * (1 + 0.6 * HNn_R_sina) / (0.3 * HNn_R_sina + 1)))


    #
