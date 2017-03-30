#!/usr/bin/python

# description: Debye's Model of Substance Electric Properties
# author:      Vladimir V. Chukin (chukin@meteolab.ru)
# version:     2017-03-03
# license:     GNU General Public License 2.0
# coding:      utf8

from abc import ABCMeta, abstractmethod, abstractproperty
from numpy import *


class Debye:

    __metaclass__ = ABCMeta

    # Speed of light (m/s)
    SOL = 299792458.0


    def __init__(self, T, f):
        """ Initialization """
        self.setup(T, f)


    def setup(self, T, f):
        """ Setup All Properties """
        # Relaxation Time (s)
        self.tau   = self.getRelaxationTime(T)
        # Low-frequency Limit of Relative Permittivity
        self.eps_s = self.getStaticPermittivity(T)
        # High-frequency Limit of Relative Permittivity
        self.eps_8 = self.getPermittivityHighFrequencyLimit(T)
        # Real Part of Relative Permittivity
        self.eps1  = self.getPermittivityRealPart(f)
        # Imaginary Part of Relative Permittivity
        self.eps2  = self.getPermittivityImaginaryPart(f)
        # Real Part of Refractive Index
        self.n     = self.getRefractiveIndexRealPart(self.eps1, self.eps2)
        # Imaginary Part of Refractive Index
        self.a     = self.getRefractiveIndexImaginaryPart(self.eps1, self.eps2)


    def getRefractiveIndexRealPart(self):
        """ Real Part of Refractive Index """
        return n


    def getRefractiveIndexImaginaryPart(self):
        """ Imaginary Part of Refractive Index """
        return a


    def getRefractiveIndexRealPart(self, eps1, eps2):
        """ Real Part of Refractive Index """
        return sqrt( (sqrt(eps1*eps1+eps2*eps2)+eps1) / 2.0 )


    def getRefractiveIndexImaginaryPart(self, eps1, eps2):
        """ Imaginary Part of Refractive Index """
        return sqrt( (sqrt(eps1*eps1+eps2*eps2)-eps1) / 2.0 )


    def getPermittivityRealPart(self):
        """ Real Part of Relative Permittivity """
        return eps1


    def getPermittivityImaginaryPart(self):
        """ Imaginary Part of Relative Permittivity """
        return eps2


    def getPermittivityRealPart(self, n, a):
        """ Real Part of Relative Permittivity """
        return n*n - a*a


    def getPermittivityImaginaryPart(self, n, a):
        """ Imaginary Part of Relative Permittivity """
        return sqrt((n*n+a*a)*(n*n+a*a) - (n*n-a*a)*(n*n-a*a))


    def getPermittivityRealPart(self, f):
        """ Real Part of Relative Permittivity according Debye Theory """
        w = 2*pi*f
        return self.eps_8 + (self.eps_s-self.eps_8) / (1.0+w*w*self.tau*self.tau)


    def getPermittivityImaginaryPart(self, f):
        """ Imaginary Part of Relative Permittivity according Debye Theory """
        w = 2*pi*f
        return (self.eps_s-self.eps_8)*w*self.tau / (1.0+w*w*self.tau*self.tau)


    def set_eps_8(self, eps1_f, eps2_f, fs, f):
        """ Real Part of High-frequency Limit of Relative Permittivity """
        self.eps_8 = (eps1_f + eps1_f*(2.0*pi*fs*tau)*(2*pi*fs*tau) - eps_s) / (2.0*pi*fs*tau)/(2*pi*fs*tau)
        self.eps1  = self.getPermittivityRealPart(f)
        self.eps_8 = eps_s - eps2_f/(2.0*pi*fs*tau) - eps2_f*2.0*pi*fs*tau
        self.eps2  = self.getPermittivityImaginaryPart(f)
        self.n     = self.getRefractiveIndexRealPart(self.eps1, self.eps2)
        self.a     = self.getRefractiveIndexImaginaryPart(self.eps1, self.eps2)


    def set_eps_81(self, eps1_f, f):
        """ Real Part of High-frequency Limit of Relative Permittivity """
        self.eps_8 = (eps1_f + eps1_f*(2.0*pi*f*tau)*(2.0*pi*f*tau) - eps_s) / (2.0*pi*f*tau)/(2.0*pi*f*tau)
        self.eps1  = self.getPermittivityRealPart(f)
        self.eps2  = self.getPermittivityImaginaryPart(f)
        self.n     = self.getRefractiveIndexRealPart(self.eps1, self.eps2)
        self.a     = self.getRefractiveIndexImaginaryPart(self.eps1, self.eps2)


    def set_eps_82(self, eps2_f, f):
        """ Imaginary Part of High-frequency Limit of Relative Permittivity """
        self.eps_8 = eps_s - eps2_f/(2.0*pi*f*tau) - eps2_f*2.0*pi*f*tau
        self.eps1  = self.eps1(f)
        self.eps2  = self.eps2(f)
        self.n     = self.getRefractiveIndexRealPart(self.eps1, self.eps2)
        self.a     = self.getRefractiveIndexImaginaryPart(self.eps1, self.eps2)

    """
    @abstractmethod
    def getRelaxationTime(self, T):
        # T - temperature (K)
        pass


    @abstractmethod
    def getStaticPermittivity(self, T):
        # T - temperature (K)
        pass


    @abstractmethod
    def getPermittivityHighFrequencyLimit(self, T):
        # T - temperature (K)
        pass
    """



class DebyeWater(Debye):


    def getRelaxationTime(self, T):
        """ Relaxation Time in Water at a given Temperature (s) """
        return 1.478886E+05 * exp(-0.267225*T+6.140507E-04*T*T-4.651129E-07*T*T*T)


    def getStaticPermittivity(self, T):
        """ Static Permittivity at a given Temperature (Hz) """
        return 87.85306*exp(-0.00456992*(T-273.15))


    def getPermittivityHighFrequencyLimit(self, T):
        """ High-frequency Limit of Permittivity at a given Temperature (Hz) """
        return 5.27137 + 0.0216474*(T-273.15) - 0.00131198*(T-273.15)*(T-273.15)




class DebyeIce(Debye):


    def getRelaxationTime(self, T):
        """ Relaxation Time in Ice at a given Temperature (s) """
        return 53258.81606 * exp(0.077697*T - 1.1187E-03*T*T + 1.9937E-06*T*T*T)


    def getStaticPermittivity(self, T):
        """ Ice Static Permittivity at a given Temperature (Hz) """
        return 1982.084329 - 22.098339*T + 0.087238*T*T - 1.1596E-04*T*T*T


    def getPermittivityHighFrequencyLimit(self, T):
        """ Ice High-frequency Limit of Permittivity at a given Temperature (Hz) """
        return 3.1
