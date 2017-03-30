#!/usr/bin/python
# coding: utf8

# description: Debye's Model of Ice Electric Properties
# author:      Vladimir V. Chukin (chukin@meteolab.ru)
# version:     2017-03-03
# license:     GNU General Public License 2.0
# coding:      utf8

from abc import ABCMeta, abstractmethod, abstractproperty
from numpy import *
import Debye


class DebyeIce(Debye):


    def __init__(self, T, f):
        """ Initialization """
        setup(T, f)


    def tau(self, T):
        """ Relaxation Time in Ice at a given Temperature (s) """
        return 53258.81606 * exp(0.077697*T - 1.1187E-03*T*T + 1.9937E-06*T*T*T)


    def eps_s(self, T):
        """ Ice Static Permittivity at a given Temperature (Hz) """
        return 1982.084329 - 22.098339*T + 0.087238*T*T - 1.1596E-04*T*T*T


    def eps_8(self, T):
        """ Ice High-frequency Limit of Permittivity at a given Temperature (Hz) """
        return 3.1
