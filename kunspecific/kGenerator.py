#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
"""
Datei: 
Beschreibung: These classes simulate generation of random step-signal with stochastic characteristics.
Autor: Luciano Auguto Kruk
Erstellt am: 01.10.2025
Version: 1.0.0
Lizenz: Please keep this header with the file.
GitHub: 
"""
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
##WWww=--  import section: --=wwWW##
import numpy as np

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>#
#                                                                                  #
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>#
class kSignalGeneratorCommonUniformNormal:
    def __init__(self, **kargs):
        assert kargs["ufn_stepduration"] > 0
        assert kargs["ufn_max"] > kargs["ufn_min"]

        self.ufn_min          = kargs["ufn_min"]
        self.ufn_max          = kargs["ufn_max"]
        self.ufn_stepduration = kargs["ufn_stepduration"]
        self.ufn_firststep    = kargs["ufn_firststep"]

        self.time_last_sample = -1.0

    def update(self, t):
        assert t > self.time_last_sample
        #print("in update({:1.02f})".format(t))

        if self.time_last_sample < 0:
            # first sample:
            if t < self.ufn_firststep:
                self.last_sample = 0.0
                self.time_next_change = self.ufn_firststep
            else:
                self.last_sample = self._sample()
                self.time_next_change = self.ufn_stepduration
            #print("time_next_change = {:f}".format(self.time_next_change))

        elif t > self.time_next_change:
            self.last_sample = self._sample()
            while self.time_next_change < t:
                self.time_next_change += self.ufn_stepduration
            #print("time_next_change = {:f}".format(self.time_next_change))
        else:
            # no new sample; no change
            pass

        self.time_last_sample = t
        return self.last_sample

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>#
#                                                                                  #
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>#
class kSignalGeneratorUniform (kSignalGeneratorCommonUniformNormal):
    def __init__(self, **kargs):
        super().__init__(**kargs)

    def _sample(self):
        return self.ufn_min + ((self.ufn_max - self.ufn_min) * np.random.rand())

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>#
#                                                                                  #
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>#
class kSignalGeneratorNormal (kSignalGeneratorCommonUniformNormal):
    def __init__(self, **kargs):
        super().__init__(**kargs)

        self.ufn_sigma = kargs["ufn_sigma"]

    def _sample(self):
        mean = 0.5 * (self.ufn_max + self.ufn_min)
        out  = mean + (self.ufn_sigma * np.random.randn())
        if (out < self.ufn_min):
            out = self.ufn_min
        elif out > self.ufn_max:
            out = self.ufn_max

        return out

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>#
#                                                                                  #
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>#
class kSignalGenerator:
    UNIFORM = 1
    NORMAL  = 2

    def __init__(self, Ts, typeGen=None, **kargs):
        """
        This class generates sampled signals with particular characteristics
        based on provided configuration.  The necessary configuration is
        defined by the type of generator.

        Ts : sample period [s]
        --------------------------------------------------------------------------------
        typeGen == self.UNIFORM:
            Generates random steps with uniform distribution among min and max values.

            ufn_min          : minimal value allowed to be sampled
            ufn_max          : maximal value allowed to be sampled
            ufn_stepduration : [s] duration of each step
            ufn_firststep    : [s] zeros are sampled up to the firststep time

        --------------------------------------------------------------------------------
        typeGen == self.NORMAL:
            Generates random step with normal distribution among min and max
            values, and around the mean of these two values.

            ufn_* as UNIFORM
            ufn_sigma       : 1x standard deviation of the sampled value
        --------------------------------------------------------------------------------
        """

        assert Ts > 0
        self.Ts        = Ts
        self.curr_time = -1.0 # negative to force a sample at t=0

        if typeGen == self.UNIFORM:
            self.obj = kSignalGeneratorUniform( **kargs )
        elif typeGen == self.NORMAL:
            self.obj = kSignalGeneratorNormal( **kargs )
        else:
            raise(NameError("this type of signal is still not supported"))

    def get_next_sample(self):
        """
        returns:
            t_sample : [s]
            sample   : respective generated sample
        """
        if self.curr_time < 0:
            self.curr_time = 0
        else:
            self.curr_time += self.Ts

        return self.curr_time, self.obj.update(self.curr_time)

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>#
#                                                                                  #
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>#
if __name__ == "__main__":

        sample_freq_Hz = 3
        max_time       = 4.1
        dt             = 1./sample_freq_Hz
        max_time       = max_time
        signal_uf      = kSignalGenerator(dt, kSignalGenerator.UNIFORM, ufn_min=1, ufn_max=1.3, ufn_stepduration=0.8, ufn_firststep=0.2)
        signal_g       = kSignalGenerator(dt, kSignalGenerator.NORMAL, ufn_min=1, ufn_max=1.3, ufn_stepduration=0.8, ufn_firststep=0.2, ufn_sigma=0.2)

        time = 0.0
        print()
        print("** UNIFORM **")
        while time < max_time:
            time += dt
            print("{:1.2f} :  {:1.02f}".format(*signal_uf.get_next_sample()))

        time = 0.0
        print()
        print("** NORMAL **")
        while time < max_time:
            time += dt
            print("{:1.2f} :  {:1.02f}".format(*signal_g.get_next_sample()))

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>#
