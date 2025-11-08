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
class kSignalTypes:
    UNIFORM = 1
    NORMAL  = 2
    STEP    = 3

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>#
#                                                                                  #
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>#
class kSignalGeneratorUpdate (kSignalTypes):
    def __init__(self, **kargs):
        self.time_last_sample = -1.0

    def update(self, t):
        assert t > self.time_last_sample
        #print("in update({:1.02f})".format(t))

        if self.time_last_sample < 0:
            # first sample:
            if t < self.ufn_firststep:
                self.last_sample = 0.0
                if self.typeGen in (self.UNIFORM, self.NORMAL, self.STEP):
                    self.time_next_change = self.ufn_firststep
                else:
                    raise(NameError("still not supported"))
            else:
                self.last_sample = self._sample()
                if self.typeGen in (self.UNIFORM, self.NORMAL):
                    self.time_next_change = self.ufn_stepduration
                elif self.typeGen == self.STEP:
                    self.time_next_change = self.ufn_highduration
            #print("time_next_change = {:f}".format(self.time_next_change))

        elif t > self.time_next_change:
            while self.time_next_change < t:
                # even if the new sample is not necessary, the call to a new
                # sample might change some internal state.
                self.last_sample = self._sample()

                if self.typeGen in (self.UNIFORM, self.NORMAL):
                    self.time_next_change += self.ufn_stepduration
                else:
                    self.time_next_change += \
                            self.ufn_highduration if self.current_state == "high" else self.ufn_lowduration
            #print("time_next_change = {:f}".format(self.time_next_change))

        else:
            # no new sample; no change
            pass

        self.time_last_sample = t
        return self.last_sample

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>#
#                                                                                  #
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>#
class kSignalGeneratorUniform (kSignalGeneratorUpdate):
    def __init__(self, **kargs):
        super().__init__(**kargs)
        self.typeGen = self.UNIFORM

        assert kargs["ufn_stepduration"] > 0
        assert kargs["ufn_max"] > kargs["ufn_min"]

        self.ufn_min          = kargs["ufn_min"]
        self.ufn_max          = kargs["ufn_max"]
        self.ufn_stepduration = kargs["ufn_stepduration"]
        self.ufn_firststep    = kargs["ufn_firststep"]


    def _sample(self):
        return self.ufn_min + ((self.ufn_max - self.ufn_min) * np.random.rand())

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>#
#                                                                                  #
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>#
class kSignalGeneratorNormal (kSignalGeneratorUpdate):
    def __init__(self, **kargs):
        super().__init__(**kargs)
        self.typeGen = self.NORMAL

        assert kargs["ufn_stepduration"] > 0
        assert kargs["ufn_max"] > kargs["ufn_min"]

        self.ufn_min          = kargs["ufn_min"]
        self.ufn_max          = kargs["ufn_max"]
        self.ufn_stepduration = kargs["ufn_stepduration"]
        self.ufn_firststep    = kargs["ufn_firststep"]
        self.ufn_sigma        = kargs["ufn_sigma"]

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
class kSignalGeneratorStep (kSignalGeneratorUpdate):
    def __init__(self, **kargs):
        super().__init__(**kargs)
        self.typeGen = self.STEP

        assert kargs["ufn_low"] < kargs["ufn_high"]
        assert kargs["ufn_highduration"] >= 0
        assert kargs["ufn_lowduration"] >= 0
        assert kargs["ufn_firststep"] >= 0

        self.ufn_low          = kargs["ufn_low"]
        self.ufn_high         = kargs["ufn_high"]
        self.ufn_highduration = kargs["ufn_highduration"]
        self.ufn_lowduration  = kargs["ufn_lowduration"]
        self.ufn_firststep    = kargs["ufn_firststep"]

        self.current_state = "low"

    def _sample(self):
        if self.current_state == "low":
            self.current_state = "high"
            out = self.ufn_high
        else:
            self.current_state = "low"
            out = self.ufn_low

        return out


#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>#
#                                                                                  #
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>#
class kSignalGenerator (kSignalTypes):

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
        typeGen == self.STEP:
            Generates a sequence of steps with same amplitude either low or high.

            ufn_low          : amplitude of low level
            ufn_high         : amplitude of high level
            ufn_highduration : [s] duration at high level
            ufn_lowduration  : [s] duration at low level
            ufn_firststep    : [s] low level is sampled up to the firststep time
        --------------------------------------------------------------------------------
        """

        assert Ts > 0
        self.Ts        = Ts
        self.curr_time = -1.0 # negative to force a sample at t=0

        # backup:
        self.typeGen = typeGen

        if typeGen == self.UNIFORM:
            self.obj = kSignalGeneratorUniform( **kargs )
        elif typeGen == self.NORMAL:
            self.obj = kSignalGeneratorNormal( **kargs )
        elif typeGen == self.STEP:
            self.obj = kSignalGeneratorStep( **kargs )
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
    import matplotlib.pylab as plt
    SMALL_SIZE = 7
    MEDIUM_SIZE = 9
    BIGGER_SIZE = 11
    plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

    sample_freq_Hz = 100
    max_time       = 20
    dt             = 1./sample_freq_Hz
    max_time       = max_time
    signal_uniform = kSignalGenerator(dt, kSignalTypes.UNIFORM, ufn_min=1, ufn_max=1.3, ufn_stepduration=0.8, ufn_firststep=0.2)
    signal_normal  = kSignalGenerator(dt, kSignalTypes.NORMAL, ufn_min=1, ufn_max=1.3, ufn_stepduration=0.8, ufn_firststep=0.2, ufn_sigma=0.2)
    signal_step    = kSignalGenerator(dt, kSignalTypes.STEP, ufn_low=-3, ufn_lowduration=1, ufn_high=3, ufn_highduration=4, ufn_firststep=1.5)

    time = 0.0
    log_time    = list()
    log_uniform = list()
    log_normal  = list()
    log_step    = list()
    while time < max_time:
        time += dt
        log_time.append(time)
        log_uniform.append(signal_uniform.get_next_sample()[1])
        log_normal.append(signal_normal.get_next_sample()[1])
        log_step.append(signal_step.get_next_sample()[1])

    plt.figure(1).clf()
    fig,ax = plt.subplots(3,1,num=1)
    fig.set_size_inches(5,4)
    for i,j in enumerate( [log_uniform, log_normal, log_step] ):
        ax[i].plot(log_time, j)
        ax[i].grid(True)
        ax[i].set_ylabel([ 'uniform', 'normal', 'step' ][i] )

    plt.show(block=False)
    plt.savefig("output_kGen.png")

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>#
