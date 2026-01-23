#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
"""
Datei: 
Beschreibung: An object to make time-measurements along any execution.
Autor: Luciano Auguto Kruk
Erstellt am: 2026.01.23
Version: 1.0.0
Lizenz: Please keep this header with the file.
GitHub: 
"""
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
##WWww=--  import section: --=wwWW##

import time # to log in [ns]

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>#
#                                                                                  #
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>#

class kTimeLogger:
    def __init__(self):
        self.logger  = list()
        self.start()

    def _get_time(self):
        return time.time_ns()

    def start(self):
        """
        Starts/Reset the timer.
        """
        self.t_start = self._get_time()

    def stamp(self):
        """
        Returns current elapsed time since start().
        : [ns] :
        """
        return self._get_time() - self.t_start


#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>#
if __name__ == "__main__":
    t = kTimeLogger()

    # some delay:
    a = list()
    for _ in range(1000):
        a.append(range(10))

    print("elapsed time = {:f}[ms]".format(t.stamp() * 1e-6))

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>#
