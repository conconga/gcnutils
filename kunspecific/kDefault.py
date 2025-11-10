#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
"""
Datei: 
Beschreibung: An object to parse the content of non-positional inputs (**kargs).
Autor: Luciano Auguto Kruk
Erstellt am: 2025.11.03
Version: 1.0.0
Lizenz: Please keep this header with the file.
GitHub: 
"""
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
##WWww=--  import section: --=wwWW##

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>#
#                                                                                  #
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>#
class kDefault:
    """
    """

    def __init__(self, **kargs):
        self.kargs = kargs
        #super().__init__(**kargs)

    def get(self, txt, default=[]):
        """
        Returns the configuration from `kargs` tagged with `txt`. If this key
        is not availble, then returns the default value.
        """
        assert default is not []

        try:
            ret = self.kargs[txt]
        except:
            ret = default

        return ret

    def __call__(self, *args):
        return self.get(*args)

    def show(self, len_header, str_header, str_value, str_unit=None):
        """
        Prints a log line with this format:

                              |- len_header
        0---------------------|---------------------------------....
        |          str_header = str_value [str_unit]
        """
        
        if str_unit is None:
            print("|{{:>{:d}s}} = {{:s}}".format(len_header-3).format(str_header, str_value))
        else:
            print("|{{:>{:d}s}} = {{:s}} [{{:s}}]".format(len_header-3).format(str_header, str_value, str_unit))

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>#
if __name__ == "__main__":
    dft = kDefault()

    print("123456789|123456789|123456789")
    dft.show(20, "var 1", str(dft))
    dft.show(20, "var 2", str(dft), "m")
    dft.show(20, "var 3", str(dft), "deg")

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>#
