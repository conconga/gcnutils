#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
"""
Datei: 
Beschreibung:
Autor: Luciano Auguto Kruk
Erstellt am: 31.10.2025
Version: 1.0.0
Lizenz: Please keep this header with the file.
GitHub: 
"""
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
##WWww=--  import section: --=wwWW##

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>#
#                                                                                  #
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>#
class kBase:
    """
    The `object.__init__()` method does not accept arguments. This class is just to
    forward the call without arguments.
    Don´t ask me why.... 
    """
    def __init__(self, **kargs):
        super().__init__()

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>#
