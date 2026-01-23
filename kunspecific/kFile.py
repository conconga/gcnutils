#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
"""
Datei: 
Beschreibung: An object to parse existing files and its name elements.
Autor: Luciano Auguto Kruk
Erstellt am: 2026.01.23
Version: 1.0.0
Lizenz: Please keep this header with the file.
GitHub: 
"""
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>
##WWww=--  import section: --=wwWW##

import os

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>#
#                                                                                  #
#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>#
class kFile:
    def __init__(self, file):
        self.path   = file
        self.isfile = True if os.path.isfile(file) else False

        if self.isfile:
            self.sizeB   = os.path.getsize(file)
            self.dirname = os.path.dirname(file)
            self.base    = os.path.splitext(os.path.basename(file))[0]
            self.ext     = os.path.splitext(os.path.basename(file))[1]
        else:
            print(f"'{self.path}' is not a file")

    def __repr__(self):
        txt = [
            f">> '{self.path:s}':",
            f"   isfile  = {self.isfile}",
            f"   size    = {self.sizeB / 1024 / 1024:1.0f} [MB]",
            f"   dirname = {self.dirname}",
            f"   base    = {self.base}",
            f"   ext     = {self.ext}",
        ]
        return "\n".join(txt)

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>#
if __name__ == "__main__":
    file = kFile("choose_an_existing_file")
    print(file)

#>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>--<<..>>#
