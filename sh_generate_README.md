#!/bin/bash

#cat README.md.before_pandoc | pandoc --from markdown+grid_tables --toc -s -f markdown-tex_math_dollars --number-sections --to gfm - > README.md
cat README.md.before_pandoc | pandoc \
--from markdown+grid_tables-latex_macros+raw_tex-tex_math_dollars-tex_math_double_backslash-tex_math_single_backslash \
--log=/tmp/pandoc.log \
--toc -s --number-sections --to gfm - > README.md

