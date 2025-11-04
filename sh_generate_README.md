#!/bin/bash

#cat README.md.before_pandoc | pandoc --from markdown+grid_tables --toc -s -f markdown-tex_math_dollars --number-sections --to gfm - > README.md

find . -type f | grep "README.md.before_pandoc" | while read nome; do
    name=`basename "${nome}"`
    dir=`dirname "${nome}"`
    echo "converting ${name} at '${dir}'..."

    #vvvvvvvvvvvvvvv#
    pushd "${dir}"

        uniquefile=`sh_uniqfile -f /tmp/`
        withouttoc="${name}"

        cat ${withouttoc} | pandoc \
        --from markdown+grid_tables-latex_macros+raw_tex-tex_math_dollars-tex_math_double_backslash-tex_math_single_backslash \
        --log=/tmp/pandoc.log \
        --toc -s --number-sections --to gfm - | sed -n '1,/^[[:blank:]]*$/p' > ${uniquefile}

        cat ${uniquefile} ${withouttoc} > README.md

    popd
    #^^^^^^^^^^^^^^^#
done

