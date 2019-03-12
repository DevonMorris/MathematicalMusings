#!/bin/zsh

setopt shwordsplit
exts="aux bbl blg brf idx ilg ind lof log lol lot out toc synctex.gz fls fdb_latexmk pdf"

for ext in $exts; do
     rm -f **/*.$ext
done
