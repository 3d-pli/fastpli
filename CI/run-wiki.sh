#!/bin/bash
set -euo pipefail

env-CI/bin/pip3 install -q -r examples/requirements.txt

for f in examples/*.ipynb; do
   echo "converting $f"
   env-CI/bin/jupyter-nbconvert --to script --stdout $f |
   sed '/get_ipython*/d' |
   sed '3 i\import itertools' |
   sed '4 i\img_counter=itertools.count()' |
   sed  "s@plt\.show()@plt\.savefig(f'{\"$f\"[:-6]}_{next(img_counter)}\.png',dpi=(150), bbox_inches='tight')@g" > $f.py
   # cat $f.py
   echo "running $f"
   env-CI/bin/python3 $f.py
   rm $f.py
   echo "done $f"
done
