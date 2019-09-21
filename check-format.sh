#!/bin/bash

script_dir=$(dirname "$0")

if hash clang-format-7 2>/dev/null; then
   clang_format="clang-format-7"
elif hash clang-format-6.0 2>/dev/null; then
   clang_format="clang-format-6.0"
elif hash clang-format-4.0 2>/dev/null; then
   clang_format="clang-format-4.0"
fi

if [ -n "$clang_format" ]; then
   if [ -z "$1" ] || [ "$1" = "c" ]; then
      find $script_dir/src -regex '.*\.\(cpp\|hpp\|cc\|cxx\|h\|cu\)' -exec $clang_format -i {} \;
      find $script_dir/tests -regex '.*\.\(cpp\|hpp\|cc\|cxx\|h\|cu\)' -exec $clang_format -i {} \;
   fi
else
   echo "no clang-format found"
fi

if [ -e .venv/bin/python3 ]; then
   python=.venv/bin/python3
elif [ -e .env/bin/python3 ] ; then
   python=.env/bin/python3
elif hash python3 2> /dev/null; then
   python=python3
else
   echo "no python found"
fi

if [ -n "$python" ]; then
   if [ -z "$1" ] || [ "$1" = "p" ]; then
      $python -m yapf -i -r -p --style google $script_dir/src
      $python -m yapf -i -r -p --style google $script_dir/tests
      $python -m yapf -i -r -p --style google $script_dir/examples
   fi
fi
