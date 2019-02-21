#!/bin/bash

script_dir=$(dirname "$0")

if hash clang-format-6.0 2>/dev/null; then
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

if hash autopep8 2>/dev/null; then
   autopep="autopep8"
   if [ -z "$1" ] || [ "$1" = "p" ]; then
      $autopep -i -a -r -j 8 $script_dir/src
      $autopep -i -a -r -j 8 $script_dir/tests
      $autopep -i -a -r -j 8 $script_dir/example
   fi
else
   echo "no autopep found"
fi
