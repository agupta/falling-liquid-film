#!/bin/bash

# Remake out/
echo Remaking out/...
rm -rf out/
mkdir out

echo Building...
qcc -Wall -O3 flow.c -o out/flow -L$BASILISK/gl -lglutils -lfb_glx -lGLU -lGLEW -lGL -lX11 -lm
cd out

echo Running...
./flow
