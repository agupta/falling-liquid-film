#!/bin/bash

# # Remake out/
# echo Remaking out/...
# rm -rf out/
# mkdir out


# Remake bin/
echo Remaking bin/...
rm -rf bin/
mkdir bin

# Build
echo Building...
qcc -Wall -O3 be_cleaned.c -o bin/be_cleaned -L$BASILISK/gl -lglutils -lfb_glx -lGLU -lGLEW -lGL -lX11 -lm

echo Running...
./bin/be_cleaned
