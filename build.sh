#!/bin/sh

set -xe

CFLAGS="-Wall -Wextra -ggdb `pkg-config --cflags raylib`"
LIBS="`pkg-config --libs raylib` -lglfw -lm -ldl -lpthread"

mkdir -p ./build/

# No hot reload
# clang $CFLAGS -o ./build/musializer ./src/plug.c ./src/musializer.c $LIBS -L./build/

# Hot reload
clang $CFLAGS -o ./build/libplug.so -fPIC -shared ./src/plug.c $LIBS
clang $CFLAGS -DHOTRELOAD -o ./build/musializer ./src/musializer.c $LIBS -L./build/
