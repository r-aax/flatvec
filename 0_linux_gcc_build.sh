#!/bin/sh

g++ \
	-std=c++17 \
	-DLINUX_GCC_BUILD \
	-I./fvconsole -I./fvcases -I./fv \
	fvconsole/*.cpp fvcases/*.cpp fv/*.cpp \
	-o linux_gcc_exe
