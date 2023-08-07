#!/bin/sh

icc \
	-std=c++17 \
	-DLINUX_ICC_BUILD \
	-I./fvconsole -I./fvcases -I./fv \
	fvconsole/*.cpp fvcases/*.cpp fv/global_stat.cpp \
	-o linux_icc_exe
