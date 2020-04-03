# Copyright 2019 SiFive, Inc #
# SPDX-License-Identifier: Apache-2.0 #

PROGRAM ?= lr

$(PROGRAM): $(wildcard *.c) $(wildcard *.h) $(wildcard *.S)

lr-x64: lr.c
	gcc -g -Wall $^ -o $@

clean:
	rm -f $(PROGRAM) $(PROGRAM).hex
