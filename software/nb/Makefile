# Copyright 2019 SiFive, Inc #
# SPDX-License-Identifier: Apache-2.0 #

PROGRAM ?= nb

$(PROGRAM): $(wildcard *.c) $(wildcard *.h) $(wildcard *.S)

nb-x64: nb.c
	gcc -Wall $^ -o $@

clean:
	rm -f $(PROGRAM) $(PROGRAM).hex
