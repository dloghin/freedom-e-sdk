# Copyright 2019 SiFive, Inc #
# SPDX-License-Identifier: Apache-2.0 #

PROGRAM ?= ct

$(PROGRAM): $(wildcard *.c) $(wildcard *.h) $(wildcard *.S)

$(PROGRAM)-x64: $(PROGRAM).c
	gcc -g -Wall $^ -o $@

clean:
	rm -f $(PROGRAM) $(PROGRAM).hex
