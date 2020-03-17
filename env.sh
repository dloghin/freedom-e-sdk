#!/bin/bash

export RISCV_PATH=$HOME/git/riscv/riscv-tools-32/riscv64-unknown-elf-gcc-8.3.0-2019.08.0-x86_64-linux-ubuntu14
export RISCV_OPENOCD_PATH=$HOME/git/riscv/riscv-tools-32/riscv-openocd-0.10.0-2019.08.2-x86_64-linux-ubuntu14
export PATH=$PATH:$RISCV_PATH/bin:$RISCV_OPENOCD_PATH/bin
