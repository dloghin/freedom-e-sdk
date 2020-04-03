#pragma once

#if !defined(__x86_64__)
#include <metal/cpu.h>
#endif

unsigned long long read_cycles() {
#if defined(__x86_64__)
  return 0;
#else
  struct metal_cpu *mycpu = metal_cpu_get(0);
  return metal_cpu_get_timer(mycpu);
#endif
}
