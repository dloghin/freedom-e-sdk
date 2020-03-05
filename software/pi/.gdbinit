set remotetimeout 300
target remote localhost:3333
load
break pi.c:177
cont
x &pi
quit
