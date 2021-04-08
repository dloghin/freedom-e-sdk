set remotetimeout 300
target remote localhost:3333
load
break sp.c:269
cont
p verified
p endc
quit

