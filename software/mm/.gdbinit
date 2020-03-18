set remotetimeout 300
target remote localhost:3333
load
break mm.c:138
cont
p endc
quit
