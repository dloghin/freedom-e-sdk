set remotetimeout 300
target remote localhost:3333
load
break euler.c:81
cont
x &e
p endc
quit
