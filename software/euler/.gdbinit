set remotetimeout 300
target remote localhost:3333
load
break euler.c:71
cont
x &e
quit
