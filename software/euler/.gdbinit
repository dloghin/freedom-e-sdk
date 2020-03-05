set remotetimeout 300
target remote localhost:3333
load
break euler.c:45
cont
x &e
x e
quit
