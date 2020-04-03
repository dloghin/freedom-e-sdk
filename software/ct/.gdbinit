set remotetimeout 300
target remote localhost:3333
load
break ct.c:322
cont
p endc
p answer
quit
