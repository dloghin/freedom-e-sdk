set remotetimeout 300
target remote localhost:3333
load
break nb.c:255
cont
p endc
p answer
quit
