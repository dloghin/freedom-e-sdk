set remotetimeout 300
target remote localhost:3333
load
break lr.c:281
cont
p endc
x &answer
quit
