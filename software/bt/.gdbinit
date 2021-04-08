set remotetimeout 300
target remote localhost:3333
load
break bt.c:269
cont
p verified
p endc
quit
