set remotetimeout 300
target remote localhost:3333
load
break mm.c:74
break mm.c:89
#cont
#p endc
#quit
