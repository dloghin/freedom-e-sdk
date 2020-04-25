set remotetimeout 300
target remote localhost:3333
load
break cifar10.c:90
break cifar10.c:101
cont
x &max
p idx
cont
p endc
quit
