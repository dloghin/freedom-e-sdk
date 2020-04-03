set remotetimeout 300
target remote localhost:3333
load
break knn.c:117
cont
p endc
p answer
quit
