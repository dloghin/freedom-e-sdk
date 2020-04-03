set remotetimeout 300
target remote localhost:3333
load
break svm.c:258
cont
p endc
x &answer
quit
