ener = readBin('ener.dat','double',100000)

length(ener)

plot(ener,type='l')

eround=round(ener,2)
table(eround)
plot(eround,type='l')


