%chk=hess.chk
%mem=32000MB
%nprocshared=8
#P PBE1PBE/6-311G(d,p) int=(ultrafine) nosymm opt=(Z-matrix, CalcAll,MaxCycle=150,MaxStep=20) pop=(hlygat) scf=(QC) scrf=(PCM,Read)

Gas-phase optimization from crystal

0 1
C
C 1 bnd2
C 2 bnd3 1 ang3
C 3 bnd4 2 ang4 1 dih4
H 4 bnd5 3 ang5 2 dih5
H 3 bnd6 2 ang6 4 dih6
H 2 bnd7 1 ang7 3 dih7
N 1 bnd8 2 ang8 3 dih8
H 8 bnd9 1 ang9 2 dih9
C 1 bnd10 2 ang10 8 dih10
O 10 bnd11 1 ang11 8 dih11
H 11 bnd12 10 ang12 1 dih12
O 10 bnd13 1 ang13 11 dih13

bnd2=1.372897
bnd3=1.399639
ang3=107.218644
bnd4=1.367404
ang4=107.261237
dih4=-0.402747
bnd5=0.980706
ang5=131.482307
dih5=-178.894319
bnd6=0.995125
ang6=125.935480
dih6=-178.809398
bnd7=0.965355
ang7=125.407167
dih7=179.989284
bnd8=1.364659
ang8=107.872888
dih8=0.114275
bnd9=0.859011
ang9=123.699634
dih9=176.799044
bnd10=1.442416
ang10=130.873774
dih10=-179.964741
bnd11=1.315339
ang11=113.826975
dih11=-179.811433
bnd12=0.954559
ang12=114.062943
dih12=-174.068759
bnd13=1.230118
ang13=123.560839
dih13=179.582320

EPS=11.0

