%chk=hess.chk
%mem=32000MB
%nprocshared=8
#P PBE1PBE/6-311G(d,p) int=(ultrafine) nosymm opt=(Z-matrix, CalcAll,MaxCycle=150,MaxStep=20) pop=(hlygat) scf=(QC) scrf=(PCM,Read)

Gas-phase optimization from crystal

0 1
C
C 1 bnd2
N 2 bnd3 1 ang3
H 3 bnd4 2 ang4 1 dih4
H 3 bnd5 2 ang5 4 dih5
O 2 bnd6 1 ang6 3 dih6
C 1 bnd7 2 ang7 3 dih7
N 7 bnd8 1 ang8 2 dih8
N 1 bnd9 2 ang9 7 dih9
O 9 bnd10 1 ang10 2 dih10
H 10 bnd11 9 ang11 1 dih11

bnd2=1.489831
bnd3=1.314650
ang3=116.364594
bnd4=0.950889
ang4=124.533985
dih4=0.376109
bnd5=1.001677
ang5=114.723885
dih5=176.505236
bnd6=1.228382
ang6=120.274045
dih6=179.052816
bnd7=1.438414
ang7=119.159000
dih7=176.299790
bnd8=1.130758
ang8=179.613474
dih8=146.572698
bnd9=1.275273
ang9=117.321734
dih9=178.863032
bnd10=1.354548
ang10=113.658383
dih10=-178.028927
bnd11=0.987588
ang11=106.394226
dih11=-176.746346

EPS=11.0

