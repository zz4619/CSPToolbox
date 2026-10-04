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
N 4 bnd5 3 ang5 2 dih5
H 5 bnd6 4 ang6 3 dih6
O 4 bnd7 3 ang7 5 dih7
H 3 bnd8 2 ang8 4 dih8
H 2 bnd9 1 ang9 3 dih9
N 1 bnd10 2 ang10 9 dih10
O 1 bnd11 2 ang11 10 dih11
H 11 bnd12 1 ang12 2 dih12

bnd2=1.417839
bnd3=1.341561
ang3=118.639475
bnd4=1.440230
ang4=120.012417
dih4=-0.406421
bnd5=1.338405
ang5=115.206923
dih5=-0.066784
bnd6=1.039759
ang6=117.066829
dih6=-178.989478
bnd7=1.261898
ang7=124.126060
dih7=-179.685828
bnd8=1.045667
ang8=124.879979
dih8=-178.545114
bnd9=1.042637
ang9=117.324061
dih9=-179.356809
bnd10=1.306261
ang10=123.587707
dih10=-178.797310
bnd11=1.335362
ang11=117.831611
dih11=179.605834
bnd12=1.031492
ang12=114.915649
dih12=-174.392749

EPS=11.0

