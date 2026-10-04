%chk=hess.chk
%mem=32000MB
%nprocshared=8
#P PBE1PBE/6-311G(d,p) int=(ultrafine) nosymm opt=(Z-matrix, CalcAll,MaxCycle=150,MaxStep=20) pop=(hlygat) scf=(QC) scrf=(PCM,Read)

Gas-phase optimization from crystal

0 1
C
C 1 bnd2
C 2 bnd3 1 ang3
N 3 bnd4 2 ang4 1 dih4
H 3 bnd5 2 ang5 4 dih5
H 2 bnd6 1 ang6 3 dih6
Cl 1 bnd7 2 ang7 3 dih7
N 1 bnd8 2 ang8 7 dih8
C 8 bnd9 1 ang9 7 dih9
Cl 9 bnd10 8 ang10 1 dih10

bnd2=1.358767
bnd3=1.382915
ang3=115.766380
bnd4=1.334709
ang4=122.192575
dih4=0.770386
bnd5=0.930202
ang5=118.860257
dih5=179.993872
bnd6=0.930053
ang6=122.108465
dih6=-179.963144
bnd7=1.722786
ang7=120.513501
dih7=179.818385
bnd8=1.326919
ang8=124.411684
dih8=178.175253
bnd9=1.327419
ang9=113.167791
dih9=-179.230418
bnd10=1.725247
ang10=114.563115
dih10=178.793637

EPS=11.0

