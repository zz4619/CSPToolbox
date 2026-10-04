%chk=hess.chk
%mem=32000MB
%nprocshared=8
#P PBE1PBE/6-311G(d,p) int=(ultrafine) nosymm opt=(Z-matrix, CalcAll,MaxCycle=150,MaxStep=20) pop=(hlygat) scf=(QC) scrf=(PCM,Read)

Gas-phase optimization from crystal

0 1
C
C 1 bnd2
C 2 bnd3 1 ang3
H 3 bnd4 2 ang4 1 dih4
F 2 bnd5 1 ang5 3 dih5
F 1 bnd6 2 ang6 3 dih6
C 1 bnd7 2 ang7 6 dih7
N 7 bnd8 1 ang8 2 dih8
H 8 bnd9 7 ang9 1 dih9
H 7 bnd10 1 ang10 8 dih10

bnd2=1.402977
bnd3=1.364369
ang3=108.358484
bnd4=0.962739
ang4=126.907762
dih4=-179.789354
bnd5=1.349657
ang5=125.131327
dih5=178.757579
bnd6=1.346856
ang6=125.518344
dih6=-179.750480
bnd7=1.370553
ang7=108.116326
dih7=179.841067
bnd8=1.367604
ang8=106.089590
dih8=-0.387078
bnd9=0.919798
ang9=120.981361
dih9=177.583728
bnd10=0.957234
ang10=127.028488
dih10=179.928167

EPS=11.0

