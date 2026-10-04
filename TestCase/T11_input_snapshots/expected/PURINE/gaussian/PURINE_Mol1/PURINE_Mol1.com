%chk=hess.chk
%mem=32000MB
%nprocshared=8
#P PBE1PBE/6-311G(d,p) int=(ultrafine) nosymm opt=(Z-matrix, CalcAll,MaxCycle=150,MaxStep=20) pop=(hlygat) scf=(QC) scrf=(PCM,Read)

Gas-phase optimization from crystal

0 1
C
C 1 bnd2
N 2 bnd3 1 ang3
C 3 bnd4 2 ang4 1 dih4
H 4 bnd5 3 ang5 2 dih5
N 2 bnd6 1 ang6 3 dih6
C 1 bnd7 2 ang7 6 dih7
N 7 bnd8 1 ang8 2 dih8
H 7 bnd9 1 ang9 8 dih9
N 1 bnd10 2 ang10 7 dih10
C 10 bnd11 1 ang11 2 dih11
H 11 bnd12 6 ang12 2 dih12
H 10 bnd13 1 ang13 11 dih13

bnd2=1.406970
bnd3=1.337469
ang3=123.869195
bnd4=1.338750
ang4=113.013674
dih4=-1.882002
bnd5=1.062508
ang5=109.848619
dih5=-169.376491
bnd6=1.379260
ang6=109.581784
dih6=178.968884
bnd7=1.393172
ang7=117.865804
dih7=179.483761
bnd8=1.329995
ang8=118.903360
dih8=1.199565
bnd9=1.000755
ang9=122.042638
dih9=-178.243167
bnd10=1.372847
ang10=105.072042
dih10=-178.820310
bnd11=1.326937
ang11=106.480880
dih11=-0.734981
bnd12=0.954115
ang12=121.245048
dih12=173.126717
bnd13=0.886144
ang13=118.423024
dih13=174.681901

EPS=11.0

