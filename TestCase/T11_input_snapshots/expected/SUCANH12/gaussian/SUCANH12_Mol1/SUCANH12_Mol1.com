%chk=hess.chk
%mem=32000MB
%nprocshared=8
#P PBE1PBE/6-311G(d,p) int=(ultrafine) nosymm opt=(Z-matrix, CalcAll,MaxCycle=150,MaxStep=20) pop=(hlygat) scf=(QC) scrf=(PCM,Read)

Gas-phase optimization from crystal

0 1
C
C 1 bnd2
C 2 bnd3 1 ang3
O 3 bnd4 2 ang4 1 dih4
O 3 bnd5 2 ang5 4 dih5
H 2 bnd6 1 ang6 3 dih6
H 2 bnd7 3 ang7 6 dih7
H 1 bnd8 2 ang8 3 dih8
H 1 bnd9 2 ang9 8 dih9
C 1 bnd10 9 ang10 8 dih10
O 10 bnd11 4 ang11 3 dih11

bnd2=1.524608
bnd3=1.502278
ang3=104.339688
bnd4=1.386549
ang4=110.238095
dih4=-1.648155
bnd5=1.198854
ang5=130.314114
dih5=-179.021414
bnd6=0.979310
ang6=113.522500
dih6=-116.869327
bnd7=0.898300
ang7=107.406874
dih7=117.437749
bnd8=0.973046
ang8=113.774416
dih8=120.628768
bnd9=0.938032
ang9=112.268135
dih9=129.060971
bnd10=1.499695
ang10=104.454276
dih10=117.325776
bnd11=1.198120
ang11=119.997013
dih11=-178.692456

EPS=11.0

