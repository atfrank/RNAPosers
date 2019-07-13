#!/bin/bash

# loop over RNAs
rnas="1aju 1akx 1am0 1byj 1DDY 1eht 1ei2 1ET4 1EVV 1f1t 1f27 1fmn 1FUF 1fyp 1j7t 1koc 1kod 1LC4 1lvj 1mwl 1nem 1nta 1ntb 1nyi 1o9m 1o15 1pbr 1q8n 1qd3 1TN1 1TN2 1tob 1u8d 1uts 1uud 1uui 1XPF 1Y26 1YLS 1yrj 1ZZ5 2au4 2b57 2be0 2bee 2CKY 2EES 2EET 2EEU 2EEV 2EEW 2esj 2et3 2et4 2ET5 2et8 2f4s 2f4t 2f4u 2fcx 2fcy 2fcz 2fd0 2g5q 2GDI 2GIS 2HO6 2HO7 2HOJ 2HOL 2HOM 2HOO 2juk 2l1v 2l94 2lwk 2m4q 2miy 2mxs 2n0j 2NPY 2O3V 2o3w 2o3x 2O3Y 2oe5 2oe8 2QWY 2tob 3B4A 3B4B 3B4C 3c44 3D0U 3D2G 3D2V 3D2X 3DIG 3DIL 3DIM 3DIX 3DIY 3DIZ 3DJ0 3DJ2 3dvv 3E5E 3E5F 3F2Q 3F2T 3F4H 3GCA 3GS8 3GX2 3GX3 3GX5 3GX6 3GX7 3IQN 3IQR 3NPQ 3RKF 3SKI 3SKL 3SLQ 3SUH 3SUX 3TD1 3TZR 3WRU 4b5r 4E8N 4E8Q 4F8U 4F8V 4FAW 4FEJ 4FEL 4FEN 4FEO 4FEP 4GPW 4GPX 4GPY 4JIY 4k32 4L81 4LVV 4LVW 4LVX 4LVY 4LVZ 4LW0 4LX5 4LX6 4NYA 4NYD 4nyg 4p3s 4P5J 4P20 4P95 4PDQ 4QLM 4QLN 4RZD 4TS2 4TZX 4TZY 4wcr 4XNR 4XW7 4XWF 4Y1I 4YAZ 4yb0 4znp 5btp 5bws 5bxk 5c45 5kx9 5lwj 5uza 6bfb"

# goto location
cd /Users/atfrank/Desktop/test_rna_poser

for rna in $rnas
do
    rnaupper=`echo "$rna" | awk '{print toupper($0)}'`
    mkdir -p renamed/${rnaupper}
    cp -v ${rna}/${rna}_lig.sd renamed/${rnaupper}/lig.sd
    cp -v ${rna}/${rna}_rdock.mol2 renamed/${rnaupper}/receptor.mol2
    cp -v ${rna}/poses/all_rmsd_sorted.sd renamed/${rnaupper}/poses.sd
    cp -v ${rna}/rmsds/all_rmsd_sorted.txt renamed/${rnaupper}/rmsd.txt
done