
rnas="2g5k 2o3w 2xnw 3fu2 3la5 3mum 3q50 3sd3 3slm 4erj 4fe5 4jf2 4lx5 4nya 4xwf 4yb0 5c7w 2b57 1f1t 2ydh 3npn 4aob 4kqy 4l81 5kpy"

for rna in ${rnas}
do

cd ${rna}/
grep -A1 "<SCORE>" poses.sd | grep -v -e "--" | grep -v "<SCORE>" > s1.txt
grep -A1 "<SCORE.INTER>" poses.sd | grep -v -e "--" | grep -v "<SCORE.INTER>" > s2.txt
grep -A1 "<SCORE.INTER.CONST>" poses.sd | grep -v -e "--" | grep -v "<SCORE.INTER.CONST>" > s3.txt
grep -A1 "<SCORE.INTER.POLAR>" poses.sd | grep -v -e "--" | grep -v "<SCORE.INTER.POLAR>" > s4.txt
grep -A1 "<SCORE.INTER.ROT>" poses.sd | grep -v -e "--" | grep -v "<SCORE.INTER.ROT>" > s5.txt
grep -A1 "<SCORE.INTER.SOLV>" poses.sd | grep -v -e "--" | grep -v "<SCORE.INTER.SOLV>" > s6.txt
grep -A1 "<SCORE.INTER.VDW>" poses.sd | grep -v -e "--" | grep -v "<SCORE.INTER.VDW>" > s7.txt
grep -A1 "<SCORE.INTER.norm>" poses.sd | grep -v -e "--" | grep -v "<SCORE.INTER.norm>" > s8.txt
grep -A1 "<SCORE.INTRA>" poses.sd | grep -v -e "--" | grep -v "<SCORE.INTRA>" > s9.txt
grep -A1 "<SCORE.INTRA.DIHEDRAL>" poses.sd | grep -v -e "--" | grep -v "<SCORE.INTRA.DIHEDRAL>" > s10.txt
grep -A1 "<SCORE.INTRA.DIHEDRAL.0>" poses.sd | grep -v -e "--" | grep -v "<SCORE.INTRA.DIHEDRAL.0>" > s11.txt
grep -A1 "<SCORE.INTRA.POLAR>" poses.sd | grep -v -e "--" | grep -v "<SCORE.INTRA.POLAR>" > s12.txt
grep -A1 "<SCORE.INTRA.POLAR.0>" poses.sd | grep -v -e "--" | grep -v "<SCORE.INTRA.POLAR.0>" > s13.txt
grep -A1 "<SCORE.INTRA.REPUL>" poses.sd | grep -v -e "--" | grep -v "<SCORE.INTRA.REPUL>" > s14.txt
grep -A1 "<SCORE.INTRA.REPUL.0>" poses.sd | grep -v -e "--" | grep -v "<SCORE.INTRA.REPUL.0>" > s15.txt
grep -A1 "<SCORE.INTRA.SOLV>" poses.sd | grep -v -e "--" | grep -v "<SCORE.INTRA.SOLV>" > s16.txt
grep -A1 "<SCORE.INTRA.SOLV.lig_0>" poses.sd | grep -v -e "--" | grep -v "<SCORE.INTRA.SOLV.lig_0>" > s17.txt
grep -A1 "<SCORE.INTRA.VDW>" poses.sd | grep -v -e "--" | grep -v "<SCORE.INTRA.VDW>" > s18.txt
grep -A1 "<SCORE.INTRA.VDW.0>" poses.sd | grep -v -e "--" | grep -v "<SCORE.INTRA.VDW.0>" > s19.txt
grep -A1 "<SCORE.INTRA.norm>" poses.sd | grep -v -e "--" | grep -v "<SCORE.INTRA.norm>" > s20.txt
grep -A1 "<SCORE.RESTR>" poses.sd | grep -v -e "--" | grep -v "<SCORE.RESTR>" > s21.txt
grep -A1 "<SCORE.RESTR.CAVITY>" poses.sd | grep -v -e "--" | grep -v "<SCORE.RESTR.CAVITY>" > s22.txt
grep -A1 "<SCORE.RESTR.norm>" poses.sd | grep -v -e "--" | grep -v "<SCORE.RESTR.norm>" > s23.txt
grep -A1 "<SCORE.SYSTEM>" poses.sd | grep -v -e "--" | grep -v "<SCORE.SYSTEM>" > s24.txt
grep -A1 "<SCORE.SYSTEM.CONST>" poses.sd | grep -v -e "--" | grep -v "<SCORE.SYSTEM.CONST>" > s25.txt
grep -A1 "<SCORE.SYSTEM.DIHEDRAL>" poses.sd | grep -v -e "--" | grep -v "<SCORE.SYSTEM.DIHEDRAL>" > s26.txt
grep -A1 "<SCORE.SYSTEM.SOLV>" poses.sd | grep -v -e "--" | grep -v "<SCORE.SYSTEM.SOLV>" > s27.txt
grep -A1 "<SCORE.SYSTEM.norm>" poses.sd | grep -v -e "--" | grep -v "<SCORE.SYSTEM.norm>" > s28.txt
grep -A1 "<SCORE.heavy>" poses.sd | grep -v -e "--" | grep -v "<SCORE.heavy>" > s29.txt
grep -A1 "<SCORE.norm>" poses.sd | grep -v -e "--" | grep -v "<SCORE.norm>" > s30.txt

paste -d " " s1.txt s2.txt s3.txt s4.txt s5.txt s6.txt s7.txt s8.txt s9.txt s10.txt s11.txt s12.txt s13.txt s14.txt s15.txt s16.txt s17.txt s18.txt s19.txt s20.txt s21.txt s22.txt s23.txt s24.txt s25.txt s26.txt s27.txt s28.txt s29.txt s30.txt > Scores.txt
rm s*
cd ../
done






