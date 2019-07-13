#!/bin/bash

# natvigate and setup
mkdir -p /home/afrankz/local_software/repo/RNAPoser/tests/input/rna-protein/ten_3drpc
cd /home/afrankz/local_software/repo/RNAPoser/tests/input/rna-protein/ten_3drpc
rm -rf 1B7F 1DFU 1JBS 1P6V 1WPU 1WSU 2ASB 2BH2 2QUX 3BX2
mkdir -p 1B7F 1DFU 1JBS 1P6V 1WPU 1WSU 2ASB 2BH2 2QUX 3BX2

# download files
wget http://biophy.hust.edu.cn/3dRPC/download/187992324/all -O 1B7F.complex.tar.gz
wget http://biophy.hust.edu.cn/3dRPC/download/193930742/all -O 1DFU.complex.tar.gz
wget http://biophy.hust.edu.cn/3dRPC/download/647193534/all -O 1JBS.complex.tar.gz
wget http://biophy.hust.edu.cn/3dRPC/download/832606522/all -O 1P6V.complex.tar.gz
wget http://biophy.hust.edu.cn/3dRPC/download/278160761/all -O 1WPU.complex.tar.gz
wget http://biophy.hust.edu.cn/3dRPC/download/500814065/all -O 1WSU.complex.tar.gz
wget http://biophy.hust.edu.cn/3dRPC/download/165464264/all -O 2ASB.complex.tar.gz
wget http://biophy.hust.edu.cn/3dRPC/download/337477182/all -O 2BH2.complex.tar.gz
wget http://biophy.hust.edu.cn/3dRPC/download/182534598/all -O 2QUX.complex.tar.gz
wget http://biophy.hust.edu.cn/3dRPC/download/776468585/all -O 3BX2.complex.tar.gz

# decompress decoys
tar -zxvf 1B7F.complex.tar.gz -C /home/afrankz/local_software/repo/RNAPoser/tests/input/rna-protein/ten_3drpc/1B7F/
tar -zxvf 1DFU.complex.tar.gz -C /home/afrankz/local_software/repo/RNAPoser/tests/input/rna-protein/ten_3drpc/1DFU/
tar -zxvf 1JBS.complex.tar.gz -C /home/afrankz/local_software/repo/RNAPoser/tests/input/rna-protein/ten_3drpc/1JBS/
tar -zxvf 1P6V.complex.tar.gz -C /home/afrankz/local_software/repo/RNAPoser/tests/input/rna-protein/ten_3drpc/1P6V/
tar -zxvf 1WPU.complex.tar.gz -C /home/afrankz/local_software/repo/RNAPoser/tests/input/rna-protein/ten_3drpc/1WPU/
tar -zxvf 1WSU.complex.tar.gz -C /home/afrankz/local_software/repo/RNAPoser/tests/input/rna-protein/ten_3drpc/1WSU/
tar -zxvf 2ASB.complex.tar.gz -C /home/afrankz/local_software/repo/RNAPoser/tests/input/rna-protein/ten_3drpc/2ASB/
tar -zxvf 2BH2.complex.tar.gz -C /home/afrankz/local_software/repo/RNAPoser/tests/input/rna-protein/ten_3drpc/2BH2/
tar -zxvf 2QUX.complex.tar.gz -C /home/afrankz/local_software/repo/RNAPoser/tests/input/rna-protein/ten_3drpc/2QUX/
tar -zxvf 3BX2.complex.tar.gz -C /home/afrankz/local_software/repo/RNAPoser/tests/input/rna-protein/ten_3drpc/3BX2/

# copy inputs
cp -v ../older-decoys/1B7F_*_original_3dRPC.pdb .
cp -v ../older-decoys/1DFU_*_original_3dRPC.pdb .
cp -v ../older-decoys/1JBS_*_original_3dRPC.pdb .
cp -v ../older-decoys/1P6V_*_original_3dRPC.pdb .
cp -v ../older-decoys/1WPU_*_original_3dRPC.pdb .
cp -v ../older-decoys/1WSU_*_original_3dRPC.pdb .
cp -v ../older-decoys/2ASB_*_original_3dRPC.pdb .
cp -v ../older-decoys/2BH2_*_original_3dRPC.pdb .
cp -v ../older-decoys/2QUX_*_original_3dRPC.pdb .
cp -v ../older-decoys/3BX2_*_original_3dRPC.pdb .

# generate files for RNAPoser
/export/apps/pymol-1.8.6/pymol/pymol -cr ../generate_rna_poser_files.py