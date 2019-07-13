#!/bin/bash
rnas="1AJU 1AKX 1AM0 1BYJ 1DDY 1EHT 1EI2 1ET4 1EVV 1F1T 1F27 1FMN 1FUF 1FYP 1J7T 1KOC 1KOD 1LC4 1LVJ 1MWL 1NEM 1NTA 1NTB 1NYI 1O9M 1O15 1PBR 1Q8N 1QD3 1TN1 1TN2 1TOB 1U8D 1UTS 1UUD 1UUI 1XPF 1Y26 1YLS 1YRJ 1ZZ5 2AU4 2B57 2BE0 2BEE 2CKY 2EES 2EET 2EEU 2EEV 2EEW 2ESJ 2ET3 2ET4 2ET5 2ET8 2F4S 2F4T 2F4U 2FCX 2FCY 2FCZ 2FD0 2G5Q 2GDI 2GIS 2HO6 2HO7 2HOJ 2HOL 2HOM 2HOO 2JUK 2L1V 2L94 2LWK 2M4Q 2MIY 2MXS 2N0J 2NPY 2O3V 2O3W 2O3X 2O3Y 2OE5 2OE8 2QWY 2TOB 3B4A 3B4B 3B4C 3C44 3D0U 3D2G 3D2V 3D2X 3DIG 3DIL 3DIM 3DIX 3DIY 3DIZ 3DJ0 3DJ2 3DVV 3E5E 3E5F 3F2Q 3F2T 3F4H 3GCA 3GS8 3GX2 3GX3 3GX5 3GX6 3GX7 3IQN 3IQR 3NPQ 3RKF 3SKI 3SKL 3SLQ 3SUH 3SUX 3TD1 3TZR 3WRU 4B5R 4E8N 4E8Q 4F8U 4F8V 4FAW 4FEJ 4FEL 4FEN 4FEO 4FEP 4GPW 4GPX 4GPY 4JIY 4K32 4L81 4LVV 4LVW 4LVX 4LVY 4LVZ 4LW0 4LX5 4LX6 4NYA 4NYD 4NYG 4P3S 4P5J 4P20 4P95 4PDQ 4QLM 4QLN 4RZD 4TS2 4TZX 4TZY 4WCR 4XNR 4XW7 4XWF 4Y1I 4YAZ 4YB0 4ZNP 5BTP 5BWS 5BXK 5C45 5KX9 5LWJ 5UZA 6BFB"

# Create Header
# Using website for getting UNICODE for special characters
echo "<!DOCTYPE html>"
echo "<html>"
echo "  <head>"
echo "    <title>RNAPoser</title>"
echo "  </head>"
echo "  <body>"
echo "    <h3 align=\"center\" >Summary of decoys used to train and test RNAPoser</h3>"
echo "    <table  align=\"center\" border=\"3\" style=\"table-layout: auto width: 250px;\">"
echo "      <tr>"
echo "          <th>Index</th>"
echo "          <th>PDB</th>"
echo "          <th># Decoys</th>"
echo "          <th>Min.  RMSD (<span>&#8491;</span>)</th>"
echo "          <th>Max.  RMSD (<span>&#8491;</span>)</th>"
echo "          <th>Mean. RMSD (<span>&#8491;</span>)</th>"
echo "          <th>Ligand</th>"
echo "          <th>Poses</th>"
echo "      </tr>"
mkdir -p images/ > /dev/null


# Loop over RNAs and populate field
i=0
for rna in $rnas
do
    i=$((i+1))
    
    # get info
    ndecoys=`awk '{print $2}' renamed/${rna}/rmsd.txt | sort -n | tail -1`
    min=`awk '{print $3}' renamed/${rna}/rmsd.txt | sort -n | head -1 | awk '{printf "%4.2f\n", $1}'`
    max=`awk '{print $3}' renamed/${rna}/rmsd.txt | sort -n | tail -1 | awk '{printf "%4.2f\n", $1}'`
    mean=`awk '{x+=$3; next} END{print x/NR}' renamed/${rna}/rmsd.txt | awk '{printf "%4.2f\n", $1}'`
    
    # Index
    echo "<tr>"
    echo "<td halign=\"center\" style=\"word-wrap: break-word;\" valign=\"top\">"
    echo "  <p> $i </p>"
    echo "</td>"
    
    # PDBID
    echo "  <td halign=\"center\" style=\"word-wrap: break-word;\" valign=\"top\">"
    echo "    <p> ${rna} </p>"
    echo "  </td>"
    echo "  </td>"
    
    # Number of decoys
    echo "  <td halign=\"center\" style=\"word-wrap: break-word;\" valign=\"top\">"
    echo "    <p> ${ndecoys} </p>"
    echo "  </td>"
    echo "  </td>"
    
    # Minimum RMSD
    echo "  <td halign=\"center\" style=\"word-wrap: break-word;\" valign=\"top\">"
    echo "    <p> ${min} </p>"
    echo "  </td>"
    echo "  </td>"
    
    # Maximum RMSD
    echo "  <td halign=\"center\" style=\"word-wrap: break-word;\" valign=\"top\">"
    echo "    <p> ${max} </p>"
    echo "  </td>"
    echo "  </td>"
    
    # Mean RMSD
    echo "  <td halign=\"center\" style=\"word-wrap: break-word;\" valign=\"top\">"
    echo "    <p> ${mean} </p>"
    echo "  </td>"
    echo "  </td>"
    
    # 2D structure of Ligand
    alt=`echo "${rna}-ligand-file"`
    echo "<td halign=\"center\" style=\"word-wrap: break-word;\" valign=\"top\">"
    echo "  <p>"
    echo "    <a href=\"images/${rna}/lig.svg\">"
    echo "      <img src=\"images/${rna}/lig.svg\" style=\"width:185px\" alt=\"${alt}\">"
    echo "    </a><br>"
    echo "  </p>"
    echo "</td>"
    
    # GIF of poses
    alt=`echo "${rna}-poses-file"`
    echo "<td halign=\"center\" style=\"word-wrap: break-word;\" valign=\"top\">"
    echo "  <p>"
    echo "    <a href=\"images/${rna}/poses.gif\">"
    echo "      <img src=\"images/${rna}/poses.gif\" style=\"width:185px\" style=\"height:199px\"alt=\"${alt}\">"
    echo "    </a><br>"
    echo "  </p>"
    echo "</td>"

    echo "</tr>"
done

echo "    </table>"
echo "  </body>"
echo "</html>"
echo ""