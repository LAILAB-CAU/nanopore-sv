In this project, the genomic features of SVs were annotated in two senarios. In one senario, TE feature was considered and the Mo17 reference genome was first annotated in the priority of TE region, exon, intron, upstream 5kbp region (Up5k), downstream 5 kbp region (Down5k), and intergenic region. In the other senario, TE region was not considered, and the Mo17 reference genome was annotated in the priority of CDS, UTR, intron, up5K, down5K, intergenic. In both senarios, two adjacent regions with the same annotation were merged into one region. SVs were annotated through two steps. 

Step 1. The partition of Mo17 reference genome: a custom script was used to annotate genome features from the downloaded gff3 file at https://download.maizegdb.org/Zm-Mo17-REFERENCE-CAU-1.0/Zm-Mo17-REFERENCE-CAU-1.0_Zm00014a.1.gff3.gz. 

In senario 1, different features had the following priority: TE region, exon, intron, upstream 5kbp region (Up5k), downstream 5 kbp region (Down5k), and intergenic region. TE annotation (file) was from the 'TE annotation' section, which conducted EDTA and annotated the Mo17 genome into 5 TE categories: Copia, Gypsy, Helitron, DNA/T, MITEs, and other TEs. The annotation results are in maize_*_uniq.

```
cd with_TE
## Generate mo17.EDTA.5.TE.type.bed from EDTA output Zm-Mo17-REFERENCE-CAU-1.0.fa.mod.EDTA.TEanno.gff3
grep -v "#" ./section4/EDTA/Zm-Mo17-REFERENCE-CAU-1.0.fa.mod.EDTA.TEanno.gff3 | awk -F"Classification=" '{print $1"\t"$2}' - | awk '{print $1"\t"$4"\t"$5"\t"$10}' -| awk -F";" '{print $1}' > ./section4/mo17.EDTA.original.TE.type.bed

## Convert TE types into Copia, Gypsy, Helitron, DNA/T, MITEs, and other TEs.
python original_type_to_5_types.py mo17.EDTA.original.TE.type.bed mo17.EDTA.5.TE.type.bed

sh combine.sh
```

In senario 2, TE region was not considered, and the Mo17 reference genome was annotated in the priority of CDS, UTR, intron, up5K, down5K, intergenic. The annotation results are in maize_*_uniq.
```
cd without_TE
sh combine.sh
```

Step 2. Assign genomic features to each SV according to the location of the SV. The location of SVs was determined as the most overlapped genomic feature performed by bedtools. The last step summarize the number of SVs in each genomic feature.
```
cd ..
sh run_annotation.sh
```