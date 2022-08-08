In this project, the genomic features of SVs were annotated in two senarios. In one senario, TE feature was considered and the Mo17 reference genome was first annotated in the priority of TE region, exon, intron, upstream 5kbp region (Up5k), downstream 5 kbp region (Down5k), and intergenic region. In the other senario, TE region was not considered, and the Mo17 reference genome was annotated in the priority of . In both senarios, two adjacent regions with the same annotation were merged into one region. SVs were annotated through two steps. 

Step 1. The partition of Mo17 reference genome: a custom script was used to annotate genome features from the downloaded gff3 file at https://download.maizegdb.org/Zm-Mo17-REFERENCE-CAU-1.0/Zm-Mo17-REFERENCE-CAU-1.0_Zm00014a.1.gff3.gz. 

In senario 1, different features had the following priority: TE region, exon, intron, upstream 5kbp region (Up5k), downstream 5 kbp region (Down5k), and intergenic region. TE annotation (file) was from the 'TE annotation' section, which conducted EDTA and annotated the Mo17 genome into 5 TE categories: Copia, Gypsy, Helitron, DNA/T, MITEs, and other TEs. 

```
sh run_step1_with_TE.sh
```

Step 2. Assign genomic features to each SV according to the location of the SV. The location of SVs was determined as the most overlapped genomic feature performed by bedtools.
```
sh run_step2.sh
```