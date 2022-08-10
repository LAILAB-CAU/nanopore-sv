The calculation of insertion time for LTR-containing TE types and non-LTR types were different. For each LTR-containg TE, such as Copia and Gypsy, sequences for its paired LTRs (coordinates deposited in from .EDTA.intact.gff3) were extracted and were subject to pairwise alignment by MUSCLE. The alignment results was used to calculate the distance between paired LTRs with distmat functions in emboss. As a result, the insertion time for this LTR-containing TE was estimated as the distance divided by the molecular clock of 1.3-e8 per site per year. 

```
perl get_the_insertion_time_of_LTR.pl ../EDTA/Zm-Mo17-REFERENCE-CAU-1.0.fa.mod.EDTA.intact.gff3 Zm-Mo17-REFERENCE-CAU-1.0.fa Mo17.final.polish.fa.intact.LTR.time
```

For non-LTR TE types, such as DNA/T, Helitron, and MITEs, sequences for 11 subtypes (EDTA annotated DNA_DTA, DNA_DTC, DNA_DTH, DNA_DTM, DNA_DTT, DNA_Helitron, MITE_DTA, MITE_DTC, MITE_DTH, MITE_DTM, MITE_DTT) were first extracted in 11 fasta files.   

```
## Generate mo17.EDTA.5.TE.type.bed from EDTA output Zm-Mo17-REFERENCE-CAU-1.0.fa.mod.EDTA.TEanno.gff3
grep -v "#" EDTA/Zm-Mo17-REFERENCE-CAU-1.0.fa.mod.EDTA.TEanno.gff3 | awk -F"Classification=" '{print $1"\t"$2}' - | awk '{print $1"\t"$4"\t"$5"\t"$10}' -| awk -F";" '{print $1}' > mo17.EDTA.original.TE.type.bed
grep 'DNA\|MITE' mo17.EDTA.original.TE.type.bed > mo17.EDTA.original.DNA.type.bed 

perl ./get_the_sequences_of_DNA_transposons.pl mo17.EDTA.original.DNA.type.bed
```

Multiple sequence alignment (MSA) were conducted within each fasta file by blastn with default settings.
```
for te in DNA_DTA DNA_DTC DNA_DTH DNA_DTM DNA_DTT DNA_Helitron MITE_DTA MITE_DTC MITE_DTH MITE_DTM MITE_DTT
do
echo $te
seqkit rmdup mo17.EDTA.original.DNA.type.bed.$te.fasta > $te.rmdup.fa 
makeblastdb -in $te.rmdup.fa -dbtype nucl -out $te.blastdb -parse_seqids 
blastn -db $te.blastdb -query $te.rmdup.fa -out ./$te.blastn.DNA_DTA.m8 -outfmt 6 -num_threads 5
done
```

In the MSA result, the insertion time for aligned pairs were again calculated by MUSCLE  alignment and followed distmat functions in emboss. If a TE aligned with multiple other TEs in MSA, only the best hit was used to calculate insertion time. 

```
for te in DNA_DTA DNA_DTC DNA_DTH DNA_DTM DNA_DTT DNA_Helitron MITE_DTA MITE_DTC MITE_DTH MITE_DTM MITE_DTT
do
echo $te
perl ../get_the_emboss_results.pl mo17.EDTA.original.DNA.type.bed.$te.fasta $te.blastn.$te.m8 ./$te.insertion_time 
done
```
At this point, the insertion time for each Mo17 TE was calculated. To calculate the insertion time of deletions, deletions were annotated as a TE category if at least 80% of the deletion overlapped with this TE category. If the deletion overlapped with multiple TEs in the assigned category, the average insertion time among those TEs were used.

```
python from_insertion_time_to_bed.py Mo17.final.polish.fa.intact.LTR.time LTR.insert_time.bed
grep Gypsy LTR.insert_time.bed | awk '{print $1"\t"$2"\t"$3"\tGypsy\t"$5}' > Gypsy.insertion_time.bed
grep Copia LTR.insert_time.bed | awk '{print $1"\t"$2"\t"$3"\tCopia\t"$5}' > Copia.insertion_time.bed


for te in DNA_DTA DNA_DTC DNA_DTH DNA_DTM DNA_DTT DNA_Helitron MITE_DTA MITE_DTC MITE_DTH MITE_DTM MITE_DTT
do
    awk -F'_' 'NR>1{print $3"\t"$4"\t"$5}' $te.insertion_time | awk -v var=$te '{print $1"\t"$2"\t"$3"\t"var"\t"$5}' > $te.insertion_time.bed
done 

cat DNA_DTA.insertion_time.bed DNA_DTC.insertion_time.bed DNA_DTH.insertion_time.bed DNA_DTM.insertion_time.bed DNA_DTT.insertion_time.bed | awk '{print $1"\t"$2"\t"$3"\tDNA\t"$5}'> DNA.insertion_time.bed

cat MITE_DTA.insertion_time.bed MITE_DTC.insertion_time.bed MITE_DTH.insertion_time.bed MITE_DTM.insertion_time.bed MITE_DTT.insertion_time.bed | awk '{print $1"\t"$2"\t"$3"\tMITE\t"$5}'> MITE.insertion_time.bed

awk '{print $1"\t"$2"\t"$3"\tHelitron\t"$5}' DNA_Helitron.insertion_time.bed > Helitron.insertion_time.bed

cat MITE.insertion_time.bed DNA.insertion_time.bed Helitron.insertion_time.bed Gypsy.insertion_time.bed Copia.insertion_time.bed > all_5_te.insertion_time.bed

# For instance, DEL.bed was used to calculate the insertion time of deletions. 

bedtools intersect -a DEL.bed -b all_5_te.insertion_time.bed -wo > DEL.overlap.5.TE.insertion_time.bed
python sv_overlap_with_each_te_with_insertion_time.py DEL.overlap.5.TE.insertion_time.bed SV_TE_overlap.sum_len.insertion_time.bed 0.8
python assign_te_label_to_SV_insertion_time.py SV_TE_overlap.sum_len.insertion_time.bed SV.TE.insertion_time.bed
```


