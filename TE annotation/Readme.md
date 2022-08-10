For the Mo17 reference genome, both structurally intact and fragmented transposable elements were annotated with EDTA (v 2.0.1). The non-redundant pan-genome TE library for NAM founder maize genomes18 (located at https://github.com/HuffordLab/NAM-genomes/tree/master/te-annotation/assets/NAM.EDTA1.8.0.MTEC02052020.TElib.clean.fa) was used as the base library. Novel TE families that were single-copy in Mo17 genome were identified by RepeatMasker and would be removed. 
```
lib=NAM.EDTA1.8.0.MTEC02052020.TElib.clean.fa
genome=Zm-Mo17-REFERENCE-CAU-1.0.fa
RepeatMasker -pa 36 -q -div 40 -lib $lib -cutoff 225 -gff $genome
```
The genic sequences were removed in the TE annotation (--cds Zm-Mo17-REFERENCE-CAU-1.0_Zm00014a.1.cds.fa). The species parameter was set to Maize (--species Maize). 

```
lib=NAM.EDTA1.8.0.MTEC02052020.TElib.clean.fa
genome=Zm-Mo17-REFERENCE-CAU-1.0.fa
cds=Zm-Mo17-REFERENCE-CAU-1.0_Zm00014a.1.cds.fa
threads=40

type=helitron
perl /public1/home/scb8190/softwares/EDTA/EDTA_raw.pl --genome $genome --species Maize --type $type -t $threads

perl /public1/home/scb8190/softwares/EDTA/EDTA.pl \
   --genome $genome \
   --species Maize \
   -t $threads \
   --anno 1 \
   --rmout $genome.out \
   --curatedlib $lib \
   --cds $cds
```
TE families were further re-grouped into six categories: Gypsy, Copia, class II DNA transposons (DNA/T), class II DNA transposon Helitron, MITEs, and Other TEs. Regions not covered by TEs were categorized as non-TE. 
```
grep -v "#" Zm-Mo17-REFERENCE-CAU-1.0.fa.mod.EDTA.TEanno.gff3 | awk -F"Classification=" '{print $1"\t"$2}' - | awk '{print $1"\t"$4"\t"$5"\t"$10}' -| awk -F";" '{print $1}' > mo17.EDTA.original.TE.type.bed

python original_type_to_5_types.py mo17.EDTA.original.TE.type.bed mo17.EDTA.5.TE.type.bed
```
With the genome-wide annotation of Mo17 genome, deletions were annotated as a TE category if at least 80% of the deletion overlapped with this TE category. DEL.bed is the bed file for all deletions. 
```
bedtools intersect -a DEL.bed -b mo17.EDTA.5.TE.type.bed -wo > DEL.overlap.5.TE.bed

python sv_overlap_with_each_te.py DEL.overlap.5.TE.bed SV_TE_overlap.sum_len.bed 0.8
python assign_te_label_to_SV.py SV_TE_overlap.sum_len.bed SV.TE.bed

awk 'BEGIN{FS=OFS="\t"} {if(NR==FNR){array[$1"\t"$2"\t"$3]=$4} else{if($0 in array) {print $0"\t"array[$0]} else {print $0"\tnon-TE"}}}' SV.TE.bed DEL.bed > DEL.TE.all.bed
```

For the annotation of insertions, the sequences of insertions were first extractd and TE families of them were identified by RepeatMasker with parameters “-q -div 40 -cutoff 225” and with pan-genome TE library. If an insertion was annotated with multiple categories, the longest category that covering at least 80% of the insertion was assigned.

```
perl extract_seq.pl INS.SViper.out.1000.vcf.uniqe.nomissing.vcf INS.fa INS_id.bed 2

awk 'BEGIN{FS=OFS="\t"} {if(NR==FNR){array[$1]=$2} else{if($5 in array) {print array[$5]"\t"$7-$6"\t"$11}}}' INS_id.bed INS.fa.out > INS.TE.temp

python prepare_INS_TE.py EDTA/INS_id.bed EDTA/INS.fa.out INS.overlap.5.TE.bed

python sv_overlap_with_each_te.py INS.overlap.5.TE.bed INS_TE_overlap.sum_len.bed 0.8
python assign_te_label_to_SV.py INS_TE_overlap.sum_len.bed INS.TE.bed

awk 'BEGIN{FS=OFS="\t"} {if(NR==FNR){array[$1"\t"$2"\t"$3]=$4} else{if($0 in array) {print $0"\t"array[$0]} else {print $0"\tnon-TE"}}}' INS.TE.bed INS.bed > INS.TE.all.bed
```
