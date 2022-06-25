Two steps to annotate where each SV is located in the genome:

Step 1. The partition of Mo17 reference genome: a custom script was used to annotate genome features from the downloaded gff3 file at https://download.maizegdb.org/Zm-Mo17-REFERENCE-CAU-1.0/. The partition of different features had the following priority: TE region, exon, intron, upstream 2 kbp region, downstream 2 kbp region, and intergenic region. Two adjacent regions with the same annotation were merged into one region. 
```
# Download Zm-Mo17-REFERENCE-CAU-1.0_mgdb.gff3 file from https://download.maizegdb.org/Zm-Mo17-REFERENCE-CAU-1.0/
sh run_step1.sh
```

Step 2. Assign genomic features to each SV according to the location of the SV. The location of SVs was determined as the most overlapped genomic feature performed by bedtools.
```
sh run_step2.sh
```