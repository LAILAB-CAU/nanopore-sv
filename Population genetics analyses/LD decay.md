## LD decay
 ```
 Input file: 
 Output file: Average_output_DEL.out.1000.vcf.unique.chromosome.vcftools.modify.ld.txt
    Distance(bp)	Average_R2
    0	0.032718450507455
    100	0.0348902858627838
    200	0.0549190959311831
 ```
 - Covert .vcf file to ped and map format 
     ```
     vcftools --vcf *.vcf --plink --out *
     ```
 - Modify .mapw file so that it can be recognized by plink1.07
 - LD calculation, obtain the pairwise LD between SNPs
    ```
    plink --file * --r2 --out * --ld-window 10000 --ld-window-kb 200 --ld-window-r2 0 --map3 --noweb
    ```
    Check plink help page at http://zzz.bwh.harvard.edu/plink/ld.shtml
 - Take an average of LD values of SNPs within 100bp