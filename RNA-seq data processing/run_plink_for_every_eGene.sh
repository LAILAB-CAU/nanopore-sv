#!/bin/sh

wanted_variation=$1  ## absolute path of SNPs and small InDels
## prepare LD files
Variation=./
plink --bfile $Variation --extract $wanted_variation --make-bed --allow-extra-chr --out $wanted_variation.plink_format
plink --bfile $wanted_variation.plink_format --r --ld-window 1000 --out $wanted_variation.plink_format
