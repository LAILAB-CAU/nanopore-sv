#!/bin/sh

wanted_variation=$1  ## absolute path of SNPs and small InDels
zscore=$2

## CAVIAR was from https://github.com/fhormoz/caviar
echo RUN_CAVIAR
source /public1/soft/modules/module.sh
module load gcc/7.3.0-wzm
module load gsl/2.5-cjj
/public1/home/sc30797/bxin/softwares/caviar-master/CAVIAR-C++/CAVIAR -l $wanted_variation.plink_format.ld.2d -z $2 -c 1 -o $wanted_variation.caviar.out
