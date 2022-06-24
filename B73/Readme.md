Four methods has been used to detect SVs in B73: Assembly-based method AnchorWave (https://github.com/baoxingsong/anchorwave), ONT-based, Pacbio-based, and short-read sequencing (SRS)-based methods. 

# 1. AnchorWave

 B73_v4 (downloaded from https://maizegdb.org/genome/assembly/Zm-B73-REFERENCE-GRAMENE-4.0) and Mo17_CAU (downloaded from https://maizegdb.org/genome/assembly/Zm-Mo17-REFERENCE-CAU-1.0) were collected. 

Step1. Extract the CDS sequence 

```
anchorwave gff2seq -r ../Zm-Mo17-REFERENCE-CAU-1.0.fa -i ../Zm-Mo17-REFERENCE-CAU-1.0_Zm00014a.1.gff3  -o Zm-Mo17_cds.fa
```

Step2: Align CDS to the ref and query genome respectively
```
minimap2 -x splice -t 6 -k 12 -a -p 0.4 -N 20 ../Zm-B73-REFERENCE-GRAMENE-4.0.fa Zm-Mo17_cds.fa > B73_to_Mo17cds.sam
minimap2 -x splice -t 6 -k 12 -a -p 0.4 -N 20 ../Zm-Mo17-REFERENCE-CAU-1.0.fa Zm-Mo17_cds.fa > Mo17_to_Mo17ref.sam
```

Step3: Genome alignment by anchorWave

```
anchorwave genoAli -i ../Zm-Mo17-REFERENCE-CAU-1.0_Zm00014a.1.gff3 -as Zm-Mo17_cds.fa -r ../Zm-Mo17-REFERENCE-CAU-1.0.fa -a B73_to_Mo17cds.sam -ar Mo17_to_Mo17ref.sam -s ../Zm-B73-REFERENCE-GRAMENE-4.0.fa -n B73_Mo17_anchors -o B73_to_Mo17anchorwave.maf -f B73_toMo17anchorwave.f.maf -w 38000 -fa3 200000 -B -6 -O1 -8 -E1 -2 -O2 -75 -E2 -1 -t 12 -IV >B73_Mo17_anchorwave.log 2>&1
```

Step4: Call SV from *.maf file, follow the pipeline described in https://bitbucket.org/bucklerlab/practicalhaplotypegraph/wiki/UserInstructions/CreatePHG_step2_MAFToGVCFPluginDetails.md

```
docker pull maizegenetics/phg
docker pull maizegenetics/phg_liquibase

WORKING_DIR=/data/bxin/nanopore_sv/maf_to_vcf
DOCKER_CONFIG_FILE=/phg/config.txt
cd $WORKING_DIR

echo "MAFToGVCFPlugin.referenceFasta=/phg/mo17.fa" > config.txt
echo "MAFToGVCFPlugin.mafFile=/phg/B73_to_Mo17anchorwave.maf" >>config.txt
echo "MAFToGVCFPlugin.gvcfOutput=/phg/B73_to_Mo17anchorwave.gvcf" >>config.txt 
echo "MAFToGVCFPlugin.sampleName=mo17_anchorwave" >> config.txt
echo "MAFToGVCFPlugin.fillGaps=false" >> config.txt
echo "MAFToGVCFPlugin.twoGvcfs=false" >> config.txt

# The the config file must contain the necessary MAFToGVCFPLugin parameters
docker run --name anchorwave_assemblies_maf --rm \
    -v ${WORKING_DIR}/:/phg/ \
    -t maizegenetics/phg:latest \
    /tassel-5-standalone/run_pipeline.pl -Xmx100G -debug -configParameters ${DOCKER_CONFIG_FILE} \
    -MAFToGVCFPlugin -endPlugin

## Seperate DEL and INS
grep -v "#" B73_to_Mo17anchorwave.gvcf | awk 'length($4)>=50||length($5)>=60{print $_}' - > B73_to_Mo17anchorwave.gvcf.temp
python3 convert_gvcf_to_sv.py B73_to_Mo17anchorwave.gvcf.temp  
rm B73_to_Mo17anchorwave.gvcf.temp

## Remove SVs with N in the sequecne 
grep -v N B73_to_Mo17anchorwave.gvcf.del > anchorwave_DEL.bed 
grep -v N B73_to_Mo17anchorwave.gvcf.ins > anchorwave_INS.bed 

## convert bed file to vcf file with SURVIVOR 
SURVIVOR bedtovcf anchorwave_DEL.bed DEL anchorwave_DEL.vcf
SURVIVOR bedtovcf anchorwave_INS.bed INS anchorwave_INS.vcf
```
