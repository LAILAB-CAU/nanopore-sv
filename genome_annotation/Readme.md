共分为两步


A: 根据gff注释文件划分基因组

     1: 获得全基因组exon与基因的坐标
   
     2：根据基因的坐标得到基因上下游2k的坐标
    
     3：合并所有exon得到merge后的 exon坐标
   
     4：根据合并后的exon得到intron的坐标
   
     5: 根据EDTA注释的TE得到TE的坐标
   
     6：按照TE>exon>intron>up2k>down2K>intergenic的优先级得到基因组上每个位点所属于的特征
   
     7：合并相邻的相同特征得到基因组划分结果
   
   

B：根据SV与基因组各个特征的overlap确定SV所在的区域

     1：将vcf格式的SV文件转化为只包含SV位置的bed 文件
   
     2：将bed文件和A种的基因组划分文件做overlap
   
     3: 根据SV与基因组各个特征的overlap的程度确定SV所在的区域
   
