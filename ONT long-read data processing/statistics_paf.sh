#!/bin/sh
paf=$1

awk 'ide=sprintf ("%.3f",$10/$11) {print $11"\t"ide}' ${paf} | perl -ane 'chomp; $i++; $t_aln+=$F[0];$t_ide+=$F[1];$ide1{$F[1]}++ if($F[1]<0.96);$ide2{$F[1]}++ if($F[1]>=0.96);$aln1{$F[0]}++ if($F[0]<8000);$aln2{$F[0]}++ if($F[0]>=8000); END {$m_ide=sprintf "%.3f", $t_ide/$i;$m_aln=sprintf "%.0f", $t_aln/$i; @ide1=sort {$b<=>$a} values %ide1; @ide2=sort {$b<=>$a} values %ide2; @aln1=sort {$b<=>$a} values %aln1; @aln2=sort {$b<=>$a} values %aln2; foreach $ide1(keys %ide1){if($ide1{$ide1}==$ide1[0]){$mode_ide1=$ide1}} foreach $ide2(keys %ide2){if($ide2{$ide2}==$ide2[0]){$mode_ide2=$ide2}} foreach $aln1(keys %aln1){if($aln1{$aln1}==$aln1[0]){$mode_aln1=$aln1}} foreach $aln2(keys %aln2){if($aln2{$aln2}==$aln2[0]){$mode_aln2=$aln2}} print "$m_ide\t$mode_ide1, $mode_ide2\t$m_aln\t$mode_aln1, $mode_aln2\n";}' > ${paf}.ide_aln

awk 'ide=sprintf ("%.3f",$10/$11) {print $11"\t"ide}' ${paf} > ${paf}.tmp

echo -e "
args=commandArgs(T)
setwd(getwd())
D=read.table(args[1], header=FALSE)
png(args[2],width=800,height=400)
par(mfrow=c(1,2))
plot(density(D\$V1),xlim=c(0,50000),main=args[1],xlab=\"Align length (bp)\")
plot(density(D\$V2),main=args[1],xlab=\"Identity\")
dev.off()
" > ${paf}.R

Rscript ${paf}.R ${paf}.tmp ${paf}.ide_aln.png

rm ${paf}.R ${paf}.tmp

awk 't=$3+1 {print $1"\t"t"\t"$4}' ${paf} | sort -k1,1 -k2,2n | /project/personal/yangm/02_software/01_package/bedtools2-2.27.1/bin/bedtools merge > ${paf}.qbed

perl -ane 'chomp; $diff=$F[2]-$F[1]+1; $hash{$F[0]}+=$diff; END{foreach $key(keys %hash){print "$key\t$hash{$key}\n"}}' ${paf}.qbed > ${paf}.qbed.cb

awk 'ARGIND==1{tmp[$1]=$2} ARGIND==2{a=sprintf("%.4f",$2/tmp[$1]);print a}' ${paf} ${paf}.qbed.cb | perl -ane 'chomp; $i++; $t+=$F[0]; $hash{$F[0]}++; END {$mean=sprintf "%.4f", $t/$i; @value=sort {$b<=>$a} values %hash; foreach $key(keys %hash){if ($hash{$key} == $value[0]){$mode=$key}} print "$mean\t$mode\n"}'> ${paf}.qbed.cb.q

awk '{print $6"\t"$8"\t"$9}' ${paf} | sort -k1,1 -k2,2n | /project/personal/yangm/02_software/01_package/bedtools2-2.27.1/bin/bedtools merge | awk 'BEGIN{t=0}{t+=$3-$2}END{a=sprintf ("%0.4f",t/2147377915); print a}' > ${paf}.gcov.txt 

perl -ane '$hash{$F[0]}++; $len=$F[3]-$F[2]; $ratio=$len/$F[1]; $i++ if ($ratio>=0.9);$j++ if ($ratio>=0.8);$k++ if ($ratio>=0.7); END {@total=keys %hash; $total=scalar @total; $r1=sprintf "%0.4f", $i/$total; $r2=sprintf "%0.4f", $j/$total; $r3=sprintf "%0.4f", $k/$total; print "$r1\t$r2\t$r3\n"}' ${paf} > ${paf}.qbed.bestalign.txt 

echo ${paf} > ${paf}.sample

perl -ane '$t_len+=$F[1]; $hash{$F[0]}++; END {@tnum=keys %hash; $tnum=scalar @tnum; print "$tnum\t$t_len\n"}' ${paf}.qbed.cb > ${paf}.treads

paste ${paf}.sample ${paf}.ide_aln ${paf}.qbed.cb.q ${paf}.gcov.txt ${paf}.qbed.bestalign.txt ${paf}.treads > ${paf}.results

rm ${paf}.sample ${paf}.qbed ${paf}.qbed.cb ${paf}.qbed.cb.q ${paf}.qbed.bestalign.txt ${paf}.gcov.txt ${paf}.ide_aln ${paf}.treads

