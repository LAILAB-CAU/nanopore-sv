open EXON, "$ARGV[0]" or die "$!";
while(<EXON>)
	{
	chomp;
	@ss=split;
	for($i=$ss[1];$i<=$ss[2];$i++)
		{
		$hash_exon{"$ss[0]\t$i"}++;
		}
	}

open INTRON, "$ARGV[1]" or die "$!";
while(<INTRON>)
        {
        chomp;
        @ss=split;
        for($i=$ss[1];$i<=$ss[2];$i++)
                {
                $hash_intron{"$ss[0]\t$i"}++;
                }
        }

open UP5K, "$ARGV[2]" or die "$!" ;
while(<UP5K>)
        {
        chomp;
        @ss=split;
        for($i=$ss[1];$i<=$ss[2];$i++)
                {
                $hash_up5K{"$ss[0]\t$i"}++;
                }
        }

open DOWN5K, "$ARGV[3]" or die "$!";
while(<DOWN5K>)
        {
        chomp;
        @ss=split;
        for($i=$ss[1];$i<=$ss[2];$i++)
                {
                $hash_down5K{"$ss[0]\t$i"}++;
                }
        }

open I ,"$ARGV[4]" or die "$!";
while(<I>)
	{
	chomp;
	@ss=split;
	for($i=$ss[1];$i<=$ss[2];$i++)
		{
		if(exists $hash_exon{"$ss[0]\t$i"})
			{
			print "$ss[0]\t$i\t$i\texon\n";
			}
		elsif(exists $hash_intron{"$ss[0]\t$i"})
			{
			print "$ss[0]\t$i\t$i\tintron\n";
			}
		elsif(exists $hash_up5K{"$ss[0]\t$i"})
			{
			print "$ss[0]\t$i\t$i\tup5K\n";
			}
		elsif(exists $hash_down5K{"$ss[0]\t$i"})
			{
			print "$ss[0]\t$i\t$i\tdown5K\n";
			}
		else
			{
			print "$ss[0]\t$i\t$i\tintergenic\n";
			}
		}
	}

