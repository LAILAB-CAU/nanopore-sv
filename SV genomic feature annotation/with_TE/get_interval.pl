open TE, "$ARGV[0]" or die "$!";
while(<TE>)
	{	
	chomp;
	@ss=split;
	for($i=$ss[1];$i<=$ss[2];$i++)
		{
		$hash_TE{"$ss[0]\t$i"}++;
		}
	}

open EXON, "$ARGV[1]" or die "$!";
while(<EXON>)
	{
	chomp;
	@ss=split;
	for($i=$ss[1];$i<=$ss[2];$i++)
		{
		$hash_exon{"$ss[0]\t$i"}++;
		}
	}

open INTRON, "$ARGV[2]" or die "$!";
while(<INTRON>)
        {
        chomp;
        @ss=split;
        for($i=$ss[1];$i<=$ss[2];$i++)
                {
                $hash_intron{"$ss[0]\t$i"}++;
                }
        }

open UP2K, "$ARGV[3]" or die "$!" ;
while(<UP2K>)
        {
        chomp;
        @ss=split;
        for($i=$ss[1];$i<=$ss[2];$i++)
                {
                $hash_up2K{"$ss[0]\t$i"}++;
                }
        }

open DOWN2K, "$ARGV[4]" or die "$!";
while(<DOWN2K>)
        {
        chomp;
        @ss=split;
        for($i=$ss[1];$i<=$ss[2];$i++)
                {
                $hash_down2K{"$ss[0]\t$i"}++;
                }
        }

open I ,"$ARGV[5]" or die "$!";
while(<I>)
	{
	chomp;
	@ss=split;
	for($i=$ss[1];$i<=$ss[2];$i++)
		{
		if(exists $hash_TE{"$ss[0]\t$i"})
			{
			print "$ss[0]\t$i\t$i\tTE\n";
			}
		elsif(exists $hash_exon{"$ss[0]\t$i"})
			{
			print "$ss[0]\t$i\t$i\texon\n";
			}
		elsif(exists $hash_intron{"$ss[0]\t$i"})
			{
			print "$ss[0]\t$i\t$i\tintron\n";
			}
		elsif(exists $hash_up2K{"$ss[0]\t$i"})
			{
			print "$ss[0]\t$i\t$i\tup2K\n";
			}
		elsif(exists $hash_down2K{"$ss[0]\t$i"})
			{
			print "$ss[0]\t$i\t$i\tdown2K\n";
			}
		else
			{
			print "$ss[0]\t$i\t$i\tintergenic\n";
			}
		}
	}

