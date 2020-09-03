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


open I ,"$ARGV[1]" or die "$!";
while(<I>)
	{
	chomp;
	@ss=split;
	for($i=$ss[1];$i<=$ss[2];$i++)
		{
		if(exists $hash_TE{"$ss[0]\t$i"})
			{
				next;
			}
		else
			{
			print "$ss[0]\t$i\t$i\tNO_TE\n";
			}
		}
	}

