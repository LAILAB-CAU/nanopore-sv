open CDS, "$ARGV[0]" or die "$!";
while(<CDS>)
	{
	chomp;
	@ss=split;
	for($i=$ss[1];$i<=$ss[2];$i++)
		{
		$hash_CDS{"$ss[0]\t$i"}++;
		}
	}
close CDS;

open UTR5, "$ARGV[1]" or die "$!";
while(<UTR5>)
        {
        chomp;
        @ss=split;
        for($i=$ss[1];$i<=$ss[2];$i++)
                {
                $hash_UTR5{"$ss[0]\t$i"}++;
                }
        }
close UTR5;
open UTR3, "$ARGV[2]" or die "$!";
while(<UTR3>)
	{
	chomp;
	@ss=split;
	for($i=$ss[1];$i<=$ss[2];$i++)
		{
		$hash_UTR3{"$ss[0]\t$i"}++;
		}
	}
close UTR3;

open INTRON, "$ARGV[3]" or die "$!";
while(<INTRON>)
        {
        chomp;
        @ss=split;
        for($i=$ss[1];$i<=$ss[2];$i++)
                {
                $hash_intron{"$ss[0]\t$i"}++;
                }
        }
close INTRON;

open UP2K, "$ARGV[4]" or die "$!" ;
while(<UP2K>)
        {
        chomp;
        @ss=split;
        for($i=$ss[1];$i<=$ss[2];$i++)
                {
                $hash_up2K{"$ss[0]\t$i"}++;
                }
        }

open DOWN2K, "$ARGV[5]" or die "$!";
while(<DOWN2K>)
        {
        chomp;
        @ss=split;
        for($i=$ss[1];$i<=$ss[2];$i++)
                {
                $hash_down2K{"$ss[0]\t$i"}++;
                }
        }

open I ,"$ARGV[6]" or die "$!";
while(<I>)
	{
	chomp;
	@ss=split;
	for($i=$ss[1];$i<=$ss[2];$i++)
		{
		if(exists $hash_CDS{"$ss[0]\t$i"})
			{
			print "$ss[0]\t$i\t$i\tcds\n";
			}
		elsif(exists $hash_UTR5{"$ss[0]\t$i"})
			{
			print "$ss[0]\t$i\t$i\t5utr\n";
			}
		elsif(exists $hash_UTR3{"$ss[0]\t$i"})
			{	
			print "$ss[0]\t$i\t$i\t3utr\n";
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

