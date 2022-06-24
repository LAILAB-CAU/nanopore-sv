#!/usr/bin/perl -w
use strict;
use warnings;

# usage:
die "perl **.pl   <input_file1> <input_file2> <out_file1> <out_file2>" if (@ARGV != 4 );

my ($input_file1,$input_file2,$out_file1,$out_file2) = @ARGV;
my @input_files = ($input_file1,$input_file2);
my @depth_files = <ABSOLUTE_PATH_OF_DEPTH_FILES/*.gz>;
my @out_files = ($out_file1,$out_file2);

#hash_depth stores the depth of each SV start position at each accession. It has two keys, chr + start pos and sample ID, its value is Y/N. SVs with depth >=2 will be assinged 0/0, otherwise ./. missing. 
my %hash_depth;
foreach my $depth_file(@depth_files){
    open DEPTH, "gzip -dc $depth_file |" or die;
    
     
    my @sample_depth_file_name = split /\//,$depth_file;
    my @sample_name = split /\./, $sample_depth_file_name[-1];
    my $depth_sample = $sample_name[0];
    

    while(my $line = <DEPTH>){
        chomp $line;
        my @array = split /\s+/, $line;

        if($array[2]>=2){
            $hash_depth{"$array[0]\t$array[1]"}{$depth_sample} = "Y";
		$hash_depth{"$array[0]\t$array[1]"}{"C_72"}= "Y";
		$hash_depth{"$array[0]\t$array[1]"}{"Z58"}= "Y";
        }
        else{
            $hash_depth{"$array[0]\t$array[1]"}{$depth_sample} = "N";
        }

    }
    close DEPTH;
}

my $jvcf_header = '';
my %hash_pos;
# Load vcf for INS and DEL, create hash with key as chr + pos, value is an array storing genotypes in all accessions
# For coding convenience, A represents DEL,C represents INS, T represent 0/0, X represents ./.
foreach my $input_file(@input_files){
    open VCF, "$input_file" or die;
    while(my $line = <VCF>){
        chomp $line;
        
        # skip # lines
        next if $line =~ /^##/;

        my @array = split /\s+/, $line;
        
        # store header to $jvcf_header
        if ($line =~ /^#CHROM/){
            if($jvcf_header eq ''){
                $line =~ s/^#//;
                
                for(my $i=0;$i<=8;$i++){
                    $jvcf_header .= "$array[$i]\t";
                }

                #Genotype in the vcf file starts from the 10th column
                for(my $i = 9; $i <= $#array; $i++){
                    # Use sample ID only
                    my @temp = split /\//, $array[$i];
                    my @sample = split /\./, $temp[-1];
                    $jvcf_header .= "$sample[0]\t";
                }
            }
            next;
        }

        # If the SV is new:
        #key:chr  pos
        #value:[GT1,GT2,...]
        if(!exists $hash_pos{"$array[0]\t$array[1]"}){
            
            my @sample_array = split /\s+/, $jvcf_header;

            for(my $i = 9; $i <= $#array; $i++){
                my @temp = split /:/,$array[$i];
                
                if($array[2]=~/DEL/ || $array[4]=~/DEL/){
                    if($temp[0] eq '1/1'){
                        $hash_pos{"$array[0]\t$array[1]"}{$sample_array[$i]} = "A";
                    }
                    else{
                        if((exists $hash_depth{"$array[0]\t$array[1]"}{$sample_array[$i]}) && ($hash_depth{"$array[0]\t$array[1]"}{$sample_array[$i]} eq "Y")){
                            $hash_pos{"$array[0]\t$array[1]"}{$sample_array[$i]} = "T";
                        }
                        else{
                            $hash_pos{"$array[0]\t$array[1]"}{$sample_array[$i]} = "X";
                        }
                    }
                }
                if($array[2]=~/INS/ || $array[4]=~/INS/){
                    if($temp[0] eq '1/1'){
                        $hash_pos{"$array[0]\t$array[1]"}{$sample_array[$i]} = "C";
                    }
                    else{
                        if((exists $hash_depth{"$array[0]\t$array[1]"}{$sample_array[$i]}) && ($hash_depth{"$array[0]\t$array[1]"}{$sample_array[$i]} eq "Y")){
                            $hash_pos{"$array[0]\t$array[1]"}{$sample_array[$i]} = "T";
                        }
                        else{
                            $hash_pos{"$array[0]\t$array[1]"}{$sample_array[$i]} = "X";
                        }
                    }
                }
            }
        }
        

        # If the SV has been stored, discard this SV. 
        else{
            delete $hash_pos{"$array[0]\t$array[1]"};
        }
    }
    close VCF;
}

my @header_array = split /\s+/, $jvcf_header;

# If a SV occurs multiple times, remove duplicates
my @sample_array;
for(my $cc =0;$cc <=$#input_files;$cc++){
    my $input_file = $input_files[$cc];
    my $out_file = $out_files[$cc];
    open VCF, "$input_file" or die;
    open UVCF, ">$out_file" or die;
    while(my $line = <VCF>){
        chomp $line;
        
        # Keep using the same annotation lines
        if($line =~ /^##/){
            print UVCF "$line\n";
            next;
        }
        my @array = split /\s+/, $line;
        
        # Use sample ID only for each accession
        
        if($line =~ /^#CHROM/){
            for(my $i = 0; $i <= 8; $i++){
                print UVCF "$array[$i]\t";
            }
            for(my $i = 9; $i <= $#array; $i++){
                my @temp = split /\//, $array[$i];
                my @sample = split /\./, $temp[-1];
                print UVCF "$sample[0]\t";
            }
            print UVCF "\n";
            next;
        }

        # keep SVs with length larger than 50bp
        if( $line =~ /SVLEN=-?(\d+);/){
            next if ($1 <50);
        }
        
        # Write information for all unique SVs
        if(exists $hash_pos{"$array[0]\t$array[1]"}){
            for(my $i = 0; $i <= 8; $i++){
                print UVCF "$array[$i]\t";
            }
            # Change ./. to 0/0 if read depth >= 2 in a certain accession
            for(my $i = 9; $i <= $#array; $i++){
                my @temp = split /:/,$array[$i];
                if($hash_depth{"$array[0]\t$array[1]"}{$header_array[$i]} eq "Y"){
                    if($temp[0] eq "1/1"){
                        print UVCF "$temp[0]\t";
                    }
                    else{
                        print UVCF "0/0\t";
                    }
                }
                else{
                    print UVCF "./.\t";
                }
            }
            print UVCF "\n";
        }
    }
    close VCF;
    close UVCF;
}