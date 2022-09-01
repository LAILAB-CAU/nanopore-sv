#!/usr/bin/perl -w
use strict;
use warnings;

#input_file1 and input_file2 are VCF files for insertion and deletions.
#output_file1 and output_file2 are output VCF files
die "perl **.pl   <input_file1> <input_file2> <out_file1> <out_file2>" if (@ARGV != 4 );

my ($input_file1,$input_file2,$out_file1,$out_file2) = @ARGV;
my @input_files = ($input_file1,$input_file2);
my @depth_files = </public1/home/sc30797/shijunpeng/depth/split_bams/*.gz>;
my @out_files = ($out_file1,$out_file2);


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
foreach my $input_file(@input_files){
    open VCF, "$input_file" or die;
    while(my $line = <VCF>){
        chomp $line;
        
        next if $line =~ /^##/;

        my @array = split /\s+/, $line;
        
        if ($line =~ /^#CHROM/){
            if($jvcf_header eq ''){
                $line =~ s/^#//;
                
                for(my $i=0;$i<=8;$i++){
                    $jvcf_header .= "$array[$i]\t";
                }

                for(my $i = 9; $i <= $#array; $i++){
                    my @temp = split /\//, $array[$i];
                    my @sample = split /\./, $temp[-1];
                    $jvcf_header .= "$sample[0]\t";
                }
            }
            next;
        }
        
        

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
        

        else{
            delete $hash_pos{"$array[0]\t$array[1]"};
        }
    }
    close VCF;
}

my @header_array = split /\s+/, $jvcf_header;


my @sample_array;
for(my $cc =0;$cc <=$#input_files;$cc++){
    my $input_file = $input_files[$cc];
    my $out_file = $out_files[$cc];
    open VCF, "$input_file" or die;
    open UVCF, ">$out_file" or die;
    while(my $line = <VCF>){
        chomp $line;
        
        if($line =~ /^##/){
            print UVCF "$line\n";
            next;
        }
        my @array = split /\s+/, $line;
        
        
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

        if( $line =~ /SVLEN=-?(\d+);/){
            next if ($1 <50);
        }
        
        if(exists $hash_pos{"$array[0]\t$array[1]"}){
            for(my $i = 0; $i <= 8; $i++){
                print UVCF "$array[$i]\t";
            }
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
