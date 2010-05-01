#!/usr/bin/perl -w
use strict;
use Getopt::Long;

my ($genome_file, $gff_file, $new_gff_file);
my ($freq_start_file, $freq_start1_file);
my ($len, $count, $cg);
my ($genome_seq, $head, $temp_genome_seq);
my ($id, $i, $j);
my ($utr, $out, $e_save, $i_save, $codon, $command, $return_seq);
my @freq;
my @freq1;

GetOptions(
           'genome=s' => \$genome_file,
           'pre=s' => \$gff_file,
           'post=s' => \$new_gff_file,
           );



my %freq_unit =('AAA',0,'AAC',1,'AAG',2,'AAT',3,'ACA',4,'ACC',5,'ACG',6,'ACT',7,
		'AGA',8,'AGC',9,'AGG',10,'AGT',11,'ATA',12,'ATC',13,'ATG',14,'ATT',15,
		'CAA',16,'CAC',17,'CAG',18,'CAT',19,'CCA',20,'CCC',21,'CCG',22,'CCT',23,
		'CGA',24,'CGC',25,'CGG',26,'CGT',27,'CTA',28,'CTC',29,'CTG',30,'CTT',31,
		'GAA',32,'GAC',33,'GAG',34,'GAT',35,'GCA',36,'GCC',37,'GCG',38,'GCT',39,
		'GGA',40,'GGC',41,'GGG',42,'GGT',43,'GTA',44,'GTC',45,'GTG',46,'GTT',47,
		'TAA',48,'TAC',49,'TAG',50,'TAT',51,'TCA',52,'TCC',53,'TCG',54,'TCT',55,
		'TGA',56,'TGC',57,'TGG',58,'TGT',59,'TTA',60,'TTC',61,'TTG',62,'TTT',63);

my $program = $0;
my $dir = substr($0, 0, length($0)-15);

get_sequence($genome_file, \$genome_seq, \$head);
$temp_genome_seq = $genome_seq;
$len = length($genome_seq);
$count = ($temp_genome_seq =~ s/C|G|c|g//g);
$temp_genome_seq="";
$cg = sprintf("%.0f", $count*100/$len);
if ($cg < 26){
    $cg=26;
}elsif($cg>69){
    $cg=69;
}

$freq_start_file = $dir."train/post_start/".$cg;
$freq_start1_file = $dir."train/post_start_r/".$cg;

# get frequency model                                                                                            
$id=0;
open(IN, $freq_start_file);
while(my $each_line=<IN>){
    chomp($each_line);
    if ($each_line!~/^\>/){
	my @temp = split(/\s+/, $each_line);
	for ($i=0; $i<=63; $i++){
	    $freq[$id][$i] = $temp[$i];
	}
	$id+=1;
    }
}
close(IN);

$id=0;
open(IN, $freq_start1_file);
while(my $each_line=<IN>){
    chomp($each_line);
    if ($each_line!~/^\>/){
	my @temp = split(/\s+/, $each_line);
	for ($i=0; $i<=63; $i++){
	    $freq1[$id][$i] = $temp[$i];
	}
	$id+=1;
    }
}
close(IN);


# find the optimal start codon with 30bp up- and downstream of start codon
open OUT, ">$new_gff_file";
open(IN, $gff_file);
while(my $each_line=<IN>){
    
    chomp($each_line);
    my @temp = split(/\s+/, $each_line);

    if ($each_line =~ /^\>/){
	print OUT $each_line."\n";
    }elsif ($temp[2] eq "+"){
	my $i=0;
	$codon = substr($genome_seq, $temp[0]-1, 3);

	while ($codon !~ /TAA|TAG|TGA/ &&$temp[0]-1-$i-35>=0  ){
	    
	    if ($codon =~ /ATG|GTG|TTG/){
		$utr = substr($genome_seq, $temp[0]-1-$i-30, 63);
	     
		my @temp1 = split(//, $utr);
		my $freq_sum = 0;
		if ($#temp1==62){
		    for($j=0; $j<=60; $j++){
			my $key = $temp1[$j].$temp1[$j+1].$temp1[$j+2];
			if (exists $freq_unit{$key}){
			    $freq_sum -= log($freq[$j][$freq_unit{$key}]);
			}else{
			    $freq_sum -= log(1/64);
			}
		    }
		}
		if ($i==0){
		    $e_save = $freq_sum;
		    $i_save = 0;
		}elsif ($freq_sum < $e_save){
		    $e_save = $freq_sum;
		    $i_save = -1*$i;
		}
	    }
	    $i += 3;
	    $codon = substr($genome_seq, $temp[0]-1-$i, 3);
	}
	print OUT eval($temp[0]+$i_save)."\t".$temp[1]."\t".$temp[2]."\t".$temp[3]."\t".$temp[4]."\t".$temp[5]."\t".$temp[6]."\n";
	
    }elsif ($temp[2] eq "-"){
	my $i=0;
	$codon = substr($genome_seq, $temp[1]-1-2, 3);
	while ($codon !~ /TTA|CTA|TCA/ && $temp[1]-2+$i+35<length($genome_seq) ){
	    
	    if ($codon =~ /CAT|CAC|CAA/){
		$utr = substr($genome_seq, $temp[1]-1-2+$i-30, 63);
		my @temp1 = split(//, $utr);
		my $freq_sum = 0;
		if ($#temp1==62){
		    for($j=0; $j<=60; $j++){
			my $key = $temp1[$j].$temp1[$j+1].$temp1[$j+2];
			if (exists $freq_unit{$key}){
			    $freq_sum -= log($freq1[$j][$freq_unit{$key}]);
			}else{
			    $freq_sum -= log(1/64);
			}
		    }
		}
		if ($i==0){
		    $e_save = $freq_sum;
		    $i_save = 0;
		}elsif ($freq_sum < $e_save){
		    $e_save = $freq_sum;
		    $i_save = $i;
		}
	    }
	    $i += 3;
	    $codon = substr($genome_seq, $temp[1]-1-2+$i, 3);
	}
	print OUT $temp[0]."\t".eval($temp[1]+$i_save)."\t".$temp[2]."\t".$temp[3]."\t".$temp[4]."\t".$temp[5]."\t".$temp[6]."\n";
    }
}
close(OUT);
close(IN);

sub get_sequence{  # file name, @variable for seq, @variable for head

    open(GENOME, $_[0])|| die("ERROR: Couldn't open genome_file $_[0]!\n");
    while( my $each_line=<GENOME>)  {

        if ($each_line =~ m/>/){
            ${$_[1]} = "";
            chomp($each_line);
            ${$_[2]} = $each_line;
        }else{
            chomp($each_line);
            ${$_[1]} .= $each_line;
        }
    }
    close(GENOME);
}



