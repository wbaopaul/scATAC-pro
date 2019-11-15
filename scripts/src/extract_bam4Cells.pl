#Run as follows:
#perl extract_bam4bcs.pl --cellbarcode_file cellbarcode_file --bam_file bam_file --output_dir output_dir --samtools_path SAMTOOLS_PATH
#bam file for given barcodes will be saved in output_dir.bam
#
#
#example:
#cellbarcode_file format:
#AAACGAAAGAAATGGG-1	
#AAACGAAAGTCTCTAG-1	
#AAACGAACAACTAGAA-1	
#AAACGAACAATCCATG-1	


use strict;


#Receive options from command line
use Getopt::Long;

GetOptions( 'cellbarcode_file=s' => \my $cellbarcode_file 
          , 'bam_file=s' => \my $bam_file  
          , 'output_dir=s' => \my $output_dir  
          , 'samtools_path=s' => \my $samtools_path  
          );



open(BARCODE, $cellbarcode_file ) or die("Cannot read $cellbarcode_file \n");
open(BAM, "$samtools_path/samtools view -h $bam_file |" ) or die("Cannot read $bam_file \n");

my $output_sam=$output_dir."/cell_barcodes.sam";
my $output_sam1=$output_dir."/non_cell_barcodes.sam";
open(OUT, ">$output_sam" ) or die("Cannot write $output_sam \n");
open(OUT1, ">$output_sam1" ) or die("Cannot write $output_sam1 \n");

my %cell_barcodes = ();

#my $header = <BARCODE>;

print("I am reading barcode information file: $cellbarcode_file \n");

my $pair_counter = 0;
while(my $line = <BARCODE>)
{
   $pair_counter++;
   chomp $line;
   my @array = split /\t/, $line;
   my $barcode = $array[0];
   $cell_barcodes{$barcode} = 1;
}

print "I have successfully read $pair_counter cell-barcode from $cellbarcode_file \n";

print "...\n";



print("I am reading BAM file: $bam_file \n");
print "...\n";
my $bam_line_counter = 0;
my $barcode_counter = 0;
my $matching_line_counter = 0;

my $btime = time;
my $barcode="XXXX";
while(my $line = <BAM>)
{
   $bam_line_counter++;
   if ($line =~ /^\@/){
        print OUT $line; ## pring header
        print OUT1 $line; ## pring header
    }
   chomp $line;
   if($bam_line_counter % 1000000 == 0)
   {
      my $etime = time;
      my $elapsed = $etime - $btime;
      my $mil =  $bam_line_counter / 1000000;
      print("I have processed $mil million reads in $elapsed seconds.\n");
   }
   $barcode_counter++;
   
   my @array0 = split /\t/, $line;

   my @array = split /:/, $array0[0];
   
   my $barcode = $array[0];


   if(exists($cell_barcodes{$barcode})) {
          $matching_line_counter++;
          my $output = $line."\n";
          print OUT $output;
   }else{
          my $output = $line."\n";
          print OUT1 $output;
   }

}#while(<BAM>)

print("I have read $bam_line_counter reads from input BAM file: $bam_file .\n");
print("I have written $matching_line_counter reads into $output_sam for selected barcodes .\n");

print "convert sam to bam file:\n";
my $output_bam=$output_dir."/cell_barcodes.bam";
my $output_bam1=$output_dir."/non_cell_barcodes.bam";
#system("$samtools_path/samtools view -@ 4 -bS $output_sam > $output_bam &");
#system("$samtools_path/samtools view -@ 4 -bS $output_sam1 > $output_bam1 &");
#system("wait");
#system("rm $output_sam");
#system("rm $output_sam1");
