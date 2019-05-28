#Run as follows:
#perl Katana.pl --cluster_file cluster_file --read_file read_file --read_length read_length --output_dir output_dir
#example:
#perl Katana.pl --cluster_file example_clusters.txt --read_file example_alignment.sam --read_length 50 --output_dir
#perl Katana.pl --cluster_file example_clusters.txt --read_file example_alignment.bed --output_dir example_output/
#cluster_file format:
#Barcode	Cluster
#AAACGAAAGAAATGGG-1	1
#AAACGAAAGTCTCTAG-1	3
#AAACGAACAACTAGAA-1	5
#AAACGAACAATCCATG-1	2

#Fragment (bed) file format:
#chrom  start	end	barcode			reads
#chr1	10079	10519	CAGTGTAGTAACCGAG-1	1
#chr1	10085	10339	CGTTCCACACAGCTTA-1	1
#chr1	10091	10320	GCATTGAGTGCAAGCA-1	1
#chr1	10151	10180	AACTGTGAGCTCGGCT-1	2

#Output file names: cluster_${cluster_name}.bed

#Output format: 
#chrom   start   end    mapping_quality_score  barcode
#chr1	 9997	10047	4	               GTCGTAAAGCTATCCA-1
#chr1	 9998	10048	0	               GCCAGACCAATGGCTT-1
#chr1	 9998	10048	12	               CATGCCTGTCCATTGA-1


use strict;


#Receive options from command line
use Getopt::Long;

GetOptions( 'cluster_file=s' => \my $cluster_file 
          , 'read_file=s' => \my $read_file  
          , 'read_length=s' => \my $read_length 
          , 'output_dir=s' => \my $output_dir  
          );

print("...\n");
print("Hello! My name is Katana.\n");
print("I will split your input read file $read_file into clusters based on $cluster_file .!\n");
print("...\n");


system("mkdir -p $output_dir");

use File::Basename;

#print "Input read file:".$read_file."\n";
(my $read_file_name,my $read_file_path,my $read_file_suffix) = fileparse($read_file,qr"\..[^.]*$");

if($read_file_suffix ne ".sam" && $read_file_suffix ne ".bed")
{
    print("!!!!!!!!!!!!!!!! INPUT ERROR : Unrecognized file type: $read_file_suffix !!!!!!!!!!!!!!!!!!!!\n");
    print("Sorry, I just work with sam or bed input files. \n");
    die("Make sure that the read file format and file extension is correct. \n");
}

if($read_file_suffix eq ".sam")
{
   print("It looks like the read file is in sam format.\n")
}
if($read_file_suffix eq ".bed")
{
   print("It looks like the read file is in bed format.\n")
}

open(CLUSTER, $cluster_file ) or die("Cannot read $cluster_file \n");
open(READ, $read_file ) or die("Cannot read $read_file \n");

my %cell_clusters = ();
my %cluster_files = ();
my %cluster_file_handles = ();

my $header = <CLUSTER>;

print("I am reading cluster information file: $cluster_file .\n");

my $pair_counter = 0;
while(my $line = <CLUSTER>)
{
   $pair_counter++;
   #print $line;
   chomp $line;
   #print "***************\n";
   #print ("xx".$line);
   my @array = split /\t/, $line;
   my $barcode = $array[0];
   my $cluster = $array[1];
   #print "********$cluster********\n";
   $cell_clusters{$barcode} = $cluster;
   my $cluster_file = $output_dir."/"."cluster_".$cluster.".sam";
   #print $cluster_file."\n";
   $cluster_files{$cluster} = $cluster_file;
}

print "I have successfully read $pair_counter cell-cluster pairs from $cluster_file \n";

print "...\n";


print "I am initializing output bed files .\n";

my $cluster_counter = 0;
foreach my $cluster (keys %cluster_files)
{
  $cluster_counter++; 
  my $cluster_file = $cluster_files{$cluster};
  $cluster_file_handles{$cluster} = IO::File->new();
  $cluster_file_handles{$cluster}->open(">$cluster_file");

}
print "I have found $cluster_counter clusters in your data and initialized one bed file for each of them.\n";

print "...\n";



print("Now I am starting to process read file: $read_file .\n");
print "...\n";
my $read_file_counter = 0;
my $barcode_counter = 0;
my $matching_line_counter = 0;

my $btime = time;







while(<READ>)
{
   $read_file_counter++;

   my $chrom = "XXXXXX";
   my $start = "XXXXXX";
   my $end = "XXXXXX";
   my $barcode = "XXXXXX";
   my $qual = "Not_Provided";

   if($read_file_counter % 1000000 == 0)
   {
      my $etime = time;
      my $elapsed = $etime - $btime;
      my $mil =  $read_file_counter / 1000000;
      print("I have processed $mil million reads in $elapsed seconds.\n");
   }
   
   if($read_file_suffix eq ".sam")
   {

	   chomp;
	   my $CBZ_index = index($_, "CB:Z:");
	   if ($CBZ_index  < 0) {
	    next;
	   } 

	   
	   my @array = split /\t/;
	   $chrom = $array[2];
	   $start = $array[3];
	   $end = $start + $read_length;
	   $qual = $array[4];


	   my $BCZ_index = index($_, "BC:Z:");
	     
	   my $barcode_start_index = $CBZ_index + 5;

	   my $barcode_length = $BCZ_index - 1 - $barcode_start_index;
	   
	   $barcode = substr($_, $barcode_start_index, $barcode_length);

   #print $barcode."\n";
   }

   if($read_file_suffix eq ".bed")
   {	      

      my @array = split /\t/;
      $chrom = $array[0];
      $start = $array[1];
      $end = $array[2];
      $barcode = $array[3];
      $qual = "Not_Provided";

   }#if($read_file_suffix eq "bed")

   $barcode_counter++;
  
   #my $output = $chrom."\t".$start."\t".$end."\t".$qual."\t".$barcode."\n";
   my $output = $_."\n";
  
   my $cluster = $cell_clusters{$barcode};
   
   if(exists($cluster_file_handles{$cluster}))
   {
      $matching_line_counter++;
      $cluster_file_handles{$cluster}->print($output)
   }

}#while(<READ>)



print("I have read $read_file_counter reads from input read file: $read_file .\n");

my $barcode_ratio = int($barcode_counter * 100 / $read_file_counter + 0.5);
my $matching_ratio = int($matching_line_counter * 100 / $read_file_counter + 0.5);

print("$barcode_counter reads (${barcode_ratio}%) have barcodes. \n");
if($barcode_ratio < 0.5)
{ 
   print("!!!!!!!!!!!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!!!!!!!");
   print("More than half of the reads do not have barcodes.\n I recommend you to check out the read file!");
}
print("$matching_line_counter (${matching_ratio}%) read barcodes match with the cluster file. \n");

if($matching_ratio < 0.5)
{ 
   print("!!!!!!!!!!!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!!!!!!!");
   print("More than half of the reads do not match with the barcodes in the cluster file.\n I recommend you to check out the cluster file!");
}
print("Output bed files are in the directory: $output_dir  . \n");

print("I was a nice business. See you next time! \n");


