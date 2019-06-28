#Run as follows:
#perl split_sam.parametric.pl --cluster_file cluster_file --sam_file sam_file --output_dir output_dir
#example:
#perl split_sam_to_bed_by_clusters.pl --cluster_file example_clusters.txt --sam_file example_alignment.sam  --output_dir example_output/
#cluster_file format:
#Barcode	Cluster
#AAACGAAAGAAATGGG-1	1
#AAACGAAAGTCTCTAG-1	3
#AAACGAACAACTAGAA-1	5
#AAACGAACAATCCATG-1	2

#Output file names: cluster_${clsuter_name}.bed

use strict;


#Receive options from command line
use Getopt::Long;

GetOptions( 'cluster_file=s' => \my $cluster_file 
          , 'sam_file=s' => \my $sam_file  
          , 'output_dir=s' => \my $output_dir  
          );


system("mkdir -p $output_dir");

open(CLUSTER, $cluster_file ) or die("Cannot read $cluster_file \n");
open(SAM, $sam_file ) or die("Cannot read $sam_file \n");

my %cell_clusters = ();
my %cluster_files = ();
my %cluster_file_handles = ();

my $header = <CLUSTER>;

print("I am reading cluster information file: $cluster_file \n");

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


print "I am initializing output sam files \n";

my $cluster_counter = 0;
foreach my $cluster (keys %cluster_files)
{
  $cluster_counter++; 
  my $cluster_file = $cluster_files{$cluster};
  $cluster_file_handles{$cluster} = IO::File->new();
  $cluster_file_handles{$cluster}->open(">$cluster_file");

}
print "I have found $cluster_counter clusters in your data and initialized one sam file for each of them.\n";

print "...\n";



print("I am reading SAM file: $sam_file \n");
print "...\n";
my $sam_line_counter = 0;
my $barcode_counter = 0;
my $matching_line_counter = 0;

my $btime = time;
my $barcode="XXXX";
while(my $line = <SAM>)
{
   $sam_line_counter++;
   chomp $line;
   if($sam_line_counter % 1000000 == 0)
   {
      my $etime = time;
      my $elapsed = $etime - $btime;
      my $mil =  $sam_line_counter / 1000000;
      print("I have processed $mil million reads in $elapsed seconds.\n");
   }
   $barcode_counter++;
   
   my @array0 = split /\t/, $line;

   my @array = split /:/, $array0[0];
   
   my $barcode = $array[0];

   #print $barcode."\n";

   if(exists($cell_clusters{$barcode})) {
       my $cluster = $cell_clusters{$barcode};
       
       if(exists($cluster_file_handles{$cluster}))
       {
          $matching_line_counter++;
          my $output = $line."\n";
          $cluster_file_handles{$cluster}->print($output)
          #print $cluster_file_handles{$cluster} "$_";
       }
   }
}#while(<SAM>)

print("I have read $sam_line_counter reads from input SAM file: $sam_file .\n");
print("I have written $sam_line_counter reads into each cluster .\n");
