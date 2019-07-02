#Input: fragment file with region position,barcode,ncounts
#perl Read_Matcher.pl --region_file example_region.bed --read_file fragments.bed  --output_file example_output.mat

use strict;

#Receive options from command line
use Getopt::Long;

my $region_file; 
my $read_file;
my $output_file;

GetOptions( 'region_file=s' => \$region_file 
          , 'read_file=s' => \$read_file  
          , 'output_file=s' => \$output_file 
          );


use File::Basename;

(my $read_file_name,my $read_file_path,my $read_file_suffix) = fileparse($read_file,qr"\..[^.]*$");

if($read_file_suffix ne ".bed" && $read_file_suffix ne ".txt")
{
    print("!!!!!!!!!!!!!!!! INPUT ERROR : Unrecognized file type: $read_file_suffix !!!!!!!!!!!!!!!!!!!!\n");
    print("Sorry, I just work with bed like input files. \n");
    die("Make sure that the read file format and file extension is correct. \n");
}


open(REGION, $region_file ) or die("Cannot read $region_file \n");
open(READ, $read_file ) or die("Cannot read $read_file \n");
open(OUT, ">$output_file" ) or die("Cannot write $output_file \n");

my %allRegions = ();

print("I am reading cluster information file: $region_file .\n");

my $region_counter = 0;
while(my $line = <REGION>)
{
   $region_counter++;
   chomp $line;
   my @array = split /\t/, $line;
   my $chrom = $array[0];
   my $start = $array[1];
   my $end = $array[2];


   my $curr_region = $chrom."_".$start."_".$end;
   if( exists($allRegions{$chrom}) )
   {
      push @{ $allRegions{$chrom} }, $curr_region;
   }else
   {
   	  my @regions = ($curr_region);
       $allRegions{$chrom} = [ @regions ];
   }  

}#while(my $line = <REGION>)


print "I have successfully read $region_counter regions from $region_file \n";

print "...\n";

print("Now I am starting to process read file: $read_file .\n");
print "...\n";
my $read_file_counter = 0;

my $btime = time;

my %overlap_matrix = ();

my $old_chrom = `head -1 $read_file | cut -f1`;
my $i = 0;

my $chrom = "XXXXXX";
my $start = "XXXXXX";
my $end = "XXXXXX";
my $barcode = "XXXXXX";
my $totalOverlap = 0;
print("Now working on $old_chrom...\n");

while(<READ>)
{
   $read_file_counter++;

	 chomp;

	 my @array = split /\t/;
	 $chrom = $array[0];
	 $start = $array[1];
	 $end = $array[2];
	 $barcode = $array[3];
     #$ndup = $array[4]; --ignore it if don't want include duplicates

    if($old_chrom ne $chrom) {
      print("Now working on $chrom...\n");
      my $i = 0;  ## set up search region as the first region every time changing chrom
    }

   if(exists($allRegions{$chrom}))
   {
      
       my @regions = [ @{ $allRegions{$chrom} } ];
         
  	   foreach my $j ($i .. $#regions)
  	   {
          my $curr_region = $regions[$j];
          my @region_coord_array = split /_/, $curr_region;
          my $region_chrom = $region_coord_array[0];
          my $region_start = $region_coord_array[1];
          my $region_end = $region_coord_array[2];

         if($start >= $region_start && $start <= $region_end) #read start is between region start and end
         {
           $totalOverlap++;
           $overlap_matrix{$curr_region}{$barcode} = $overlap_matrix{$curr_region}{$barcode} + 1; 
           $i = $j;  ## if region and reads are sorted, next time search from the current region
           last;  ## suppose there are no overlap between regions

         }elsif($end >= $region_start && $end <= $region_end) #read end is between region start and end
         {
           $totalOverlap++;
           $overlap_matrix{$curr_region}{$barcode} = $overlap_matrix{$curr_region}{$barcode} + 1; 
           $i = $j;
           last;           
         }
               

  	   }#foreach

   } 
  $old_chrom = $chrom;
}#while(<READ>)

my $etime = time;
my $elapsed = $etime - $btime;


print("I have read $read_file_counter reads from input read file: $read_file .\n");
print("I have found $totalOverlap reads overlapped with $region_file .\n");

print("It takes $elapsed seconds.\n");
print("Now I am writing the output to $output_file .\n");


$btime = time;
foreach my $region (keys %overlap_matrix){
    foreach my $cell (keys $overlap_matrix{$region})
    {
        print OUT $region."\t".$cell."\t".$overlap_matrix{$region}{$cell}."\n";
    }
}


my $etime = time;
my $elapsed = $etime - $btime;

print("It takes another $elapsed seconds for writing results.\n");

print("Output matrix is in the file $output_file  . \n");






