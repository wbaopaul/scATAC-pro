#Sam file input

#perl get_region_barcode_mat.pl --region_file example_regions.bed --read_file example_reads.sam --read_length 100 --output_mat_file example_sam_output.bed.txt --output_barcode_file barcodes.txt 



use strict;


#Receive options from command line
use Getopt::Long;

GetOptions( 'region_file=s' => \my $region_file 
          , 'read_file=s' => \my $read_file  
          , 'read_length=s' => \my $read_length 
          , 'output_mat_file=s' => \my $output_mat_file  
          , 'output_barcode_file=s' => \my $output_barcode_file  
          );




#system("mkdir -p $output_dir");

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

open(REGION, $region_file ) or die("Cannot read $region_file \n");
open(READ, $read_file ) or die("Cannot read $read_file \n");
open(OUT, ">$output_mat_file" ) or die("Cannot write $output_mat_file \n");
open(OUT_barcode, ">$output_barcode_file" ) or die("Cannot write $output_barcode_file \n");

my %regions = ();


#my $header = <REGION>;

print("I am reading region information file: $region_file .\n");

my $region_counter = 0;
while(my $line = <REGION>)
{
   $region_counter++;
   #print $line;
   chomp $line;
   #print "***************\n";
   #print ("xx".$line);
   my @array = split /\t/, $line;
   my $chrom = $array[0];
   my $start = $array[1];
   my $end = $array[2];
   my $name = $array[3];

   my $region_id = $chrom."_".$start."_".$end;
   $regions{$region_id} = $name;

}

print "I have successfully read $region_counter regions from $region_file \n";

print "...\n";

print("Now I am starting to process reads file: $read_file .\n");
print "...\n";
my $read_file_counter = 0;
my $barcode_counter = 0;
my $matching_line_counter = 0;

my $btime = time;

my %cells = ();
my %overlap_matrix = ();

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
   
   my $read_length_half;
   if($read_file_suffix eq ".sam")
   {

	   chomp;
	   
	   my @array = split /\t/;

	   $chrom = $array[2];
	   $start = $array[3];
	   $end = $start + $read_length;
	   $qual = $array[4];

	   my $tmp = $array[0];
	   my @tmp_array = split /:/, $tmp;
           $barcode = $tmp_array[0];   ## the barcode was indicated in the name of each read in the demuliplexed fastq file, thus bam file
           

	   $read_length_half = ($read_length+1) / 2;

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
      $read_length = $end - $start;
      $read_length_half = ($read_length+1) / 2;

   }#if($read_file_suffix eq "bed")

   $cells{$barcode} =  $cells{$barcode} + 1;

   foreach my $region_id (keys %regions)
   {

      my @region_coord_array = split /_/, $region_id;
      my $region_chrom = $region_coord_array[0];
      my $region_start = $region_coord_array[1];
      my $region_end = $region_coord_array[2];

      if($chrom eq $region_chrom)
      {
         if($start >= $region_start && $start <= $region_end) #Read start is between region start and end
         { 
           #if($region_end - $start >= $read_length_half) {$overlap_matrix{$region_id}{$barcode} = $overlap_matrix{$region_id}{$barcode} + 1; }
           #use following less stringent one
           $overlap_matrix{$region_id}{$barcode} = $overlap_matrix{$region_id}{$barcode} + 1; 
           $matching_line_counter++;

         }elsif($end >= $region_start && $end <= $region_end) #Read end is between region start and end
         {
           #if($end - $region_start >= $read_length_half) {$overlap_matrix{$region_id}{$barcode} = $overlap_matrix{$region_id}{$barcode} + 1; }           
           #less stringent:
           $overlap_matrix{$region_id}{$barcode} = $overlap_matrix{$region_id}{$barcode} + 1;           
           $matching_line_counter++;
         }
         
      }#if($chrom eq $region_chrom)

   }#foreach my $region_id (keys %regions)
    
   $barcode_counter++;
  
}#while(<READ>)



print("I have read $read_file_counter reads from input read file: $read_file .\n");
print("I have read $matching_line_counter reads are overlapped .\n");


print("Now I am writing the output to $output_mat_file .\n");

print OUT "region_position";

#print OUT_region "\n";


foreach my $cell (keys %cells)
{
   print OUT "\t";      
   print OUT "$cell";      
   print OUT_barcode "$cell \t $cells{$cell}\n";      
}

print OUT "\n";      

foreach my $region_id (keys %regions)
{
   print OUT $region_id;
   foreach my $cell (keys %cells)
   {
      print OUT "\t";
      if(exists($overlap_matrix{$region_id}{$cell}))
      {
          print OUT $overlap_matrix{$region_id}{$cell};
      }else
      {
          print OUT "0";
      }
   }
   print OUT "\n";
}


#my $barcode_ratio = int($barcode_counter * 100 / $read_file_counter + 0.5);
#my $matching_ratio = int($matching_line_counter * 100 / $read_file_counter + 0.5);

#print("$barcode_counter reads (${barcode_ratio}%) have barcodes. \n");
#if($barcode_ratio < 0.5)
#{ 
#   print("!!!!!!!!!!!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!!!!!!!");
#   print("More than half of the reads do not have barcodes.\n I recommend you to check out the read file!");
#}

#print("$matching_line_counter (${matching_ratio}%) read barcodes overlap with the region file. \n");

#if($matching_ratio == 0)
#{ 
#   print("!!!!!!!!!!!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!!!!!!!");
#   print("Matching rate is 0. I recommmend you to check your input files!");
#}
print("Output matrix is in the file $output_mat_file  . \n");
print("Output barcodes are in the file $output_barcode_file  . \n");






