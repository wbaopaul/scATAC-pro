#Sam file input

## given an input sam and region, output a vector recodes the total # of reads in Mitocondrial 

#perl cal_frac_mito.pl --read_file example_reads.sam --output_file output.txt 



use strict;


#Receive options from command line
use Getopt::Long;

GetOptions(  'read_file=s' => \my $read_file  
          , 'output_file=s' => \my $output_file  
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

open(READ, $read_file ) or die("Cannot read $read_file \n");
open(OUT, ">$output_file" ) or die("Cannot write $output_file \n");



print("Now I am starting to process reads file: $read_file .\n");
print "...\n";
my $read_file_counter = 0;
my $matching_line_counter = 0;

my $btime = time;

my %cells = ();
my %overlap_vec = ();

while(<READ>)
{
   $read_file_counter++;

   my $chrom = "XXXXXX";
   my $barcode = "XXXXXX";

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

	   my $tmp = $array[0];
	   my @tmp_array = split /:/, $tmp;
           $barcode = $tmp_array[0];   ## the barcode was indicated in the name of each read in the demuliplexed fastq file, thus bam file
           

       #print $barcode."\n";
   }


   $cells{$barcode} =  $cells{$barcode} + 1;



      if($chrom eq 'chrM')
      {
           $overlap_vec{$barcode} = $overlap_vec{$barcode} + 1; 
           $matching_line_counter++;

         
      }

    
  
}#while(<READ>)

close READ;

print("I have read $read_file_counter reads from input read file: $read_file .\n");
print("And $matching_line_counter reads from input reads are overlapped in given region.\n");


print("Now I am writing the output to $output_file .\n");



print OUT "Barcode \t Total_reads \t Total_Mito_read\n";      
foreach my $cell (keys %cells)
{
   print OUT "$cell \t $cells{$cell} \t"; 
   if(exists($overlap_vec{$cell}))
   {
     print OUT "$overlap_vec{$cell}\n";
   }else
   {
     print OUT "0\n";
   }      
}



print("Output is in the file $output_file  . \n");

close OUT;




