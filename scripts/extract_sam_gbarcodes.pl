#Sam file input
#perl extract_sam_gbarcodes.pl --barcode_file barcode_list.bed --read_file example_reads.sam  --output_file example_sam_output.sam

# should using sam file with header

use strict;


#Receive options from command line
use Getopt::Long;


GetOptions( 'barcode_file=s' => \my $barcode_file 
          , 'read_file=s' => \my $read_file  
          , 'output_file=s' => \my $output_file 
          );


use File::Basename;

(my $read_file_name,my $read_file_path,my $read_file_suffix) = fileparse($read_file,qr"\..[^.]*$");

if($read_file_suffix ne ".sam")
{
    print("!!!!!!!!!!!!!!!! INPUT ERROR : Unrecognized file type: $read_file_suffix !!!!!!!!!!!!!!!!!!!!\n");
    print("Sorry, I just work with sam  input files. \n");
    die("Make sure that the read file format and file extension is correct. \n");
}

if($read_file_suffix eq ".sam")
{
   print("It looks like the read file is in sam format.\n")
}

open(BARCODES, $barcode_file ) or die("Cannot read $barcode_file \n");
open(READ, $read_file ) or die("Cannot read $read_file \n");
open(OUT, ">$output_file" ) or die("Cannot write $output_file \n");


print("I am reading given barcodes information: $barcode_file .\n");

my $barcode_counter = 0;
my %barcode_list = ();
while(<BARCODES>)
{
   chomp;
   $barcode_counter++;
   $barcode_list{$_} = 0; 
   #print $_;

}#while(my $line = <BARCODES>)


print "I have successfully read $barcode_counter regions from $barcode_file \n";

print "...\n";

print("Now I am starting to process read file: $read_file .\n");
print "...\n";
my $read_file_counter = 0;
my $matching_line_counter = 0;

my $btime = time;


my $barcode = 'XXX';

while(<READ>)
{
     print OUT "$_";
     if(/\@PG/){
        last;
     }
}


while(<READ>)
{


     $read_file_counter++;
     #chomp $line;
     my @array = split /\t/, $_;	

     my @tmp_array = split /:/, $array[0];

	   
	$barcode = $tmp_array[0];



       if( exists($barcode_list{$barcode}) )
       {
          $matching_line_counter++;
          print OUT "$_";
       }

    
  
}#while(<READ>)



print("Output sam for given barcodes is in the file $output_file  . \n");
print("There are $matching_line_counter reads in selected barcodes . \n");






