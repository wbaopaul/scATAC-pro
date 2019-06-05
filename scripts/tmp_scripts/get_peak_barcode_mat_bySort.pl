#Sam file input
#perl Read_Matcher.pl --region_file example_regions.bed --read_file example_reads.sam  --output_file example_sam_output.bed.txt

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

open(REGION, $region_file ) or die("Cannot read $region_file \n");
open(READ, $read_file ) or die("Cannot read $read_file \n");
open(OUT, ">$output_file" ) or die("Cannot write $output_file \n");

my %chroms = ();

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


   my $region_id = $chrom."_".$start."_".$end;
   if( exists($chroms{$chrom}) )
   {
      push @{ $chroms{$chrom} }, $region_id;
   }else
   {
   	  my @regions = ($region_id);
      @{ $chroms{$chrom} } = [ @regions ];
   }  

}#while(my $line = <REGION>)


print "I have successfully read $region_counter regions from $region_file \n";

print "...\n";

print("Now I am starting to process read file: $read_file .\n");
print "...\n";
my $read_file_counter = 0;

my $btime = time;

my %cells = ();
my %overlap_matrix = ();

my $old_chrom = `head -1 $read_file | cut -f3`;
my $i = 0;

while(<READ>)
{
   $read_file_counter++;

	 chomp;
	 my $chrom = "XXXXXX";
     my $start = "XXXXXX";
     my $end = "XXXXXX";
     my $barcode = "XXXXXX";

	   my @array = split /\t/;
	   $chrom = $array[2];
	   $start = $array[3];
	   $end = $array[7];
	  
    ## only keep one read per pair
     if($end < $start){
      next;  
     }

     if($old_chrom ne $chrom) {
     	my $i = 0;  ## set up search region as the first region every time changing chrom
     }

     my @tmp_array = split /:/, $array[0];

	   
	 $barcode = $tmp_array[0];
    

   if(exists($chroms{$chrom}))
   {
       $cells{$barcode} = 1;
       my @regions = @{ $chroms{$chrom} };
       
	   foreach my $j ($i .. $#regions)
	   {

          my $region_id = $regions[$j];
		  my @region_coord_array = split /_/, $region_id;
		  my $region_chrom = $region_coord_array[0];
		  my $region_start = $region_coord_array[1];
		  my $region_end = $region_coord_array[2];

		     if($start >= $region_start && $start <= $region_end) #Read start is between region start and end
		     {
		       $overlap_matrix{$region_id}{$barcode} = $overlap_matrix{$region_id}{$barcode} + 1; 
		       $i = $j;  ## if region and reads are sorted, next time search from the current region
               last;  ## suppose there are no overlap between regions

		     }elsif($end >= $region_start && $end <= $region_end) #Read end is between region start and end
		     {
		       $overlap_matrix{$region_id}{$barcode} = $overlap_matrix{$region_id}{$barcode} + 1; 
		       $i = $j;
               last;           
		     }
		     

	   }#foreach

   } 
  $old_chrom = $chrom;
}#while(<READ>)



print("I have read $read_file_counter reads from input read file: $read_file .\n");

my $etime = time;
my $elapsed = $etime - $btime;

print("It takes $elapsed seconds.\n");

$btime = time;

print("Now I am writing the output to $output_file .\n");

print OUT "region_position";

foreach my $cell (keys %cells)
{
   print OUT "\t".$cell;      
}

print OUT "\n";

foreach my $chrom (keys %chroms)
{

    my @regions = @{ $chroms{$chrom} };

	foreach my $region_id (@regions)
	{
	   print OUT $region_id;
	   foreach my $cell (keys %cells)
	   {
		  print OUT "\t";
		  if(exists($overlap_matrix{$region_id}{$cell}) )
		  {
		      print OUT $overlap_matrix{$region_id}{$cell};
		  }else
		  {
		      print OUT "0";
		  }
	   }
	   print OUT "\n";
	}

}

my $etime = time;
my $elapsed = $etime - $btime;

print("It takes another $elapsed seconds for writing results.\n");

print("Output matrix is in the file $output_file  . \n");






