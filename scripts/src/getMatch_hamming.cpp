#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

int hammingDist(std::string str1, std::string str2)
{
  int count = 0;
  for(int i=0; i < str1.size(); i++)
  {
    if (str1[i] != str2[i])
      count++;
  }
  return count;
}


// Given a few barcodes, return the first white barcodes within
// a fixed hamming distance, return None if not found any matched
// [[Rcpp::export]]
std::vector< std::string >  barcode_match(std::vector< std::string >  bcs, 
                                          std::vector< std::string >  white_bcs, 
                                          int hdist_cutoff = 1){
  int bcs_size = bcs.size();
  int wbcs_size = white_bcs.size();
  std::vector< std::string >  match_str(bcs_size);
  std::string nstr("None");
  int k = 0;
  for(int i = 0; i < bcs_size; i++){
    match_str[i] = nstr;
    for(int j = k; j < wbcs_size; j++){
      if(hammingDist(bcs[i], white_bcs[j]) <= hdist_cutoff){
        match_str[i] = white_bcs[j];
        k = j; // supporse both bcs are sorted
        break;
      }
    }
  }
  return(match_str);
}