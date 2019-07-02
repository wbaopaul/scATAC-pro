sourceCpp(code='
  #include <Rcpp.h>
  using namespace Rcpp;
  // for each read, return 1 if it overlap with any region, otherwise 0;
  // should sort both reads and regions
  // [[Rcpp::export]]
  IntegerVector getOverlaps_tss2read(DataFrame reads, DataFrame regions) {
      NumericVector start1 = reads["start"];
      NumericVector end1 = reads["end"];
      NumericVector start2 = regions["start"];
      NumericVector end2 = regions["end"];
          
      int n1 = start1.size(), n2 = start2.size();
      NumericVector midP1(n1), len1(n1), len2(n2), midP2(n2);
      IntegerVector over1(n1);
          
      len1 = (end1 - start1)/2;
      midP1 = (end1 + start1)/2;
          
      len2 = (end2 - start2)/2;
      midP2 = (end2 + start2)/2;
      int k1 = 0, k2=0;
      for(int i=0; i<n1; i++){
          over1[i] = 0;
          for(int j=k1; j<n2; j++){
            if((fabs(midP1[i] - midP2[j]) <= (len1[i]+len2[j]))){
              over1[i]++;
              k1 = j;
              k2 = j;
              k2++;
              while(k2 < n2 && (fabs(midP1[i] - midP2[k2]) <= (len1[i]+len2[k2]))){
                over1[i] ++;
                k2++;
              }
              break;
            }
              
            
          }
          

      }
          
         return(over1);
    }'
)