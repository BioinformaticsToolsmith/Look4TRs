
#include "RepeatAlignE.h"

namespace motif{

  RepeatAlignE::RepeatAlignE(std::string input){
    if(input.size() > 200000){
      std::cout << input << std::endl;
    }
    seq = std::string(" ");
    seq.append(input);	
    alignRepeatWithBackTracking();		
  }
  
  RepeatAlignE::~RepeatAlignE(){
    delete[] cost;
  }
  
  void RepeatAlignE::alignRepeatWithBackTracking(){
    std::transform(seq.begin(), seq.end(),seq.begin(), ::toupper);
    std::string b = seq;
    int * rowA = (int*) calloc(1, sizeof(int));
    int * rowB;
    rowA[0] = 0;
    for(int i = 1; i < seq.size(); i++){
      rowB = (int*) calloc(i + 1, sizeof(int));
      rowB[0] = 0;
      for(int j = 1; j < i; j++){
       int s;
       if(seq[i] == b[j]){
         s = 1;
       } else {					
         s = -1;
       }
       
       rowB[j] = std::max(0, std::max(rowA[j] - 2, std::max(rowA[j - 1] + s, rowB[j - 1] - 2)) );
     }
     delete[] rowA;
     rowA = rowB;
   }
   cost = rowB;
 }
 
  /*
   *  This method will return the last row of the 2d matrix.
   */
   std::vector<int> RepeatAlignE::getLastRow(){
    
    std::vector<int> lastRow;
    for (int i = 0; i < seq.size(); i++) {
      lastRow.push_back(cost[i]);
    }
    return lastRow;
  }
}

