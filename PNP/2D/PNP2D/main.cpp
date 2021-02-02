/***********************************
 **@file main.cpp
 **@author Peimeng Yin   
 **date 2013 3 22
 **
 **brief   main function
 **
 **********************************/
#include "BH2d.h"

int main(int argc,char *argv[])
{
  try{
	      clock_t start,finish;
          double totaltime;
          start=clock();
          BH2D the_app(argv[1]);
          the_app.run();
          finish=clock();
          totaltime=(double)(finish-start)/CLOCKS_PER_SEC;
          std::cout<<"CPU execute time T= "<<totaltime<<" seconds"<<std::endl;
  }
  catch(std::exception& e){
    std::cerr<<"Exception caughted:"<<std::endl
             <<e.what()
             <<std::endl;
  }
  catch(...){
    std::cerr<<"Exception caughted:"<<std::endl
             <<"unknown exception caughted."
             <<std::endl;
  }
  return 0;
}
