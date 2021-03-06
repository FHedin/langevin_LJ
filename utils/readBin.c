/**
 * \file readBin.c
 *
 * \brief Basic file for reading content of energy binary file generated by the program
 *
 * \authors Florent Hédin (École des Ponts - ParisTech) \n
 *          Tony Lelièvre (École des Ponts - ParisTech)
 *
 * \copyright Copyright (c) 2016-2017, Florent Hédin, Tony Lelièvre, and École des Ponts - ParisTech \n
 *            All rights reserved. \n
 *            The 3-clause BSD license is applied to this software. \n
 *            See LICENSE.txt
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>

int main(int argc, char** argv)
{
  
  FILE* fi = fopen(argv[1],"rb");
  
  uint64_t saved;
  fread(&saved,sizeof(uint64_t),1,fi);
  
  printf("Number of frames saved : %lu\n",saved);
  
  double time=0.;
  double epot=0.,ekin=0.,etot=0.;
  for(uint64_t i=0;i<saved;i++)
  {
    //read time and energy terms
    fread(&time,sizeof(double),1,fi);
    fread(&epot,sizeof(double),1,fi);
    fread(&ekin,sizeof(double),1,fi);
    fread(&etot,sizeof(double),1,fi);
    
    printf("time (ps) \t %lf \t epot (kj/mol) \t %lf \t ekin (kj/mol) \t %lf \t etot (kj/mol) \t %lf\n",time,epot,ekin,etot);
  }
  
  fclose(fi);
  
  return 0;
  
}
