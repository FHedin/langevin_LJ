#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("One argument must be supplied (energy file name).n", call.=FALSE)
}

fname <- args[1]

# read one 64 unsigned integer corresponding to number of steps
nsteps <- readBin(fname,'integer',size=8,signed=FALSE)

print( paste("Number of steps in energy file",fname,"is :",nsteps) )

data <- matrix(NA,nrow=nsteps,ncol=4)

for(i in 1:nsteps)
{
  # read 4 double precision
  data[i,1] = readBin(fname,'double',n=1)
  data[i,2] = readBin(fname,'double',n=1)
  data[i,3] = readBin(fname,'double',n=1)
  data[i,4] = readBin(fname,'double',n=1)
}

