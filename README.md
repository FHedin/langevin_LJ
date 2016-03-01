----------------------------------------------
# langevin_LJ
----------------------------------------------

C program for Langevin MD simulations applied to Lennard Jones clusters only.

----------------------------------------------
## COMPILE & INSTALL
----------------------------------------------
A C compiler compatible with the C99 standard or newer is required.

You will need to Download or Compile the OpenMM library, a required dependancy : 

see https://simtk.org/home/openmm
and https://github.com/pandegroup/openmm

Be sure to have CMAKE installed (http://www.cmake.org/), available on most repositories.

Tested compilers:
	* gcc 4.6.3 and newer versions
	* clang 3.4 and newer versions
	* icc 15.0.1

Create a build directory and move to that directory: 
  * mkdir build && cd build

For buildind a debug or release (default) version : 
  * cmake -DCMAKE_BUILD_TYPE=Debug ..
  * cmake -DCMAKE_BUILD_TYPE=Release ..

Debug builds are slower but useful when debugging with gdb or valgrind.

Then, once cmake built a Makefile, just execute :
  * make

For a verbose make, use : 
  * make VERBOSE=1

Please never edit this autogenerated Makefile, edit the CMakeLists.txt instead.

For specifying another compiler on linux (for example clang or Intel icc): 
  * CC=clang cmake ..
  * CC=icc cmake ..
    
----------------------------------------------
## DOCUMENTATION
----------------------------------------------
Check the doc directory that contains programming documentation.

For generating the documentation, run doxygen in the current directory (http://www.doxygen.org)

----------------------------------------------
## LICENSING (all files excepted subdirectory dSFMT)
----------------------------------------------
Copyright (c) 2011-2016, Florent Hédin, Markus Meuwly, and the University of Basel.
All rights reserved.

The 3-clause BSD license is applied to this software.

See LICENSE.txt

----------------------------------------------
## NOTE CONCERNING dSFMT
----------------------------------------------
This project uses as random numbers generator the dSFMT code : 
http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/SFMT/

All files in the subdirectory dSFMT are licensed with a 3-clause BSD license and copyrighted as following :
Copyright (C) 2007,2008 Mutsuo Saito, Makoto Matsumoto and Hiroshima
University. All rights reserved.

See dSFMT/LICENSE_dSFMT.txt

----------------------------------------------
## OpenMM platform
----------------------------------------------

The software should automatically detect the fastest OpenMM Platform available on your computer (i.e. CUDA, OpenCL, ...)
If not it will run with the slow Reference platform : it is most probably because the OpenMM directory with the 'plugins' library is not found.
You may need to export the following OPENMM_PLUGIN_DIR environment variable to solve the problem : 

For example for a custom installation in /home/$USER/bin/openmm.

  * export OPENMM_PLUGIN_DIR=$HOME/bin/openmm/lib/plugins


