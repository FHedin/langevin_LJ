/**
 * \file global.h
 *
 * \brief Header file included everywhere and containing variables / macros / structures of common use
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

#ifndef GLOBAL_H_INCLUDED
#define GLOBAL_H_INCLUDED

#include <stdint.h>
#include <inttypes.h>

#ifndef STDRAND
#include "dSFMT.h"
#endif

#ifndef FILENAME_MAX
#define FILENAME_MAX    4096
#endif

// macro for powers
/**
 * \def X2(a)
 * \brief A macro that squares a number
 */
#define X2(a)   (a)*(a)
/**
 * \def X3(a)
 * \brief A macro that provides the power of 3 of a number
 */
#define X3(a)   X2(a)*(a)
/**
 * \def X4(a)
 * \brief A macro that provides the power of 4 of a number
 */
#define X4(a)   X2(a)*X2(a)
/**
 * \def X6(a)
 * \brief A macro that provides the power of 6 of a number
 */
#define X6(a)   X4(a)*X2(a)
/**
 * \def X12(a)
 * \brief A macro that provides the power of 12 of a number
 */
#define X12(a)  X6(a)*X6(a)

//define where is the null file
#ifdef __unix__
#define NULLFILE "/dev/null"
#else //for MS windows
#define NULLFILE "nul"
#endif

/// if stdout has been redirected to a file from command line call ( -o option)
extern uint32_t is_stdout_redirected;

/**
 * @brief This structure holds useful variables used across the simulations,
 * it is almost always passed by pointer from one function to another one .
 */
typedef struct
{
  uint32_t natom ;    ///< Number of atoms
  
  uint8_t  platform;  ///< The OpenMM platform desired by the user ; by default fastest chosen by openMM itself
  
  uint8_t  method;    ///< MD integration method string : 'LANGEVIN' or 'BROWNIAN' (case insensitive)
  
  uint64_t nsteps ;   ///< Number of steps as a 64 bits integer to allow really long simulations (i.e. more than 2 billions)

  double inid ;       ///< An initial distance term used when randomly assigning coordinates to atoms when generating a cluster
  
  double T ;          ///< Temperature : in Kelvin
  uint8_t integrator; ///< The type on integrator used : Langevin (0) or Brownian (1)
  double friction ;   ///< Friction for Langevin/Brownian integrator : in ps^-1
  double timestep;    ///< Timestep for Langevin/Brownian integrator : in ps
  
  double cuton;       ///< cuton value for non-bonded  interactions
  double cutoff;      ///< cutoff value for non-bonded interactions

#ifndef STDRAND
  dsfmt_t dsfmt;      ///< A structure used by the dSFMT random numbers generator
  uint32_t *seeds;    ///< An array of seeds used for intialising the dSFMT random numbers generator
#endif
  uint32_t nrn ;      ///< a counter to know how many random numbers from the rn array we have used
  double *rn ;        ///< to avoid calling too often dSFMT, numbers are "cached" i.e. stored in an array ; see rand.c and rand.h
} DATA;

/**
 * @brief A structure holding the mass charge and Lennard-Jones parameters for a given atom type
 * See http://www.sklogwiki.org/SklogWiki/index.php/Lennard-Jones_model
 */
typedef struct
{
  char sym[4];    ///< atomic symbol
  double mass;    ///< atomic mass
  double charge;  ///< atomic charge ; unused in current code
  double sig ;    ///< L-J sigma parameter
  double eps ;    ///< L-J epsilon parameter
} PARAMS;

/**
 * @brief A structure representing the center of mass of a system, simply a point in
 * a 3-Dim space
 */
typedef struct
{
  double cx,cy,cz;
} CM;

/**
 * @brief A structure representing an atom
 */
typedef struct
{
  /// trick for openMM : anonymous union for either accesing coordinates atomically or directly as a structure
  union
  {
    struct {double x,y,z;};         ///< X,Y,Z coordinates in anonymous struct so accessible atomically
    struct {double x,y,z;}xyz;      ///< X,Y,Z coordinates but wrapped as a named structure : for mapping it later to OpenMM Vec3 type
  };
  char sym[4] ;   ///< atomic symbol
  PARAMS pars;   ///< substructure containing FF parameters
} ATOM;

/**
 * @brief A structure holding energy terms at a curent step
 * Using anonymous struct and unions the energy can be accessed either by term of as an array
 */
typedef struct
{
  union
  {
    struct {double epot,ekin,etot;};
    double ene[3];
  };
} ENERGIES;

#endif // GLOBAL_H_INCLUDED
