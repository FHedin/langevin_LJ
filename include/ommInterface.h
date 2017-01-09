/**
 * \file ommInterface.h
 *
 * \brief Header file for ommInterface.c
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


#ifndef OMMINTERFACE_H_INCLUDED
#define OMMINTERFACE_H_INCLUDED

#include "global.h"
#include "OpenMMCWrapper.h"

typedef struct {
  OpenMM_System*      system;
  OpenMM_Context*     context;
  OpenMM_Integrator*  integrator;
  const char*         platformName;
} MyOpenMMData;

typedef enum
{
  AUTO= -1, //< AUTO : let OpenMM find the fastest platform
  REF = 0,  //< REF  : on cpu, not optimised, nor parallellised
  CPU = 1,  //< CPU  : on cpu, optimised, parallellised with OpenMP
  CUDA= 2,  //< OCL  : on cpu or gpu or any accelerating device available
  OCL = 3   //< CUDA : on nvidia gpu only probably the fastest
} PLATFORMS;

extern const char* ommPlatformName[4];

typedef enum
{
  LANGEVIN = 0,     //< code will use a Langevin integrator from OpenMM
  BROWNIAN = 1      //< code will use a Brownian integrator (i.e. overdamped Langevin) from OpenMM
} INTEGRATORS;

extern const char* integratorsName[2];

MyOpenMMData* init_omm(ATOM atoms[], DATA* dat);

void doNsteps_omm(MyOpenMMData* omm, int numSteps);

void getState_omm(MyOpenMMData* omm, int wantEnergy, 
                  double* timeInPs, ENERGIES* energies, double* currentTemperature,
                  ATOM atoms[], DATA* dat);

void infos_omm(const MyOpenMMData* omm);

void terminate_omm(MyOpenMMData* omm);

#endif // OMMINTERFACE_H_INCLUDED
