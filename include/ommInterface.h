/**
 * \file ommInterface.h
 *
 * \brief Header file for ommInterface.c
 *
 * \authors Florent Hédin
 *
 * \copyright Copyright (c) 2011-2016, Florent Hédin, Markus Meuwly, and the University of Basel. \n
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
  LANGEVIN = 0,     //< code will use a Langevin integrator from OpenMM
  BROWNIAN = 1      //< code will use a Brownian integrator from OpenMM
} INTEGRATORS;

MyOpenMMData* init_omm(ATOM atoms[], DATA* dat);
void doNsteps_omm(MyOpenMMData* omm, int numSteps);
void getState_omm(MyOpenMMData* omm, int wantEnergy, 
                  double* timeInPs, double* energyInKcal,
                  ATOM atoms[], DATA* dat);
void terminate_omm(MyOpenMMData* omm);

#endif // OMMINTERFACE_H_INCLUDED
