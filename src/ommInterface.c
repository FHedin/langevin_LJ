/**
 * \file ommInterface.c
 *
 * \brief Functions for interfacing this code with the openMM library
 *
 * \authors Florent Hédin
 *
 * \copyright Copyright (c) 2016, Florent Hedin.\n
 *            All rights reserved. \n
 *            The 3-clause BSD license is applied to this software. \n
 *            See LICENSE.txt
 *
 */

#include <stdlib.h>

#include "logger.h"
#include "ommInterface.h"

/*
 * modification of omm example file HelloSodiumChlorideInC.c
 */

/* --------------------------------------------------------------------------
 *                      INITIALIZE OpenMM DATA STRUCTURES
 * --------------------------------------------------------------------------
 */
MyOpenMMData* init_omm(ATOM atoms[], DATA* dat)
{

  /* Allocate space to hold OpenMM objects while we're using them. */
  MyOpenMMData* omm = (MyOpenMMData*)malloc(sizeof(MyOpenMMData));
  
  /* These are temporary OpenMM objects used and discarded here. */
  OpenMM_Vec3Array*       initialPosInNm;
  OpenMM_StringArray*     pluginList;
  OpenMM_NonbondedForce*  nonbond;
  OpenMM_Platform*        platform;
  OpenMM_Integrator*      lintegrator;
  
  OpenMM_Vec3 a = {2.0,0.0,0.0};
  OpenMM_Vec3 b = {0.0,2.0,0.0};
  OpenMM_Vec3 c = {0.0,0.0,2.0};
  
  /* Load all available OpenMM plugins from their default location. */
  pluginList = OpenMM_Platform_loadPluginsFromDirectory(
      OpenMM_Platform_getDefaultPluginsDirectory());
  OpenMM_StringArray_destroy(pluginList);

  /* Create a System and Force objects within the System. Retain a reference
    * to each force object so we can fill in the forces. Note: the OpenMM
    * System takes ownership of the force objects; don't delete them yourself. */
  omm->system = OpenMM_System_create();
  nonbond     = OpenMM_NonbondedForce_create();
//   OpenMM_NonbondedForce_setNonbondedMethod(nonbond,OpenMM_NonbondedForce_CutoffPeriodic);
//   OpenMM_NonbondedForce_setCutoffDistance(nonbond,0.8);
  
  OpenMM_System_addForce(omm->system, (OpenMM_Force*)nonbond);
  
//   OpenMM_System_setDefaultPeriodicBoxVectors(omm->system,&a,&b,&c);

  /* Specify the atoms and their properties:
  *  (1) System needs to know the masses.
  *  (2) NonbondedForce needs charges,van der Waals properties (in MD units!).
  *  (3) Collect default positions for initializing the simulation later.
  */
  initialPosInNm = OpenMM_Vec3Array_create(0);
  for(uint32_t n=0; n < dat->natom; n++)
  {
    OpenMM_System_addParticle(omm->system, atoms[n].pars.mass);
    
    // add particle
    OpenMM_NonbondedForce_addParticle(nonbond,atoms[n].pars.charge,
                                      atoms[n].pars.sig,
                                      atoms[n].pars.eps
    );
    
    // open mm expects nanometers for the unit of distance, so we need to convert from Angstroems to nm 
    OpenMM_Vec3 posInAng,posInNm;
    posInAng.x = atoms[n].x;
    posInAng.y = atoms[n].y;
    posInAng.z = atoms[n].z;
    posInNm = OpenMM_Vec3_scale(posInAng,OpenMM_NmPerAngstrom);
    
    // add coordinates to openmm vector structure
    OpenMM_Vec3Array_append(initialPosInNm, posInNm);
  }
  
  /* Choose an Integrator for advancing time, and a Context connecting the
  * System with the Integrator for simulation. Let the Context choose the
  * best available Platform. Initialize the configuration from the default
  * positions we collected above. Initial velocities will be zero but could
  * have been set here. */
  INTEGRATORS integType = (INTEGRATORS) dat->integrator;
  switch(integType)
  {
    case LANGEVIN:
      lintegrator = (OpenMM_Integrator*)OpenMM_LangevinIntegrator_create(
                                          dat->T, dat->friction, 
                                          dat->timestep);
      break;

    case BROWNIAN:
      lintegrator = (OpenMM_Integrator*)OpenMM_BrownianIntegrator_create(
                                          dat->T, dat->friction, 
                                          dat->timestep);
      break;

    default:
      LOG_PRINT(LOG_ERROR,"Error : invalid integrator type %d\n",integType);
      exit(-1);
      break;
  }
  
  omm->integrator = lintegrator;
  
  omm->context = OpenMM_Context_create(omm->system, omm->integrator);
  OpenMM_Context_setPositions(omm->context, initialPosInNm);
//   OpenMM_Context_setPeriodicBoxVectors(omm->context,&a,&b,&c);

  platform = OpenMM_Context_getPlatform(omm->context);
  omm->platformName = OpenMM_Platform_getName(platform);
  
  // set velocities to initial temperature
  OpenMM_Context_setVelocitiesToTemperature(omm->context,dat->T,dat->seeds[0]);
    
  return omm;
}

// -----------------------------------------------------------------------------
//                     TAKE MULTIPLE STEPS USING OpenMM 
// -----------------------------------------------------------------------------
void doNsteps_omm(MyOpenMMData* omm, int numSteps)
{
  OpenMM_Integrator_step(omm->integrator, numSteps);
}

/* --------------------------------------------------------------------------
 *                    COPY STATE BACK TO CPU FROM OPENMM
 * -------------------------------------------------------------------------- */
void getState_omm(MyOpenMMData* omm, int wantEnergy, 
                  double* timeInPs, double* energyInKcal,
                  ATOM atoms[], DATA* dat)
{
  OpenMM_State*           state;
  const OpenMM_Vec3Array* posArrayInNm;
  int                     infoMask;
  
  infoMask = OpenMM_State_Positions;
  if (wantEnergy) {
      infoMask += OpenMM_State_Velocities; /*for kinetic energy (cheap)*/
      infoMask += OpenMM_State_Energy;     /*for pot. energy (expensive)*/
  }
  /* Forces are also available (and cheap). */
  
  /* State object is created here and must be explicitly destroyed below. */
  state = OpenMM_Context_getState(omm->context, infoMask, 0);
  *timeInPs = OpenMM_State_getTime(state); /* OpenMM time is in ps already. */

  /* Positions are maintained as a Vec3Array inside the State. This will give
    * us access, but don't destroy it yourself -- it will go away with the State. */
  posArrayInNm = OpenMM_State_getPositions(state);
  for (uint32_t n=0; n < dat->natom; n++)
  {
    OpenMM_Vec3 posInAng = OpenMM_Vec3_scale(*OpenMM_Vec3Array_get(posArrayInNm,n),OpenMM_AngstromsPerNm);
    atoms[n].x = posInAng.x;
    atoms[n].y = posInAng.y;
    atoms[n].z = posInAng.z;
  }

    /* If energy has been requested, obtain it and convert from kJ to kcal. */
    *energyInKcal = 0;
    if (wantEnergy) 
        *energyInKcal = (   OpenMM_State_getPotentialEnergy(state) 
                          + OpenMM_State_getKineticEnergy(state));

    OpenMM_State_destroy(state);
  
}

// -----------------------------------------------------------------------------
//                     DEALLOCATE OpenMM OBJECTS
// -----------------------------------------------------------------------------
void terminate_omm(MyOpenMMData* omm) {
    /* Clean up top-level heap allocated objects that we're done with now. */
    OpenMM_Context_destroy(omm->context);
    OpenMM_Integrator_destroy(omm->integrator);
    OpenMM_System_destroy(omm->system);
    free(omm);
}

