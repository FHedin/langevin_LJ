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

#include "global.h"

#include "OpenMMCWrapper.h"

// typedef struct {
//     OpenMM_System*      system;
//     OpenMM_Context*     context;
//     OpenMM_Integrator*  integrator;
// } MyOpenMMData;

void init_omm(DATA* dat)
{
  
  OpenMM_System*         system;
  OpenMM_Integrator*     integrator;
  OpenMM_Context*        context;
  OpenMM_Platform*       platform;
  OpenMM_NonbondedForce* nonbond; 
  OpenMM_Vec3Array*      initPosInNm;
  OpenMM_StringArray*    pluginList;
  
  /* Load any shared libraries containing GPU implementations. */
  pluginList = OpenMM_Platform_loadPluginsFromDirectory(
      OpenMM_Platform_getDefaultPluginsDirectory());
  OpenMM_StringArray_destroy(pluginList);

  /* Create a system with nonbonded forces. System takes ownership
      of Force; don't destroy it yourself. */
  system  = OpenMM_System_create();
  nonbond = OpenMM_NonbondedForce_create(); 
  OpenMM_System_addForce(system, (OpenMM_Force*)nonbond);

  /* Create N atoms. */
  initPosInNm = OpenMM_Vec3Array_create(dat->natom);
  
}
