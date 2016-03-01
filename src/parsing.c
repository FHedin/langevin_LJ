/**
 * \file parsing.c
 *
 * \brief File containing functions for parsing the input file
 *
 * \authors Florent Hédin (University of Basel, Switzerland) \n
 *          Markus Meuwly (University of Basel, Switzerland)
 *
 * \copyright Copyright (c) 2011-2016, Florent Hédin, Markus Meuwly, and the University of Basel. \n
 *            All rights reserved. \n
 *            The 3-clause BSD license is applied to this software. \n
 *            See LICENSE.txt
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <math.h>

#include "global.h"
#include "io.h"
#include "tools.h"
#include "logger.h"

///the array of params size
static uint32_t pars_size = 0 ;

/**
 * @brief his function parses the input file, fills fields of the DATA structure,
 * and allocates the ATOM list.
 * 
 * @note strcasecmp(...) is used so the input file is case insensitive
 * 
 * @param fname Path to the input file to open
 * @param dat Common data
 * @param at The atom list
 */
void parse_from_file(char fname[], DATA *dat, ATOM **at)
{
    PARAMS *pars=NULL;

    char buff1[FILENAME_MAX]="", *buff2=NULL, *buff3=NULL ;

    FILE *ifile=NULL;
    ifile=fopen(fname,"r");

    if (ifile==NULL)
    {
        LOG_PRINT(LOG_ERROR,"Error while opening the file '%s'\n",fname);
        exit(-1);
    }

    ///iterate over each line of the text file
    while(fgets(buff1,FILENAME_MAX,ifile)!=NULL)
    {
        // skip comment line, but print it to LOG_INFO , may be useful for debugging a bad input file
        if (buff1[0]=='#')
        {
            LOG_PRINT(LOG_INFO,"Skipping line %s",buff1);
            continue;
        }

        buff2=strtok(buff1," \n\t");

        while (buff2 != NULL)
        {
            buff3=strtok(NULL," \n\t");

            ///to know which MD method we use
            if (!strcasecmp(buff2,"METHOD"))
            {
                if (!strcasecmp(buff3,"LANGEVIN"))
                  sprintf(dat->method,"%s",buff3);
                else if (!strcasecmp(buff3,"BROWNIAN"))
                  sprintf(dat->method,"%s",buff3);
                else
                {
                    LOG_PRINT(LOG_ERROR,"%s %s is unknown. Should be LANGEVIN or BROWNIAN.\n",buff2,buff3);
                    exit(-1);
                }
                
                char *friction=NULL , *tstep=NULL;
                friction = strtok(NULL," \n\t");
                friction = strtok(NULL," \n\t");
                tstep = strtok(NULL," \n\t");
                tstep = strtok(NULL," \n\t");
                dat->friction = atof(friction);
                dat->timestep = atof(tstep);
            }
            /// get the nonbonded parameters
            if (!strcasecmp(buff2,"NONBOND"))
            {
              // error if something else than nopbc given 
              if (strcasecmp(buff3,"NOPBC"))
              {
                LOG_PRINT(LOG_ERROR,"%s is not a valid keyword for PBC. Should be NOPBC and nothing else with current code.\n",buff3);
                exit(-1);
              }
              
              char *cuton=NULL , *cutoff=NULL;
              cuton = strtok(NULL," \n\t");
              
              // nocutoff desired
              if (!strcasecmp(cuton,"NOCUT"))
              {
                dat->cuton  = INFINITY;
                dat->cutoff = INFINITY;
              }
              else
              {
                cuton = strtok(NULL," \n\t");
                cutoff = strtok(NULL," \n\t");
                cutoff = strtok(NULL," \n\t");
                dat->cuton = atof(cuton);
                dat->cutoff = atof(cutoff);
              }
                
            }
            /// section where saving of energy, coordinates and trajectory is handled
            else if (!strcasecmp(buff2,"SAVE"))
            {
                ///energy saving
                if (!strcasecmp(buff3,"ENER"))
                {
                    char *title=NULL , *each=NULL;
                    title = strtok(NULL," \n\t\'");
                    sprintf(io.etitle,"%s",title);
                    each = strtok(NULL," \n\t");
                    each = strtok(NULL," \n\t");
                    io.esave = (uint32_t) atoi(each);
                }
                ///coordinates saving
                else if (!strcasecmp(buff3,"COOR"))
                {
                    char *what=NULL , *type=NULL , *title=NULL ;
                    what = strtok(NULL," \n\t");
                    ///for saving initial coordinates i.e. before simulation starts
                    if (!strcasecmp(what,"FIRST"))
                    {
                        type = strtok(NULL," \n\t");
                        title = strtok(NULL," \n\t\'");
                        sprintf(io.crdtitle_first,"%s",title);
                    }
                    ///for saving last coordinates i.e. after simulation ends
                    else if (!strcasecmp(what,"LAST"))
                    {
                        type = strtok(NULL," \n\t");
                        title = strtok(NULL," \n\t\'");
                        sprintf(io.crdtitle_last,"%s",title);
                    }
                    else if (!strcasecmp(what,"TRAJ"))
                    {
                        char *each=NULL;
                        ///only dcd type for the moment so useless variable
                        type = strtok(NULL," \n\t");
                        title = strtok(NULL," \n\t\'");
                        sprintf(io.trajtitle,"%s",title);
                        each = strtok(NULL," \n\t"); //junk
                        each = strtok(NULL," \n\t");
                        io.trsave = (uint32_t) atoi(each);
                    }
                }
            }
            /// define number of atoms 
            else if (!strcasecmp(buff2,"NATOMS"))
            {
                dat->natom = (uint32_t) atoi(buff3);
                *at = malloc(dat->natom*sizeof(ATOM));
                build_cluster(*at,dat,0,dat->natom,-1);	///< initialise the cluster with atoms at infinity initially
            }
            /// define temperature
            else if (!strcasecmp(buff2,"TEMP"))
                dat->T = atof(buff3);
            /// define number of steps as 64 bits integer
            else if (!strcasecmp(buff2,"NSTEPS"))
                dat->nsteps = (uint64_t) strtoull(buff3,NULL,0);    //atol(buff3);
            /// define a list of LJ parameters
            else if (!strcasecmp(buff2,"PARAMS"))
            {
                char *type=buff3 , *mass=NULL, *epsi=NULL , *sigma=NULL ;

                /// the array grows each time a LJPARAMS section is detected
                pars=(PARAMS*)realloc(pars,(pars_size+1)*sizeof(PARAMS));

                sprintf(pars[pars_size].sym,"%s",type);

                mass=strtok(NULL," \n\t");
                mass=strtok(NULL," \n\t");
                pars[pars_size].mass=atof(mass);
                
                pars[pars_size].charge=0.0;
                
                epsi=strtok(NULL," \n\t");
                epsi=strtok(NULL," \n\t");
                pars[pars_size].eps=atof(epsi);

                sigma=strtok(NULL," \n\t");
                sigma=strtok(NULL," \n\t");
                pars[pars_size].sig=atof(sigma);

                pars_size++;
            }
            /// for building atom list manually
            else if (!strcasecmp(buff2,"ATOM"))
            {
                uint32_t i=0,j=0,k=0,l=0;
                char *from=buff3, *to=NULL, *type=NULL, *coor=NULL;
                
                to=strtok(NULL," \n\t");
                to=strtok(NULL," \n\t");
                type=strtok(NULL," \n\t");
                coor=strtok(NULL," \n\t");
                coor=strtok(NULL," \n\t");
                
                j = (uint32_t) atoi(from) - 1;
                if (!strcasecmp(to,"END"))
                    k = dat->natom;
                else
                    k = (uint32_t) atoi(to);
                
                if (k<=0 || k>dat->natom)
                    k=dat->natom;

                LOG_PRINT(LOG_INFO,"Building an atomic list from index %d to %d and of type  %s.\n",j,k-1,type);
                
                for(i=j; i<k; i++)
                {
                    sprintf((*at)[i].sym,"%s",type);
                    for(l=0; l<pars_size; l++)
                    {
                        if (!strcasecmp(pars[l].sym,(*at)[i].sym))
                        {
                            (*at)[i].pars.mass=pars[l].mass;
                            (*at)[i].pars.charge=pars[l].charge;
                            (*at)[i].pars.eps=pars[l].eps;
                            (*at)[i].pars.sig=pars[l].sig;
                            break;
                        }
                    }
                }
                ///randomly distribute atoms
                if(!strcasecmp(coor,"RANDOM"))
                {
                    build_cluster(*at,dat,j,k,1);
                }
                ///put all atoms at origin
                else if (!strcasecmp(coor,"ZERO"))
                {
                    build_cluster(*at,dat,j,k,0);
                }
                ///read starting cnfiguration from file
                else if (!strcasecmp(coor,"FILE"))
                {
                    // initial structure read from an xyz file
                    char *fi=NULL;
                    fi=strtok(NULL," \n\t'");

                    FILE *start=NULL;
                    start = fopen(fi,"r");
                    if (start==NULL)
                    {
                        LOG_PRINT(LOG_ERROR,"Error while opening initial structure file %s\n",fi);
                        exit(-1);
                    }

                    read_xyz(*at,dat,start);

                    fclose(start);
                }
                else
                {
                    ///random is default
                    build_cluster(*at,dat,j,k,1);
                }
            }

            buff2=strtok(NULL," \n\t");
        }
    }

    fclose(ifile);
    free(pars);
}
