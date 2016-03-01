/**
 * \file main.c
 *
 * \brief C program for MD simulation of Lennard Jones clusters.
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
#include <time.h>
#include <math.h>

#include "global.h"
#include "tools.h"
#include "rand.h"
#include "io.h"
#include "parsing.h"
#include "logger.h"

// -----------------------------------------------------------------------------------------

/*
 * Some global variables initialisation here
 */

/*
 * Initialise the io structure containing file names,
 * by default everything is discarded to NULLFILE
 */
IODAT io = {NULLFILE,NULLFILE,NULLFILE,NULLFILE,1000,1000};

// pointer to FILE for trajectory, coordinates, energy
FILE *traj=NULL;
FILE *crdfile=NULL;
FILE *efile=NULL;

/*
 * boolean like values
 * is the stdout redirected ?
 */
uint32_t is_stdout_redirected=0;

/*
 * Errors, warning, etc ... --> logging.
 * Default Level is LOG_WARNING, which means that everything which is at least
 * a warning is printed (it includes Errors also).
 * See logger.h for other possibilities.
 */
LOG_LEVELS LOG_SEVERITY = LOG_WARNING;

/*
 *  End of global variables initialisation
 */

// -----------------------------------------------------------------------------------------

//prototypes of functions written in this main.c
void run_md(DATA *dat, ATOM at[]);
void help(char **argv);

// -----------------------------------------------------------------------------------------
/**
 * \brief   Main entry of the program
 * \details The 2 variables \b argc and \b argv are used for extracting
 *          command line parameters.
 * \param   argc Number of arguments, at least one as \b argv contains at least the program's name.
 * \param   argv[] Array of character strings containing the arguments.
 * \return  On exit returns EXIT_SUCCESS, EXIT_FAILURE otherwise.
 */
int main(int argc, char** argv)
{
    /* arguments parsing, we need at least "prog_name -i an_input_file"
     * prints some more instructions if needed
     */
    if (argc < 3)
    {
        fprintf(stdout,"[Info] No input file ! \n");
        help(argv);
        return EXIT_SUCCESS;
    }

    uint32_t i;
    char seed[128] = "";
    char inpf[FILENAME_MAX] = "";

    DATA dat ;
    ATOM *at = NULL;

    // default trajectory mode is binary dcd
    write_traj= &(write_dcd);

    // arguments parsing
    for (i=1; i<(uint32_t)argc; i++)
    {
        // get name of input file
        if (!strcasecmp(argv[i],"-i"))
        {
            sprintf(inpf,"%s",argv[++i]);
        }
        // get user specified seed, 128 characters max, keep it as a string for the moment
        else if (!strcasecmp(argv[i],"-seed"))
        {
            sprintf(seed,"%s",argv[++i]);
        }
        // reopen stdout to user specified file
        else if (!strcasecmp(argv[i],"-o"))
        {
            freopen(argv[++i],"w",stdout);
            is_stdout_redirected = 1 ;
        }
        // specify the logging level
        else if (!strcasecmp(argv[i],"-log"))
        {
            if (!strcasecmp(argv[++i],"no"))
            {
                LOG_SEVERITY = LOG_NOTHING;
            }
            else if (!strcasecmp(argv[i],"err"))
            {
                LOG_SEVERITY = LOG_ERROR;
            }
            else if (!strcasecmp(argv[i],"warn"))
            {
                LOG_SEVERITY = LOG_WARNING;
            }
            else if (!strcasecmp(argv[i],"info"))
            {
                LOG_SEVERITY = LOG_INFO;
            }
            else if (!strcasecmp(argv[i],"dbg"))
            {
                LOG_SEVERITY = LOG_DEBUG;
            }
            else
                fprintf(stdout,"[Warning] Unknown log level '%s' : default value used.\n\n",argv[i]);
        }
        // print help and proper exit
        else if ( !strcasecmp(argv[i],"-h") || !strcasecmp(argv[i],"-help") || !strcasecmp(argv[i],"--help") )
        {
            help(argv);
            return EXIT_SUCCESS;
        }
        // error if unknown command line option
        else
        {
            fprintf(stdout,"[Error] Argument '%s' is unknown.\n",argv[i]);
            help(argv);
            exit(-2);
        }
    }

    //prepare log files if necessary
    init_logfiles();

    // Print date and some env. variables
    fprintf(stdout,"Welcome to %s ! Command line arguments succesfully parsed, now intialising parameters...\n\n",argv[0]);
    fprintf(stdout,"Logging level is : %s : see the documentation to see which .log files are generated, and what they contain.\n\n",get_loglevel_string());
    fprintf(stdout,"Now printing some local informations : \n");
    fprintf(stdout,"DATE : %s\n",get_time());
    fprintf(stdout,"HOSTNAME : %s\n",getenv("HOSTNAME"));
    fprintf(stdout,"USER : %s\n",getenv("USER"));
    fprintf(stdout,"PWD : %s\n",getenv("PWD"));

    /*
     * Random numbers can be generated by using the standard functions from the C library (no guarantee on the quality)
     * or by using the dSFMT generator (default, extremely stable, no redudancy )
     * if STDRAND is defined we use the C library.
     *
     * Random numbers are stored in a double array dat.rn, of size dat.nrn
     *
     * If no string seed was passed by command line we generate one by using the unix timestamp
     *  -For STDRAND, this is directly used for srand()
     *  -For dSFMT, an array of integers generated by using the string seed is sent to dsfmt_init_by_array
     *
     */
    if (!strlen(seed))
        sprintf(seed,"%d",(uint32_t)time(NULL)) ;
    LOG_PRINT(LOG_INFO,"seed = %s \n",seed);
    dat.nrn = 2048 ;
    dat.rn = calloc(dat.nrn,sizeof dat.rn);
#ifdef STDRAND
    srand( (uint32_t)labs(atol(seed)) );
#else
    dat.seeds = calloc(strlen(seed),sizeof dat.seeds);
    for (i=0; i<(uint32_t)strlen(seed); i++)
        dat.seeds[i] = (uint32_t) seed[i] << 8;
    for (i=0; i<(uint32_t)strlen(seed); i++)
        dat.seeds[i] *= (dat.seeds[strlen(seed)-1]+i+1);

    dsfmt_init_by_array(&(dat.dsfmt),dat.seeds,(int32_t)strlen(seed));

    for (i=0; i<(uint32_t)strlen(seed); i++)
    {
        LOG_PRINT(LOG_INFO,"dat.seeds[%d] = %d \n",strlen(seed)-1-i,dat.seeds[strlen(seed)-1-i]);
    }
#endif
    
    // parse input file, initialise atom list
    parse_from_file(inpf,&dat,&at);

    // summary of parameters to output file
    fprintf(stdout,"\nStarting program in sequential mode\n\n");

    fprintf(stdout,"Seed   = %s \n\n",seed);

    fprintf(stdout,"Using OpenMM toolkit for energy and integration\n");

    fprintf(stdout,"Energy      saved each %d  steps in file %s\n",io.esave,io.etitle);
    fprintf(stdout,"Trajectory  saved each %d  steps in file %s\n",io.trsave,io.trajtitle);
    fprintf(stdout,"Initial configuration saved in file %s\n",io.crdtitle_first);
    fprintf(stdout,"Final   configuration saved in file %s\n\n",io.crdtitle_last);

    // again print parameters
    fprintf(stdout,"method   = %s\n",dat.method);
    fprintf(stdout,"natom    = %d\n",dat.natom);
    fprintf(stdout,"nsteps   = %"PRIu64"\n",dat.nsteps);
    fprintf(stdout,"T        = %lf \n",dat.T);
    fprintf(stdout,"friction = %lf \n",dat.friction);
    fprintf(stdout,"tstep    = %lf \n",dat.timestep);
    fprintf(stdout,"nb cuton    = %lf \n",dat.cuton);
    fprintf(stdout,"nb cutoff   = %lf \n\n",dat.cutoff);
    
    run_md(&dat,at);

    fprintf(stdout,"End of program\n");

    // free memory and exit properly
    free(dat.rn);
#ifndef STDRAND
    free(dat.seeds);
#endif
    free(at);

    // closing log files is the last thing to do as errors may occur at the end
    close_logfiles();

    return EXIT_SUCCESS;
}

// -----------------------------------------------------------------------------------------
/**
 * \brief   This function simply prints a basic help message.
 *
 * \details If any of \b -h or \b -help or \b --help are provided on the command line this help message is printed.\n
 *          If no command line parameter is present this message is also printed.\n
 *          If an unknown command line parameter is present this message is also printed.
 *
 * \param   argv Simply the same array of command line parameters, from function \b #main.
 */
void help(char **argv)
{
  fprintf(stdout,"Need at least one argument : %s -i an_input_file\n",argv[0]);
  fprintf(stdout,"optional args : -seed [a_rnd_seed] -o [output_file] -log [logging level, one of { no | err | warn | info | dbg }] \n");
  fprintf(stdout,"Example : \n %s -i input_file -seed 1330445520 -o out.txt -log info \n\n",argv[0]);
  fprintf(stdout,"The default logging level is 'warn' \n");
}

// -----------------------------------------------------------------------------------------
// OPEN MM data structures and code only included after this point : more modularity
// -----------------------------------------------------------------------------------------

#include "ommInterface.h"

/**
 * \brief   This function starts a Langevin or Brownian MD simulation using OpenMM
 *
 * \details This function is first in charge of opening all the output (coordinates, trajectory and energy) files.\n
 *          Then it runs the MD using openMM code.\n
 *          In the end it prints results, close the files and goes back to the function \b #main.
 *
 * \param   dat is a structure containing control parameters common to all simulations.
 * \param   at[] is an array of structures ATOM containing coordinates and other variables.
 */
void run_md(DATA *dat, ATOM at[])
{
  
  LOG_PRINT(LOG_INFO,"Forcing energy save frequency to be the same than trajectory save frequency.");
  io.esave = (io.esave == io.trsave) ? io.esave : io.trsave;
  
  if(!strcasecmp(dat->method,"LANGEVIN"))
  {
    dat->integrator = LANGEVIN;
  }
  else if(!strcasecmp(dat->method,"BROWNIAN"))
  {
    dat->integrator = BROWNIAN;
  }
  
  // initialise openMM code : fastest platform (usually cuda) will be selected automatically
  MyOpenMMData* omm = init_omm(at,dat);
  
  fprintf(stdout,"OpenMM automatically initialised with fastest platform : %s\n\n",omm->platformName);
  
  // print to info log file more infos concerning platform selected
  infos_omm(omm);
  
  //open required output files
  crdfile=fopen(io.crdtitle_first,"wt");
  efile=fopen(io.etitle,"wb");
  traj=fopen(io.trajtitle,"wb");

  //write initial coordinates at step 0
  write_xyz(at,dat,0,crdfile);
  fclose(crdfile);
  
  // TODO : code calling openMM for performing MD
  double time;
  double energy;
  
  // get initial energy
  getState_omm(omm,1,&time,&energy,at,dat);
  fprintf(stdout,"time (ps) \t %lf \t energy (kj/mol) %lf\n",time,energy);
  

  
  uint64_t steps = 0;
  do
  {
  // do some steps
  doNsteps_omm(omm,io.trsave);
  
  //get time energy and coordinates
  getState_omm(omm,1,&time,&energy,at,dat);
  fprintf(stdout,"time (ps) \t %lf \t energy (kj/mol) %lf\n",time,energy);
  
  steps += io.trsave;
  
  //write trajectory
  write_traj(at,dat,steps);
  
  //write energy
  fwrite(&energy,sizeof(double),1,efile);
  
  }while(steps < dat->nsteps);
  
  // END TODO
    
  terminate_omm(omm);
  
  crdfile=fopen(io.crdtitle_last,"wt");
  //write lqst coordinates
  write_xyz(at,dat,steps,crdfile);
  fclose(crdfile);
}
