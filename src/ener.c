/**
 * \file ener.c
 *
 * \brief Functions for getting hard coded Potential and Gradient (Lennard-Jones and Aziz)
 *          For used defined lua potentials see instead plugins.lua.c
 *
 * \authors Florent Hedin (University of Basel, Switzerland) \n
 *          Markus Meuwly (University of Basel, Switzerland)
 *
 * \copyright Copyright (c) 2011-2015, Florent Hedin, Markus Meuwly, and the University of Basel. \n
 *            All rights reserved. \n
 *            The 3-clause BSD license is applied to this software. \n
 *            See LICENSE.txt
 *
 */

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "global.h"
#include "tools.h"
#include "ener.h"

#ifndef K_CONSTRAINT
#define K_CONSTRAINT    4.00
#endif

/* How to call this function :
 *
 *  get_LJV(at,&dat,-1) is for total energy of the whole system.
 *
 *  get_LJV(at,&dat,candidate_atom_number) is for energy evaluation of candidate_atom_number only.
 *
 */
double get_LJ_V(ATOM at[], DATA *dat, int32_t candidate)
{

    uint32_t i,j;
    double dx1,dy1,dz1;
    double dx2,dy2,dz2;
    double dcm;
    double d2, epsi_g, sig_g;
    double energy = 0.0;

    dat->E_constr = 0.0;
    CM cm = getCM(at,dat);

    if (candidate==-1)
    {
        for (i=0; i<(dat->natom); i++)
        {
            dx1=at[i].x;
            dy1=at[i].y;
            dz1=at[i].z;

            dcm = X2(cm.cx-dx1) +  X2(cm.cy-dy1) + X2(cm.cz-dz1) ;
            dat->E_constr += getExtraPot(dcm,at[i].ljp.sig,at[i].ljp.eps);

            for (j=i+1; j<(dat->natom); j++)
            {
                dx2=at[j].x;
                dy2=at[j].y;
                dz2=at[j].z;

                d2 = X2(dx2-dx1) +  X2(dy2-dy1) + X2(dz2-dz1) ;
                epsi_g = sqrt( at[i].ljp.eps * at[j].ljp.eps );
                sig_g  = 0.5*( at[i].ljp.sig + at[j].ljp.sig );

                energy += 4.0 * epsi_g *( (X12(sig_g))/(X6(d2)) - (X6(sig_g))/(X3(d2)) );
            }
        }
    }
    else
    {
        i = (uint32_t) candidate;

        dx1=at[i].x;
        dy1=at[i].y;
        dz1=at[i].z;

        dcm = X2(cm.cx-dx1) +  X2(cm.cy-dy1) + X2(cm.cz-dz1) ;
        dat->E_constr += getExtraPot(dcm,at[i].ljp.sig,at[i].ljp.eps);

        for (j=0; j<(dat->natom); j++)
        {
            if (j!=i)
            {
                dx2=at[j].x;
                dy2=at[j].y;
                dz2=at[j].z;

                d2 = X2(dx2-dx1) +  X2(dy2-dy1) + X2(dz2-dz1) ;
                epsi_g = sqrt( at[i].ljp.eps * at[j].ljp.eps );
                sig_g  = 0.5*( at[i].ljp.sig + at[j].ljp.sig );

                energy += 4.0 * epsi_g *( (X12(sig_g))/(X6(d2)) - (X6(sig_g))/(X3(d2)) );
            }
        }
    }

    return energy;
}

void get_LJ_DV(ATOM at[], DATA *dat, double fx[], double fy[], double fz[])
{
    uint32_t i=0 , j=0 ;
    double epsi_g=0.0 , sig_g=0.0;
    double dx=0.0 , dy=0.0 , dz=0.0 , d2=0.0 ;
    double de=0.0 ;

    for (i=0 ; i < dat->natom ; i++ )
    {
        fx[i] = 0.0 ;
        fy[i] = 0.0 ;
        fz[i] = 0.0 ;
        for (j=0 ; j < dat->natom ; j++ )
        {
            if (i==j) continue ;
            epsi_g = sqrt( at[i].ljp.eps * at[j].ljp.eps );
            sig_g  = 0.5*( at[i].ljp.sig + at[j].ljp.sig );
            dx = at[i].x - at[j].x ;
            dy = at[i].y - at[j].y ;
            dz = at[i].z - at[j].z ;
            d2  = dx*dx + dy*dy + dz*dz ;
            de = -24.0*epsi_g*( 2.0*(X12(sig_g))/(X6(d2)) - (X6(sig_g))/(X3(d2)) ) / d2 ;
            fx[i] += de*dx;
            fy[i] += de*dy;
            fz[i] += de*dz;
        }
    }
}

double getExtraPot(double d2, double sig, double eps)
{
    double vc = d2/(X2(K_CONSTRAINT*sig));
    vc=pow(vc,10.0);
    vc*=eps;

    return vc;
}
