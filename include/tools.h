/**
 * \file tools.h
 *
 * \brief Header file for tools.c
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

#ifndef TOOLS_H_INCLUDED
#define TOOLS_H_INCLUDED

/// fill a coordinates vector with one two or three random numbers
void get_vector(DATA *dat,int32_t mv_direction, double vec[3]);

///build an initial cluster of atoms for starting simulation
void build_cluster(ATOM at[], DATA *dat, uint32_t from, uint32_t to, int32_t mode);
///check if some atoms are too close from each other
int32_t  no_conflict(ATOM at[],uint32_t i);

///get centre of mass of the system
CM getCM(ATOM at[],DATA *dat);
///recentre the system to origin
void recentre(ATOM at[], DATA *dat);

#endif // TOOLS_H_INCLUDED
