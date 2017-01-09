/**
 * \file rand.h
 *
 * \brief Header file of rand.c file
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

#ifndef RAND_H_INCLUDED
#define RAND_H_INCLUDED

/// get a uniformly distributed random number 
double get_next(DATA *dat);

/// get a normally distributed random number
// double get_BoxMuller(DATA *dat);

/// if we want to test the random numbers generators
// void test_norm_distrib(DATA *dat, uint32_t n);

#endif // RAND_H_INCLUDED
