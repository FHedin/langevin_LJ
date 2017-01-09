/**
 * \file memory.h
 *
 * \brief Header file for memory.c
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

#ifndef MEMORY_H_INCLUDED
#define MEMORY_H_INCLUDED

#include <stdint.h>
#include <inttypes.h>

void** calloc_2D(uint32_t dim1, uint32_t dim2, size_t si);
void free_2D(uint32_t dim1, ...);

void*** calloc_3D(uint32_t dim1, uint32_t dim2, uint32_t dim3, size_t si);
void free_3D(uint32_t dim1, uint32_t dim2, ...);


#endif // MEMORY_H_INCLUDED
