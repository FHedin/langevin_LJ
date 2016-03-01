/**
 * \file parsing.h
 *
 * \brief Header File of parsing.c
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

#ifndef PARSING_H_INCLUDED
#define PARSING_H_INCLUDED

/// will parse the input file
void parse_from_file(char fname[], DATA *dat, ATOM **at);

#endif // PARSING_H_INCLUDED
