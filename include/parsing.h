/**
 * \file parsing.h
 *
 * \brief Header File of parsing.c
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

#ifndef PARSING_H_INCLUDED
#define PARSING_H_INCLUDED

/// will parse the input file
void parse_from_file(char fname[], DATA *dat, ATOM **at);

#endif // PARSING_H_INCLUDED
