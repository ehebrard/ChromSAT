/******************************************************************************/
//
//  ReactPartialCol, PartialCol, ReactTabucol and Tabucol graph coloring
//  heuristics. Reference code for the paper
//  "A Reactive Tabu Search Using Partial Solutions for the
//  Graph Coloring Problem" by Ivo Bloechliger and Nicolas Zuffery.
//
//  Copyright (C) 2003 - Ecole Poyltechnique Federale de Lausanne - EPFL, 
//  Laboratory of Operations Research South Est-ROSE, 1015 Lausanne, Switzerland
//  Written by Ivo Bloechliger, Ivo.Bloechliger@epfl.ch
//  http://rose.epfl.ch/~bloechli/coloring/
//
/******************************************************************************/
//
// This program is distributed under the terms of the GNU General Public License
// as published by the Free Software Foundation. In paticular, this program is 
// distributed WITHOUT ANY WARRANTY; without even the implied warranty of 
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. The EPFL shall in no 
// case be liable for any damage of any kind in connection with the use of this
// program.  See the GNU General Public License for more details 
// (http://www.gnu.org/copyleft/gpl.html#SEC1).
//
/******************************************************************************/

#ifndef TERMCOLORS_INCLUDED
#define TERMCOLORS_INCLUDED

#include <iomanip>

#define RED "\e[1;37;41m"
#define REDBLINK "\e[5;37;41m"
#define YELLOW "\e[1;37;42m"
#define WHITE "\e[0;37;40m\e[1;37;40m"
#define MANGENTA "\e[1;35;40m"
#define CYAN "\e[1;36;40m"


#define NODE(X) RED << setw(4) << X << WHITE
#define ERROR(X) REDBLINK << X << WHITE
#define COLOR(X) YELLOW << setw(3) << X << WHITE
#define NUMBER(X) MANGENTA << setw(9) << X << WHITE
#define NUM(X) CYAN << setw(4) << X << WHITE



#endif
