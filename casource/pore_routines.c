
/****************************************************************/
/*      Copyright (c) 1993 Peter D Lee                          */
/*      Copyright (c) 1998 Dept. of Materials, ICSTM            */
/*      All Rights Reserved                                     */
/*      The copyright notice above does not evidence any        */
/*      actual or intended publication of such source code,     */
/*      and is an unpublished work by Dept. of Materials, ICSTM.*/
/*      continuing D Phil work from University of Oxford        */
/****************************************************************/
/* This code is part of the umats routines developed at in the  */
/* Materials Processing Group, Dept. of Materials, ICSTM.       */
/*      email p.d.lee or r.atwood @imperial.ac.uk for details   */
/****************************************************************/

/********************************************************************************/
/*  This version is distributed under a BSD style public license, as follows:   */
/*                                                                              */
/*  Copyright (c) 2007, Dept. of Materials, Imperial College London             */
/*  All rights reserved.                                                        */
/*  Redistribution and use in source and binary forms, with or without          */
/*  modification, are permitted provided that the following conditions          */
/*  are met:                                                                    */
/*                                                                              */
/*  * Redistributions of source code must retain the above copyright            */
/*  notice, this list of conditions and the following disclaimer.               */
/*                                                                              */
/*  * Redistributions in binary form must reproduce the above                   */
/*  copyright notice, this list of conditions and the following                 */
/*  disclaimer in the documentation and/or other materials provided             */
/*  with the distribution.                                                      */
/*                                                                              */
/*  * Neither the name of the Dept. of Materials, Imperial College London, nor  */
/*  the names of its contributors may be used to endorse or promote products    */
/*  derived from this software without specific prior written permission.       */
/*                                                                              */
/*  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS         */
/*  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT           */
/*  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR       */
/*  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT        */
/*  OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,       */
/*  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED    */
/*  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR      */
/*  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF      */
/*  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING        */
/*  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS          */
/*  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                */
/********************************************************************************/
/*END of LICENSE NOTICE*/

/* $Id: pore_routines.c 1342 2008-07-23 15:45:00Z  $*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "machine.h"
#include "blocks.h"
#include "pore.h"
#include "pore_routines.h"
extern int getxyz (int cellnum, int *nc, int *Cell);

/*****************************************************/
/* add all appropriate cells at the boundary of a pore */

void add_all_cells (BB_struct * bp, CA_FLOAT * fs, PORE_str * c_p)
{
  int *nni, *nnip, *nniend, min_cell, cellnum;
  int check_flag, check_cell;

  CA_FLOAT *fsp, min_fs, this_fs;
  p_c_list *bdy;
  p_c_node *node, *checknode, *newnode;

  nni = nnip = bp->nbhd.nnq;
  nniend = nni + 6;
  bdy = c_p->boundary;
  node = bdy->first;

  while (node != NULL) {

    cellnum = node->cellnum;
    fsp = fs + cellnum;

    for (nnip = nni; nnip < nniend; nnip++) {
      this_fs = *(fsp + *nnip);
      if ((this_fs >= LIQUID) && (this_fs < SOLID)) {
        check_cell = cellnum + *nnip;
        check_flag = 1;
        /* check if it is already in the pore */
        for (checknode = bdy->first; checknode != NULL; checknode = checknode->next) {
          if (check_cell == checknode->cellnum) {
            check_flag = 0;
            break;
          }
        }                       /* end loop to check if already in pore */
        /* add to the list if not already in pore */
        if (check_flag) {
          add_to_list (bdy, new_node (check_cell));
          c_p->ncells++;
        }
      }                         /* end check if liquid/ partly solid cell */
    }                           /* end of neigbour chek loop */
    node = node->next;
  }                             /* end traverse bdy list */
  return;
}                               /* end of add_all_cells */

/*****************************************************/
/* Find the gradient of the hydrogen solute field */
/* this is used in the next subroutine find_new_cell */
/* THIS ROUTINE is specific for 6 Neigbour central FD scheme */
/*************************************************************/
/* still needs to handle NOT_CASTING somehow */
CA_FLOAT get_grad (CA_FLOAT * solp, int *nni, int ncsb, int cellnum)
{
  CA_FLOAT gradx, grady, gradz;
  CA_FLOAT gradsq;

/**there is a BUG here  \todo  Fix this bug rewrite so as to use PADDED sol array!  -- general -- hard -- maybe obsolete by xly curvature?*/
#ifdef GRAD
  gradx = *(solp + nni[0]) - *(solp + nni[1]);
  grady = *(solp + nni[2]) - *(solp + nni[3]);
#ifdef EXTERNAL_3D
  if ((cellnum + nni[4]) < 0) {
    gradz = 2 * (*(solp + nni[5]) - *solp);
  } else if ((cellnum + nni[5]) > ncsb) {
    gradz = 2 * (*solp - *(solp + nni[4]));
  } else {
    gradz = *(solp + nni[4]) - *(solp + nni[5]);
  }
#else
  gradz = 0.0;
#endif /* EXTERNAL_3D */
  gradsq = gradx * gradx + grady * grady + gradz * gradz;
#else
  gradsq = 1;
#endif /* GRAD */
  return gradsq;
}

/*****************************************************/
/* find one cell which has minimum of fraction solid */
/* neighbouring the boundary of the pore */
/* also from these fine the one with minium grad H   */
/* and choose randombly from among multiple choices   */
/*****************************************************/
int find_new_cell (BB_struct * bp, CA_FLOAT * fs, CA_FLOAT * sol, PORE_str * c_p)
{

  int *nni, *nnip, *nniend, min_cell, cellnum;
  int check_flag = 0, check_cell;
  int *choose_list, nchoose = 0, choice;

  CA_FLOAT *fsp, *solp, min_grad, min_fs, this_fs, this_grad, gradsq;
  p_c_list *bdy;
  p_c_node *node, *checknode;

  min_grad = LARGEGRAD;
  choose_list = (int *) calloc (6 * c_p->ncells, sizeof (int));

  nni = nnip = bp->nbhd.nnq;
  nniend = nni + 6;
  bdy = c_p->boundary;

  node = bdy->first;
  nchoose = 1;
  min_fs = 1.0;
  min_cell = -1;

  /* check the boundary cells */
  for (node = bdy->first; node != NULL; node = node->next) {
    cellnum = node->cellnum;
    fsp = fs + cellnum;
    solp = sol + cellnum;

    for (nnip = nni; nnip < nniend; nnip++) {
      check_cell = cellnum + *nnip;
      /* restart loop if  out of bounds */
      if (check_cell < 0 || check_cell >= bp->ncsb) {
        continue;
      }

      this_fs = *(fsp + *nnip);

      /* or not in casting */
      if (this_fs == NOT_CASTING)
        continue;
      /* or greater than previous miniumum */
      if (this_fs > min_fs)
        continue;

      check_flag = 1;

      /* check if it is already in the pore */
      for (checknode = bdy->first; checknode != NULL; checknode = checknode->next) {
        if (check_cell == checknode->cellnum) {
          check_flag = 0;
          break;
        }
      }

      if (check_flag) {

        /*randomize if equal fraction solid and gradient */
        /* don't need to calculate the gradient if fs is greater */
        if (this_fs == min_fs) {
          /* find the gradient of hydrogen in the target cell */
          this_grad = get_grad (solp, nni, bp->ncsb, cellnum);
          if (this_grad == min_grad) {
            /*equal fs and grad, so add to choice list */
            choose_list[nchoose++] = check_cell;
          } else if (this_grad < min_grad) {
            /* new miniumum grad so reset choice list */
            min_cell = check_cell;
            min_grad = this_grad;
            nchoose = 1;
            choose_list[0] = check_cell;
          }
          /* otherwise, this_grad > min_grad so do nothing */
        } else {
          /* new miniumum fs so reset choice list */
          min_cell = check_cell;
          min_fs = this_fs;
          min_grad = this_grad;
          nchoose = 1;
          choose_list[0] = check_cell;
        }                       /*end add/reset choice list (choose_list) */
      }                         /*end if check_flag */
    }                           /* end of neighbour chek loop */
  }                             /* end traverse bdy list */
  /* select a random cell if there is more than one candidate */
  if (nchoose > 1) {
    choice = (int) (floor (((double) (nchoose)) * drand48 ()));
    min_cell = choose_list[choice];
  }
  free (choose_list);
  return (min_cell);
}                               /* end of find_new_cell */

/*****************************************************/
/* remove a node from the list -- does not destroy it */
/*****************************************************/
void remove_from_list (p_c_list * list, p_c_node * node)
{
  if (list->first == node) {
    /* the cell is first on the list */
    list->first = node->next;
    if (list->first == NULL) {
      /* the cell is the only one in the list */
      list->last = NULL;
      node->next = node->previous = NULL;
      return;
    }
    list->first->previous = NULL;
    node->next = node->previous = NULL;
    return;
  }

  if (list->last == node) {
    /* the cell is last on the list */
    list->last = node->previous;
    if (list->last == NULL) {
      /* the cell is the only one in the list */
      list->first = NULL;
      return;
    }
    list->last->next = NULL;
    return;
  }

  /* ordinary deletion */
  node->previous->next = node->next;
  node->next->previous = node->previous;
  return;
}                               /* end of remove_from_list */

/* add a node to the list */
void add_to_list (p_c_list * list, p_c_node * node)
{

  if (list->last == NULL) {
    /* the list is empty */
    list->first = node;
    node->previous = node->next = NULL;
    list->last = node;
    return;
  }

  /* add to the end */

  node->previous = list->last;
  list->last->next = node;
  list->last = node;
  node->next = NULL;
  return;
}                               /* end of add_to_list */

/* create a new node and initialise the cellnumber */
p_c_node *new_node (int cellnum)
{
  p_c_node *node;

  node = calloc (1, sizeof (p_c_node));
#ifdef DBM_VERBOSE
  fprintf (stderr, "new_node:, cell, %i, mem, %x \n", cellnum, node);
#endif /*DBM_VERBOSE */
  node->cellnum = cellnum;
  return (node);
}

/* find a node by cellnum */
p_c_node *find_node (p_c_list * list, int cellnum)
{
  p_c_node *node;

  node = list->first;
  while (node != NULL) {
    if (node->cellnum == cellnum)
      return (node);
    node = node->next;
  };
  return (NULL);
}

/* free a node -- should be removed first! */
void delete_node (p_c_node * node)
{
  if ((node->next != NULL) || (node->previous != NULL)) {
    fprintf (stderr, "ERROR:delete_node: cannot remove node, as it is attatched to something!\n");
    return;
  }
#ifdef DBM
  fprintf (stderr, "delete_node:, cell, %i, mem, %x \n", node->cellnum, node);
#endif /*DBM*/
    free (node);
}

/* move a node from boundary to body */
void move_to_body (PORE_str * this_pore, p_c_node * node)
{
  remove_from_list (this_pore->boundary, node);
  add_to_list (this_pore->body, node);
}                               /* end of move_to_body */

/* traverse the list and print out cellnum */
void print_list (p_c_list * list, FILE * outfile)
{
  p_c_node *node;

  node = list->first;
  fprintf (outfile, "Cellnum list:");
  while (node != NULL) {
    fprintf (outfile, ",%i", node->cellnum);
    node = node->next;
  }
  fprintf (outfile, "\n");
}

/*************************************************/
/* traverse the list and find extent of the pore */
/*************************************************/
void find_minmax (p_c_list * list, int cellminmax[])
{
  /* 0 - min x */
  /* 1 - min y */
  /* 2 - min z */
  /* 3 - max x */
  /* 4 - max y */
  /* 5 - max z */
  int i;
  p_c_node *node;

  for (i = 0; i < 3; i++) {
    cellminmax[i] = 100000;
    cellminmax[i + 3] = -100000;
  }

  node = list->first;
  while (node != NULL) {
    for (i = 0; i < 3; i++) {
      if (node->Cell[i] < cellminmax[i]) {
        cellminmax[i] = node->Cell[i];
      }
      if (node->Cell[i] > cellminmax[i + 3]) {
        cellminmax[i + 3] = node->Cell[i];
      }
    }
    node = node->next;
  }
}

/*************************************************/
/* traverse the list and find extent of the pore */
/*************************************************/
void free_porelist (p_c_list * list)
{
  int i;
  p_c_node *node;

  node = list->first;
  while (node != NULL) {
    remove_from_list (list, node);
    delete_node (node);
    node = list->first;
  }
}                               /* end free_porelist */

#ifdef TEST_PORE_LIST
void main ()
{
  PORE_str pore;
  p_c_node *newnode;
  int i;

  pore.body = calloc (1, sizeof (p_c_list));
  pore.boundary = calloc (1, sizeof (p_c_list));

  for (i = 0; i < 20; i++) {
    newnode = new_node (i);
    add_to_list (pore.boundary, newnode);
  }
  print_list (pore.boundary, stdout);

  newnode = find_node (pore.boundary, 12);
  remove_from_list (pore.boundary, newnode);
  print_list (pore.boundary, stdout);
  delete_node (newnode);

  newnode = pore.boundary->first;
  while (newnode != NULL) {
    remove_from_list (pore.boundary, newnode);
    delete_node (newnode);
    print_list (pore.boundary, stdout);
    newnode = pore.boundary->first;
  }

}
#endif /*TEST_PORE_LIST */

/* Little subroutine to get rcs id into the object code */
/* so you can use ident on the compiled program  */
/* also you can call this to print out or include the rcs id in a file */
char const *rcs_id_pore_routines_c ()
{
  static char const rcsid[] = "$Id: pore_routines.c 1342 2008-07-23 15:45:00Z  $";

  return (rcsid);
}

/* end of rcs_id_subroutine */

/*
*/
