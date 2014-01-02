/*    Copyright (C) David M. Rogers, 2014
 *    
 *    David M. Rogers <predictivestatmech@gmail.com>
 *    Nonequilibrium Stat. Mech. Research Group
 *    Department of Chemistry
 *    University of South Florida
 *
 *    This file is part of rbtree.
 *
 *    rbtree is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    rbtree is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with rbtree.  If not, see <http://www.gnu.org/licenses/>.
 */

/* The cmp function operates between nodes (void *N)-s.
 * These must store L, R (void *)-s at N + coff.
 * The black (0) / red (1) bit is used at
 * the masked bit of N+boff.
 */
typedef struct {
    int (*cmp)(const void *, const void *);
    unsigned int coff, boff;
    unsigned char mask; // contains a one where red/black bit is set.
    void *nil;
} rbop_t;

void new_tree(void *N, const rbop_t *o);
void *add_node(void **N, void *A, const rbop_t *o);
void *del_node(void **N, const void *A, const rbop_t *o);
void *lookup_node(void *N, const void *A, const rbop_t *o);

// returns mask or 0
unsigned char get_mask(const void *N, const rbop_t *o);
