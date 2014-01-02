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

// stdio is only really needed for printing error messages
#include <stdio.h>
#include "rbtree.h"

/*****************  Red/Black Trees (in your data str) **************/

static void color_red(void *N, const rbop_t *o) {
    unsigned char *u = N + o->boff;
    *u |= o->mask;
}
static void color_black(void *N, const rbop_t *o) {
    unsigned char *u = N + o->boff;
    *u &= ~(o->mask);
}
// this only understands m zero (black) or nonzero (red)
static void set_mask(void *N, const unsigned char m, const rbop_t *o) {
    if(m) color_red(N, o);
    else color_black(N, o);
}
// returns mask or 0
unsigned char get_mask(const void *N, const rbop_t *o) {
    const unsigned char *u = N + o->boff;
    return *u & o->mask;
}

static void set_left(void *N, void *x, const rbop_t *o) {
    void **u = N + o->coff;
    *u = x;
}
static void set_right(void *N, void *x, const rbop_t *o) {
    void **u = N + o->coff + sizeof(void *);
    *u = x;
}
static void *get_left(void *N, const rbop_t *o) {
    void **u = N + o->coff;
    return *u;
}
static void *get_right(void *N, const rbop_t *o) {
    void **u = N + o->coff + sizeof(void *);
    return *u;
}


/* Magic internal data structure. */
typedef struct rbtrav_s rbtrav_t;
struct rbtrav_s {
    rbtrav_t *up;
    void *N;
    int d; // direction taken from up to N
};

static int recurse_tree(void **rep, void *A, void *G, void *P, void *C,
                    int dg, int dp, const rbop_t *o); // add
static int pop_extreme(void **ret, const void *A, int dir, rbtrav_t *self,
                    const rbop_t *o); // del


void *lookup_node(void *N, const void *A, const rbop_t *o) {
    void *C = N;
    int d;
    
    while(C != o->nil) {
        N = C;
        d = o->cmp(A, C);
        if(d < 0) C = get_left(C, o);
        else if(d > 0) C = get_right(C, o);
        else break;
    }
    return C;
}

// Now the serious stuff.
void *add_node(void **N, void *A, const rbop_t *o) {
    void *R = o->nil;
    int st, d;
    
    if(*N == o->nil) {
        new_tree(A, o);
        *N = A;
        return o->nil;
    }

    /*d = o->cmp(A, *N);
    if(d < 0) R = get_left(*N, o);
    else if(d > 0) R = get_right(*N, o);
    else { // replaces root
        R = *N;
        color_black(A, o);
        set_left(A, get_left(R, o), o);
        set_right(A, get_right(R, o), o);
        *N = A;
        return R;
    }*/
    st = recurse_tree(&R, A, o->nil, (void *)N - o->coff, *N, 0, -1, o);

    if(st) {
        if(st == 1) {
            color_black(*N, o); // Have been reddened, turn back!
        } else {
            fprintf(stderr, "Unknown error in add_node "
                        "-- is this a red-black tree?\n");
        }
    }
    return R;
}


/*  Red-black binary trees are implemented using a recursive,
 *  finite state machine that traverses the tree but
 *  back-tracks where necessary to achieve re-coloring.
 *  Keeping the last three nodes:
 *  G -dg-> P -dp-> C [current] (-d-> A[proposed])
 *  in memory keeps the number of states low.
 *
 * Starting state:
 *   0. initial traversal to addition point
 *   1. resuming a re-coloring
 *
 * Return value:
 *   0. done - return nil to caller
 *   1. re-color
 *   2. re-color, starting at parent
 *
 */
static int recurse_tree(void **rep, void *A, void *G, void *P, void *C,
                 int dg, int dp, const rbop_t *o) {
    void *N, *U;
    int d, st;

    d = o->cmp(A, C);
    if(d < 0) N = get_left(C, o);
    else if(d > 0) N = get_right(C, o);
    else { // replacement case
        if(dp < 0) set_left(P, A, o);
        else      set_right(P, A, o);
        set_mask(A, get_mask(C, o), o);
        set_left(A, get_left(C, o), o);
        set_right(A, get_right(C, o), o);
        *rep = C;
        return 0;
    }
    if(N != o->nil) {
        st = recurse_tree(rep, A, P, C, N, dp, d, o);
        switch(st) { // return cases
            case 0:
                return st;
            case 1:
                goto lbl1;
            case 2:
                return 1; // popping the stack once to get to lbl1
        }
    }
    // begin actual insertion, replacing a black leaf-node somewhere.
    color_red(A, o);
    set_left(A, o->nil, o);
    set_right(A, o->nil, o);
    if(d < 0) set_left(C, A, o);
    if(d > 0) set_right(C, A, o);
    N = A; // re-label child (don't ref. A below)

lbl1: // just colored/inserted C -d-> N (red), check for violations
    if(get_mask(C, o) == 0) { // black
        return 0; // done.
    }
    if(dp < 0) {
        U = get_right(P, o);
        if(U != o->nil && get_mask(U, o)) { // red
            //show_tree("case1", (Ast *)P, 0);
            color_red(P, o);
            color_black(C, o); color_black(U, o);
            //show_tree("case1-fin", (Ast *)P, 1);
            return 2; // popping the stack twice to get to lbl1
        }
        // have an inward-leaning chain P -dp-> C -d-> N of red nodes
        if(d > 0) { // move N up and rotate left
            //show_tree("case2", (Ast *)P, 0);
            set_left(P, N, o);
            set_right(C, get_left(N, o), o);
            set_left(N, C, o);
            U = C; C = N; N = U; // swap(C, N);
            //show_tree("case2-fin", (Ast *)P, 1);
            //d  = -d;
        }
        // have an outward-leaning chain G -dg-> P -dp-> C -d-> N of red nodes
        //show_tree("case3", (Ast *)P, 0);
        color_black(C, o); // right-rotate G - P - C
        color_red(P, o);
        if(dg < 0) set_left(G, C, o);
        else      set_right(G, C, o);
        set_left(P, get_right(C, o), o);
        set_right(C, P, o);
        //show_tree("case3-fin", (Ast *)C, 1);
    } else { // as above, but now right <-> left
        U = get_left(P, o);
        if(U != o->nil && get_mask(U, o)) { // red
            //show_tree("case1", (Ast *)P, 0);
            color_red(P, o);
            color_black(C, o); color_black(U, o);
            //show_tree("case1-fin", (Ast *)P, 1);
            return 2; // popping the stack twice to get to lbl1
        }
        // have an inward-leaning chain P -dp-> C -d-> N of red nodes
        if(d < 0) { // rotate right
            //show_tree("case2", (Ast *)P, 0);
            set_right(P, N, o);
            set_left(C, get_right(N, o), o);
            set_right(N, C, o);
            U = C; C = N; N = U; // swap(C, N);
            //show_tree("case2-fin", (Ast *)P, 1);
            // d = -d;
        } // left-rotate G - P - C
        // have an outward-leaning chain G -dg-> P -dp-> C -d-> N of red nodes
        //show_tree("case3", (Ast *)P, 0);
        color_black(C, o);
        color_red(P, o);
        if(dg < 0) set_left(G, C, o);
        else      set_right(G, C, o);
        set_right(P, get_left(C, o), o);
        set_left(C, P, o);
        //show_tree("case3-fin", (Ast *)C, 1);
    }
    return 0;
}

// Setup node as root of a new tree.
void new_tree(void *N, const rbop_t *o) {
    color_black(N, o);
    set_left(N, o->nil, o);
    set_right(N, o->nil, o);
}

// Store N in ret and walk up the stack with the replacement.
// This ensures the node is consistently replaced everywhere.
void replace_node(rbtrav_t *ret, void *N, const rbop_t *o) {
    //printf("Replacing %s with %s\n", ((Ast *)ret->N)->name,
    //                    ((Ast *)N)->name);

    // replace link with N
    if(ret->d < 0) set_left(ret->up->N, N, o);
    else          set_right(ret->up->N, N, o);

    set_mask(N, get_mask(ret->N, o), o); // put N in ret's place
    set_left(N, get_left(ret->N, o), o);
    set_right(N, get_right(ret->N, o), o);
    ret->N = N;
}

/*
void print_stack(rbtrav_t *self, const rbop_t *o) {
    for(; self->up != NULL; self=self->up) {
        if(self->N == o->nil) {
          printf("nil");
        } else {
          printf(" %s%c(%d)", ((Ast *)self->N)->name,
                ' '+(get_mask(self->N, o)>0)*('*'-' '), self->d);
        }
    }
    printf("\n");
}*/

/* deletion helper routine for node with <2 children
 * Note that these maintain the red/black property, and
 * so are not useful for simple/final listings.
 *
 * modes are:
 *  a) ret is location to put removed node
 *     A is comparison node key
 *     dir = 0 / undefined
 *  b) ret is location of node to replace with leaf
 *     A is nil
 *     dir = -1 (find smallest leaf) or +1 (find largest)
 *  
 *  for both modes, self points to the input stack
 *  {current node, direction taken to get there, and parent stack}
 *  and o, of course, defines the data structure.
 */
// smallest <= dir < 0; largest <= dir > 0
#define PARENT (self->up->N)
#define DP (self->up->d)
#define GRAND (self->up->up->N)
static int pop_extreme(void **ret, const void *A, int dir,
        rbtrav_t *self, const rbop_t *o) {
    rbtrav_t next = { .up = self };
    void *C, *S, *SL, *SR;
    int d;

    if(ret == NULL) goto recall;

    if(A != o->nil) { // initial descent to A
        if(self->N == o->nil) { // not found
            *ret = o->nil;
            return 0;
        }
        d = o->cmp(A, self->N);
        if(d != 0) {
            if(d < 0) next.N = get_left(self->N, o);
            else      next.N = get_right(self->N, o);
            next.d = d;
            if(pop_extreme(ret, A, 0, &next, o))
                goto chk;
            return 0;
        }
        //printf("Start replacement mode: ");
        //print_stack(self, o);
        // d = 0 case falls through to replacement mode start
        *ret = self->N; // Return the deleted node in ret.

        // Replace mode will replace P's link to N.
        if(self->d < 0) {
            ret = (void **)self;
        } else {
            ret = (void **)self;
        }
        SL = get_left(self->N, o);
        SR = get_right(self->N, o);
        // Replace nil (if present), then base try order on dp (mostly random).
        if(SL == o->nil) {
            d = get_mask(self->N, o);
            C = SR;
            goto replaceme;
        } else if(SR == o->nil) {
            d = get_mask(self->N, o);
            C = SL;
            goto replaceme;
        } else {
            dir = self->d;
            if(dir < 0) next.N = get_right(self->N, o);
            else        next.N = get_left(self->N, o);
            next.d = -dir; // step opposite min/max direction dir.
            if(pop_extreme(ret, o->nil, dir, &next, o))
                goto chk;
            return 0;
        }
        A = o->nil;
    }

    // A = nil: descend to replacement leaf (dir has been set)
    if(dir < 0) next.N = get_left(self->N, o);
    else       next.N = get_right(self->N, o);

    if(next.N != o->nil) { // not there yet
        next.d = dir;
        if(pop_extreme(ret, o->nil, dir, &next, o))
            goto chk;
        return 0;
    }

    // there: N = extreme element to unlink (no child in direction = dir)
    //  and sits in G -dp-> P -(self->d)-> N (-!dir-> C)?
    if(dir < 0) C = get_right(self->N, o);
    else        C = get_left(self->N, o);

    d = get_mask(self->N, o); // save N's old color
    replace_node((rbtrav_t *)ret, self->N, o); // replace ret
replaceme:
    self->N = C; // replace N with child to continue algo.
    if(self->d < 0) set_left(PARENT, C, o); // unlink N from P
    else           set_right(PARENT, C, o);
    //printf("Found replacement: ");
    //print_stack(self, o);

    if(d) // replaced red node
        return 0;
    if(C != o->nil && get_mask(C, o)) { // replaced black with red node
        color_black(C, o);
        return 0;
    }
    // N isn't used past this point (except as parent from recursive call...)
chk:
    //printf("Checking: ");
    //print_stack(self, o);
    if(self->up->up == NULL) { // N is root.
        return 0;
    }
    if(self->d < 0) S = get_right(PARENT, o);
    else           S =  get_left(PARENT, o);
    if(S == o->nil) {
        fprintf(stderr, "Algorithm error: sibling = nil\n");
        return 0;
    }
    if(S != o->nil && get_mask(S, o)) { // S (red)
        //show_tree("case2.dot", (Ast *)PARENT, 0);
        color_red(PARENT, o); // implies no recursive checking.
        color_black(S, o);
        if(self->d < 0) {
            if(DP < 0) set_left(GRAND, S, o);
            else      set_right(GRAND, S, o);
            SL = get_left(S, o); // new S (black)
            set_left(S, PARENT, o);
            set_right(PARENT, SL, o);
        }
        else {
            if(DP < 0) set_left(GRAND, S, o);
            else      set_right(GRAND, S, o);
            SL = get_right(S, o); // new S (black)
            set_right(S, PARENT, o);
            set_left(PARENT, SL, o);
        }
        //printf("Handled case 2.\n");
        //show_tree("case2-fin.dot", (Ast *)S, 0);
        // adds visit to S in stack by re-writing traversal
        // and calling myself again
        next.N = self->N; // delegate self's work to next level
        // not set? since N isn't used anyway
        next.d = self->d;
        self->N = PARENT; // take over P's work
        // self->d = self->d; (parent is in d direction from S)
        // self->up->d = self->up->d; (S is in parent's direction from G)
        PARENT = S; // visit S instead of G (parent stack still G)

        if(pop_extreme(NULL, o->nil, 0, &next, o)) goto chk;
        return 0;
//recall:
        //printf("Extended stack: ");
        //print_stack(self, o);
        //printf("Resuming (%s)\n", ((Ast *)self->N)->name);
    }
recall: // still checking.  ret, A and dir don't matter at this point.
    if(self->d < 0) S = get_right(PARENT, o); // find new sibling
    else             S = get_left(PARENT, o);

    if(get_mask(S, o)) {
        fprintf(stderr, "Red-black algo. error: S is not black!\n");
        show_tree("err.dot", PARENT, 1);
        return 0;
    }
    if(self->d < 0) {
        SL = get_left(S, o);
        SR = get_right(S, o);
    } else {
        SL = get_right(S, o);
        SR = get_left(S, o);
    }
    if(SR == o->nil || get_mask(SR, o) == 0) {
        color_red(S, o);
        if(SL == o->nil || get_mask(SL, o) == 0) { // S{L,R} (black)
            // S can be red.
            if(get_mask(PARENT, o) == 0) { // P (black)
                //printf("Handled case 3\n");
                return 1;
            } else {
                color_black(PARENT, o); // P (red) -> P (black) done
                //printf("Handled case 4\n");
                return 0;
            }
        } else { // S_L (red) S_R (black)
            //show_tree("case5.dot", (Ast *)PARENT, 0);
            color_black(SL, o);
            if(self->d < 0) {
                set_right(PARENT, SL, o);
                set_left(S, get_right(SL, o), o);
                set_right(SL, S, o);
                SR = S;
                S = SL;
                SL = get_left(S, o);
            } else {
                set_left(PARENT, SL, o);
                set_right(S, get_left(SL, o), o);
                set_left(SL, S, o);
                SR = S;
                S = SL;
                SL = get_right(S, o);
            }
            //show_tree("case5-fin.dot", (Ast *)PARENT, 0);
            //printf("Handled case 5\n");
        }
    }
    if(get_mask(SR, o) == 0) {
        fprintf(stderr, "Algorithm error: SR is not red!\n");
        return 0;
    }
    //show_tree("case6.dot", (Ast *)PARENT, 0);
    //printf("Stack at start of case6: ");
    //print_stack(self, o);
    // rotate left at P
    if(self->d < 0) {
        if(DP < 0) set_left(GRAND, S, o);
        else       set_right(GRAND, S, o);
        set_right(PARENT, SL, o);
        set_left(S, PARENT, o);
    } else {
        if(DP < 0) set_left(GRAND, S, o);
        else       set_right(GRAND, S, o);
        set_left(PARENT, SL, o);
        set_right(S, PARENT, o);
    }
    set_mask(S, get_mask(PARENT, o), o);
    color_black(PARENT, o);
    color_black(SR, o);
    //printf("Handled case 6\n");
    //show_tree("case6-fin.dot", (Ast *)S, 0);

    return 0;
}

/* Returns the node if deleted,
 * nil if not present
 */
void *del_node(void **N, const void *A, const rbop_t *o) {
    rbtrav_t prev, next;
    void *R;
   
    if(*N == o->nil || *N == NULL) return o->nil;

    // Install a sentinel, hacked to look like root is its left child.
    prev.up = NULL;
    prev.N = (void *)N - o->coff;
    prev.d = 0;

    next.up = &prev;
    next.N = *N;
    next.d = -1;

    pop_extreme(&R, A, 0, &next, o);
    return R;
}

