/*
This file is part of ``jkdtree'', a library for working with jkd-trees.
Copyright (C) 2007-2011 John Tsiombikas <nuclear@member.fsf.org>

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
3. The name of the author may not be used to endorse or promote products
   derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO
EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY
OF SUCH DAMAGE.
*/
typedef float pix_t;

#ifndef _JKDTREE_H_
#define _JKDTREE_H_

#ifdef __cplusplus
extern "C" {
#endif

struct jkdtree;
struct jkdres;


/* create a jkd-tree for "k"-dimensional data */
struct jkdtree *jkd_create(int k);
/* free the struct jkdtree */
void jkd_free(struct jkdtree *tree);

/* remove all the elements from the tree */
void jkd_clear(struct jkdtree *tree);

/* if called with non-null 2nd argument, the function provided
 * will be called on data pointers (see jkd_insert) when nodes
 * are to be removed from the tree.
 */
void jkd_data_destructor(struct jkdtree *tree, void (*destr)(void*));

/* insert a node, specifying its position, and optional data */
int jkd_insert(struct jkdtree *tree, const double *pos, pix_t data);
int jkd_insertf(struct jkdtree *tree, const float *pos, pix_t data);
int jkd_insert3(struct jkdtree *tree, double x, double y, double z, pix_t data);
int jkd_insert3f(struct jkdtree *tree, float x, float y, float z, pix_t data);

/* Find the nearest node from a given point.
 *
 * This function returns a pointer to a result set with at most one element.
 */
struct jkdres *jkd_nearest(struct jkdtree *tree, const double *pos);
struct jkdres *jkd_nearestf(struct jkdtree *tree, const float *pos);
struct jkdres *jkd_nearest3(struct jkdtree *tree, double x, double y, double z);
struct jkdres *jkd_nearest3f(struct jkdtree *tree, float x, float y, float z);

/* Find the N nearest nodes from a given point.
 *
 * This function returns a pointer to a result set, with at most N elements,
 * which can be manipulated with the jkd_res_* functions.
 * The returned pointer can be null as an indication of an error. Otherwise
 * a valid result set is always returned which may contain 0 or more elements.
 * The result set must be deallocated with jkd_res_free after use.
 */
/*
struct jkdres *jkd_nearest_n(struct jkdtree *tree, const double *pos, int num);
struct jkdres *jkd_nearest_nf(struct jkdtree *tree, const float *pos, int num);
struct jkdres *jkd_nearest_n3(struct jkdtree *tree, double x, double y, double z);
struct jkdres *jkd_nearest_n3f(struct jkdtree *tree, float x, float y, float z);
*/

/* Find any nearest nodes from a given point within a range.
 *
 * This function returns a pointer to a result set, which can be manipulated
 * by the jkd_res_* functions.
 * The returned pointer can be null as an indication of an error. Otherwise
 * a valid result set is always returned which may contain 0 or more elements.
 * The result set must be deallocated with jkd_res_free after use.
 */
struct jkdres *jkd_nearest_range(struct jkdtree *tree, const double *pos, double range);
struct jkdres *jkd_nearest_rangef(struct jkdtree *tree, const float *pos, float range);
struct jkdres *jkd_nearest_range3(struct jkdtree *tree, double x, double y, double z, double range);
struct jkdres *jkd_nearest_range3f(struct jkdtree *tree, float x, float y, float z, float range);

/* frees a result set returned by jkd_nearest_range() */
void jkd_res_free(struct jkdres *set);

/* returns the size of the result set (in elements) */
int jkd_res_size(struct jkdres *set);

/* rewinds the result set iterator */
void jkd_res_rewind(struct jkdres *set);

/* returns non-zero if the set iterator reached the end after the last element */
int jkd_res_end(struct jkdres *set);

/* advances the result set iterator, returns non-zero on success, zero if
 * there are no more elements in the result set.
 */
int jkd_res_next(struct jkdres *set);

/* returns the data pointer (can be null) of the current result set item
 * and optionally sets its position to the pointers(s) if not null.
 */
pix_t jkd_res_item(struct jkdres *set, double *pos);
pix_t jkd_res_itemf(struct jkdres *set, float *pos);
void *jkd_res_item3(struct jkdres *set, double *x, double *y, double *z);
void *jkd_res_item3f(struct jkdres *set, float *x, float *y, float *z);

/* equivalent to jkd_res_item(set, 0) */
pix_t jkd_res_item_data(struct jkdres *set);


#ifdef __cplusplus
}
#endif

#endif	/* _JKDTREE_H_ */
