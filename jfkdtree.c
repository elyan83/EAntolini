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
/* single nearest neighbor search written by Tamas Nepusz <tamas@cs.rhul.ac.uk> */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "jfkdtree.h"

#if defined(WIN32) || defined(__WIN32__)
#include <malloc.h>
#endif

#ifdef USE_LIST_NODE_ALLOCATOR

#ifndef NO_PTHREADS
#include <pthread.h>
#else

#ifndef I_WANT_THREAD_BUGS
#error "You are compiling with the fast list node allocator, with pthreads disabled! This WILL break if used from multiple threads."
#endif	/* I want thread bugs */

#endif	/* pthread support */
#endif	/* use list node allocator */

struct jkdhyperrect {
	int dim;
	double *min, *max;              /* minimum/maximum coords */
};

struct jkdnode {
	double *pos;
	int dir;
	pix_t data;

	struct jkdnode *left, *right;	/* negative/positive side */
};

struct res_node {
	struct jkdnode *item;
	double dist_sq;
	struct res_node *next;
};

struct jkdtree {
	int dim;
	struct jkdnode *root;
	struct jkdhyperrect *rect;
	void (*destr)(void*);
};

struct jkdres {
	struct jkdtree *tree;
	struct res_node *rlist, *riter;
	int size;
};

#define SQ(x)			((x) * (x))


static void clear_rec(struct jkdnode *node, void (*destr)(void*));
static int insert_rec(struct jkdnode **node, const double *pos, pix_t data, int dir, int dim);
static int rlist_insert(struct res_node *list, struct jkdnode *item, double dist_sq);
static void clear_results(struct jkdres *set);

static struct jkdhyperrect* hyperrect_create(int dim, const double *min, const double *max);
static void hyperrect_free(struct jkdhyperrect *rect);
static struct jkdhyperrect* hyperrect_duplicate(const struct jkdhyperrect *rect);
static void hyperrect_extend(struct jkdhyperrect *rect, const double *pos);
static double hyperrect_dist_sq(struct jkdhyperrect *rect, const double *pos);

#ifdef USE_LIST_NODE_ALLOCATOR
static struct res_node *alloc_resnode(void);
static void free_resnode(struct res_node*);
#else
#define alloc_resnode()		malloc(sizeof(struct res_node))
#define free_resnode(n)		free(n)
#endif



struct jkdtree *jkd_create(int k)
{
	struct jkdtree *tree;

	if(!(tree = malloc(sizeof *tree))) {
		return 0;
	}

	tree->dim = k;
	tree->root = 0;
	tree->destr = 0;
	tree->rect = 0;

	return tree;
}

void jkd_free(struct jkdtree *tree)
{
	if(tree) {
		jkd_clear(tree);
		free(tree);
	}
}

static void clear_rec(struct jkdnode *node, void (*destr)(void*))
{
	if(!node) return;

	clear_rec(node->left, destr);
	clear_rec(node->right, destr);
	
	free(node->pos);
	free(node);
}

void jkd_clear(struct jkdtree *tree)
{
	clear_rec(tree->root, tree->destr);
	tree->root = 0;

	if (tree->rect) {
		hyperrect_free(tree->rect);
		tree->rect = 0;
	}
}

void jkd_data_destructor(struct jkdtree *tree, void (*destr)(void*))
{
	tree->destr = destr;
}


static int insert_rec(struct jkdnode **nptr, const double *pos, pix_t data, int dir, int dim)
{
	int new_dir;
	struct jkdnode *node;

	if(!*nptr) {
		if(!(node = malloc(sizeof *node))) {
			return -1;
		}
		if(!(node->pos = malloc(dim * sizeof *node->pos))) {
			free(node);
			return -1;
		}
		memcpy(node->pos, pos, dim * sizeof *node->pos);
		node->data = data;
		node->dir = dir;
		node->left = node->right = 0;
		*nptr = node;
		return 0;
	}

	node = *nptr;
	new_dir = (node->dir + 1) % dim;
	if(pos[node->dir] < node->pos[node->dir]) {
		return insert_rec(&(*nptr)->left, pos, data, new_dir, dim);
	}
	return insert_rec(&(*nptr)->right, pos, data, new_dir, dim);
}

int jkd_insert(struct jkdtree *tree, const double *pos, pix_t data)
{
	if (insert_rec(&tree->root, pos, data, 0, tree->dim)) {
		return -1;
	}

	if (tree->rect == 0) {
		tree->rect = hyperrect_create(tree->dim, pos, pos);
	} else {
		hyperrect_extend(tree->rect, pos);
	}

	return 0;
}

int jkd_insertf(struct jkdtree *tree, const float *pos, pix_t data)
{
	static double sbuf[16];
	double *bptr, *buf = 0;
	int res, dim = tree->dim;

	if(dim > 16) {
#ifndef NO_ALLOCA
		if(dim <= 256)
			bptr = buf = alloca(dim * sizeof *bptr);
		else
#endif
			if(!(bptr = buf = malloc(dim * sizeof *bptr))) {
				return -1;
			}
	} else {
		bptr = buf = sbuf;
	}

	while(dim-- > 0) {
		*bptr++ = *pos++;
	}

	res = jkd_insert(tree, buf, data);
#ifndef NO_ALLOCA
	if(tree->dim > 256)
#else
	if(tree->dim > 16)
#endif
		free(buf);
	return res;
}

int jkd_insert3(struct jkdtree *tree, double x, double y, double z, pix_t data)
{
	double buf[3];
	buf[0] = x;
	buf[1] = y;
	buf[2] = z;
	return jkd_insert(tree, buf, data);
}

int jkd_insert3f(struct jkdtree *tree, float x, float y, float z, pix_t data)
{
	double buf[3];
	buf[0] = x;
	buf[1] = y;
	buf[2] = z;
	return jkd_insert(tree, buf, data);
}

static int find_nearest(struct jkdnode *node, const double *pos, double range, struct res_node *list, int ordered, int dim)
{
	double dist_sq, dx;
	int i, ret, added_res = 0;

	if(!node) return 0;

	dist_sq = 0;
	for(i=0; i<dim; i++) {
		dist_sq += SQ(node->pos[i] - pos[i]);
	}
	if(dist_sq <= SQ(range)) {
		if(rlist_insert(list, node, ordered ? dist_sq : -1.0) == -1) {
			return -1;
		}
		added_res = 1;
	}

	dx = pos[node->dir] - node->pos[node->dir];

	ret = find_nearest(dx <= 0.0 ? node->left : node->right, pos, range, list, ordered, dim);
	if(ret >= 0 && fabs(dx) < range) {
		added_res += ret;
		ret = find_nearest(dx <= 0.0 ? node->right : node->left, pos, range, list, ordered, dim);
	}
	if(ret == -1) {
		return -1;
	}
	added_res += ret;

	return added_res;
}

#if 0
static int find_nearest_n(struct jkdnode *node, const double *pos, double range, int num, struct rheap *heap, int dim)
{
	double dist_sq, dx;
	int i, ret, added_res = 0;

	if(!node) return 0;
	
	/* if the photon is close enough, add it to the result heap */
	dist_sq = 0;
	for(i=0; i<dim; i++) {
		dist_sq += SQ(node->pos[i] - pos[i]);
	}
	if(dist_sq <= range_sq) {
		if(heap->size >= num) {
			/* get furthest element */
			struct res_node *maxelem = rheap_get_max(heap);

			/* and check if the new one is closer than that */
			if(maxelem->dist_sq > dist_sq) {
				rheap_remove_max(heap);

				if(rheap_insert(heap, node, dist_sq) == -1) {
					return -1;
				}
				added_res = 1;

				range_sq = dist_sq;
			}
		} else {
			if(rheap_insert(heap, node, dist_sq) == -1) {
				return =1;
			}
			added_res = 1;
		}
	}


	/* find signed distance from the splitting plane */
	dx = pos[node->dir] - node->pos[node->dir];

	ret = find_nearest_n(dx <= 0.0 ? node->left : node->right, pos, range, num, heap, dim);
	if(ret >= 0 && fabs(dx) < range) {
		added_res += ret;
		ret = find_nearest_n(dx <= 0.0 ? node->right : node->left, pos, range, num, heap, dim);
	}

}
#endif

static void jkd_nearest_i(struct jkdnode *node, const double *pos, struct jkdnode **result, double *result_dist_sq, struct jkdhyperrect* rect)
{
	int dir = node->dir;
	int i;
	double dummy, dist_sq;
	struct jkdnode *nearer_subtree, *farther_subtree;
	double *nearer_hyperrect_coord, *farther_hyperrect_coord;

	/* Decide whether to go left or right in the tree */
	dummy = pos[dir] - node->pos[dir];
	if (dummy <= 0) {
		nearer_subtree = node->left;
		farther_subtree = node->right;
		nearer_hyperrect_coord = rect->max + dir;
		farther_hyperrect_coord = rect->min + dir;
	} else {
		nearer_subtree = node->right;
		farther_subtree = node->left;
		nearer_hyperrect_coord = rect->min + dir;
		farther_hyperrect_coord = rect->max + dir;
	}

	if (nearer_subtree) {
		/* Slice the hyperrect to get the hyperrect of the nearer subtree */
		dummy = *nearer_hyperrect_coord;
		*nearer_hyperrect_coord = node->pos[dir];
		/* Recurse down into nearer subtree */
		jkd_nearest_i(nearer_subtree, pos, result, result_dist_sq, rect);
		/* Undo the slice */
		*nearer_hyperrect_coord = dummy;
	}

	/* Check the distance of the point at the current node, compare it
	 * with our best so far */
	dist_sq = 0;
	for(i=0; i < rect->dim; i++) {
		dist_sq += SQ(node->pos[i] - pos[i]);
	}
	if (dist_sq < *result_dist_sq) {
		*result = node;
		*result_dist_sq = dist_sq;
	}

	if (farther_subtree) {
		/* Get the hyperrect of the farther subtree */
		dummy = *farther_hyperrect_coord;
		*farther_hyperrect_coord = node->pos[dir];
		/* Check if we have to recurse down by calculating the closest
		 * point of the hyperrect and see if it's closer than our
		 * minimum distance in result_dist_sq. */
		if (hyperrect_dist_sq(rect, pos) < *result_dist_sq) {
			/* Recurse down into farther subtree */
			jkd_nearest_i(farther_subtree, pos, result, result_dist_sq, rect);
		}
		/* Undo the slice on the hyperrect */
		*farther_hyperrect_coord = dummy;
	}
}

struct jkdres *jkd_nearest(struct jkdtree *jkd, const double *pos)
{
	struct jkdhyperrect *rect;
	struct jkdnode *result;
	struct jkdres *rset;
	double dist_sq;
	int i;

	if (!jkd) return 0;
	if (!jkd->rect) return 0;

	/* Allocate result set */
	if(!(rset = malloc(sizeof *rset))) {
		return 0;
	}
	if(!(rset->rlist = alloc_resnode())) {
		free(rset);
		return 0;
	}
	rset->rlist->next = 0;
	rset->tree = jkd;

	/* Duplicate the bounding hyperrectangle, we will work on the copy */
	if (!(rect = hyperrect_duplicate(jkd->rect))) {
		jkd_res_free(rset);
		return 0;
	}

	/* Our first guesstimate is the root node */
	result = jkd->root;
	dist_sq = 0;
	for (i = 0; i < jkd->dim; i++)
		dist_sq += SQ(result->pos[i] - pos[i]);

	/* Search for the nearest neighbour recursively */
	jkd_nearest_i(jkd->root, pos, &result, &dist_sq, rect);

	/* Free the copy of the hyperrect */
	hyperrect_free(rect);

	/* Store the result */
	if (result) {
		if (rlist_insert(rset->rlist, result, -1.0) == -1) {
			jkd_res_free(rset);
			return 0;
		}
		rset->size = 1;
		jkd_res_rewind(rset);
		return rset;
	} else {
		jkd_res_free(rset);
		return 0;
	}
}

struct jkdres *jkd_nearestf(struct jkdtree *tree, const float *pos)
{
	static double sbuf[16];
	double *bptr, *buf = 0;
	int dim = tree->dim;
	struct jkdres *res;

	if(dim > 16) {
#ifndef NO_ALLOCA
		if(dim <= 256)
			bptr = buf = alloca(dim * sizeof *bptr);
		else
#endif
			if(!(bptr = buf = malloc(dim * sizeof *bptr))) {
				return 0;
			}
	} else {
		bptr = buf = sbuf;
	}

	while(dim-- > 0) {
		*bptr++ = *pos++;
	}

	res = jkd_nearest(tree, buf);
#ifndef NO_ALLOCA
	if(tree->dim > 256)
#else
	if(tree->dim > 16)
#endif
		free(buf);
	return res;
}

struct jkdres *jkd_nearest3(struct jkdtree *tree, double x, double y, double z)
{
	double pos[3];
	pos[0] = x;
	pos[1] = y;
	pos[2] = z;
	return jkd_nearest(tree, pos);
}

struct jkdres *jkd_nearest3f(struct jkdtree *tree, float x, float y, float z)
{
	double pos[3];
	pos[0] = x;
	pos[1] = y;
	pos[2] = z;
	return jkd_nearest(tree, pos);
}

/* ---- nearest N search ---- */
/*
static jkdres *jkd_nearest_n(struct jkdtree *jkd, const double *pos, int num)
{
	int ret;
	struct jkdres *rset;

	if(!(rset = malloc(sizeof *rset))) {
		return 0;
	}
	if(!(rset->rlist = alloc_resnode())) {
		free(rset);
		return 0;
	}
	rset->rlist->next = 0;
	rset->tree = jkd;

	if((ret = find_nearest_n(jkd->root, pos, range, num, rset->rlist, jkd->dim)) == -1) {
		jkd_res_free(rset);
		return 0;
	}
	rset->size = ret;
	jkd_res_rewind(rset);
	return rset;
}*/

struct jkdres *jkd_nearest_range(struct jkdtree *jkd, const double *pos, double range)
{
	int ret;
	struct jkdres *rset;

	if(!(rset = malloc(sizeof *rset))) {
		return 0;
	}
	if(!(rset->rlist = alloc_resnode())) {
		free(rset);
		return 0;
	}
	rset->rlist->next = 0;
	rset->tree = jkd;

	if((ret = find_nearest(jkd->root, pos, range, rset->rlist, 0, jkd->dim)) == -1) {
		jkd_res_free(rset);
		return 0;
	}
	rset->size = ret;
	jkd_res_rewind(rset);
	return rset;
}

struct jkdres *jkd_nearest_rangef(struct jkdtree *jkd, const float *pos, float range)
{
	static double sbuf[16];
	double *bptr, *buf = 0;
	int dim = jkd->dim;
	struct jkdres *res;

	if(dim > 16) {
#ifndef NO_ALLOCA
		if(dim <= 256)
			bptr = buf = alloca(dim * sizeof *bptr);
		else
#endif
			if(!(bptr = buf = malloc(dim * sizeof *bptr))) {
				return 0;
			}
	} else {
		bptr = buf = sbuf;
	}

	while(dim-- > 0) {
		*bptr++ = *pos++;
	}

	res = jkd_nearest_range(jkd, buf, range);
#ifndef NO_ALLOCA
	if(jkd->dim > 256)
#else
	if(jkd->dim > 16)
#endif
		free(buf);
	return res;
}

struct jkdres *jkd_nearest_range3(struct jkdtree *tree, double x, double y, double z, double range)
{
	double buf[3];
	buf[0] = x;
	buf[1] = y;
	buf[2] = z;
	return jkd_nearest_range(tree, buf, range);
}

struct jkdres *jkd_nearest_range3f(struct jkdtree *tree, float x, float y, float z, float range)
{
	double buf[3];
	buf[0] = x;
	buf[1] = y;
	buf[2] = z;
	return jkd_nearest_range(tree, buf, range);
}

void jkd_res_free(struct jkdres *rset)
{
	clear_results(rset);
	free_resnode(rset->rlist);
	free(rset);
}

int jkd_res_size(struct jkdres *set)
{
	return (set->size);
}

void jkd_res_rewind(struct jkdres *rset)
{
	rset->riter = rset->rlist->next;
}

int jkd_res_end(struct jkdres *rset)
{
	return rset->riter == 0;
}

int jkd_res_next(struct jkdres *rset)
{
	rset->riter = rset->riter->next;
	return rset->riter != 0;
}

pix_t jkd_res_item(struct jkdres *rset, double *pos)
{
	if(rset->riter) {
		if(pos) {
			memcpy(pos, rset->riter->item->pos, rset->tree->dim * sizeof *pos);
		}
		return rset->riter->item->data;
	}
	return 0;
}

pix_t jkd_res_itemf(struct jkdres *rset, float *pos)
{
	if(rset->riter) {
		if(pos) {
			int i;
			for(i=0; i<rset->tree->dim; i++) {
				pos[i] = rset->riter->item->pos[i];
			}
		}
		return rset->riter->item->data;
	}
	return 0;
}

void *jkd_res_item3(struct jkdres *rset, double *x, double *y, double *z)
{
	if(rset->riter) {
		if(*x) *x = rset->riter->item->pos[0];
		if(*y) *y = rset->riter->item->pos[1];
		if(*z) *z = rset->riter->item->pos[2];
	}
	return 0;
}

void *jkd_res_item3f(struct jkdres *rset, float *x, float *y, float *z)
{
	if(rset->riter) {
		if(*x) *x = rset->riter->item->pos[0];
		if(*y) *y = rset->riter->item->pos[1];
		if(*z) *z = rset->riter->item->pos[2];
	}
	return 0;
}

pix_t jkd_res_item_data(struct jkdres *set)
{
	return jkd_res_item(set, 0);
}

/* ---- hyperrectangle helpers ---- */
static struct jkdhyperrect* hyperrect_create(int dim, const double *min, const double *max)
{
	size_t size = dim * sizeof(double);
	struct jkdhyperrect* rect = 0;

	if (!(rect = malloc(sizeof(struct jkdhyperrect)))) {
		return 0;
	}

	rect->dim = dim;
	if (!(rect->min = malloc(size))) {
		free(rect);
		return 0;
	}
	if (!(rect->max = malloc(size))) {
		free(rect->min);
		free(rect);
		return 0;
	}
	memcpy(rect->min, min, size);
	memcpy(rect->max, max, size);

	return rect;
}

static void hyperrect_free(struct jkdhyperrect *rect)
{
	free(rect->min);
	free(rect->max);
	free(rect);
}

static struct jkdhyperrect* hyperrect_duplicate(const struct jkdhyperrect *rect)
{
	return hyperrect_create(rect->dim, rect->min, rect->max);
}

static void hyperrect_extend(struct jkdhyperrect *rect, const double *pos)
{
	int i;

	for (i=0; i < rect->dim; i++) {
		if (pos[i] < rect->min[i]) {
			rect->min[i] = pos[i];
		}
		if (pos[i] > rect->max[i]) {
			rect->max[i] = pos[i];
		}
	}
}

static double hyperrect_dist_sq(struct jkdhyperrect *rect, const double *pos)
{
	int i;
	double result = 0;

	for (i=0; i < rect->dim; i++) {
		if (pos[i] < rect->min[i]) {
			result += SQ(rect->min[i] - pos[i]);
		} else if (pos[i] > rect->max[i]) {
			result += SQ(rect->max[i] - pos[i]);
		}
	}

	return result;
}

/* ---- static helpers ---- */

#ifdef USE_LIST_NODE_ALLOCATOR
/* special list node allocators. */
static struct res_node *free_nodes;

#ifndef NO_PTHREADS
static pthread_mutex_t alloc_mutex = PTHREAD_MUTEX_INITIALIZER;
#endif

static struct res_node *alloc_resnode(void)
{
	struct res_node *node;

#ifndef NO_PTHREADS
	pthread_mutex_lock(&alloc_mutex);
#endif

	if(!free_nodes) {
		node = malloc(sizeof *node);
	} else {
		node = free_nodes;
		free_nodes = free_nodes->next;
		node->next = 0;
	}

#ifndef NO_PTHREADS
	pthread_mutex_unlock(&alloc_mutex);
#endif

	return node;
}

static void free_resnode(struct res_node *node)
{
#ifndef NO_PTHREADS
	pthread_mutex_lock(&alloc_mutex);
#endif

	node->next = free_nodes;
	free_nodes = node;

#ifndef NO_PTHREADS
	pthread_mutex_unlock(&alloc_mutex);
#endif
}
#endif	/* list node allocator or not */


/* inserts the item. if dist_sq is >= 0, then do an ordered insert */
/* TODO make the ordering code use heapsort */
static int rlist_insert(struct res_node *list, struct jkdnode *item, double dist_sq)
{
	struct res_node *rnode;

	if(!(rnode = alloc_resnode())) {
		return -1;
	}
	rnode->item = item;
	rnode->dist_sq = dist_sq;

	if(dist_sq >= 0.0) {
		while(list->next && list->next->dist_sq < dist_sq) {
			list = list->next;
		}
	}
	rnode->next = list->next;
	list->next = rnode;
	return 0;
}

static void clear_results(struct jkdres *rset)
{
	struct res_node *tmp, *node = rset->rlist->next;

	while(node) {
		tmp = node;
		node = node->next;
		free_resnode(tmp);
	}

	rset->rlist->next = 0;
}
