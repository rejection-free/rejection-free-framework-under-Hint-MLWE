#include "lazer.h"

/* poly.c */
void
poly_addscale2 (poly_t r, const int_t a, poly_t b, int crt)
{
  polyring_srcptr ring = poly_get_ring (r);
  poly_t tmp;

  poly_alloc (tmp, ring);

  poly_scale (tmp, a, b);
  poly_add (r, r, tmp, crt);

  poly_free (tmp);
}

void
poly_subscale2 (poly_t r, const int_t a, poly_t b, int crt)
{
  polyring_srcptr ring = poly_get_ring (r);
  poly_t tmp;

  poly_alloc (tmp, ring);

  poly_scale (tmp, a, b);
  poly_sub (r, r, tmp, crt);

  poly_free (tmp);
}

/* spolyvec.c */

void
spolyvec_get_subvec (spolyvec_t subvec, spolyvec_t vec, unsigned int elem_init,
                    unsigned int nelems)
{
  poly_ptr poly, poly2;
  unsigned int i, k =0;

  subvec->nelems = 0;
  for (i = 0; i < nelems; i++) {
    poly = spolyvec_get_elem (vec, elem_init+i);

    if (poly != NULL)
    {
      poly2 = spolyvec_insert_elem (subvec, i);
      poly_set (poly2, poly);
      k++;
    }
  }
  subvec->sorted = 1;
  subvec->nelems = k;
}

/* spolymat.c */

/* r != b */
void
spolymat_neg (spolymat_t r, spolymat_t b)
{
  poly_ptr ri, bi;
  unsigned int brow, bcol, i;

  ASSERT_ERR (r->nelems_max >= b->nelems_max);

  r->nelems = 0;
  _SMAT_FOREACH_ELEM (b, i)
  {
    bi = spolymat_get_elem (b, i);
    brow = spolymat_get_row (b, i);
    bcol = spolymat_get_col (b, i);

    ri = spolymat_insert_elem (r, brow, bcol);
    poly_neg(ri, bi);
  }
  r->nelems = b->nelems;
  r->sorted = b->sorted;
}

void
spolymat_get_submat (spolymat_ptr r, spolymat_ptr a, const unsigned int nrows_start, const unsigned int ncols_start)
{
  poly_ptr poly, poly2;
  unsigned int i, arow, acol,k =0;

  ASSERT_ERR (r->nelems_max >= nrows_output * ncols_output);
  ASSERT_ERR (r->ring == a->ring);
  ASSERT_ERR (r->nrows + nrows_start <= a->nrows);
  ASSERT_ERR (r->ncols + ncols_start <= a->ncols);
  ASSERT_ERR (a->sorted);

  r->nelems = 0;
  _SMAT_FOREACH_ELEM (a, i)
  {
    poly = spolymat_get_elem (a, i);
    arow = spolymat_get_row (a, i);
    acol = spolymat_get_col (a, i);
    if (arow >= nrows_start && acol >= ncols_start)
    {
      poly2 = spolymat_insert_elem (r, arow - nrows_start, acol - ncols_start);
      poly_set (poly2, poly);
      k++;
    }
  }

  r->sorted = 1;
  r->nelems = k;
}

void
spolymat_get_submat_upperdiag (spolymat_ptr r, spolymat_ptr a, const unsigned int nrows_start, const unsigned int ncols_start)
{
  poly_ptr poly, poly2;
  unsigned int i, arow, acol, k=0;

  ASSERT_ERR (r->ring == a->ring);
  ASSERT_ERR (r->nrows + nrows_start <= a->nrows);
  ASSERT_ERR (r->ncols + ncols_start <= a->ncols);
  ASSERT_ERR (a->sorted);

  r->nelems = 0;
  _SMAT_FOREACH_ELEM (a, i)
  {
    poly = spolymat_get_elem (a, i);
    arow = spolymat_get_row (a, i);
    acol = spolymat_get_col (a, i);
    if (arow >= nrows_start && acol >= ncols_start && acol >= arow)
    {
      poly2 = spolymat_insert_elem (r, arow - nrows_start, acol - ncols_start);
      poly_set (poly2, poly);
      k++;
    }
  }
  r->sorted = 1;
  r->nelems = k;
}

void
spolymat_set_squared_block_matrix (spolymat_ptr r, spolymat_ptr a1, spolymat_ptr a2, spolymat_ptr a3, spolymat_ptr a4)
{
  poly_ptr poly, poly2;
  unsigned int i, row, col;

  ASSERT_ERR (r->nelems_max >= a1->nelems + a2->nelems + a3->nelems + a4->nelems);
  ASSERT_ERR (r->ring == a->ring);
  ASSERT_ERR (r->nrows == a1->nrows + a3->nrows);
  ASSERT_ERR (r->ncols == a1->ncols + a2->ncols);
  ASSERT_ERR (a1->ncols ==  a3->ncols);
  ASSERT_ERR (a2->ncols ==  a4->ncols);
  ASSERT_ERR (a1->nrows ==  a2->nrows);
  ASSERT_ERR (a4->nrows ==  a3->nrows);
  ASSERT_ERR (a1->sorted);
  ASSERT_ERR (a2->sorted);
  ASSERT_ERR (a3->sorted);
  ASSERT_ERR (a4->sorted);

  r->nelems = 0;
  _SMAT_FOREACH_ELEM (a1, i)
  {
    poly = spolymat_get_elem (a1, i);
    row = spolymat_get_row (a1, i);
    col = spolymat_get_col (a1, i);

    poly2 = spolymat_insert_elem (r, row, col);
    poly_set (poly2, poly);
  }

  _SMAT_FOREACH_ELEM (a2, i)
  {
    poly = spolymat_get_elem (a2, i);
    row = spolymat_get_row (a2, i);
    col = spolymat_get_col (a2, i);

    poly2 = spolymat_insert_elem (r, row, a1->ncols + col);
    poly_set (poly2, poly);
  }

  _SMAT_FOREACH_ELEM (a3, i)
  {
    poly = spolymat_get_elem (a3, i);
    row = spolymat_get_row (a3, i);
    col = spolymat_get_col (a3, i);

    poly2 = spolymat_insert_elem (r, a1->nrows+row, col);
    poly_set (poly2, poly);
  }

  _SMAT_FOREACH_ELEM (a4, i)
  {
    poly = spolymat_get_elem (a4, i);
    row = spolymat_get_row (a4, i);
    col = spolymat_get_col (a4, i);

    poly2 = spolymat_insert_elem (r, a1->nrows+row, a1->ncols + col);
    poly_set (poly2, poly);
  }
  r->sorted = 1;
  r->nelems = a1->nelems + a2->nelems + a3->nelems + a4->nelems;
}

void spolymat_set_upper_block_matrix (spolymat_ptr r, spolymat_ptr a1, spolymat_ptr a2, spolymat_ptr a3){
  poly_ptr poly, poly2;
  unsigned int i, row, col;

  ASSERT_ERR (r->nelems_max >= a1->nelems + a2->nelems + a3->nelems);
  ASSERT_ERR (r->ring == a->ring);
  ASSERT_ERR (r->nrows == a1->nrows + a3->nrows);
  ASSERT_ERR (r->ncols == a1->ncols + a2->ncols);
  ASSERT_ERR (a2->ncols ==  a3->ncols);
  ASSERT_ERR (a1->nrows ==  a2->nrows);
  ASSERT_ERR (a1->sorted);
  ASSERT_ERR (a2->sorted);
  ASSERT_ERR (a3->sorted);

  r->nelems = 0;
  _SMAT_FOREACH_ELEM (a1, i)
  {
    poly = spolymat_get_elem (a1, i);
    row = spolymat_get_row (a1, i);
    col = spolymat_get_col (a1, i);

    poly2 = spolymat_insert_elem (r, row, col);
    poly_set (poly2, poly);
  }

  _SMAT_FOREACH_ELEM (a2, i)
  {
    poly = spolymat_get_elem (a2, i);
    row = spolymat_get_row (a2, i);
    col = spolymat_get_col (a2, i);

    poly2 = spolymat_insert_elem (r, row, a1->ncols+col);
    poly_set (poly2, poly);
  }

  _SMAT_FOREACH_ELEM (a3, i)
  {
    poly = spolymat_get_elem (a3, i);
    row = spolymat_get_row (a3, i);
    col = spolymat_get_col (a3, i);

    poly2 = spolymat_insert_elem (r, a1->nrows+row, a1->ncols+col);
    poly_set (poly2, poly);
  }

  r->sorted = 1;
  r->nelems = a1->nelems + a2->nelems + a3->nelems;
}

void spolymat_set_block_matrix (spolymat_ptr r, spolymat_ptr a1, spolymat_ptr a2){
  poly_ptr poly, poly2;
  unsigned int i, row, col;

  ASSERT_ERR (r->nelems_max >= a1->nelems + a2->nelems);
  ASSERT_ERR (r->ring == a->ring);
  ASSERT_ERR (r->nrows == a1->nrows + a2->nrows);
  ASSERT_ERR (r->ncols == a1->ncols);
  ASSERT_ERR (a1->ncols ==  a2->ncols);
  ASSERT_ERR (a1->sorted);
  ASSERT_ERR (a2->sorted);

  r->nelems = 0;
  _SMAT_FOREACH_ELEM (a1, i)
  {
    poly = spolymat_get_elem (a1, i);
    row = spolymat_get_row (a1, i);
    col = spolymat_get_col (a1, i);

    poly2 = spolymat_insert_elem (r, row, col);
    poly_set (poly2, poly);
  }

  _SMAT_FOREACH_ELEM (a2, i)
  {
    poly = spolymat_get_elem (a2, i);
    row = spolymat_get_row (a2, i);
    col = spolymat_get_col (a2, i);

    poly2 = spolymat_insert_elem (r, a1->nrows +row, col);
    poly_set (poly2, poly);
  }

  r->sorted = 1;
  r->nelems = a1->nelems + a2->nelems;
}

/* stopwatch.c */
STOPWATCH_T (stopwatch_modified_quad_many_prove);
STOPWATCH_T (stopwatch_modified_quad_many_verify);

STOPWATCH_T (stopwatch_modified_quad_prove);
STOPWATCH_T (stopwatch_modified_quad_verify);

STOPWATCH_T (stopwatch_modified_quad_eval_prove);
STOPWATCH_T (stopwatch_modified_quad_eval_prove_compute_h);
STOPWATCH_T (stopwatch_modified_quad_eval_verify);
STOPWATCH_T (stopwatch_modified_quad_eval_schwartz_zippel_quad);
STOPWATCH_T (stopwatch_modified_quad_eval_schwartz_zippel_lin);
STOPWATCH_T (stopwatch_modified_quad_eval_schwartz_zippel_const);

void
print_stopwatch_modified_quad_many_prove (unsigned int indent)
{
  STOPWATCH_PRINT (stopwatch_modified_quad_many_prove, STOPWATCH_MSEC, indent);
  print_stopwatch_modified_quad_prove (indent + INCINDENT);
}

void
print_stopwatch_modified_quad_many_verify (unsigned int indent)
{
  STOPWATCH_PRINT (stopwatch_modified_quad_many_verify, STOPWATCH_MSEC, indent);
  print_stopwatch_modified_quad_verify (indent + INCINDENT);
}

void
print_stopwatch_modified_quad_prove (UNUSED unsigned int indent)
{
  STOPWATCH_PRINT (stopwatch_modified_quad_prove, STOPWATCH_MSEC, indent);
}

void
print_stopwatch_modified_quad_verify (UNUSED unsigned int indent)
{
  STOPWATCH_PRINT (stopwatch_modified_quad_verify, STOPWATCH_MSEC, indent);
}

void
print_stopwatch_modified_quad_eval_prove (unsigned int indent)
{
  STOPWATCH_PRINT (stopwatch_modified_quad_eval_prove, STOPWATCH_MSEC, indent);
  STOPWATCH_PRINT (stopwatch_modified_quad_eval_prove_compute_h, STOPWATCH_MSEC,
                   indent + INCINDENT);
  STOPWATCH_PRINT (stopwatch_modified_quad_eval_schwartz_zippel_quad,
                   STOPWATCH_MSEC, indent + INCINDENT);
  STOPWATCH_PRINT (stopwatch_modified_quad_eval_schwartz_zippel_lin, STOPWATCH_MSEC,
                   indent + INCINDENT);
  STOPWATCH_PRINT (stopwatch_modified_quad_eval_schwartz_zippel_const,
                   STOPWATCH_MSEC, indent + INCINDENT);
  print_stopwatch_modified_quad_many_prove (indent + INCINDENT);
}


void
print_stopwatch_modified_quad_eval_verify (unsigned int indent)
{
  STOPWATCH_PRINT (stopwatch_modified_quad_eval_verify, STOPWATCH_MSEC, indent);
  STOPWATCH_PRINT (stopwatch_modified_quad_eval_schwartz_zippel_quad,
                   STOPWATCH_MSEC, indent + INCINDENT);
  STOPWATCH_PRINT (stopwatch_modified_quad_eval_schwartz_zippel_lin, STOPWATCH_MSEC,
                   indent + INCINDENT);
  STOPWATCH_PRINT (stopwatch_modified_quad_eval_schwartz_zippel_const,
                   STOPWATCH_MSEC, indent + INCINDENT);
  print_stopwatch_modified_quad_many_verify (indent + INCINDENT);
}