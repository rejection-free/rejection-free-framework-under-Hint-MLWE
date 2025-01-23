#include "lazer.h"
#include "stopwatch.h"
#include <mpfr.h>

static void
_schwartz_zippel_poly_2 (spolymat_t R2, spolyvec_t r1, poly_t r0,
                       uint8_t hash[32], spolymat_ptr R2i[],
                       spolyvec_ptr r1i[], poly_ptr r0i[], unsigned int N,
                       const modified_abdlop_params_t params)
{
  polyring_srcptr Rq = params->ring;
  int_srcptr q = polyring_get_mod (Rq);
  const unsigned int log2q = polyring_get_log2q (Rq);
  const unsigned int m1 = params->m1;
  const unsigned int l = params->l;
  spolyvec_t r1tmp, r1tmp2;
  spolymat_t R2tmp, R2tmp2;
  const unsigned int nelems = 2 * (m1 + l);
  poly_ptr mui;
  polyvec_t mu;
  unsigned int i;

  polyvec_alloc (mu, Rq, N);
  spolyvec_alloc (r1tmp, Rq, nelems, nelems);
  spolyvec_alloc (r1tmp2, Rq, nelems, nelems);
  spolymat_alloc (R2tmp, Rq, nelems, nelems,
                  (nelems * nelems - nelems) / 2 + nelems);
  spolymat_alloc (R2tmp2, Rq, nelems, nelems,
                  (nelems * nelems - nelems) / 2 + nelems);
  
  if (r0 != NULL)
    poly_set_zero (r0);

  polyvec_urandom (mu, q, log2q, hash, 0);

  _VEC_FOREACH_ELEM (mu, i)
  {
    mui = polyvec_get_elem (mu, i);

    if (R2i[i] != NULL)
      {
        spolymat_scale2 (R2tmp2, mui, R2i[i]);
        spolymat_add (R2tmp, R2, R2tmp2, 0);
        spolymat_set (R2, R2tmp);
      }
    if (r1i[i] != NULL)
      {
        spolyvec_scale2 (r1tmp2, mui, r1i[i]);
        spolyvec_add (r1tmp, r1, r1tmp2, 0);
        spolyvec_set (r1, r1tmp);
      }

    if (r0 == NULL)
      continue;

    if (!(r0i[i] == NULL))
      poly_addmul (r0, mui, r0i[i], 0);
  }
  spolymat_fromcrt (R2);
  spolyvec_fromcrt (r1);
  if (r0 != NULL)
    poly_fromcrt (r0);

  spolyvec_free (r1tmp);
  spolyvec_free (r1tmp2);
  spolymat_free (R2tmp);
  spolymat_free (R2tmp2);
  polyvec_free (mu);

  ASSERT_ERR (spolymat_is_upperdiag (R2));
}

/*
 * hash hash of tA1, tB.
 * t must be a subvector of (tB,tBext).
 */
void
modified_quad_many_prove (uint8_t hash[32], polyvec_t tB, poly_t c, polyvec_t z1,
                     polyvec_t z21, polyvec_t h, polyvec_t randencs1, polyvec_t m,
                     polyvec_t s2, polyvec_t tA2, polymat_t A1,
                     polymat_t A2prime, polymat_t Bprime, spolymat_ptr R2i[],
                     spolyvec_ptr r1i[], unsigned int N,
                     const uint8_t seed[32], const modified_abdlop_params_t params)
{
#if ASSERT == ASSERT_ENABLED
  unsigned int i;
#endif
  polyring_srcptr Rq = params->ring;
  const unsigned int m1 = params->m1;
  const unsigned int l = params->l;
  const unsigned int nelems = 2 * (m1 + l);
  spolymat_t R2;
  spolyvec_t r1;

  DEBUG_PRINTF (DEBUG_PRINT_FUNCTION_ENTRY, "%s", "modified_quad_many_prove");
  STOPWATCH_START (stopwatch_modified_quad_many_prove, "modified_quad_many_prove begin");

  ASSERT_ERR (N > 0); /* use quad if only one quad eq is needed. */
  ASSERT_ERR (params->lext == 1);
  ASSERT_ERR (poly_get_ring (c) == Rq);
  ASSERT_ERR (polyvec_get_ring (z1) == Rq);
  ASSERT_ERR (polyvec_get_nelems (z1) == 2*m1);
  ASSERT_ERR (polyvec_get_ring (z21) == Rq);
  ASSERT_ERR (polyvec_get_nelems (z21) == params->m2 - params->kmsis);
  ASSERT_ERR (polyvec_get_ring (h) == Rq);
  ASSERT_ERR (polyvec_get_nelems (h) == params->kmsis);
  ASSERT_ERR (polyvec_get_ring (tA2) == Rq);
  ASSERT_ERR (polyvec_get_nelems (tA2) == params->kmsis);
  ASSERT_ERR (polyvec_get_ring (tB) == Rq);
  ASSERT_ERR (polyvec_get_nelems (tB) == l + params->lext);
  ASSERT_ERR (polyvec_get_ring (randencs1) == Rq);
  ASSERT_ERR (polyvec_get_nelems (randencs1) == 2*m1);
  ASSERT_ERR (polyvec_get_ring (s2) == Rq);
  ASSERT_ERR (polyvec_get_nelems (s2) == params->m2);
  ASSERT_ERR (polymat_get_ring (A1) == Rq);
  ASSERT_ERR (polymat_get_nrows (A1) == params->kmsis);
  ASSERT_ERR (polymat_get_ncols (A1) == 2*m1);
  ASSERT_ERR (polymat_get_ring (A2prime) == Rq);
  ASSERT_ERR (polymat_get_nrows (A2prime) == params->kmsis);
  ASSERT_ERR (polymat_get_ncols (A2prime) == params->m2 - params->kmsis);
  ASSERT_ERR (polymat_get_ring (Bprime) == Rq);
  ASSERT_ERR (polymat_get_nrows (Bprime) == l + params->lext);
  ASSERT_ERR (polymat_get_ncols (Bprime) == params->m2 - params->kmsis);
#if ASSERT == ASSERT_ENABLED
  for (i = 0; i < N; i++)
    {
      ASSERT_ERR (spolymat_is_upperdiag (R2i[i]));
      ASSERT_ERR (R2i[i] == NULL
                  || spolymat_get_nrows (R2i[i]) == 2 * (m1 + l));
      ASSERT_ERR (R2i[i] == NULL
                  || spolymat_get_ncols (R2i[i]) == 2 * (m1 + l));
      ASSERT_ERR (r1i[i] == NULL || r1i[i]->nelems_max == 2 * (m1 + l));
    }
#endif

  spolyvec_alloc (r1, Rq, 2 * (m1 + l), 2 * (m1 + l));
  spolymat_alloc (R2, Rq, nelems, nelems,
                  (nelems * nelems - nelems) / 2 + nelems);

  _schwartz_zippel_poly_2 (R2, r1, NULL, hash, R2i, r1i, NULL, N, params);

  /* seed can be passed directly to quad sub-protocol since it is not used. */
  modified_quad_prove (hash, tB, c, z1, z21, h, randencs1, m, s2, tA2, A1, A2prime, Bprime,
                  R2, r1, seed, params);

  spolymat_free (R2);
  spolyvec_free (r1);

  STOPWATCH_STOP (stopwatch_modified_quad_many_prove);
  DEBUG_PRINTF (DEBUG_PRINT_FUNCTION_RETURN, "%s", "lnp_modified_many_prove end");
}

int
modified_quad_many_verify (uint8_t hash[32], poly_t c, polyvec_t z1, polyvec_t z21,
                      polyvec_t h, polyvec_t tA1, polyvec_t tB, polymat_t A1,
                      polymat_t A2prime, polymat_t Bprime, spolymat_ptr R2i[],
                      spolyvec_ptr r1i[], poly_ptr r0i[], unsigned int N,
                      const modified_abdlop_params_t params)
{
#if ASSERT == ASSERT_ENABLED
  unsigned int i;
#endif
  polyring_srcptr Rq = params->ring;
  const unsigned int m1 = params->m1;
  const unsigned int l = params->l;
  const unsigned int nelems = 2 * (m1 + l);
  spolymat_t R2;
  spolyvec_t r1;
  poly_t r0;
  int b;

  STOPWATCH_START (stopwatch_modified_quad_many_verify, "modified_quad_many_verify");

  ASSERT_ERR (params->lext == 1);
  ASSERT_ERR (poly_get_ring (c) == Rq);
  ASSERT_ERR (polyvec_get_ring (z1) == Rq);
  ASSERT_ERR (polyvec_get_nelems (z1) == 2*m1);
  ASSERT_ERR (polyvec_get_ring (z21) == Rq);
  ASSERT_ERR (polyvec_get_nelems (z21) == params->m2 - params->kmsis);
  ASSERT_ERR (polyvec_get_ring (h) == Rq);
  ASSERT_ERR (polyvec_get_nelems (h) == params->kmsis);
  ASSERT_ERR (polyvec_get_ring (tA1) == Rq);
  ASSERT_ERR (polyvec_get_nelems (tA1) == params->kmsis);
  ASSERT_ERR (polyvec_get_ring (tB) == Rq);
  ASSERT_ERR (polyvec_get_nelems (tB) == l + params->lext);
  ASSERT_ERR (polymat_get_ring (A1) == Rq);
  ASSERT_ERR (polymat_get_nrows (A1) == params->kmsis);
  ASSERT_ERR (polymat_get_ncols (A1) == 2*m1);
  ASSERT_ERR (polymat_get_ring (A2prime) == Rq);
  ASSERT_ERR (polymat_get_nrows (A2prime) == params->kmsis);
  ASSERT_ERR (polymat_get_ncols (A2prime) == params->m2 - params->kmsis);
  ASSERT_ERR (polymat_get_ring (Bprime) == Rq);
  ASSERT_ERR (polymat_get_nrows (Bprime) == l + params->lext);
  ASSERT_ERR (polymat_get_ncols (Bprime) == params->m2 - params->kmsis);
#if ASSERT == ASSERT_ENABLED
  for (i = 0; i < N; i++)
    {
      ASSERT_ERR (spolymat_is_upperdiag (R2i[i]));
      ASSERT_ERR (R2i[i] == NULL
                  || spolymat_get_nrows (R2i[i]) == 2 * (m1 + l));
      ASSERT_ERR (R2i[i] == NULL
                  || spolymat_get_ncols (R2i[i]) == 2 * (m1 + l));
      ASSERT_ERR (r1i[i] == NULL || r1i[i]->nelems_max == 2 * (m1 + l));
    }
#endif

  spolymat_alloc (R2, Rq, nelems, nelems,
                  (nelems * nelems - nelems) / 2 + nelems);
  spolyvec_alloc (r1, Rq, 2 * (m1 + l), 2 * (m1 + l));
  poly_alloc (r0, Rq);

  _schwartz_zippel_poly_2 (R2, r1, r0, hash, R2i, r1i, r0i, N, params);

  b = modified_quad_verify (hash, c, z1, z21, h, tA1, tB, A1, A2prime, Bprime, R2, r1, r0, params);

  spolymat_free (R2); 
  spolyvec_free (r1);
  poly_free (r0);

  STOPWATCH_STOP (stopwatch_modified_quad_many_verify);
  return b;
}
