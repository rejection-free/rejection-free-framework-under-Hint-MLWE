#include "lazer.h"
#include "stopwatch.h"

/*
 * hash hash of tA1, tB,
 * tB = (tB_,t, Bprime=(Bprime_,bextprime)
 */
void
rf_quad_prove (uint8_t hash[32], polyvec_t tB, poly_t c, polyvec_t z1,
                polyvec_t z21, polyvec_t h, polyvec_t randencs1, polyvec_t m,
                polyvec_t s2, polyvec_t tA2, polymat_t A1, polymat_t A2prime,
                polymat_t Bprime, spolymat_t R2, spolyvec_t r1,
                const uint8_t seed[32], const rf_abdlop_params_t params)
{
  polyring_srcptr Rq = params->ring;
  const unsigned int kmsis = params->kmsis;
  const unsigned int m1 = params->m1;
  const unsigned int m2 = params->m2;
  const unsigned int l = params->l;
  const unsigned int lext = params->lext;
  const unsigned int log2q = polyring_get_log2q (Rq);
  const unsigned int d = polyring_get_deg (Rq);
  int_srcptr q = polyring_get_mod (Rq);
  const unsigned int nlimbs = int_get_nlimbs (q);
  dcompress_params_srcptr dcomp_param = params->dcompress;
  int_srcptr gamma = dcompress_get_gamma (dcomp_param);
  int_srcptr m_ = dcompress_get_m (dcomp_param);
  const unsigned int log2m = dcompress_get_log2m (dcomp_param);
  polyvec_t tsub, asub, bsub, asub_auto, bsub_auto, y_rand, y21, y22, s21, s22, t, subv, R2y,
      tmp, tmp1, tmp2, y, s, w, w1, w0, y1, y2, cs1, cs2;
  polymat_t Bprime_, bextprime;

  spolymat_t R2prime, R2_ss_prime, R2_sm_prime, R2_ss, R2_sm, R2_mm, bR2_sm, nR2_ss, bR2_ss;
  spolyvec_t r1prime, r1_s, r1_m;
  
  

  INTVEC_T (z1coeffs, d * 2 * m1, nlimbs);
  INTVEC_T (z2coeffs, d * m2, nlimbs);
  INTVEC_T (cs1coeffs, d * 2 * m1, nlimbs);
  INTVEC_T (cs2coeffs, d * m2, nlimbs);
  INT_T (norm, nlimbs * 2);
  shake128_state_t hstate;
  coder_state_t cstate;
  rng_state_t rngstate;
  uint32_t dom;
  poly_t g1, g0;
  /* buff for encoding of t,v,w1 */
  uint8_t out[CEIL (2 * log2q * d + log2m * d * kmsis, 8) + 1];
  unsigned int outlen;
  uint8_t cseed[32]; /* challenge seed */
  uint8_t yseed[32]; /* mask seed */
  int rej = 1;

  DEBUG_PRINTF (DEBUG_PRINT_FUNCTION_ENTRY, "%s", "rf_quad_prove begin");
  STOPWATCH_START (stopwatch_rf_quad_prove, "rf_quad_prove");

  ASSERT_ERR (lext == 1);
  ASSERT_ERR (spolymat_is_upperdiag (R2));
  ASSERT_ERR (spolymat_get_ring (R2) == Rq);
  ASSERT_ERR (spolymat_get_nrows (R2) == 2 * (m1 + l));
  ASSERT_ERR (spolymat_get_ncols (R2) == 2 * (m1 + l));
  ASSERT_ERR (spolyvec_get_ring (r1) == Rq);
  ASSERT_ERR (r1->nelems_max == 2 * (m1 + l));
  ASSERT_ERR (poly_get_ring (c) == Rq);
  ASSERT_ERR (polyvec_get_ring (z1) == Rq);
  ASSERT_ERR (polyvec_get_nelems (z1) == m1);
  ASSERT_ERR (polyvec_get_ring (z21) == Rq);
  ASSERT_ERR (polyvec_get_nelems (z21) == m2 - kmsis);
  ASSERT_ERR (polyvec_get_ring (h) == Rq);
  ASSERT_ERR (polyvec_get_nelems (h) == kmsis);
  ASSERT_ERR (polyvec_get_ring (tA2) == Rq);
  ASSERT_ERR (polyvec_get_nelems (tA2) == kmsis);
  ASSERT_ERR (polyvec_get_ring (tB) == Rq);
  ASSERT_ERR (polyvec_get_nelems (tB) == l + lext);
  ASSERT_ERR (polyvec_get_ring (randencs1) == Rq);
  ASSERT_ERR (polyvec_get_nelems (randencs1) == 2*m1);
  ASSERT_ERR (polyvec_get_ring (m) == Rq);
  ASSERT_ERR (polyvec_get_nelems (m) == l + lext);
  ASSERT_ERR (polyvec_get_ring (s2) == Rq);
  ASSERT_ERR (polyvec_get_nelems (s2) == m2);
  ASSERT_ERR (polymat_get_ring (A1) == Rq);
  ASSERT_ERR (polymat_get_nrows (A1) == kmsis);
  ASSERT_ERR (polymat_get_ncols (A1) == 2 * m1);
  ASSERT_ERR (polymat_get_ring (A2prime) == Rq);
  ASSERT_ERR (polymat_get_nrows (A2prime) == kmsis);
  ASSERT_ERR (polymat_get_ncols (A2prime) == m2 - kmsis);
  ASSERT_ERR (polymat_get_ring (Bprime) == Rq);
  ASSERT_ERR (polymat_get_nrows (Bprime) == l + lext);
  ASSERT_ERR (polymat_get_ncols (Bprime) == m2 - kmsis);

  polyvec_alloc (y_rand, Rq, m1);
  polyvec_alloc (y1, Rq, 2 * m1);
  polyvec_alloc (y2, Rq, m2);
  polyvec_alloc (cs1, Rq, 2 * m1);
  polyvec_alloc (cs2, Rq, m2);

  polyvec_alloc (R2y, Rq, 2 * (2 * m1 + l));
  polyvec_alloc (tmp, Rq, 2 * (2 * m1 + l));
  if (l > 0)
    polyvec_alloc (tmp2, Rq, l);
    
  polyvec_alloc (y, Rq, 2 * (2 * m1 + l));
  polyvec_alloc (s, Rq, 2 * (2 * m1 + l));
  polyvec_alloc (w, Rq, kmsis);
  polyvec_alloc (w1, Rq, kmsis);
  polyvec_alloc (w0, Rq, kmsis);
  poly_alloc (g1, Rq);
  poly_alloc (g0, Rq);

  spolymat_alloc (R2prime, Rq, 2 * (2*m1+l), 2 * (2*m1+l),  2 * (2 * m1 * 2 * m1 + 2 * m1) + 8 *  m1 * l + (2*l * 2*l + 2*l)/2  );
  spolymat_alloc (R2_ss, Rq, 2 * m1, 2 * m1, (2 * m1 * 2 * m1 + 2 * m1)/2);
  spolymat_alloc (bR2_ss, Rq, 2 * m1, 2 * m1, (2 * m1 * 2 * m1 + 2 * m1)/2);
  spolymat_alloc (nR2_ss, Rq, 2 * m1, 2 * m1, (2 * m1 * 2 * m1 + 2 * m1)/2);
  spolymat_alloc (R2_ss_prime, Rq, 4 * m1, 4 * m1, 2*(2 * m1 * 2 * m1 + 2 * m1));
  spolymat_alloc (R2_sm, Rq, 2 * m1, 2 * l, 4 *  m1 * l);
  spolymat_alloc (bR2_sm, Rq, 2 * m1, 2 * l, 4 *  m1 * l);
  spolymat_alloc (R2_sm_prime, Rq, 4 * m1, 2 * l, 8 *  m1 * l);
  spolymat_alloc (R2_mm, Rq, 2 * l, 2 * l, (2*l * 2*l + 2*l)/2 );
  spolyvec_alloc (r1prime, Rq, 2*m1, 2*m1);
  spolyvec_alloc (r1_s, Rq, 2*m1, 2*m1);
  spolyvec_alloc (r1_m, Rq, 2*l, 2*l);

  polyvec_get_subvec (t, tB, l, lext, 1);
  polymat_get_submat (bextprime, Bprime, l, 0, lext, m2 - kmsis, 1, 1);
  if (l > 0)
    polymat_get_submat (Bprime_, Bprime, 0, 0, l, m2 - kmsis, 1, 1);

  /* reuse y1, y2=(y21,y22) as z1, z2=(z21,z22) later */

  polyvec_set_coeffvec2 (y1, z1coeffs);
  polyvec_set_coeffvec2 (y2, z2coeffs);
  polyvec_set_coeffvec2 (cs1, cs1coeffs);
  polyvec_set_coeffvec2 (cs2, cs2coeffs);

  polyvec_get_subvec (s21, s2, 0, m2 - kmsis, 1);
  polyvec_get_subvec (s22, s2, m2 - kmsis, kmsis, 1);

  polyvec_get_subvec (y21, y2, 0, m2 - kmsis, 1);
  polyvec_get_subvec (y22, y2, m2 - kmsis, kmsis, 1);

  /* s = (<s1>,<m>) */

  polyvec_get_subvec (asub, s, 0, m1, 2);
  polyvec_get_subvec (asub_auto, s, 1, m1, 2);
  polyvec_get_subvec (tsub, randencs1, 0, m1, 2);
  polyvec_set (asub, tsub);
  polyvec_auto (asub_auto, tsub);

  polyvec_get_subvec (asub, s, 2*m1, m1, 2);
  polyvec_get_subvec (asub_auto, s, 2*m1+1, m1, 2);
  polyvec_get_subvec (tsub, randencs1, 1, m1, 2);
  polyvec_set (asub, tsub);
  polyvec_auto (asub_auto, tsub);

  if (l > 0)
    {
      polyvec_get_subvec (bsub, s, m1 * 4, l, 2);
      polyvec_get_subvec (bsub_auto, s, m1 * 4 + 1, l, 2);
      polyvec_get_subvec (subv, m, 0, l, 1);
      polyvec_set (bsub, subv);
      polyvec_auto (bsub_auto, subv);
    }

  if (l > 0)
    {
      polyvec_get_subvec (bsub, y, 2*y1->nelems, l, 2);
      polyvec_get_subvec (bsub_auto, y, 2*y1->nelems + 1, l, 2);
    }

  rng_init (rngstate, seed, 0);
  rng_urandom (rngstate, yseed, 32);

  /* Computation of the new r1prime and R2 prime */
  spolyvec_get_subvec(r1_s,r1,0,2*m1);
  spolyvec_get_subvec(r1_m,r1,2*m1,2*l);
  spolyvec_scale(r1prime,params->base,r1_s);
     
  spolymat_get_submat_upperdiag(R2_ss, R2, 0, 0);
  spolymat_get_submat_upperdiag(R2_mm, R2, 2*m1, 2*m1);
  spolymat_get_submat(R2_sm, R2, 0, 2*m1);

  spolymat_scale(bR2_sm,params->base,R2_sm);
  spolymat_scale(bR2_ss,params->base,R2_ss);
  spolymat_neg(nR2_ss,R2_ss);

  spolymat_set_squared_block_matrix (R2_ss_prime, R2_ss, bR2_ss, bR2_ss, nR2_ss);
  spolymat_set_block_matrix (R2_sm_prime, R2_sm, bR2_sm);
  spolymat_set_upper_block_matrix (R2prime, R2_ss_prime, R2_sm_prime, R2_mm);

  dom = 0;
  while (1)
    {
      /* y1, y2 */

      polyvec_grandom (y_rand, params->log2sigma1, yseed, dom++);
      randomized_encoding(y1, y_rand, seed, params->log2sigma1, params->base);

      polyvec_grandom (y2, params->log2sigma2, yseed, dom++);

      /* w */

      polyvec_set (w, y22);
      polyvec_addmul (w, A1, y1, 1);
      polyvec_addmul (w, A2prime, y21, 1); // XXX w correct, 7 rejs

      polyvec_dcompress_decompose (
          w1, w0, w, dcomp_param); // w1*gamma+w0 == w XXX correct

      /* y */
      /* y = (<y1>,-<By2>) */

      polyvec_get_subvec (asub, y, 0, m1, 2);
      polyvec_get_subvec (asub_auto, y, 1, m1, 2);
      polyvec_get_subvec (tsub, y1, 0, m1, 2);
      polyvec_set (asub, tsub);
      polyvec_auto (asub_auto, tsub);
      polyvec_get_subvec (asub, y, 2*m1, m1, 2);
      polyvec_get_subvec (asub_auto, y, 2*m1+1, m1, 2);
      polyvec_get_subvec (tsub, y1, 1, m1, 2);
      polyvec_set (asub, tsub);
      polyvec_auto (asub_auto, tsub);

      if (l > 0)
        {
          polyvec_mul (tmp2, Bprime_, y21);
          polyvec_set (bsub, tmp2);
          polyvec_auto (bsub_auto, tmp2);
          polyvec_neg_self (bsub);
          polyvec_neg_self (bsub_auto);
        }

      polyvec_fromcrt (y);

      /* g_1 */
      polyvec_get_subvec (tmp1, y, 2*m1, 2*m1, 1);
      polyvec_dot2 (g1, r1prime, tmp1);
      polyvec_get_subvec (tmp1, y, 0, 2*m1, 1);
      poly_adddot2 (g1, r1_s, tmp1, 0);
      polyvec_get_subvec (tmp1, y, 4*m1, 2*l, 1);
      poly_adddot2 (g1, r1_m, tmp1, 0);
      poly_fromcrt (g1);
      
      polyvec_mulsparse (R2y, R2prime, s);
      polyvec_fromcrt (R2y); // reduce XXX
      poly_adddot (g1, y, R2y, 0);
      
      polyvec_mulsparse (R2y, R2prime, y);
      polyvec_fromcrt (R2y);
      poly_adddot (g1, s, R2y, 0);
      poly_fromcrt (g1);
       /* t */

      poly_set (polyvec_get_elem (t, 0), g1);
      polyvec_addmul (t, bextprime, s21, 0);
      polyvec_fromcrt (t);

      /* g0 */
      
      polyvec_dot (g0, y, R2y);
      poly_addmul2 (g0, bextprime, y21, 0);
      poly_fromcrt (g0);

      /* encode */

      polyvec_mod (t, t);
      polyvec_redp (t, t);

      poly_mod (g0, g0);
      poly_redp (g0, g0);

      coder_enc_begin (cstate, out);
      coder_enc_urandom3 (cstate, t, q, log2q);
      coder_enc_urandom2 (cstate, g0, q, log2q);
      coder_enc_urandom3 (cstate, w1, m_, log2m);
      coder_enc_end (cstate);

      outlen = coder_get_offset (cstate);
      ASSERT_ERR (outlen % 8 == 0);
      ASSERT_ERR (outlen / 8
                  <= CEIL (2 * log2q * d + log2m * d * kmsis, 8) + 1);
      outlen >>= 3; /* nbits to nbytes */

      shake128_init (hstate);
      shake128_absorb (hstate, hash, 32);
      shake128_absorb (hstate, out, outlen);
      shake128_squeeze (hstate, cseed, 32);

      poly_urandom_autostable (c, params->omega, params->log2omega, cseed, 0);

      polyvec_scale2 (cs1, c, randencs1);
      polyvec_scale2 (cs2, c, s2); // XXX correct
      polyvec_add (y1, y1, cs1, 0);
      polyvec_add (y2, y2, cs2, 0); // XXX correct, XXX sample sign with bimodal ?

      polyvec_fromcrt (y1);
      polyvec_fromcrt (cs1);
      
      polyvec_fromcrt (y2);
      polyvec_fromcrt (cs2);

      polyvec_subscale2 (y22, c, tA2, 0);
      polyvec_sub (y22, y22, w0, 0); // XXX correct

      polyvec_l2sqr (norm, y2);

      rej = int_gt (norm, params->Bsqr);
      if (rej)
        {
          DEBUG_PRINTF (DEBUG_PRINT_REJ, "%s", "reject on s1");
          continue;
        }

      break;
    }
  /* update fiat-shamir hash */
  memcpy (hash, cseed, 32);

  polyvec_scale (w1, gamma, w1);
  polyvec_sub (w1, w1, y22, 0); // correct XXX

  /* output proof (z1,z21,h,c) */
  polyvec_set (z1, y1);
  polyvec_set (z21, y21);
  
  polyvec_mod (w1, w1);
  polyvec_redc (w1, w1); // XXX
  polyvec_dcompress_make_ghint (h, y22, w1, dcomp_param);

  /* cleanup */
  shake128_clear (hstate);
  rng_clear (rngstate);

  polyvec_free (y_rand);
  polyvec_free (y1);
  polyvec_free (y2);
  polyvec_free (cs1);
  polyvec_free (cs2);

  polyvec_free (R2y);
  polyvec_free (tmp);
  if (l > 0)
    polyvec_free (tmp2);
  polyvec_free (y);
  polyvec_free (s);
  polyvec_free (w);
  polyvec_free (w1);
  polyvec_free (w0);
  poly_free (g1);
  poly_free (g0);
  
  spolymat_free (R2_ss);
  spolymat_free (bR2_ss);
  spolymat_free (nR2_ss);
  spolymat_free (R2_sm);
  spolymat_free (bR2_sm);
  spolymat_free (R2_mm);
  spolymat_free (R2prime);
  spolymat_free(R2_sm_prime);
  spolymat_free(R2_ss_prime);

  spolyvec_free (r1prime);
  spolyvec_free (r1_s);
  spolyvec_free (r1_m);

  STOPWATCH_STOP (stopwatch_rf_quad_prove);
  DEBUG_PRINTF (DEBUG_PRINT_FUNCTION_RETURN, "%s", "rf_quad_prove end");
}

int
rf_quad_verify (uint8_t hash[32], poly_t c, polyvec_t z1, polyvec_t z21,
                 polyvec_t h, polyvec_t tA1, polyvec_t tB, polymat_t A1,
                 polymat_t A2prime, polymat_t Bprime, spolymat_t R2,
                 spolyvec_t r1, poly_t r0, const rf_abdlop_params_t params)
{
  polyring_srcptr Rq = params->ring;
  const unsigned int kmsis = params->kmsis;
  const unsigned int m1 = params->m1;
  const unsigned int m2 = params->m2;
  const unsigned int l = params->l;
  const unsigned int lext = params->lext;
  const unsigned int log2q = polyring_get_log2q (Rq);
  const unsigned int d = polyring_get_deg (Rq);
  int_srcptr q = polyring_get_mod (Rq);
  const unsigned int nlimbs = int_get_nlimbs (q);
  dcompress_params_srcptr dcomp_param = params->dcompress;
  
  int_srcptr m = dcompress_get_m (dcomp_param);
  const unsigned int D = dcompress_get_d (dcomp_param);
  const unsigned int log2m = dcompress_get_log2m (dcomp_param);
  int_srcptr m_ = dcompress_get_m (dcomp_param);
  INT_T (l2sqr, nlimbs * 2);
  INT_T (bnd, nlimbs * 2);
  INT_T (linf, nlimbs);
  INT_T (tmp, 1);
  unsigned int outlen;
  shake128_state_t hstate;
  polyvec_t tsub, suba, suba_auto, subb, subb_auto, tB_, t, z, tmp4, tmp3, w1, tmp1, f;
  polymat_t Bprime_, bextprime;
  spolymat_t R2prime, R2_ss_prime, R2_sm_prime, R2_ss, R2_sm, R2_mm, bR2_sm, nR2_ss, bR2_ss;
  spolyvec_t r1prime, r1_s, r1_m;
  poly_t tmp2, c2, v;
  coder_state_t cstate;
  /* buff for encoding of t,v,w1 */
  uint8_t out[CEIL (2 * log2q * d + log2m * d * kmsis, 8) + 1];
  uint8_t cseed[32];
  int skip, accept = 0;

  STOPWATCH_START (stopwatch_rf_quad_verify, "rf_quad_verify");

  ASSERT_ERR (lext == 1);
  ASSERT_ERR (poly_get_ring (c) == Rq);
  ASSERT_ERR (polyvec_get_ring (z1) == Rq);
  ASSERT_ERR (polyvec_get_nelems (z1) == m1);
  ASSERT_ERR (polyvec_get_ring (z21) == Rq);
  ASSERT_ERR (polyvec_get_nelems (z21) == m2 - kmsis);
  ASSERT_ERR (polyvec_get_ring (h) == Rq);
  ASSERT_ERR (polyvec_get_nelems (h) == kmsis);
  ASSERT_ERR (polyvec_get_ring (tA1) == Rq);
  ASSERT_ERR (polyvec_get_nelems (tA1) == kmsis);
  ASSERT_ERR (polyvec_get_ring (tB) == Rq);
  ASSERT_ERR (polyvec_get_nelems (tB) == l + lext);
  ASSERT_ERR (polymat_get_ring (A1) == Rq);
  ASSERT_ERR (polymat_get_nrows (A1) == kmsis);
  ASSERT_ERR (polymat_get_ncols (A1) == m1);
  ASSERT_ERR (polymat_get_ring (A2prime) == Rq);
  ASSERT_ERR (polymat_get_nrows (A2prime) == kmsis);
  ASSERT_ERR (polymat_get_ncols (A2prime) == m2 - kmsis);
  ASSERT_ERR (polymat_get_ring (Bprime) == Rq);
  ASSERT_ERR (polymat_get_nrows (Bprime) == l + lext);
  ASSERT_ERR (polymat_get_ncols (Bprime) == m2 - kmsis);
  ASSERT_ERR (spolymat_is_upperdiag (R2));
  ASSERT_ERR (spolymat_get_ring (R2) == Rq);
  ASSERT_ERR (spolymat_get_nrows (R2) == 2 * (m1 + l));
  ASSERT_ERR (spolymat_get_ncols (R2) == 2 * (m1 + l));
  ASSERT_ERR (spolyvec_get_ring (r1) == Rq);
  ASSERT_ERR (r1->nelems_max == 2 * (m1 + l));

  polyvec_alloc (z, Rq, 2 * (2 * m1 + l));
  polyvec_alloc (tmp3, Rq, 2 * (2 * m1 + l));
  polyvec_alloc (w1, Rq, kmsis);
  polyvec_alloc (tmp1, Rq, kmsis);
  polyvec_alloc (f, Rq, 1);
  poly_alloc (tmp2, Rq);
  poly_alloc (c2, Rq);
  poly_alloc (v, Rq);


  spolymat_alloc (R2prime, Rq, 2 * (2*m1+l), 2 * (2*m1+l),  2 * (2 * m1 * 2 * m1 + 2 * m1) + 8 *  m1 * l + (2*l * 2*l + 2*l)/2  );
  spolymat_alloc (R2_ss, Rq, 2 * m1, 2 * m1, (2 * m1 * 2 * m1 + 2 * m1)/2);
  spolymat_alloc (bR2_ss, Rq, 2 * m1, 2 * m1, (2 * m1 * 2 * m1 + 2 * m1)/2);
  spolymat_alloc (nR2_ss, Rq, 2 * m1, 2 * m1, (2 * m1 * 2 * m1 + 2 * m1)/2);
  spolymat_alloc (R2_ss_prime, Rq, 4 * m1, 4 * m1, 2*(2 * m1 * 2 * m1 + 2 * m1));
  spolymat_alloc (R2_sm, Rq, 2 * m1, 2 * l, 4 *  m1 * l);
  spolymat_alloc (bR2_sm, Rq, 2 * m1, 2 * l, 4 *  m1 * l);
  spolymat_alloc (R2_sm_prime, Rq, 4 * m1, 2 * l, 8 *  m1 * l);
  spolymat_alloc (R2_mm, Rq, 2 * l, 2 * l, (2*l * 2*l + 2*l)/2 );
  spolyvec_alloc (r1prime, Rq, 2*m1, 2*m1);
  spolyvec_alloc (r1_s, Rq, 2*m1, 2*m1);
  spolyvec_alloc (r1_m, Rq, 2*l, 2*l);

  polyvec_get_subvec (t, tB, l, lext, 1);
  polymat_get_submat (bextprime, Bprime, l, 0, lext, m2 - kmsis, 1, 1);
  if (l > 0)
    {
      polyvec_get_subvec (tB_, tB, 0, l, 1);
      polymat_get_submat (Bprime_, Bprime, 0, 0, l, m2 - kmsis, 1, 1);
    }

  /* Computation of the new r1prime and R2 prime */
  spolyvec_get_subvec(r1_s,r1,0,2*m1);
  spolyvec_get_subvec(r1_m,r1,2*m1,2*l);
  spolyvec_scale(r1prime,params->base,r1_s);
      
  spolymat_get_submat_upperdiag(R2_ss, R2, 0, 0);
  spolymat_get_submat_upperdiag(R2_mm, R2, 2*m1, 2*m1);
  spolymat_get_submat(R2_sm, R2, 0, 2*m1);

  spolymat_scale(bR2_sm,params->base,R2_sm);
  spolymat_scale(bR2_ss,params->base,R2_ss);
  spolymat_neg(nR2_ss,R2_ss);

  spolymat_set_squared_block_matrix (R2_ss_prime, R2_ss, bR2_ss, bR2_ss, nR2_ss);
  spolymat_set_block_matrix (R2_sm_prime, R2_sm, bR2_sm);
  spolymat_set_upper_block_matrix (R2prime, R2_ss_prime, R2_sm_prime, R2_mm);

  /* recover w1 */

  polyvec_mul (tmp1, A1, z1);
  polyvec_addmul (tmp1, A2prime, z21, 0);
  poly_lshift (tmp2, c, D);
  polyvec_subscale2 (tmp1, tmp2, tA1, 0);

  polyvec_mod (tmp1, tmp1); // XXX
  polyvec_dcompress_use_ghint (w1, h, tmp1, dcomp_param);

  /* recover v */

  polyvec_scale2 (f, c, t);
  polyvec_submul (f, bextprime, z21, 0);

  polyvec_get_subvec (suba, z, 0, m1, 2);
  polyvec_get_subvec (suba_auto, z, 1, m1, 2);
  polyvec_get_subvec (tsub, z1, 0, m1, 2);
  polyvec_set (suba, tsub);
  polyvec_auto (suba_auto, tsub);

  polyvec_get_subvec (suba, z, 2*m1, m1, 2);
  polyvec_get_subvec (suba_auto, z, 2*m1+1, m1, 2);
  polyvec_get_subvec (tsub, z1, 1, m1, 2);
  polyvec_set (suba, tsub);
  polyvec_auto (suba_auto, tsub);

  if (l > 0)
    {
      polyvec_get_subvec (subb, z, 2 * 2 * m1, l, 2);
      polyvec_get_subvec (subb_auto, z, 2 * 2 * m1 + 1, l, 2);
      polyvec_scale2 (subb, c, tB_);
      polyvec_submul (subb, Bprime_, z21, 0);
      polyvec_auto (subb_auto, subb);
    }
  polyvec_fromcrt (z);
  
  poly_set (v, r0);                            /* r0 */
  poly_mul (v, c, v);                          /* c * r0 */
  
  polyvec_get_subvec (tmp4, z, 0, 2*m1, 1);
  poly_adddot2 (v, r1_s, tmp4, 0);
  
  polyvec_get_subvec (tmp4, z, 2*m1, 2*m1, 1);
  poly_adddot2 (v, r1prime, tmp4, 0);
  polyvec_get_subvec (tmp4, z, 4*m1, 2*l, 1);
  poly_adddot2 (v, r1_m, tmp4, 0);

  poly_fromcrt (v);                            // XXX reduce
  poly_mul (v, c, v);                          /* c*r1*z + c^2*r0 */
  poly_sub (v, v, polyvec_get_elem (f, 0), 0); /* c*r1*z + c^2*r0 - f */
  polyvec_mulsparse (tmp3, R2prime, z);             /* R2*z */
  polyvec_fromcrt (tmp3);                      // XXX reduce
  poly_adddot (v, z, tmp3, 0); /* z*R2*z + c*r1*z + c^2*r0 - f*/
  poly_fromcrt (v);

  /* recover challenge from t, w1, v */
  polyvec_mod (t, t);
  polyvec_redp (t, t);
  poly_mod (v, v);
  poly_redp (v, v);
  
  coder_enc_begin (cstate, out);
  coder_enc_urandom3 (cstate, t, q, log2q);
  coder_enc_urandom2 (cstate, v, q, log2q);
  coder_enc_urandom3 (cstate, w1, m_, log2m);
  coder_enc_end (cstate);

  outlen = coder_get_offset (cstate);
  ASSERT_ERR (outlen % 8 == 0);
  ASSERT_ERR (outlen / 8 <= CEIL (2 * log2q * d + log2m * d * kmsis, 8) + 1);
  outlen >>= 3; /* nbits to nbytes */

  shake128_init (hstate);
  shake128_absorb (hstate, hash, 32);
  shake128_absorb (hstate, out, outlen);
  shake128_squeeze (hstate, cseed, 32);

  poly_urandom_autostable (c2, params->omega, params->log2omega, cseed, 0);
  skip = poly_eq (c, c2);
  if (!skip){
    goto ret;
  }
  
  /* check bounds */

  /* 2*linf(h) <= m */
  polyvec_linf (linf, h);
  skip = int_le (linf, m);
  if (!skip)
    goto ret;

  /* update fiat-shamir hash */
  memcpy (hash, cseed, 32);
  accept = 1;
ret:
  /* cleanup */
  shake128_clear (hstate);
  polyvec_free (z);
  polyvec_free (tmp3);
  polyvec_free (w1);
  polyvec_free (tmp1);
  polyvec_free (f);
  poly_free (tmp2);
  poly_free (c2);
  poly_free (v);
  spolymat_free (R2_ss);
  spolymat_free (bR2_ss);
  spolymat_free (nR2_ss);
  spolymat_free (R2_sm);
  spolymat_free (bR2_sm);
  spolymat_free (R2_mm);
  spolymat_free (R2prime);
  spolymat_free(R2_sm_prime);
  spolymat_free(R2_ss_prime);

  spolyvec_free (r1prime);
  spolyvec_free (r1_s);
  spolyvec_free (r1_m);

  STOPWATCH_STOP (stopwatch_rf_quad_verify);
  return accept;
}
