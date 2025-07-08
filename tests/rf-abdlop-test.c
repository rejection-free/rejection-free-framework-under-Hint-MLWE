#include "lazer.h"
#include "rf-abdlop-params1.h"
#include "rf-abdlop-params2.h"
#include "rf-abdlop-params3.h"
#include "rf-abdlop-params4.h"
#include "rf-abdlop-params5.h"
#include "test.h"

static void test_mabdlop (uint8_t seed[32], const rf_abdlop_params_t params);

int
main (void)
{
  unsigned int i;
  uint8_t seed[32];

  lazer_init();

  for (i = 0; i < 1; i++)
    {
      bytes_urandom (seed, sizeof (seed));
      test_mabdlop (seed, rf_params1);
    }
    
  for (i = 0; i < 1; i++)
    {
      bytes_urandom (seed, sizeof (seed));
      test_mabdlop (seed, rf_params2);
    }
    
  for (i = 0; i < 1; i++)
    {
      bytes_urandom (seed, sizeof (seed));
      test_mabdlop (seed, rf_params3);
    }
    
  for (i = 0; i < 1; i++)
    {
      bytes_urandom (seed, sizeof (seed));
      test_mabdlop (seed, rf_params4);
    }
  for (i = 0; i < 1; i++)
    {
      bytes_urandom (seed, sizeof (seed));
      test_mabdlop (seed, rf_params5);
    }  
  TEST_PASS ();
}

void
test_mabdlop(uint8_t seed[32], const rf_abdlop_params_t params){
  uint8_t hashp[32] = { 0 };
  uint8_t hashv[32] = { 0 };
  uint8_t hashcomm[32] = { 0 };
  polyring_srcptr Rq = params->ring;
  INT_T (lo, Rq->q->nlimbs);
  INT_T (hi, Rq->q->nlimbs);
  uint8_t buf[2];
  uint32_t dom;
  unsigned int i;
  int validity;
  polymat_t A1, A2prime, Bprime, A1err, A2primeerr;
  polyvec_t z1err, z21err, herr, tA1err, s1, randencs1, s2, m, tA1, tA2, tB, z1, z21, h;
  poly_t c;
  const unsigned int l = params->l + params->lext;
  
  poly_alloc (c, Rq);
  polyvec_alloc (z1err, Rq, 2*params->m1);
  polyvec_alloc (z21err, Rq, params->m2 - params->kmsis);
  polyvec_alloc (herr, Rq, params->kmsis);
  polyvec_alloc (tA1err, Rq, params->kmsis);
  polyvec_alloc (s1, Rq, params->m1);
  polyvec_alloc (randencs1, Rq, 2*params->m1);
  polyvec_alloc (s2, Rq, params->m2);
  if (l > 0)
    polyvec_alloc (m, Rq, params->l + params->lext);
  polyvec_alloc (tA1, Rq, params->kmsis);
  polyvec_alloc (tA2, Rq, params->kmsis);
  if (l > 0)
    polyvec_alloc (tB, Rq, params->l + params->lext);
  polyvec_alloc (z1, Rq, 2*params->m1);
  polyvec_alloc (z21, Rq, params->m2 - params->kmsis);
  polyvec_alloc (h, Rq, params->kmsis);
  polymat_alloc (A1, Rq, params->kmsis, 2*params->m1);
  polymat_alloc (A2prime, Rq, params->kmsis, params->m2 - params->kmsis);
  if (l > 0)
    polymat_alloc (Bprime, Rq, params->l + params->lext,
                   params->m2 - params->kmsis);
  polymat_alloc (A1err, Rq, params->kmsis, 2*params->m1);
  polymat_alloc (A2primeerr, Rq, params->kmsis, params->m2 - params->kmsis);

  dom = 0;
  int_set_i64 (lo, -1);
  int_set_i64 (hi, 1);
  polyvec_urandom_bnd (s1, lo, hi, seed, dom++);
  
  polyvec_grandom(s2, params->log2sigma2, seed, dom++);
  if (l > 0)
    polyvec_urandom (m, Rq->q, Rq->log2q, seed, dom++);

  /* generate public parameters */
  rf_abdlop_keygen (A1, A2prime, Bprime, seed, params);
  memcpy (hashp, seed, 32);
  
  rf_abdlop_commit (tA1, tA2, tB, s1, randencs1, m, s2, A1, A2prime, Bprime, seed, params);
  rf_abdlop_hashcomm (hashp, tA1, tB, params);
  rf_abdlop_prove (hashp, c, z1, z21, h, tA2, randencs1, s2, A1, A2prime, seed, params);
  /* expect successful verification */
  
  memcpy (hashv, seed, 32);
  rf_abdlop_hashcomm (hashv, tA1, tB, params);
  validity = rf_abdlop_verify (hashv, c, z1, z21, h, tA1, A1, A2prime, params);
  TEST_EXPECT (validity == 1);
  TEST_EXPECT (memcmp (hashp, hashv, 32) == 0);

  memcpy (hashcomm, seed, 32);
  rf_abdlop_hashcomm (hashcomm, tA1, tB, params);
  for (i = 0; i < 1; i++)
    {
      /* expect verification failures */
      bytes_urandom (buf, sizeof (buf));
      memcpy (hashv, hashcomm, 32);
      hashv[buf[0] % 32] ^= (1 << (buf[1] % 8));
      validity = rf_abdlop_verify (hashv, c, z1, z21, h, tA1, A1, A2prime, params);
      TEST_EXPECT (validity == 0);
     

      polyvec_brandom (z1err, 1, seed, dom++);
      polyvec_add (z1err, z1err, z1, 0);
      memcpy (hashv, hashcomm, 32);
      validity = rf_abdlop_verify (hashv, c, z1err, z21, h, tA1, A1, A2prime, params);
      TEST_EXPECT (validity == 0);
      
      polyvec_brandom (z21err, 1, seed, dom++);
      polyvec_add (z21err, z21err, z21, 0);
      memcpy (hashv, hashcomm, 32);
      validity = rf_abdlop_verify (hashv, c, z1, z21err, h, tA1, A1, A2prime, params);
      TEST_EXPECT (validity == 0);
      
      polyvec_brandom (herr, 1, seed, dom++);
      polyvec_add (herr, herr, h, 0);
      memcpy (hashv, hashcomm, 32);
      validity = rf_abdlop_verify (hashv, c, z1, z21, herr, tA1, A1, A2prime, params);
      TEST_EXPECT (validity == 0);
      
      polyvec_brandom (tA1err, 1, seed, dom++); /* sometimes fails XXX */
      polyvec_add (tA1err, tA1err, tA1, 0);
      memcpy (hashv, hashcomm, 32);
      validity = rf_abdlop_verify (hashv, c, z1, z21, h, tA1err, A1, A2prime, params);
      TEST_EXPECT (validity == 0);
      
      polymat_urandom (A1err, Rq->q, Rq->log2q, seed, dom++);
      polymat_add (A1err, A1err, A1, 0);
      memcpy (hashv, hashcomm, 32);
      validity = rf_abdlop_verify (hashv, c, z1, z21, h, tA1, A1err, A2prime, params);
      TEST_EXPECT (validity == 0);
      
      polymat_urandom (A2primeerr, Rq->q, Rq->log2q, seed, dom++);
      polymat_add (A2primeerr, A2primeerr, A2prime, 0);
      memcpy (hashv, hashcomm, 32);
      validity = rf_abdlop_verify (hashv, c, z1, z21, h, tA1, A1, A2primeerr, params);
      TEST_EXPECT (validity == 0);
    }
  poly_free (c);
  polyvec_free (z1err);
  polyvec_free (z21err);
  polyvec_free (herr);
  polyvec_free (tA1err);
  polyvec_free (s1);
  polyvec_free (randencs1);
  polyvec_free (s2);
  if (l > 0)
    polyvec_free (m);
  polyvec_free (tA1);
  polyvec_free (tA2);
  if (l > 0)
    polyvec_free (tB);
  polyvec_free (z1);
  polyvec_free (z21);
  polyvec_free (h);
  polymat_free (A1);
  polymat_free (A2prime);
  if (l > 0)
    polymat_free (Bprime);
  polymat_free (A1err);
  polymat_free (A2primeerr);
}