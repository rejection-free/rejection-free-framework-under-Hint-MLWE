typedef struct
{
  const polyring_srcptr ring;
  const dcompress_params_srcptr dcompress;
  const int_srcptr base;
  /* dimensions  */
  const unsigned int m1;   /* length of "short" message s1 */
  const unsigned int m2;   /* length of randomness s2 */
  const unsigned int l;    /* length of "large" message m */
  const unsigned int lext; /* length of extension of m */
  const unsigned int kmsis;
  /* norms */
  const int_srcptr Bsqr; /* floor (B^2) */
  const int64_t omega;   /* challenges uniform in [-omega,omega], o(c)=c */
  const unsigned int log2omega;
  const uint64_t eta; /* sqrt(l1(o(c)*c)) <= eta XXX sqrt? */
  /* standard deviations  */
  const int_srcptr sigma1; /* standard deviation of s1*/
  const unsigned int log2sigma1; /* sigma1 = 2^log2sigma1 */
  const int_srcptr sigma2; /* standard deviation of s2*/
  const unsigned int log2sigma2; /* sigma1 = 2^log2sigma2 */
  const int_srcptr fraks1; /* sqrt(2)fraks1 is the standard deviation of the mask y1*/
  const unsigned int log2stdev1; /* fraks1 = 2^log2stdev1 */
  const int_srcptr fraks2; /* sqrt(2)fraks2 is the standard deviation of the mask y2*/
  const unsigned int log2stdev2; /* fraks2 = 2^log2stdev2 */
} modified_abdlop_params_struct;
typedef modified_abdlop_params_struct modified_abdlop_params_t[1];
typedef modified_abdlop_params_struct *modified_abdlop_params_ptr;
typedef const modified_abdlop_params_struct *modified_abdlop_params_srcptr;

typedef struct
{
  const modified_abdlop_params_srcptr quad_eval;
  const modified_abdlop_params_srcptr quad_many;
  const unsigned int lambda;

} modified_quad_eval_params_struct;
typedef modified_quad_eval_params_struct modified_quad_eval_params_t[1];
typedef modified_quad_eval_params_struct *modified_quad_eval_params_ptr;
typedef const modified_quad_eval_params_struct *modified_quad_eval_params_srcptr;

void poly_addscale2 (poly_t r, const int_t a, poly_t b, int crt);
void poly_subscale2 (poly_t r, const int_t a, poly_t b, int crt);

void spolyvec_get_subvec (spolyvec_t subvec, spolyvec_t vec, unsigned int elem_init, unsigned int nelems);

void spolymat_get_submat (spolymat_ptr r, spolymat_ptr a, const unsigned int nrows_start, const unsigned int ncols_start);
void spolymat_get_submat_upperdiag (spolymat_ptr r, spolymat_ptr a, const unsigned int nrows_start, const unsigned int ncols_start);
void spolymat_set_squared_block_matrix (spolymat_ptr r, spolymat_ptr a1, spolymat_ptr a2, spolymat_ptr a3, spolymat_ptr a4);
void spolymat_set_upper_block_matrix (spolymat_ptr r, spolymat_ptr a1, spolymat_ptr a2, spolymat_ptr a3);
void spolymat_set_block_matrix (spolymat_ptr r, spolymat_ptr a1, spolymat_ptr a2);
void spolymat_neg (spolymat_t r, spolymat_t b);

void modified_abdlop_keygen (polymat_t A1, polymat_t A2prime, polymat_t Bprime,
               const uint8_t seed[32], const modified_abdlop_params_t params);

void modified_abdlop_commit (polyvec_t tA1, polyvec_t tA2, polyvec_t tB, polyvec_t s1, polyvec_t rand_enc_s1,
               polyvec_t m, polyvec_t s2, polymat_t A1, polymat_t A2prime,
               polymat_t Bprime, const uint8_t seed[32], const modified_abdlop_params_t params);

void modified_abdlop_enccomm (uint8_t *buf, size_t *buflen, polyvec_t tA1, polyvec_t tB,
                const modified_abdlop_params_t params);

void modified_abdlop_hashcomm (uint8_t hash[32], polyvec_t tA1, polyvec_t tB,
                 const modified_abdlop_params_t params);

void modified_abdlop_prove (uint8_t hash[32], poly_t c, polyvec_t z1, polyvec_t z21,
              polyvec_t h, polyvec_t tA2, polyvec_t rand_enc_s1, polyvec_t s2,
              polymat_t A1, polymat_t A2prime, const uint8_t seed[32],
              const modified_abdlop_params_t params);

int modified_abdlop_verify (uint8_t hash[32], poly_t c, polyvec_t z1, polyvec_t z21,
              polyvec_t h, polyvec_t tA1, polymat_t A1, polymat_t A2prime,
              const modified_abdlop_params_t params);

void modified_quad_prove (uint8_t hash[32], polyvec_t tB, poly_t c, polyvec_t z1,
              polyvec_t z21, polyvec_t h, polyvec_t randencs1, polyvec_t m,
              polyvec_t s2, polyvec_t tA2, polymat_t A1,
              polymat_t A2prime, polymat_t Bprime, spolymat_t R2,
              spolyvec_t r1, const uint8_t seed[32],
              const modified_abdlop_params_t params);
int modified_quad_verify (uint8_t hash[32], poly_t c, polyvec_t z1, polyvec_t z21,
                     polyvec_t h, polyvec_t tA1, polyvec_t tB, polymat_t A1,
                     polymat_t A2prime, polymat_t Bprime, spolymat_t R2,
                     spolyvec_t r1, poly_t r0, const modified_abdlop_params_t params);

void modified_quad_many_prove (uint8_t hash[32], polyvec_t tB, poly_t c,
                          polyvec_t z1, polyvec_t z21, polyvec_t h,
                          polyvec_t randencs1, polyvec_t m, polyvec_t s2,
                          polyvec_t tA2, polymat_t A1, polymat_t A2prime,
                          polymat_t Bprime, spolymat_ptr R2i[],
                          spolyvec_ptr r1i[], unsigned int N,
                          const uint8_t seed[32],
                          const modified_abdlop_params_t params);
int modified_quad_many_verify (uint8_t hash[32], poly_t c, polyvec_t z1,
                          polyvec_t z21, polyvec_t h, polyvec_t tA1,
                          polyvec_t tB, polymat_t A1, polymat_t A2prime,
                          polymat_t Bprime, spolymat_ptr R2i[],
                          spolyvec_ptr r1i[], poly_ptr r0i[], unsigned int N,
                          const modified_abdlop_params_t params);

void modified_quad_eval_prove (uint8_t hash[32], polyvec_t tB, polyvec_t h,
                          poly_t c, polyvec_t z1, polyvec_t z21,
                          polyvec_t hint, polyvec_t randencs1, polyvec_t m,
                          polyvec_t s2, polyvec_t tA2, polymat_t A1,
                          polymat_t A2prime, polymat_t Bprime,
                          spolymat_ptr R2i[], spolyvec_ptr r1i[],
                          unsigned int N, spolymat_ptr Rprime2i[],
                          spolyvec_ptr rprime1i[], poly_ptr rprime0i[],
                          unsigned int M, const uint8_t seed[32],
                          const modified_quad_eval_params_t params);
int modified_quad_eval_verify (uint8_t hash[32], polyvec_t h, poly_t c,
                          polyvec_t z1, polyvec_t z21, polyvec_t hint,
                          polyvec_t tA1, polyvec_t tB, polymat_t A1,
                          polymat_t A2prime, polymat_t Bprime,
                          spolymat_ptr R2i[], spolyvec_ptr r1i[],
                          poly_ptr r0i[], unsigned int N,
                          spolymat_ptr Rprime2i[], spolyvec_ptr rprime1i[],
                          poly_ptr rprime0i[], unsigned int M,
                          const modified_quad_eval_params_t params);

void print_stopwatch_modified_quad_eval_prove (unsigned int indent);
void print_stopwatch_modified_quad_eval_verify (unsigned int indent);
void print_stopwatch_modified_quad_many_prove (unsigned int indent);
void print_stopwatch_modified_quad_many_verify (unsigned int indent);
void print_stopwatch_modified_quad_prove (unsigned int indent);
void print_stopwatch_modified_quad_verify (unsigned int indent);

void decompose(polyvec_t output, poly_t input, const int_t base, const int_t degree);
void regular_encoding(polyvec_t output, polyvec_t input, const int_t base);
void randomized_encoding(polyvec_t output, polyvec_t input,  const uint8_t seed[32], const unsigned int log2sigma, const int_t base);

__END_DECLS
#endif