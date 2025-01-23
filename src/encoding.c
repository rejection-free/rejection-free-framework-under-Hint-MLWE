#include "lazer.h"
#include "stopwatch.h"

/*
Encode a single element input in Rd into output in Rd^2 with b-base decomposition
*/
void decompose(polyvec_t output, poly_t input, const int_t base, const int_t degree)
{
    size_t i;
    intvec_ptr i32;
    polyring_srcptr rq = polyvec_get_ring(output);
    const unsigned int d = degree->limbs[0];
    INTVEC_T (o32,2*d,rq->q->nlimbs);
    INTVEC_T (c32,2*d,rq->q->nlimbs);
    INTVEC_T (output_coeff,2*d,rq->q->nlimbs);
    int_t check;
    int_alloc(check,rq->q->nlimbs);
    INT_T (temp0, rq->q->nlimbs); 
    INT_T (temp1, rq->q->nlimbs); 
    int_ptr v;
    int_set_i64 (check, -1); 
    intvec_set_zero(o32);
    intvec_set_zero(c32);
    intvec_set_zero(output_coeff);
    i32 = poly_get_coeffvec(input);

    for (i = 0; i < d; i++)
      {
        v = intvec_get_elem (i32, i);
        int_set_zero(temp0);
        int_set_one(temp1);
        int_set_i64 (check, -1);
        if (int_eq (v, check) == 1)
          {
            int_set(temp1, base);
          }
        else
        {
            if (v->neg == 1)
            {
              int_add(v ,v , rq->q);
            }
            int_div(temp1,temp0,v,base);
        }
        
        intvec_set_elem(o32,2*i, temp0);
        intvec_set_elem(o32,2*i+1, temp1);
      }
    
    int_set_i64 (temp0, 2);  
    int_div(check, temp1,base,temp0);
    for (i = 0; i < 2 * d; i++)
      {
        v = intvec_get_elem (o32, i);
        int_gt(v,check);
        if (int_gt(v,check) == 1)
          {
            int_sub(v,v,base);
            v = intvec_get_elem (c32, i);
            int_set_i64( v, 1);
          }
      }
    int_free(check);
    for (i = 0; i < 2*d; i++)
      {
        v = intvec_get_elem (output_coeff, i);
        if (i < d)
        {
          int_sub(v,intvec_get_elem(o32, 2*i), intvec_get_elem (c32, 2*i+1));
        } else {
          int_add(v,intvec_get_elem(o32, 2*(i-d)+1), intvec_get_elem (c32, 2*(i-d)));
        }
        
    }

    polyvec_set_coeffvec(output, output_coeff);
    
}

/*
Regular encoding that takes in input a vector s1 of size m1, and output a vector output of size 2*m1 such that the output is the gadget decomposition in base b of the 
*/
void regular_encoding(polyvec_t output, polyvec_t input, const int_t base)
{
    polyring_srcptr rq = polyvec_get_ring(input);
    INT_T (d, 1);
    int_set_i64 (d, rq->d);

    const unsigned int m1 = polyvec_get_nelems(input);
    
    POLYVEC_T (temp_out, rq, 2);
    
    poly_ptr temp_in;
    for (size_t i = 0; i < m1; i++)
    {
        polyvec_get_subvec (temp_out, output, 2*i, 2, 1);
        temp_in = polyvec_get_elem(input, i);
        decompose(temp_out, temp_in, base,d);
    }
}

/*
Randomized encoding that takes in input a vector s1 of size m1, and output a vector output of size 2*m1 with standard deviation sigma
*/
void randomized_encoding(polyvec_t output, polyvec_t input,  const uint8_t seed[32], const unsigned int log2sigma, const int_t base)
{
    polyring_srcptr rq = polyvec_get_ring(input);
    uint32_t dom;
    rng_state_t rngstate;
    uint8_t eseed[32]; /* error seed */
    const unsigned int m1 = polyvec_get_nelems(input);
    POLYVEC_T (error, rq, 2*m1);
    POLYVEC_T (multiplied_error, rq, 2*m1);
    POLYVEC_T (temp_error, rq, 2);
    poly_ptr tempe;
    POLY_T (temp0, rq);
    POLY_T (temp1, rq);
    rng_init (rngstate, seed, 0);
    rng_urandom (rngstate, eseed, 32);
    regular_encoding(output, input, base);
    dom = 0;
    polyvec_grandom(error, log2sigma, eseed, dom);
    for (size_t i = 0; i < m1; i++)
    {
        poly_set(temp0, polyvec_get_elem_src(error,2*i));
        poly_set(temp1, polyvec_get_elem_src(error,2*i+1));
        tempe = polyvec_get_elem(multiplied_error,2*i);
        poly_set(tempe, temp1);
        poly_addscale2(tempe, base, temp0,0);
        poly_neg_self(tempe);

        tempe = polyvec_get_elem(multiplied_error,2*i+1);
        poly_set(tempe, temp0);
        poly_subscale2(tempe, base, temp1,0);
    }
    polyvec_add(output, output, multiplied_error, 0);
}

/*
Decoding function : TODO
*/
// void decoding(polyvec_t output, polyvec_t input, const modified_abdlop_params_t params)
// {
//     polyring_srcptr rq = params->ring;
//     int_srcptr b = params->base;
//     const unsigned int m1 = params->m1;
// }