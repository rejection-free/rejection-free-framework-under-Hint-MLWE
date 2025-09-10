import mpmath as mp
from mpmath import mpf, nstr
import sys
import numpy
import time
from sage.stats.distributions.discrete_gaussian_polynomial import DiscreteGaussianDistributionPolynomialSampler
from estimator import *

# We do not claim the knowledge of this script and only modify the existant ones from the script subdirectory from the LaZer Library.

mp.mp.prec = 512 # precision for mp
prec = 8  # precision for nstr

load("common_code.sage")
blockPrint()
loaded = 1
verbose = 1
code = 1
# alpha_1 and alpha_2 (bounds for the binding property)
def alpha_1():
    global sigma_1
    global m_1
    global d
    global b
    global k
    return (b + 1) * sigma_1 * mp.sqrt(k * m_1 * d)

def alpha_2():
    global sigma_2
    global m_2
    global d
    global b
    return (b + 1) * sigma_2 * mp.sqrt(m_2 * d)

# beta_1 and beta_2 (bounds for the knowledge soudness)
def beta_1():
    global sigma_1
    global m_1
    global d
    global b
    global k
    global eta
    global frak_s_1
    return (b + 1) * (eta * sigma_1 + mp.sqrt(2) * frak_s_1) * mp.sqrt(k * m_1 * d)

def beta_2():
    global sigma_2
    global m_2
    global d
    global b
    global eta
    global frak_s_2

    return (eta * sigma_2 + frak_s_2) * mp.sqrt(m_2 * d) + mpf(eta) * 2 ** (D-1) * mp.sqrt(kmsis*d) + (gamma * mp.sqrt(kmsis * d))/mpf(2)

# varsigma_1 and varsigma_2 (simulatability)

def varsigma_1():
    global sigma_1
    global eta
    global frak_s_1
    return ceil(mp.sqrt(1/2 * (1/(sigma_1^2) + (eta^2)/(mpf(2)*frak_s_1^2) )^(-1)))

def varsigma_2():
    global sigma_2
    global eta
    global frak_s_2
    return ceil(mp.sqrt(1/2 * (1/(sigma_2^2) + (eta^2)/(frak_s_2^2) )^(-1)))

# bound on the extracted MSIS solution (Binding)
def bound_binding():
    return 2 * mp.sqrt(alpha_1() ** 2 + alpha_2() ** 2)

# bound B on the extracted MSIS solution (Know. Soundness)
def bound_know_soundness():
    global eta
    return 8 * mpf(eta) * mp.sqrt(beta_1() ** 2 + beta_2() ** 2)

# Early checks to launch the scripts

assert len(sys.argv) == 2
params_file = sys.argv[1]

### security parameters
LAMBDA = 128 
ROOT_HERMITE_128 = 1.0044 # 128-bit security 
param_sec = mpf(2^(-LAMBDA))

load(params_file)
name = "rf_"+ name

m_1 = m_1 + l # We keep l as the script of parameters can be used for the lnp framework from the LaZer Library
            # We just show here that we do not need to separate bounded and unbounded message: setting m_1 to (m_1 + l) and l to 0
            # Third point: we enhance k * as we define the commitment scheme over the randomized encoding of size k * m_1
            # and need to provide security for this length and not for m_1 only      
l = 0
### choices
k = 2
b = 2^(log2q//k)
p = b ** k + 1
lext = 0


if d not in [64, 128]:
    err("d not in [64,128]")
log2d = mp.log(d, 2)

D = 0       # dropping low-order bits of t_A
gamma = 0   # dropping low-order bits of w

# number of irreducible factors of X^d + 1 modulo each q_i,  q_i = 2l+1 (mod 4l)
L = 2
NADDS = 128  # chose P big enough for this many additions

# challenge space
if d == 64 and L == 2 and log2q >= 4:
    omega = 8
    eta = 140
    Csize = 2 ** 129
elif d == 128 and L == 2 and log2q >= 4:
    omega = 2
    eta = 59
    Csize = 2 ** 147
else:
    err("challenge space undefined")

# sample from [-omega,omega] <=> sample from [0,2*omega] - omega
omega_bits = ceil(mp.log(2*omega+1, 2))

# Relations parameters : standard deviations and bounds
lambda_classic = 1
lambda_ideal_p = mp.sqrt(p)
# the distribution of the message s_1 is not restricted

### PARAMETERS FOR S1

sigma_1 = smoothing_param_bound(k * m_1 * d, lambda_classic, param_sec) # SD of the randomized encoding function

frak_s_1 = max(smoothing_param_bound(k * m_1 * d, lambda_classic, param_sec)/mp.sqrt(2),mp.sqrt(2)/(b-1) * smoothing_param_bound(m_1 * d, lambda_ideal_p,param_sec))        # SD of the randomized encoding function of y_1

while mp.sqrt(2) * (b-1) * varsigma_1() < smoothing_param_bound(k * m_1 * d, lambda_ideal_p,param_sec):
    frak_s_1 *= mp.sqrt(2)

### PARAMETERS FOR S2

kmlwe = 0           # MLWE dim, to be determined
easy_mlwe_dim = 0   # lower bound for MLWE dim
hard_mlwe_dim = 64  # guess for upper bound for MLWE dim

sigma_2 = mpf(2)*smoothing_param_bound(hard_mlwe_dim*d, lambda_classic, param_sec)                            # SD of the randomness s_2 : later (depends on length of randomness s2)
frak_s_2 = mp.sqrt(2)*smoothing_param_bound(hard_mlwe_dim*d, lambda_classic, param_sec)                           # SD of y_2 : later (depends on length of randomness s2)

# find upper actual bound (and possibly improve lower bound)
while True:
    sigma_2 = mpf(2)*smoothing_param_bound(hard_mlwe_dim*d, lambda_classic, param_sec)                      
    frak_s_2 = mp.sqrt(2)*smoothing_param_bound(hard_mlwe_dim*d, lambda_classic, param_sec)
    while varsigma_2() < mp.sqrt(2)*smoothing_param_bound(hard_mlwe_dim*d, lambda_classic, param_sec):
        frak_s_2 *= mp.sqrt(2)                       
    delta_mlwe = max(findMLWEdelta(hard_mlwe_dim, d, 2 ** log2q, sigma_2),findMLWEdelta(hard_mlwe_dim, d, 2 ** log2q, varsigma_2()))

    if delta_mlwe <= ROOT_HERMITE_128:
        print(f"MLWE dim {kmlwe} : hard")
        break
    print(f"MLWE dim {kmlwe} : easy")
    easy_mlwe_dim = hard_mlwe_dim
    hard_mlwe_dim *= 2

# binary search for smallest MLWE dimension that is still hard
while True:
    kmlwe = (easy_mlwe_dim + hard_mlwe_dim) / 2
    sigma_2 = mpf(2)*smoothing_param_bound(kmlwe*d, lambda_classic, param_sec)
    frak_s_2 = mp.sqrt(2)*smoothing_param_bound(kmlwe*d, lambda_classic, param_sec)
    while varsigma_2() < mp.sqrt(2)*smoothing_param_bound(kmlwe*d, lambda_classic, param_sec):
        frak_s_2 *= mp.sqrt(2)                       
    delta_mlwe = max(findMLWEdelta(kmlwe, d, 2 ** log2q, sigma_2),findMLWEdelta(kmlwe, d, 2 ** log2q, varsigma_2()))
    
    if delta_mlwe <= ROOT_HERMITE_128:
        print(f"MLWE dim {kmlwe} : hard")
        hard_mlwe_dim = kmlwe
    else:
        print(f"MLWE dim {kmlwe} : easy")
        easy_mlwe_dim = kmlwe
    if hard_mlwe_dim == easy_mlwe_dim + 1:
        kmlwe = hard_mlwe_dim
        print(f"found MLWE dim : {kmlwe}")
        break

# Find an appropriate Module-SIS dimension
kmsis = 0   # dimension of the MSIS problem
while True:
    kmsis += 1
    m_2 = kmlwe + kmsis
    
    print(f"d {d}")
    print(f"2^log2q {2^log2q}")
    print(f"kmsis {kmsis}")
    print(f"Bound {bound_know_soundness()}")
    print(f"delta {get_delta_msis(bound_know_soundness(), kmsis, d, 2 ** log2q)}")
    print(f"delta {get_delta_msis(bound_binding(), kmsis, d, 2 ** log2q)}")
    if get_delta_msis(bound_know_soundness(), kmsis, d, 2 ** log2q) < ROOT_HERMITE_128 and bound_know_soundness() < 2 ** log2q and get_delta_msis(bound_binding(), kmsis, d, 2 ** log2q) < ROOT_HERMITE_128 and bound_binding() < 2 ** log2q:
        break

# Find the largest possible gamma which makes the MSIS solution still small.
gamma = 2 ** log2q 
while True:       # searching for right gamma
    gamma /= 2
    if get_delta_msis(bound_know_soundness(), kmsis, d, 2 ** log2q) < ROOT_HERMITE_128 and bound_know_soundness() < 2 ** log2q and get_delta_msis(bound_binding(), kmsis, d, 2 ** log2q) < ROOT_HERMITE_128 and bound_binding() < 2 ** log2q:
        break

# Finding exact values for q, b and gamma:
true_gamma_found = false    # Boolean for finding correct gamma
b = 2 ** (log2q//k) - 1             
while true_gamma_found == false:
    b = b+1
    p = b^2 + 1  
    if is_prime(p):     # we need p to be prime
        if p%8 == 5:    # we need p to be congruent to 5 modulo 8
            div_b = divisors(p-1)     # consider divisors of b^2
            for i in div_b:                
                if gamma*4/5 < i and i <= gamma and is_even(i): # find a divisor which is close to gamma
                    gamma = i # we found a good candidate for gamma
                    true_gamma_found = true

m = (p-1) / gamma

# Find the largest possible D which makes the MSIS solution small
D = log2q
while D != 0:
    D -= 1
    if get_delta_msis(bound_know_soundness(), kmsis, d, p) < ROOT_HERMITE_128 and bound_know_soundness() < p and get_delta_msis(bound_binding(), kmsis, d, p) < ROOT_HERMITE_128 and bound_binding() < p and 2 ** (D-1)*omega*d < gamma:
        break

# assert all the conditions

assert sigma_1 >= smoothing_param_bound(k*m_1*d, lambda_classic, param_sec)
assert mp.sqrt(2) * frak_s_1 >= smoothing_param_bound(k*m_1*d, lambda_classic, param_sec)
assert mpf(b-1) * frak_s_1 >= mp.sqrt(2) *smoothing_param_bound(k*m_1*d, lambda_ideal_p, param_sec)
assert mp.sqrt(2) * (b-1) * varsigma_1() >=  smoothing_param_bound(k*m_1*d, lambda_ideal_p, param_sec)

assert sigma_2 >= mpf(2)*smoothing_param_bound(kmlwe*d, lambda_classic, param_sec)
assert frak_s_2 >= mp.sqrt(2)*smoothing_param_bound(kmlwe*d, lambda_classic, param_sec)
assert varsigma_2() >= mp.sqrt(2)*smoothing_param_bound(kmlwe*d, lambda_classic, param_sec)

# update MLWE root hermite factor with exact p
delta_mlwe = max(findMLWEdelta(kmlwe, d, p, varsigma_2()),findMLWEdelta(kmlwe, d, p, sigma_2))

# computation of the proof size 

logp = ceil(mp.log(p,2))

full_size = kmsis * d * (logp - D) + lext * d * logp 

hint = 2.25 * kmsis * d

challenge = ceil(mp.log(2*omega+1,2)) * d 

length_z_1 = (k * m_1) * d * (ceil(mp.log( eta * (b+1) * sigma_1 + mp.sqrt(2)*frak_s_1,2)) + 2.57)

length_z_2 = m_2 * d * (ceil(mp.log( eta * sigma_2 + frak_s_2,2)) + 2.57)

enablePrint()
printv(f"auto-generated by rf-abdlop-codegen.sage from {params_file}.")
printv(f"")

if not (kmlwe >= 0 and kmlwe == m_2 - kmsis):
    err("protocol not simulatable because of the parameters")

printv(
    f"the commitment scheme is binding under MSIS({kmsis},{k*m_1 + m_2}) with bound={nstr(bound_binding(), prec)})") 

printv(
    f"the commitment scheme is hiding under MLWE({kmsis},{kmlwe}) with sd={nstr(sigma_2, prec)})") 

printv(
    f"protocol is simulatable under Hint-MLWE implied by MLWE({kmsis},{kmlwe}) with sd={nstr(varsigma_2(), prec)})") 

eknow = mpf(1)/mpf(Csize)
printv(
    f"protocol is knowledge-sound with knowledge error <= 2^({nstr(mp.ceil(mp.log(eknow,2)),prec)}) under MSIS({kmsis},{k*m_1 + m_2}) with bound={nstr(bound_know_soundness(), prec)})")

# print params
printv(f"")
printv(f"Ring")
printv(f"degree d = {d}")
printv(f"modulus p = {p}, log(q) ~ {nstr(mp.log(p,2),prec)}")
printv(f"base decomposition b = {b}")
printv(f"")
printv(f"Compression")
printv(f"D = {D}")
printv(f"gamma = {gamma}, log(gamma) ~ {nstr(mp.log(gamma,2),prec)}")
printv(f"m = (q-1)/gamma = {m}, log(m) ~ {nstr(mp.log(m,2),prec)}")
printv(f"")
printv(f"Dimensions of secrets")
printv(f"s1: m_1 = {(m_1)}")
printv(f"s2: m_2 = {m_2}")
printv(f"")
printv(f"Size of secrets")
printv(f"s1 unbounded")
printv(f"s2 a gaussian element distributed with standard deviation {nstr(sigma_2, prec)}")
printv(f"")
printv(f"Challenge space")
printv(
    f"c uniform in [-omega,omega] = [{-omega},{omega}], o(c)=c, sqrt(l1(o(c)*c)) <= eta = {eta}")
printv(f"")
printv(f"")
printv(f"Security")
printv(f"MSIS dimension: {kmsis}")
printv(
    f"MSIS root hermite factor: {nstr(get_delta_msis(bound_know_soundness(), kmsis, d, p), prec)}")
printv(f"MLWE dimension: {kmlwe}")
printv(f"MLWE root hermite factor: {nstr(mpf(delta_mlwe), prec)}")
printv(f"")
printv(f"Proof size of the rf abdlop opening ")

printv(f"Total proof size in KB:  { round((full_size + challenge + length_z_1 + length_z_2 + hint)/(2^13) , 2)}")
printv(f"full-sized polynomials in KB: {round(full_size/(2^13) , 2)}")
printv(f"challenge c in KB: {round(challenge/(2^13) , 2)}")
printv(f"short-sized polynomials in KB: {round((length_z_1 + length_z_2 + hint)/(2^13) , 2)}")
printv(f"")

q_nlimbs = int2limbs(p, -1)[1]
if m % 2 == 0:
    mby2 = m / 2
else:
    mby2 = 0

minP = min_P(d, p, NADDS)
moduli = moduli_list(nbit, d, minP)[0]
P = prod(moduli)
assert P >= minP
nmoduli = len(moduli)
Pmodq = redc(Mod(P, p), p)
Ppmodq = []
Ppmodq_str = []
for i in range(len(moduli)):
    Ppmodq_str += [f"{name}_Ppmodq_{i}"]
    Ppmodq += [redc(Mod(P/moduli[i], p), p)]
Ppmodq_array = strlist2ptrarray(Ppmodq_str)

out = ""
out += f"""
#include "lazer.h"
{int_t(f"{name}_q", p)}
{int_t(f"{name}_qminus1", p - 1)}
{int_t(f"{name}_b", b, q_nlimbs)}
{int_t(f"{name}_m", m, q_nlimbs)}
{int_t(f"{name}_mby2", mby2, q_nlimbs)}
{int_t(f"{name}_gamma", gamma, q_nlimbs)}
{int_t(f"{name}_gammaby2", gamma / 2, q_nlimbs)}
{int_t(f"{name}_pow2D", 2^D, q_nlimbs)}
{int_t(f"{name}_pow2Dby2", 2^D / 2, q_nlimbs)}
{int_t(f"{name}_Bsq", floor(bound_know_soundness()^2), 2*q_nlimbs)}
{int_t(f"{name}_sigma_1", int(mp.nint(sigma_1^2)), 2*q_nlimbs)}
{int_t(f"{name}_sigma_2", int(mp.nint(sigma_2^2)), 2*q_nlimbs)}
{int_t(f"{name}_frak_s1", int(mp.nint(frak_s_1^2)), 2*q_nlimbs)}
{int_t(f"{name}_frak_s2", int(mp.nint(frak_s_2^2)), 2*q_nlimbs)}
{int_t(f"{name}_inv2", redc(1/2 % p, p))}
{int_t(f"{name}_Pmodq", Pmodq, q_nlimbs)}
"""
for i in range(len(moduli)):
    out += int_t(f"{name}_Ppmodq_{i}", Ppmodq[i], q_nlimbs) + f"\n"
out += f"""
static const int_srcptr {name}_Ppmodq[] = {Ppmodq_array};
static const polyring_t {name}_ring = {{{{{name}_q, {d}, {ceil(mp.log(p-1,2))}, {log2d}, moduli_d{d}, {nmoduli}, {name}_Pmodq, {name}_Ppmodq, {name}_inv2}}}};
static const dcompress_params_t {name}_dcomp = {{{{ {name}_q, {name}_qminus1, {name}_m, {name}_mby2, {name}_gamma, {name}_gammaby2, {name}_pow2D, {name}_pow2Dby2, {D}, {m % 2}, {ceil(mp.log(m,2))} }}}};
static const rf_abdlop_params_t {name} = {{{{ {name}_ring, {name}_dcomp, {name}_b, {m_1}, {m_2},{0}, {lext}, {kmsis}, {name}_Bsq, {omega}, {omega_bits}, {eta}, {name}_sigma_1, {ceil(mp.log(sigma_1,2))}, {name}_sigma_2, {ceil(log(sigma_2,2))}, {name}_frak_s1, {ceil(mp.log(frak_s_1,2))}, {name}_frak_s2, {ceil(mp.log(frak_s_2,2))}}}}};
"""

printc(out)

sys.exit(int(0))

# Original scripts and proofs are coming from :
# [1] Lattice-Based Zero-Knowledge Proofs Under a Few Dozen Kilobytes
# https://doi.org/10.3929/ethz-b-000574844
