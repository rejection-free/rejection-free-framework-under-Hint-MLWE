# Implementation of a Rejection-free Zero-Knowledge framework of proof

Modification of the [LaZer library](https://github.com/lazer-crypto/lazer) in order to add our implementation 
for the paper **Rejection free Zero-Knowledge framework of proof under Hint-MLWE**. 
We do not claim any existing structures from the original library related to the [LSS24](https://eprint.iacr.org/2024/1846.pdf). 
Our implementation only concerns the use of the library's functions 
to modify existing frameworks accordingly.

The repository contains the implementation of our proof framework 
using the modified commitment scheme. It exploits the encoding protocol extracted from [HSS24](https://eprint.iacr.org/2024/306.pdf). The repository also contains several parameter sets in the `tests` sub-directory, defined in the files `params*.sage`. 

Using the `scripts` files, the user can manually generate the complete list of concrete parameters associated with these inputs for different uses: 
- a single instance of the commitment scheme,
- one instance of a quadratic function
- an instance of several quadratic functions
- an instance of several evaluations (i.e. a constant coefficient equal to zero).

Clone the repository
--------------------

As our code exploits the actual structure of the code from the [LaZer library](https://github.com/lazer-crypto/lazer) and from the [LWE-estimator](https://github.com/malb/lattice-estimator).
we manage to use and link the github repository from the Library directly in ours.
Then, in order to have the entire list of files from our code with the libraries, run:

`git clone --recursive git@github.com:rejection-free/rejection-free-framework-under-Hint-MLWE.git`

Be aware of the dependencies
----------------------------
We provide the list of our hardware and software used to build and run our code, and their version.

- ubuntu 22.04
- kernel version 6.8.0-51-generic
- gcc version 11.4.0
- make 4.3
- cmake version 3.22.1
- sagemath version 9.5
- python3 3.10,

If you want to fully build the entire [LaZer library](https://github.com/lazer-crypto/lazer), keep in mind that 
the entire list of requirements are differents and the library must be built independently from this code.

Compile the code
----------------
To build the C library along with the required LaZer parts, 
from the base directory, run:

`make all`

To build the list of parameters specified in the `params` subfolder, from the base directory, run:

`make params`

or specifically `make params-abdlop`, `make params-quad`, `make params-eval`.

Keep in mind that this will overwrite the actual associated files `modified-__-__-params%.h` in the `tests` subdirectory.

To build the C list of tests, from the base directory, run:

`make check`

Optional: use make's `-j` option to speed up compilation.

In order to clean your repository after compilation, run:

`make clean`

Concrete modification of the library
------------------------------------

In order to give a simpler library, we have chosen to 
select and modify the list of files compiled in the [LaZer library](https://github.com/lazer-crypto/lazer), as it was required to have access to a CPU that can  use the AVX512 set of instructions. 

In general, we base our implementation in order to be joined to
the actual [LaZer library](https://github.com/lazer-crypto/lazer), reusing the global structure, types, 
and semantics.   

We also use and modify in consequence the scripts in `scripts` to provide the 
security of our different protocols. 

LWE Estimator
-------------

We imported the [LWE estimator](https://github.com/malb/lattice-estimator) as the subfolder `scripts/lattice-estimator` from 

    Martin R. Albrecht, Rachel Player and Sam Scott. On the concrete hardness of Learning with Errors.
    Journal of Mathematical Cryptology. Volume 9, Issue 3, Pages 169–203, ISSN (Online) 1862-2984,
    ISSN (Print) 1862-2976 DOI: 10.1515/jmc-2015-0016, October 2015

The first time you launch `make params` it will extract the subfolder `estimator` for this repository providing the python files 
to the sage scripts. 

Parameters
----------

We provide in the `params` subfolder our 5 sets of parameters. In the `tests` subfolder, you will have the entire list 
of parameters `modified-_-_-params%.h` for each possible instance of our framework. 

Using the `scripts` files, the user can manually generate the complete list of concrete parameters associated with these inputs for different uses: 
- an instance of the commitment scheme for the 5 parameters sets,
- an instance of the commitment scheme + a single quadratic relations for the 5 parameters sets,
- an instance of the commitment scheme + 3 quadratic relations for the 5 parameters sets,
- an instance of the commitment scheme + 3 quadratic relations + 3 evaluations for the 5 parameters sets.