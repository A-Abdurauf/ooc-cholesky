/**
 *
 * @file cpltmg.c
 *
 *  PLASMA computational routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @generated c Fri Apr  1 11:02:55 2016
 *
 **/
#include "common.h"

/***************************************************************************/
/**
 *
 * @ingroup PLASMA_Complex32_t
 *
 *  PLASMA_cpltmg - Generate a test matrix by tiles.
 *
 *******************************************************************************
 *
 * @param[in] mtxtype
 ***** @arg PlasmaMatrixCauchy:
 *     See http://www.mathworks.fr/fr/help/matlab/ref/gallery.html#f84-1019317
 *
 *     Cauchy matrix
 *
 *     Returns an n-by-n matrix C such that, C(i,j) = 1/(i + j).
 *
 ***** @arg PlasmaMatrixChebvand:
 *     See http://www.mathworks.fr/fr/help/matlab/ref/gallery.html#f84-999859
 *
 *     Vandermonde-like matrix for the Chebyshev polynomials
 *
 *     Produces the (primal) Chebyshev Vandermonde matrix based on the vector of
 *     points p, which define where the Chebyshev polynomial is calculated.
 *
 *     If seed != 0, C(i,j) = Ti – 1(p(j)) where Ti – 1 is the Chebyshev
 *     polynomial of degree i – 1, and p is a vector of N equally spaced points on
 *     the interval [0,1].
 *
 ***** @arg PlasmaMatrixCircul:
 *     See http://www.mathworks.fr/fr/help/matlab/ref/gallery.html#f84-999880
 *
 *     Circulant matrix
 *
 *     A circulant matrix has the property that each row is obtained from the
 *     previous one by cyclically permuting the entries one step forward. It is
 *     a special Toeplitz matrix in which the diagonals "wrap around."
 *
 *     The eigensystem of C (n-by-n) is known explicitly: If t is an nth root of
 *     unity, then the inner product of v and w = [1 t t2 ... t(n – 1)] is an
 *     eigenvalue of C and w(n:-1:1) is an eigenvector, where v is the first
 *     column of C.
 *
 ***** @arg PlasmaMatrixCompan:
 *     See http://www.mathworks.fr/fr/help/matlab/ref/compan.html
 *
 *     Companion matrix
 *
 *     A = compan(u) returns the corresponding companion matrix whose first row is
 *     -u(2:n)/u(1), where u is a vector of polynomial coefficients. The
 *     eigenvalues of compan(u) are the roots of the polynomial.
 *
 ***** @arg PlasmaMatrixCondex:
 *     See http://www.mathworks.fr/fr/help/matlab/ref/gallery.html#f84-999898
 *     gallery('condex',n,4,100)
 *
 *     Returns a "counter-example" matrix to a condition estimator. It has order n
 *     and scalar parameter theta (default 100).
 *
 *     LAPACK (RCOND): It is the inverse of this matrix that is a counter-example.
 *
 ***** @arg PlasmaMatrixDemmel:
 *     See [1] J. Demmel, Applied Numerical Linear Algebra, SIAM,
 *             Philadelphia, 1997
 *
 *     Returns a matrix defined by:
 *        A = D * ( I + 1e-7* rand(n)), where D = diag(10^{14*(0:n-1)/n})
 *
 ***** @arg PlasmaMatrixDorr:
 *     See http://www.mathworks.fr/fr/help/matlab/ref/gallery.html#f84-999936
 *
 *     Diagonally dominant, ill-conditioned, tridiagonal matrix
 *
 *     Returns the n-by-n matrix, row diagonally dominant, tridiagonal
 *     matrix that is ill-conditioned for small nonnegative values of
 *     theta. The default value of theta is 0.01. The Dorr matrix
 *     itself is the same as gallery('tridiag',c,d,e).
 *
 ***** @arg PlasmaMatrixFiedler:
 *     See http://www.mathworks.fr/fr/help/matlab/ref/gallery.html#f84-999960
 *
 *     Fiedler matrix of size n-by-n is defined throug a random vector c
 *     of size n, such that each element is equal to abs(n(i)-n(j)).
 *
 *     Matrix A has a dominant positive eigenvalue and all the other
 *     eigenvalues are negative.
 *
 *     Explicit formulas for inv(A) and det(A) are given in
 *     [Todd, J., Basic Numerical Mathematics, Vol. 2: Numerical Algebra,
 *     Birkhauser, Basel, and Academic Press, New York, 1977, p. 159] and
 *     attributed to Fiedler. These indicate that inv(A) is tridiagonal
 *     except for nonzero (1,n) and (n,1) elements.
 *
 ***** @arg PlasmaMatrixFoster:
 *     See [1] L. V. Foster, Gaussian Elimination with Partial
 *             Pivoting Can Fail in Practice, SIAM J. Matrix
 *             Anal. Appl., 15 (1994), pp. 1354-1362.
 *
 *         [2] L. V. Foster, The growth factor and efficiency of
 *             Gaussian elimination with rook pivoting,
 *             J. Comput. Appl. Math., 86 (1997), pp. 177-194
 *
 *     A pathological case for LU with gaussian elimination.
 *
 ***** @arg PlasmaMatrixHadamard:
 *     See http://www.mathworks.fr/fr/help/matlab/ref/hadamard.html
 *
 *     Initialize the tile A to create the Hadamard matrix of order gN.
 *
 *     Hadamard matrices are matrices of 1's and -1's whose columns are orthogonal,
 *
 *       H'*H = gN*I
 *
 *       where [gN gN]=size(H) and I = eye(gN,gN) ,.
 *
 *      They have applications in several different areas, including
 *      combinatorics, signal processing, and numerical analysis.
 *
 *      An n-by-n Hadamard matrix with n > 2 exists only if rem(n,4) =
 *      0. This function handles only the cases where n is a power of
 *      2.
 *
 ***** @arg PlasmaMatrixHankel:
 *     See http://en.wikipedia.org/wiki/Hankel_matrix
 *
 *     Hankel matrix
 *
 *     In linear algebra, a Hankel matrix (or catalecticant matrix), named after
 *     Hermann Hankel, is a square matrix with constant skew-diagonals (positive
 *     sloping diagonals), e.g.:
 *
 *     \f[ \begin{bmatrix}
 *     a & b & c & d & e \\
 *     b & c & d & e & f \\
 *     c & d & e & f & g \\
 *     d & e & f & g & h \\
 *     e & f & g & h & i \\
 *     \end{bmatrix}
 *     \f].
 *
 *     A(i,j) = A(i-1,j+1)
 *
 ***** @arg PlasmaMatrixHilb:
 *     See http://www.mathworks.fr/fr/help/matlab/ref/hilb.html
 *
 *     Hilbert Matrix
 *
 *     The Hilbert matrix is a notable example of a poorly conditioned
 *     matrix. The elements of the Hilbert matrices are:
 *       H(i,j) = 1/(i * + j – 1).
 *
 ***** @arg PlasmaMatrixHouse:
 *     See http://www.mathworks.fr/fr/help/matlab/ref/gallery.html#f84-999993
 *
 *     Householder matrix
 *
 *     Generates a random column vector of size M, and returns the housholder matrix
 *     H = eye(n,n) - beta*v*v' that satisfies the relationship
 *
 *     H*x = -sign(x(1))*norm(x)*e1
 *
 *     where e1 is the first column of eye(n,n). Note that if x is complex, then
 *     sign(x) exp(i*arg(x)) (which equals x./abs(x) when x is nonzero).
 *
 ***** @arg PlasmaMatrixInvhess:
 *     See http://www.mathworks.fr/fr/help/matlab/ref/gallery.html#f84-1000000
 *
 *     Inverse of an upper Hessenberg matrix
 *
 *     A = gallery('invhess',x,y), where x is a length n vector and y
 *     is a length n-1 vector, returns the matrix whose lower triangle
 *     agrees with that of ones(n,1)*x' and whose strict upper
 *     triangle agrees with that of [1 y]*ones(1,n).
 *
 *     The matrix is nonsingular if x(1) ~= 0 and x(i+1) ~= y(i) for
 *     all i, and its inverse is an upper Hessenberg matrix. Argument
 *     y defaults to -x(1:n-1).
 *
 *     If x is a scalar, invhess(x) is the same as invhess(1:x).
 *
 *     Here: gallery('invhess', gM)
 *
 ***** @arg PlasmaMatrixKms:
 *     See http://www.mathworks.fr/fr/help/matlab/ref/gallery.html#f84-1000026
 *
 *     Kac-Murdock-Szego Toeplitz matrix
 *
 *     Returns the n-by-n Kac-Murdock-Szego Toeplitz matrix such that
 *     A(i,j) = rho^(abs(i-j)), for real rho.
 *
 *     For complex rho, the same formula holds except that elements
 *     below the diagonal are conjfugated. rho defaults to 0.5.
 *
 *     The KMS matrix A has these properties:
 *         - An LDL' factorization with L inv(gallery('triw',n,-rho,1))',
 *           and D(i,i) (1-abs(rho)^2)*eye(n), except D(1,1) = 1.
 *         - Positive definite if and only if 0 < abs(rho) < 1.
 *         - The inverse inv(A) is tridiagonal.Symmetric Hankel matrix
 *
 *     In this function, rho is set to 0.5 and cannot be changed.
 *
 ***** @arg PlasmaMatrixLangou:
 *     Generates a pathological case for LU with gaussian elimination.
 *
 *     Returns a random matrix on which, the columns from N/4 to N/2
 *     are scaled down by eps.
 *     These matrices fails on LU with partial pivoting, but Hybrid
 *     LU-QR algorithms manage to recover the scaled down columns.
 *
 ***** @arg PlasmaMatrixLehmer:
 *     See http://www.mathworks.fr/fr/help/matlab/ref/gallery.html#f84-1000049
 *
 *     Symmetric positive definite matrix
 *
 *     Returns the symmetric positive definite n-by-n matrix such that
 *     A(i,j) = i/j for j >= i.
 *
 *     The Lehmer matrix A has these properties:
 *         - A is totally nonnegative.
 *         - The inverse inv(A) is tridiagonal and explicitly known.
 *         - The order n <= cond(A) <= 4*n*n.Matrix associated with the
 *           Riemann hypothesis
 *
 ***** @arg PlasmaMatrixLotkin:
 *     See http://www.mathworks.fr/fr/help/matlab/ref/gallery.html#f84-1000062
 *
 *     Lotkin matrix
 *
 *     Returns the Hilbert matrix with its first row altered to all
 *     ones. The Lotkin matrix A is nonsymmetric, ill-conditioned, and
 *     has many negative eigenvalues of small magnitude. Its inverse
 *     has integer entries and is known explicitly.
 *
 ***** @arg PlasmaMatrixMinij:
 *     See http://www.mathworks.fr/fr/help/matlab/ref/gallery.html#f84-1000066
 *
 *     Symmetric positive definite matrix
 *
 *     Returns the n-by-n symmetric positive definite matrix with
 *     A(i,j) = min(i,j).
 *
 *     The minij matrix has these properties:
 *         - The inverse inv(A) is tridiagonal and equal to -1 times the
 *           second difference matrix, except its (n,n) element is 1.
 *         - Givens' matrix, 2*A-ones(size(A)), has tridiagonal inverse
 *           and eigenvalues 0.5*sec((2*r-1)*pi/(4*n))^2, where r=1:n.
 *         - (n+1)*ones(size(A))-A has elements that are max(i,j) and a
 *           tridiagonal inverse.
 *
 ***** @arg PlasmaMatrixMoler:
 *     See http://www.mathworks.fr/fr/help/matlab/ref/gallery.html#f84-1000074
 *
 *     Symmetric positive definite matrix
 *
 *     Returns the symmetric positive definite n-by-n matrix U'*U,
 *     where U = gallery('triw',n,alpha).
 *
 *     For the default alpha = -1, A(i,j) = min(i,j)-2, and A(i,i) =
 *     i. One of the eigenvalues of A is small.
 *
 ***** @arg PlasmaMatrixOrthog:
 *     See http://www.mathworks.fr/fr/help/matlab/ref/gallery.html#f84-1000083
 *
 *     Orthogonal and nearly orthogonal matrices
 *
 *     Returns the matrix Q of order n, such that:
 *        Q(i,j) = sqrt(2/(n+1)) * sin(i*j*pi/(n+1))
 *
 *     Symmetric eigenvector matrix for second difference matrix.
 *
 ***** @arg PlasmaMatrixParter:
 *     See http://www.mathworks.fr/fr/help/matlab/ref/gallery.html#f84-1000116
 *
 *     Toeplitz matrix with singular values near pi.
 *     Returns the tile A, such that the elment of the matrix are 1/(i-j+0.5).
 *
 *     C is a Cauchy matrix and a Toeplitz matrix. Most of the
 *     singular values of C are very close to pi.
 *
 ***** @arg PlasmaMatrixRandom:
 *     Random matrix with values between -0.5 and 0.5.
 *
 ***** @arg PlasmaMatrixRiemann:
 *     See http://www.mathworks.fr/fr/help/matlab/ref/gallery.html#f84-1000232
 *
 *     Matrix associated with the Riemann hypothesis
 *
 *     Returns an n-by-n matrix for which the Riemann hypothesis is
 *     true if and only if for every eps > 0.
 *
 *     The Riemann matrix is defined by:
 *
 *        A = B(2:n+1,2:n+1)
 *
 *        where B(i,j) = i-1 if i divides j, and B(i,j) = -1 otherwise.
 *
 *     The Riemann matrix has these properties:
 *         - Each eigenvalue e(i) satisfies abs(e(i)) <= m-1/m, where m = n+1.
 *         - i <= e(i) <= i+1 with at most m-sqrt(m) exceptions.
 *         - All integers in the interval (m/3, m/2] are eigenvalues.
 *
 ***** @arg PlasmaMatrixRis:
 *     See http://www.mathworks.fr/fr/help/matlab/ref/gallery.html#f84-1000243
 *
 *     Symmetric Hankel matrix
 *     Returns a symmetric gN-by-gN Hankel matrix with elements.
 *
 *         A(i,j) = 0.5/(n-i-j+1.5)
 *
 *     The eigenvalues of A cluster around π/2 and –π/2. This matrix
 *     was invented by F.N. Ris.
 *
 ***** @arg PlasmaMatrixToeppd:
 *     See http://www.mathworks.fr/fr/help/matlab/ref/gallery.html#f84-1000272
 *
 *     A toeppd matrix is an n-by-n symmetric, positive semi-definite (SPD)
 *     Toeplitz matrix composed of the sum of m rank 2 (or, for certain theta,
 *     rank 1) SPD Toeplitz matrices. Specifically,
 *
 *     T = w(1)*T(theta(1)) + ... + w(m)*T(theta(m))
 *
 *     where T(theta(k)) has (i,j) element cos(2*pi*theta(k)*(i-j)).
 *
 *     In this matrix generation: w = rand(m,1), and theta = rand(m,1).
 *
 ***** @arg PlasmaMatrixWilkinson:
 *     See http://www.mathworks.fr/fr/help/matlab/ref/wilkinson.html
 *
 *     Wilkinson's eigenvalue test matrix
 *
 *     Returns one of J. H. Wilkinson's eigenvalue test matrices. It
 *     is a symmetric, tridiagonal matrix with pairs of nearly, but
 *     not exactly, equal eigenvalues.
 *
 ***** @arg PlasmaMatrixWright:
 *     See [3] S. J. Wright, A collection of problems for which
 *             Gaussian elimination with partial pivoting is unstable,
 *             SIAM J. SCI. STATIST. COMPUT., 14 (1993), pp. 231-238.
 *
 *     A pathological case for LU with gaussian elimination.
 *
 * @param[in] M
 *          The number of rows of A.
 *
 * @param[in] N
 *          The order of the matrix A. N >= 0.
 *
 * @param[out] A
 *          On exit, The matrix A generated.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,M).
 *
 * @param[in] seed
 *          The seed used in the random generation.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 *******************************************************************************
 *
 * @sa PLASMA_cpltmg_Tile
 * @sa PLASMA_cpltmg_Tile_Async
 * @sa PLASMA_cpltmg
 * @sa PLASMA_dpltmg
 * @sa PLASMA_spltmg
 *
 ******************************************************************************/
int PLASMA_cpltmg( PLASMA_enum mtxtype, int M, int N,
                    PLASMA_Complex32_t *A, int LDA,
                    unsigned long long int seed )
{
    int NB;
    int status;
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    PLASMA_desc descA;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_cpltmg", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    /* Check input arguments */
    if (M < 0) {
        plasma_error("PLASMA_cpltmg", "illegal value of M");
        return -1;
    }
    if (N < 0) {
        plasma_error("PLASMA_cpltmg", "illegal value of N");
        return -2;
    }
    if (LDA < max(1, M)) {
        plasma_error("PLASMA_cpltmg", "illegal value of LDA");
        return -4;
    }
    /* Quick return */
    if (min(M, N) == 0)
        return PLASMA_SUCCESS;

    /* Tune NB depending on M, N & NRHS; Set NBNB */
    status = plasma_tune(PLASMA_FUNC_CGEMM, M, N, 0);
    if (status != PLASMA_SUCCESS) {
        plasma_error("PLASMA_cpltmg", "plasma_tune() failed");
        return status;
    }

    /* Set NT */
    NB = PLASMA_NB;
    plasma_sequence_create(plasma, &sequence);
    descA = plasma_desc_init(
        PlasmaComplexFloat, NB, NB, NB*NB,
        LDA, N, 0, 0, M, N);
    descA.mat = A;

    /* Call the tile interface */
    PLASMA_cpltmg_Tile_Async( mtxtype, &descA, seed, sequence, &request );

    plasma_ciptile2lap( descA, A, NB, NB, LDA, N,  sequence, &request);
    plasma_dynamic_sync();

    status = sequence->status;
    plasma_sequence_destroy(plasma, sequence);

    return status;
}

/***************************************************************************//**
 *
 * @ingroup PLASMA_Complex32_t_Tile
 *
 *  PLASMA_cpltmg_Tile - Generate a random matrix by tiles.
 *  Tile equivalent of PLASMA_cpltmg().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] mtxtype
 *         See PLASMA_cpltmg()
 *
 * @param[in] A
 *          On exit, The random matrix A generated.
 *
 * @param[in] seed
 *          The seed used in the random generation.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa PLASMA_cpltmg
 * @sa PLASMA_cpltmg_Tile_Async
 * @sa PLASMA_cpltmg_Tile
 * @sa PLASMA_dpltmg_Tile
 * @sa PLASMA_spltmg_Tile
 *
 ******************************************************************************/
int PLASMA_cpltmg_Tile( PLASMA_enum mtxtype, PLASMA_desc *A,
                         unsigned long long int seed )
{
    plasma_context_t *plasma;
    PLASMA_sequence *sequence = NULL;
    PLASMA_request request = PLASMA_REQUEST_INITIALIZER;
    int status;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_cpltmg_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    plasma_sequence_create(plasma, &sequence);
    PLASMA_cpltmg_Tile_Async( mtxtype, A, seed, sequence, &request );
    plasma_dynamic_sync();
    status = sequence->status;
    plasma_sequence_destroy(plasma, sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup PLASMA_Complex32_t_Tile_Async
 *
 *  PLASMA_cpltmg_Tile_Async - Generate a random matrix by tiles.
 *  Non-blocking equivalent of PLASMA_cpltmg_Tile().
 *  May return before the computation is finished.
 *  Allows for pipelining of operations at runtime.
 *
 *******************************************************************************
 *
 * @param[in] mtxtype
 *         See PLASMA_cpltmg()
 *
 * @param[in] A
 *          Descriptor of the matrix A to generate.
 *          On exit, The random matrix A generated.
 *
 * @param[in] sequence
 *          Identifies the sequence of function calls that this call belongs to
 *          (for completion checks and exception handling purposes).
 *
 * @param[out] request
 *          Identifies this function call (for exception handling purposes).
 *
 *******************************************************************************
 *
 * @sa PLASMA_cpltmg
 * @sa PLASMA_cpltmg_Tile
 * @sa PLASMA_cpltmg_Tile_Async
 * @sa PLASMA_dpltmg_Tile_Async
 * @sa PLASMA_spltmg_Tile_Async
 *
 ******************************************************************************/
int PLASMA_cpltmg_Tile_Async( PLASMA_enum mtxtype, PLASMA_desc *A,
                               unsigned long long int seed,
                               PLASMA_sequence *sequence,
                               PLASMA_request  *request)
{
    PLASMA_desc descA;
    plasma_context_t *plasma;

    plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA_cpltmg_Tile", "PLASMA not initialized");
        return PLASMA_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        plasma_fatal_error("PLASMA_cpltmg_Tile", "NULL sequence");
        return PLASMA_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        plasma_fatal_error("PLASMA_cpltmg_Tile", "NULL request");
        return PLASMA_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == PLASMA_SUCCESS)
        request->status = PLASMA_SUCCESS;
    else
        return plasma_request_fail(sequence, request, PLASMA_ERR_SEQUENCE_FLUSHED);

    /* Check descriptors for correctness */
    if (plasma_desc_check(A) != PLASMA_SUCCESS) {
        plasma_error("PLASMA_cpltmg_Tile", "invalid descriptor");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    } else {
        descA = *A;
    }
    /* Check input arguments */
    if (descA.nb != descA.mb) {
        plasma_error("PLASMA_cpltmg_Tile", "only square tiles supported");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }

    /* Quick return */
    if (min( descA.m, descA.n ) == 0)
        return PLASMA_SUCCESS;

    switch( mtxtype ) {
    case PlasmaMatrixCircul:
        plasma_dynamic_call_4(plasma_pcpltmg_circul,
            PLASMA_desc,            descA,
            unsigned long long int, seed,
            PLASMA_sequence*,       sequence,
            PLASMA_request*,        request);
        break;

    case PlasmaMatrixChebvand:
        plasma_dynamic_call_3(plasma_pcpltmg_chebvand,
            PLASMA_desc,            descA,
            PLASMA_sequence*,       sequence,
            PLASMA_request*,        request);
        break;

    case PlasmaMatrixCondex:
        plasma_cpltmg_condex( descA, sequence, request );
        break;

    case PlasmaMatrixFiedler:
        if (descA.m != descA.n) {
            plasma_error("PLASMA_cpltmg_Tile_Async", "Only square fiedler matrices can be generated");
            return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
        }
        plasma_dynamic_call_4(plasma_pcpltmg_fiedler,
            PLASMA_desc,            descA,
            unsigned long long int, seed,
            PLASMA_sequence*,       sequence,
            PLASMA_request*,        request);
        break;

    case PlasmaMatrixHankel:
        plasma_dynamic_call_4(plasma_pcpltmg_hankel,
            PLASMA_desc,            descA,
            unsigned long long int, seed,
            PLASMA_sequence*,       sequence,
            PLASMA_request*,        request);
        break;

    case PlasmaMatrixHouse:
        plasma_cpltmg_house( descA, seed, sequence, request );
        break;

    case PlasmaMatrixToeppd:
        plasma_dynamic_call_4(plasma_pcpltmg_toeppd,
            PLASMA_desc,            descA,
            unsigned long long int, seed,
            PLASMA_sequence*,       sequence,
            PLASMA_request*,        request);
        break;

    case PlasmaMatrixCauchy:
    case PlasmaMatrixCompan:
    case PlasmaMatrixDemmel:
    case PlasmaMatrixDorr:
    case PlasmaMatrixFoster:
    case PlasmaMatrixHadamard:
    case PlasmaMatrixHilb:
    case PlasmaMatrixInvhess:
    case PlasmaMatrixKms:
    case PlasmaMatrixLangou:
    case PlasmaMatrixLehmer:
    case PlasmaMatrixLotkin:
    case PlasmaMatrixMinij:
    case PlasmaMatrixMoler:
    case PlasmaMatrixOrthog:
    case PlasmaMatrixParter:
    case PlasmaMatrixRandom:
    case PlasmaMatrixRiemann:
    case PlasmaMatrixRis:
    case PlasmaMatrixWilkinson:
    case PlasmaMatrixWright:
        plasma_parallel_call_5(plasma_pcpltmg,
            PLASMA_enum,            mtxtype,
            PLASMA_desc,            descA,
            unsigned long long int, seed,
            PLASMA_sequence*,       sequence,
            PLASMA_request*,        request);
        break;

    default:
        plasma_error("PLASMA_cpltmg_Tile", "Illegal value of mtxtype");
        return plasma_request_fail(sequence, request, PLASMA_ERR_ILLEGAL_VALUE);
    }

    return PLASMA_SUCCESS;
}
