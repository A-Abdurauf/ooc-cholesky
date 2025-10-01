/**
 * @file core_dlapst.c
 *
 *  PLASMA computational routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.8.0
 * @author Gregoire Pichon
 * @author Azzam Haidar
 * @date 2014-07
 * @precisions normal d -> s
 *
 **/
#include "common.h"
#include <math.h>

/**
 *******************************************************************************
 *
 * @ingroup CORE_double
 *
 *  CORE_dlapst - It is a modified version of the LAPACK routine D inspired from
 *  the ScaLAPACK dlapst.
 *
 *  Define a permutation INDX that sorts the numbers in D
 *  in increasing order (if ID = 'I') or
 *  in decreasing order (if ID = 'D' ).
 *
 *  Use Quick Sort, reverting to Insertion sort on arrays of
 *  size <= 20. Dimension of STACK limits N to about 2**32.
 *
 *******************************************************************************
 *
 * @param[in] type
 *          The type of sorting
 *          'PlasmaIncreasingOrder' for increasing order
 *          'PlasmaDecreasingOrder' for decreasing order
 *
 * @param[in] n
 *          n specifies the dimension of the array D
 *
 * @param[in] D
 *          D is the array to be sorted of dimension n.
 *
 * @param[out] INDX
 *          On exit, contains the permutation which sorts the array D. INDX is
 *          of dimension n.
 *
 *******************************************************************************/
#define SELECT 20

int
CORE_dlapst( PLASMA_enum type, int n, const double *D, int *INDX )
{
    int i, j;
    int stack[32][2];
    int stkpnt;

    int start, endd;
    int itmp;

    double d1, d2, d3;
    double dmnmx;

    if ( (type != PlasmaIncreasingOrder) &&
         (type != PlasmaDecreasingOrder) ){
        coreblas_error( 1, "Sorting Type unknown\n" );
        return -1;
    }

    /* Initialize the permutation array */
    for (i=0; i<n; i++){
        INDX[i] = i;
    }

    stkpnt = 0;
    stack[0][0] = 0;
    stack[0][1] = n-1;

    while ( stkpnt >= 0 ){
        start = stack[stkpnt][0];
        endd  = stack[stkpnt][1];
        stkpnt--;

        /* Set size smaller than SELECT */
        if ( (endd-start <= SELECT)  &&
             (endd-start >  0 ) )
        {
            /*
             * Do insertion sort on D( start:endd )
             */
            if (type == PlasmaDecreasingOrder) {
                /* Sort into decreasing order */
                for (i=start; i<=endd; i++){
                    for (j=i; j>start; j--){
                        if ( D[INDX[j]] > D[INDX[j-1]] ){
                            itmp      = INDX[j];
                            INDX[j]   = INDX[j-1];
                            INDX[j-1] = itmp;
                        }
                        else{
                            break;
                        }
                    }
                }
            }
            else {
                /* Sort into increasing order */
                for (i=start; i<=endd; i++){
                    for (j=i; j>start; j--){
                        if ( D[INDX[j]] < D[INDX[j-1]] ){
                            itmp      = INDX[j];
                            INDX[j]   = INDX[j-1];
                            INDX[j-1] = itmp;
                        }
                        else{
                            break;
                        }
                    }
                }
            }
        }
        /* Set size larger than SELEC */
        else if ( endd-start > SELECT ){
            /*
             * Partition D( START:ENDD ) and stack parts, largest one first
             */

            /* Choose partition entry as median of 3 */
            d1 = D[INDX[start]];
            d2 = D[INDX[endd ]];
            i = (start + endd) / 2;
            d3 = D[INDX[i]];
            if ( d1 < d2 ){
                if ( d3 < d1 ){
                    dmnmx = d1;
                }
                else if ( d3 < d2 ){
                    dmnmx = d3;
                }
                else{
                    dmnmx = d2;
                }
            }
            else{
                if ( d3 < d2 ){
                    dmnmx = d2;
                }
                else if ( d3 < d1 ){
                    dmnmx = d3;
                }
                else{
                    dmnmx = d1;
                }
            }

            if (type == PlasmaDecreasingOrder) {
                /* Sort into decreasing order */
                i = start-1;
                j = endd+1;

              start_split_dec:
                j--;
                while ( D[INDX[j]] < dmnmx ){
                    j--;
                }

                i++;
                while ( D[INDX[i]] > dmnmx ){
                    i++;
                }

                if ( i < j ) {
                    itmp = INDX[i];
                    INDX[i] = INDX[j];
                    INDX[j] = itmp;
                }
                else {
                    goto start_split_dec;
                }

                if ( j-start > endd-j-1 ){
                    stkpnt++;
                    stack[stkpnt][0] = start;
                    stack[stkpnt][1] = j;
                    stkpnt++;
                    stack[stkpnt][0] = j+1;
                    stack[stkpnt][1] = endd;
                }
                else{
                    stkpnt++;
                    stack[stkpnt][0] = j+1;
                    stack[stkpnt][1] = endd;
                    stkpnt++;
                    stack[stkpnt][0] = start;
                    stack[stkpnt][1] = j;
                }
            }
            else {
                /* Sort into increasing order */

                i = start-1;
                j = endd+1;

              start_split_inc:
                j--;
                while ( D[INDX[j]] > dmnmx ){
                    j--;
                }

                i++;
                while ( D[INDX[i]] < dmnmx ){
                    i++;
                }

                if ( i < j ){
                    itmp = INDX[i];
                    INDX[i] = INDX[j];
                    INDX[j] = itmp;
                    goto start_split_inc;
                }

                if ( j-start > endd-j-1 ){
                    stkpnt++;
                    stack[stkpnt][0] = start;
                    stack[stkpnt][1] = j;
                    stkpnt++;
                    stack[stkpnt][0] = j+1;
                    stack[stkpnt][1] = endd;
                }
                else{
                    stkpnt++;
                    stack[stkpnt][0] = j+1;
                    stack[stkpnt][1] = endd;
                    stkpnt++;
                    stack[stkpnt][0] = start;
                    stack[stkpnt][1] = j;
                }
            }
        } /* end solving one subproblem of size SELECT */
    } /* end sorting subproblems of size SELECT */

    return 0;
}
