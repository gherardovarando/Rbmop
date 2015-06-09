/************************************************************************
*                                                                       *
*       AUTHOR : GHERARDO VARANDO, LUIS RODRIGUEZ-LUJAN, DONALD DUCK    *
*       AIM : Header of deboor_eval.c  					            *
*       DATA : 30 APRIL 2015.                                           *
*                                                                       *
*************************************************************************/  

#include <R.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "deboor_eval.h"

#define MAX_DBOR(a,b) (a)>(b)?(a):(b)

/**  
 * Evaluate MoPs using Deboor's algorithm
 *
 * @param [in] t 
 * @param [in] k
 * @param [in] knots_len Knots vector length
 * @param [in] knots
 * @param [in] ctr_len ctr vector length
 * @param [in] ctr
 * @param [in] min
 * @param [out] retval Evaluation value
 */
void deboor_eval(double* t, int* k, int* knots_len, double* knots, int* ctr_len, double* ctr, double* min, double *retval){

    // Iterators 
    int i,j;

    // Index and multiplicity
    int idx,multiplicity;

    // Ctrl aux vector
    double alpha;
    double *aux ;


    if( *t<knots[0] || *t>knots[ (*knots_len)-1 ]) { *retval = *min; return; }
    else if (*t == knots[(*knots_len)-1]) { *retval = MAX_DBOR((*min),ctr[(*ctr_len)-1]); return ;}
    else if (*t == knots[0]) { *retval = MAX_DBOR((*min),ctr[0]); return; }
    else{
        locate_and_multiplicity(*knots_len,knots,*t,&idx,&multiplicity);
        if(*k == 1 ){ *retval = ctr[idx] ; return ;}
        else{
            aux = malloc(sizeof(double)* (*ctr_len));
            memcpy(aux,ctr,sizeof(double)*(*ctr_len));
            for(j=1;j<( (*k) - multiplicity+1); j++){
                for(i=idx-*k+j+1;i<(idx-multiplicity+1);i++){
                    alpha = (*t-knots[i])/(knots[i+*k-j] - knots[i]);
                    aux[i] = (1-alpha)*ctr[i-1]+alpha*ctr[i];
                }
                memcpy(ctr,aux,sizeof(double)*(*ctr_len));
            }
            free(aux);
            *retval = MAX_DBOR(ctr[idx-multiplicity],(*min));
            return ;
        }
    }
}

/**
 * Returns the position of the given element and its multiplicity in the vecotr
 *
 * @param [in] len ordered vector length
 * @param [in] orderedVector
 * @param [in] value value to locate
 * @param [out] pos last position of value in orderedVector
 * @param [out] multiplicity number of occurences of value in ordered vector
 */
inline void locate_and_multiplicity(int len, double* orderedVector, double value, int* pos, int* multiplicity){
    for(*multiplicity=0,*pos=0;*pos<len && orderedVector[*pos]<=value; (*pos)++){
        if(orderedVector[*pos]==value) (*multiplicity)++;
    }
    // Correct to point to the last position
    (*pos)--;
}
