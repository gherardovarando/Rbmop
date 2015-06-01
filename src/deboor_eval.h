
/************************************************************************
*                                                                       *
*       AUTHOR : GHERARDO VARANDO, LUIS RODRIGUEZ-LUJAN, DONALD DUCK    *
*       AIM : Header of deboor_eval.c						            *
*       DATA : 30 APRIL 2015.                                           *
*                                                                       *
*************************************************************************/  

#ifndef _DEBOOR_EVAL_H_
#define _DEBOOR_EVAL_H_
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
void deboor_eval(double* t, int* k, int* knots_len, double* knots, int* ctr_len, double* ctr, double* min, double *retval);

/**
 * Returns the position of the given element and its multiplicity in the vecotr
 *
 * @param [in] len ordered vector length
 * @param [in] orderedVector
 * @param [in] value value to locate
 * @param [out] pos last position of value in orderedVector
 * @param [out] multiplicity number of occurences of value in ordered vector
 */
inline void locate_and_multiplicity(int len, double *orderedVector, double value, int* pos, int* multiplicity);
#endif
