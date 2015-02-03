/**
 * @file wrapperCPP.h
 *
 * @date Jan 12, 2015
 * @author meyerx
 * @brief
 */
#ifndef WRAPPERCPP_H_
#define WRAPPERCPP_H_

#include "tree.h"

#ifdef __cplusplus
extern "C" {
#endif



/* Model */
typedef struct model_cpp model_cpp;

model_cpp* new_Model_CPP(const uint aNRComb);
void delete_Model_CPP(model_cpp* model);
double Model_CPP_executeCond(model_cpp* model, struct node* n, double *vec);
void Model_CPP_setQ(model_cpp* model, double **Qm);
void Model_CPP_setQNull(model_cpp* model, double **Qm);
void Model_CPP_matInverse(model_cpp* model,double **Qm,double **QmInverse);
#ifdef __cplusplus
}
#endif



#endif /* WRAPPERCPP_H_ */
