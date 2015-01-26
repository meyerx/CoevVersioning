/**
 * @file ModelCPP.h
 *
 * @date Jan 13, 2015
 * @author meyerx
 * @brief
 */
#ifndef MODELCPP_H_
#define MODELCPP_H_

#include "tree.h"

#include <iostream>
#include <Core>
#include <Dense>
#include <Eigenvalues>

#include "wrapperCPP.h"

using namespace Eigen;
using namespace std;

struct model_cpp {};

class Model_CPP : public model_cpp {
public:
	Model_CPP(const uint aNrComb);
	~Model_CPP();

	double executeCond(struct node* n, double *vec);
	void setQ(double **Qm);

private:
	unsigned int nrComb;

	VectorXd D;
	MatrixXd V, invV;

	void condLikeFunction(double *condLike, const double *gi, const double *gj, const double ti, const double tj) const;
	void computeBranchConditionalVector(const double t, const double *g, VectorXd &h) const;

};

inline Model_CPP* real(model_cpp* model) { return static_cast<Model_CPP*>(model); }

model_cpp* new_Model_CPP(const uint aNRComb) { return new Model_CPP(aNRComb); }
void delete_Model_CPP(model_cpp* model) { delete real(model); }
double Model_CPP_executeCond(model_cpp* model, struct node* n, double *vec) {
	return real(model)->executeCond(n, vec);
}
void Model_CPP_setQ(model_cpp* model, double **Qm) {
	real(model)->setQ(Qm);
}


#endif /* MODELCPP_H_ */
