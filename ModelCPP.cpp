/**
 * @file ModelCPP.cpp
 *
 * @date Jan 13, 2015
 * @author meyerx
 * @brief
 */
#include "ModelCPP.hpp"

Model_CPP::Model_CPP(const uint aNRComb) {
	notInvertible =  false;
	nrComb = aNRComb;
}

Model_CPP::~Model_CPP() {
}


void Model_CPP::computeBranchConditionalVector(const double t, const double *g, VectorXd &h) const {

	ArrayXd dt;
	VectorXd dt2;
	VectorXd vecG = Map< const VectorXd >(g, nrComb, 1);

	dt = D.array() * t;
	dt2 = dt.exp().matrix();
	h = (V * dt2.asDiagonal() * invV) * vecG;
	/*cout << " H : " << (V * dt2.asDiagonal() * invV);
	cout << endl;*/
}

void Model_CPP::condLikeFunction(double *condLike, const double *gi, const double *gj, const double ti, const double tj) const {

	VectorXd hi, hj;
	computeBranchConditionalVector(ti, gi, hi);
	computeBranchConditionalVector(tj, gj, hj);

	for(uint i=0; i< nrComb; ++i){
		condLike[i] = hi[i] * hj[i];
	}
}


double Model_CPP::executeCond(struct node* n, double *vec) {

    if(notInvertible) {
	    uint assign;
    	for(assign=0;assign<nrComb;assign++){
        vec[assign]= -std::numeric_limits<double>::infinity();
    	}
  	  return -std::numeric_limits<double>::infinity();
    }
    if (n != NULL) {
        if((n->left == NULL) && (n->right == NULL))// is leaf
        {
        	uint assign;
            for (assign=0;assign<nrComb;assign++){
            	vec[assign]=n->probVector[assign];
            }

            return n->dist;
        }else{
            double l_probVector[nrComb];
            double r_probVector[nrComb];
            double l_dist=executeCond(n->left, l_probVector);
            double r_dist=executeCond(n->right, r_probVector);
            condLikeFunction(vec, l_probVector, r_probVector, l_dist, r_dist);

            return 	n->dist;
        }
    }else{
     cout << "DEBUG : errorrr!!!\n" << endl;
     return 0.;
    }
}

void Model_CPP::setQ(double **Qm) {
    MatrixXd Q(nrComb, nrComb);
    uint myI, myJ=0;
    for (myI = 0; myI < nrComb; myI++){
    	for (myJ = 0; myJ < nrComb; myJ++){
    		Q(myI, myJ) = Qm[myI][myJ];
    	}
    }

	EigenSolver<MatrixXd > es(Q, true);
	D = es.eigenvalues().real().eval();
	V = es.eigenvectors().real().eval();
	//invV = V.inverse().eval();

	Eigen::FullPivLU<MatrixXd> lu(V);
	if(!lu.isInvertible()){
		notInvertible = true;
		//cout << "not inversible" << endl;
	} else {
		notInvertible = false;
		invV = lu.inverse();
	}

	/*cout << "Q : " << endl << Q; cout << endl;
	cout << "D : " << endl << D; cout << endl;
	cout << "V : " << endl << V; cout << endl;
	cout << "invV : " << endl << invV; cout << endl;*/
}

void Model_CPP::matInverse(double **Qm,double **QmInverse) {
            MatrixXd QmNew(nrComb, nrComb);
            MatrixXd QmInverseNew(nrComb, nrComb);
            uint myI, myJ=0;
            for (myI = 0; myI < nrComb; myI++){
    	            for (myJ = 0; myJ < nrComb+1; myJ++){
    		QmNew(myI, myJ) = Qm[myI][myJ];
    	            }
            }
	QmInverseNew = QmNew.inverse();
	for (myI = 0; myI < nrComb; myI++){
    	            for (myJ = 0; myJ < nrComb; myJ++){
    		QmInverse[myI][myJ] = QmInverseNew(myI, myJ) ;
    	            }
            }
	
   /*MatrixXf A(3,3);
   A << 1, 2, 1,
        2, 1, 0,
        -1, 1, 2;
   MatrixXf AInv(3,3)=A.inverse() ;
   cout << "Here is the matrix A:\n" << A << endl;
   cout << "The determinant of A is " << A.determinant() << endl;
   cout << "The inverse of A is:\n" << A.inverse() << endl;
   
   
   for (myI = 0; myI < 3; myI++){
    for (myJ = 0; myJ < 3; myJ++){
    cout << " " <<AInv(myI,myJ);
   }
   cout << endl;
   }
   */
   
}

