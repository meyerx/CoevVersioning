/*
 ** ml.c: Apply the maximum likelihood method.
 **
 ** Linda Dib & Wim Hordijk   Last modified: 11 december 2014
 */

#include "def.h"
#include "model.h"
#include "wrapperCPP.h"

//%%%%%%%%%%%%%%%%%%%%%%%%%%% ML FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%

double logLikelihood(int *VectCoevComb, double* AS, double* AD, double* AW1,
		double* AW2) {
	double *a, **Q;
	int loop_j = 0, i;
	double ALikelihood = 0.0;

	/*
	 ** Allocate memory.
	 */
	a = (double *) malloc(nrComb * sizeof(double));
	Q = (double **) malloc(nrComb * sizeof(double *));
	//Qtransposed = (double **)malloc (nrComb * sizeof (double *));
	for (i = 0; i < nrComb; i++) {
		Q[i] = (double *) malloc(nrComb * sizeof(double));
		//Qtransposed[i] = (double *)malloc (nrComb * sizeof (double));
	}

	setQ(Q, *AS, *AD, *AW1, *AW2, VectCoevComb);

	model_cpp *modelCPP = new_Model_CPP(nrComb);
	Model_CPP_setQ(modelCPP, Q);
	Model_CPP_executeCond(modelCPP, &root, a);
	ALikelihood = 0.0;
	for (loop_j = 0; loop_j < nrComb; loop_j++) {
		ALikelihood = ALikelihood + a[loop_j];
	}
	delete_Model_CPP(modelCPP);

	//transposeMatrix(Q,Qtransposed);
	//double *MatrixA=(double *)malloc(nrComb*nrComb*sizeof(double));//INPUT
	//double *MatrixB=(double *)malloc(nrComb*nrComb*sizeof(double));//INPUT
	//double *MatrixC=(double *)malloc(nrComb*nrComb*sizeof(double));//OUTPUT
	//for (loop_j=0;loop_j<nrComb;loop_j++){ a[loop_j]=0;}
	//executeCond(Q, Qtransposed, &root, a, MatrixA, MatrixB, MatrixC);
	//ALikelihood=0.0;
	//for (loop_j=0;loop_j<nrComb;loop_j++){
	//  ALikelihood=ALikelihood+a[loop_j];
	// }

	//free(MatrixA);
	//free(MatrixB);
	//free(MatrixC);
	free(a);
	for (i = 0; i < nrComb; i++) {
		free(Q[i]);
		//free (Qtransposed[i]);
	}
	free(Q);
	//free (Qtransposed);

	double ALogLikelihood = log(ALikelihood);
	if (ALogLikelihood != ALogLikelihood) {
		ALogLikelihood = -INFINITY;
	}
	return ALogLikelihood;
}

double logLikelihoodNull(double* AW1, double* AW2) {
	double *a, **Q;
	int loop_j = 0, i;
	double ALikelihood = 0.0;

	/*
	 ** Allocate memory.
	 */
	a = (double *) malloc(nrComb * sizeof(double));
	Q = (double **) malloc(nrComb * sizeof(double *));
	//Qtransposed = (double **)malloc (nrComb * sizeof (double *));
	for (i = 0; i < nrComb; i++) {
		Q[i] = (double *) malloc(nrComb * sizeof(double));
		//Qtransposed[i] = (double *)malloc (nrComb * sizeof (double));
	}

	setQNull(Q, *AW1, *AW2);

	model_cpp *modelCPP = new_Model_CPP(nrComb);
	Model_CPP_setQ(modelCPP, Q);
	Model_CPP_executeCond(modelCPP, &root, a);
	ALikelihood = 0.0;
	//printf("A :\t");
	for (loop_j = 0; loop_j < nrComb; loop_j++) {
		// if error during eig or inv, set nan to -inf
		ALikelihood = ALikelihood + a[loop_j];
		//printf("a[%d] = %f\t", loop_j, log(a[loop_j]));
	}
	//printf("\n");
	delete_Model_CPP(modelCPP);

	/* transposeMatrix(Q,Qtransposed);
	 double *MatrixA=(double *)malloc(nrComb*nrComb*sizeof(double));//INPUT
	 double *MatrixB=(double *)malloc(nrComb*nrComb*sizeof(double));//INPUT
	 double *MatrixC=(double *)malloc(nrComb*nrComb*sizeof(double));//OUTPUT
	 for (loop_j=0;loop_j<nrComb;loop_j++){ a[loop_j]=0;}
	 executeCond(Q, Qtransposed, &root, a, MatrixA, MatrixB, MatrixC);
	 ALikelihood=0.0;
	 for (loop_j=0;loop_j<nrComb;loop_j++){
	 ALikelihood=ALikelihood+a[loop_j];
	 }

	 free(MatrixA);
	 free(MatrixB);
	 free(MatrixC);*/
	free(a);
	for (i = 0; i < nrComb; i++) {
		free(Q[i]);
		//free (Qtransposed[i]);
	}
	free(Q);
	//free (Qtransposed);

	double ALogLikelihood = log(ALikelihood);
	if (ALogLikelihood != ALogLikelihood) {
		ALogLikelihood = -INFINITY;
	}
	return ALogLikelihood;
}

double optimize(unsigned n, const double* x, double *grad, void *my_func_data) {
	int * func_data = (int *) my_func_data;
	double *AS = (double *) malloc(sizeof(double));
	double *AD = (double *) malloc(sizeof(double));
	double *AW1 = (double *) malloc(sizeof(double));
	double *AW2 = (double *) malloc(sizeof(double));
	// printf("optimize %e, %e, %e\n", x[0],x[1],x[2]);
	*AS = x[0]; //s
	*AD = x[1]; //d
	*AW1 = x[2]; //w
	*AW2 = x[3]; //w
	int loop_i = 0;
	int *AVectorComb = (int *) malloc(nrComb * sizeof(int));
	for (loop_i = 0; loop_i < nrComb; loop_i++) {
		AVectorComb[loop_i] = func_data[loop_i];
	}
	//call function: Likelihood
	double ALogLikelihood = logLikelihood(AVectorComb, AS, AD, AW1, AW2);
	free(AVectorComb);
	free(AS);
	free(AD);
	free(AW1);
	free(AW2);
	return ALogLikelihood;

}

double optimize_null(unsigned n, const double* x, double *grad,
		void *my_func_data) {
	double *AW1 = (double *) malloc(sizeof(double));
	double *AW2 = (double *) malloc(sizeof(double));

	*AW1 = x[0];
	*AW2 = x[1];
	double ALogLikelihood = logLikelihoodNull(AW1, AW2);
	//printf("AW1 = %f - AW2 = %f - LogLik = %f\n", *AW1, *AW2, ALogLikelihood);
	//getchar();
	free(AW1);
	free(AW2);

	return ALogLikelihood;

}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MAIN FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%
int ml() {
	int status, loop_j;
	time_t start, end;
	volatile long unsigned counter;
	status = 0;

	FILE *fp;
	fp = fopen(outFile, "w+");

	int AComb = 0;
	int *AVectorComb = (int *) malloc(nrComb * sizeof(int));
	printf("nrProfiles: %d\n", nrProfiles);
	while (AComb < nrProfiles) {

		start = time(NULL);
		for (loop_j = 0; loop_j < nrComb; loop_j++) {
			AVectorComb[loop_j] = profiles[AComb][loop_j];
			fprintf(fp, "%d", profiles[AComb][loop_j]);
			printf("%d", profiles[AComb][loop_j]);
		}
		fprintf(fp, "\t");
		printf("\n");
		//OPT NULL ////////////////////////////////////////////////////

		if (opt == 1) {
			double lb_null[2] = { 1e-9, 01e-9 }; // lower bounds
			double ub_null[2] = { 100, 100 }; // upper bounds
			nlopt_opt opt_null;

			//opt_null = nlopt_create(NLOPT_LD_LBFGS, 2); // algorithm and dimensionality
			//opt_null = nlopt_create(NLOPT_LN_COBYLA, 2);
			//opt_null = nlopt_create(NLOPT_LN_BOBYQA, 2);
			opt_null = nlopt_create(NLOPT_LN_NELDERMEAD, 2);
			nlopt_set_lower_bounds(opt_null, lb_null);
			nlopt_set_upper_bounds(opt_null, ub_null);

			nlopt_set_max_objective(opt_null, optimize_null, NULL);
			nlopt_set_xtol_rel(opt_null, 1e-10);

			double x_null[2] = { r1, r2 }; // some initial guess
			double max_nullf; // the maximum objective value, upon return
			nlopt_set_maxtime(opt_null, 10.0);
			if (nlopt_optimize(opt_null, x_null, &max_nullf) < 0) {
				fprintf(fp, "nlopt_failed\tX\tX\t");

			} else {
				fprintf(fp, "%e\t%e\t%e\t", x_null[0], x_null[1], max_nullf);
				//printf("%e\t%e\t%e\t", x_null[0],x_null[1], max_nullf);
			}
			nlopt_destroy(opt_null);

			//OPT COEV  ////////////////////////////////////////////////
			double lb[4] = { 1e-9, 1e-9, 1e-9, 1e-9 }; // lower bounds
			double ub[4] = { 100, 100, 100, 100 }; // upper bounds
			nlopt_opt opt;

			//opt = nlopt_create(NLOPT_LD_LBFGS, 4); // algorithm and dimensionality
			//opt = nlopt_create(NLOPT_LN_COBYLA, 4);
			//opt = nlopt_create(NLOPT_LN_BOBYQA, 4);
			opt = nlopt_create(NLOPT_LN_NELDERMEAD, 4);
			nlopt_set_lower_bounds(opt, lb);
			nlopt_set_upper_bounds(opt, ub);

			nlopt_set_max_objective(opt, optimize, AVectorComb);
			nlopt_set_xtol_rel(opt, 1e-10);
			nlopt_set_maxtime(opt, 10.0);
			double x[4] = { s, d, r1, r2 }; // some initial guess
			double maxf; // the maximum objective value, upon return

			if (nlopt_optimize(opt, x, &maxf) < 0) {
				fprintf(fp, "nlopt_failed\tX\tX\tX\tX\t");
			} else {
				fprintf(fp, "%e\t%e\t%e\t%e\t%e\t\n", x[0], x[1], x[2], x[3],
						maxf);
				//printf("%e\t%e\t%e\t%e\t%e\n", x[0], x[1], x[2], x[3], maxf);
			}
			nlopt_destroy(opt);
			end = time(NULL);
			//printf("The profile used %f seconds.\n", difftime(end, start));

		} else {

			double *AS = (double *) malloc(sizeof(double));
			double *AD = (double *) malloc(sizeof(double));
			double *AR1 = (double *) malloc(sizeof(double));
			double *AR2 = (double *) malloc(sizeof(double));

			*AS = s;
			*AD = d;
			*AR1 = r1;
			*AR2 = r2;

			//call function: Likelihood
			double ALogLikelihood = 0.0, ALogLikelihoodNull = 0.0;
			ALogLikelihoodNull = logLikelihoodNull(AR1, AR2);
			ALogLikelihood = logLikelihood(AVectorComb, AS, AD, AR1, AR2);

			printf("The ALogLikelihoodNull: %g.\nThe ALogLikelihood: %g.\n\n",
					ALogLikelihoodNull, ALogLikelihood);

			free(AS);
			free(AD);
			free(AR1);
			free(AR2);

		}
		AComb++;
	}
	fprintf(fp, "\n");
	fclose(fp);
	free(AVectorComb);

	End_of_Routine:
	/*
	 ** Return the status.
	 */
	return (status);
}
