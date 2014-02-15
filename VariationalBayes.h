#ifndef __VARIATIONALBAYES
#define __VARIATIONALBAYES

#include "LdaBase.h"
#include "util.h"

class VariationalBayes : public LdaBase {
	private:
		// variables
		double *p_log_phi;								// log(phi) in original LDA model
		double *p_theta;								// theta in original LDA model, get from gamma
		double *var_gamma;								// gamma in variational distribution
		double *var_old_gamma;							// gamma in variational distribution(last iteration)
		double *var_digamma;							// digamma	
		double *var_phi;								// phi in variational distribution

		int *doc_len;									// tokens in each document
		
		int DOC_MAX_LEN;

		double EM_CONVERGED;
		double VAR_CONVERGED;

		// functions 
		void init();	
		void mle();
		double variational_inference(int d);
	public:
		VariationalBayes(string path, int k, int t, double alpha, double beta, double em_converged, double var_converged);
		~VariationalBayes();
		void LearnTopics();
		double* getPhi();
		double* getTheta();
};

#endif