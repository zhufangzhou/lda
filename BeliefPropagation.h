#ifndef __BELIEFPROPAGATION
#define __BELIEFPROPAGATION

#include "util.h"
#include "LdaBase.h"

/* synchronous Belief Propagation */
class sBP : public LdaBase {
	private:
		
		/* functions */
		void init();
	public:
		sBP(string path, int k, int t, double alpha, double beta);
		~sBP();
		void LearnTopics();
		double* getPhi();
		double* getTheta();
};

/* asynchronous Belief Propagation */
class aBP : public LdaBase {
	private:

		/* functions */
		void init();
	public:
		aBP(string path, int k, int t, double alpha, double beta);
		~aBP();
		void LearnTopics();
		double* getPhi();
		double* getTheta();
};

/* residual Belief Propagation ---- residual accumulate by documents */
class RBP_doc : public LdaBase {
	private:
		double *residual;
		double *mu_new;
		int *seq;

		/* functions */
		void init();
	public:
		RBP_doc(string path, int k, int t, double alpha, double beta);
		~RBP_doc();
		void LearnTopics();
		double* getPhi();
		double* getTheta();
};
#endif
