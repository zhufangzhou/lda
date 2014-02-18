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
};

/* residual Belief Propagation ---- residual accumulated by documents */
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
};

/* residual Belief Propagation ---- residual accumulated by vocabulary */
class RBP_voc : public LdaBase {
	private:
		double *residual;
		double *mu_new;
		int *seq;

		/* functions */
		void init();
	public:
		RBP_voc(string path, int k, int t, double alpha, double beta);
		~RBP_voc();
		void LearnTopics();
};
#endif
