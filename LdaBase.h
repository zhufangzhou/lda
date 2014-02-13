#ifndef __LDABASE
#define __LDABASE

#include "util.h"
class LdaBase {
	protected:
		// Variables
		int W;									// vocabulary size
		int D;									// number of documents
		int NNZ;								// number of not zero
		int K;									// number of topics
		int T;									// number of iterations
		double ALPHA;							// hyper-parameter
		double KALPHA;							// ALPHA*K
		double BETA;							// hyper-parameter
		double WBETA;							// BETA*W
		/* input sparse matrix WD(CSC compressed) ------ first dimension is document, second dimension is word */
		int *ir;								
		int *jc;
		double *pr;

		double *mu;
		// sufficient statistics
		double *theta, *thetatot;
		double *phi, *phitot;

		string path;

		// Functions
		void ReadData(string path);
	public:
		LdaBase();
		LdaBase(string path, int K, int T, double alpha, double beta);
		~LdaBase();
		virtual void LearnTopics() = 0;
};

#endif
