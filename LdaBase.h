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

		int tokens;								// number of tokens in the dataset

		double *mu;
		// sufficient statistics
		double *theta, *thetatot;
		double *phi, *phitot;

		bool ParameterSet;						// default false, set true after call learnTopics()

		string path;

		// Functions
		void ReadData(string path, bool document_major);			// read dataset which is document major
	public:
		LdaBase(string path, int k, int t, double alpha, double beta, bool document_major);
		~LdaBase();
		virtual void LearnTopics() = 0;
		virtual double* getTheta() = 0;
		virtual double* getPhi() = 0;
};

#endif
