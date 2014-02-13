#ifndef __GIBBSSAMPLING
#define __GIBBSSAMPLING

#include "LdaBase.h"
#include "util.h"

class GibbsSampling: public LdaBase {
	private:
		int BURN_IN;								// burn-in
		int SAMPLE_LAG;								// leave an internal of (SAMPLE_LAG) iteration between subsequent read-outs to obtain decorrelated states of the Markov chain

		// after read in the dataset, restore to this format-----z is topic of each token, wd is word index of each token, doc is document index of each token 
		int *z;
		int *wd;
		int *doc;

		int tokens;									// number of tokens in the dataset
		double *p;
		void init();								// initialize Markov chain
		int sampleTopic(int token);					// sample a topic with other fixed
	public:
		GibbsSampling();
		GibbsSampling(string path, int k, int t, int burn_in, int sample_lag, double alpha, double beta);
		~GibbsSampling();
		void LearnTopics();	
		double* getPhi();
		double* getTheta();
};


#endif
