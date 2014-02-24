#ifndef __GIBBSSAMPLING
#define __GIBBSSAMPLING

#include "LdaBase.h"
#include "util.h"

class GibbsSampling: public LdaBase {
	private:
		int BURN_IN;								// burn-in

		// after read in the dataset, restore to this format-----z is topic of each token, wd is word index of each token, doc is document index of each token 
		int *z;
		int *wd;
		int *doc;

		double *p;
		void init();								// initialize Markov chain
		int sampleTopic(int token);					// sample a topic with other fixed
	public:
		GibbsSampling(string path, int k, int t, int burn_in, double alpha, double beta);
		~GibbsSampling();
		void LearnTopics();	
};

/* Fast GibbsSampling ----> 'Fast Collapsed Gibbs Sampling For Latent Dirichlet Allocation'*/
class fastGibbsSampling : public LdaBase {
	private:
		int BURN_IN;
		int *z;
		int *wd;
		int *doc;

		int topic_new;

		double *sum_p;
		double *phi_norm;
		double *theta_norm;
		double *p;
		int *theta_idx;
		int *theta_ridx;
		int *d_last;
		int *w_last;
		double min_phitot;

		/* functions */
		void init();
		int sampleTopic(int token, int iter);	
		void update_sort(int n, double *value, int *idx, int *ridx, int ip, int im);
	public:
		fastGibbsSampling(string path, int k, int t, int burn_in, double alpha, double beta);
		~fastGibbsSampling();
		void LearnTopics();
};

#endif
