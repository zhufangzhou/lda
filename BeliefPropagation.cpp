#include "util.h"
#include "BeliefPropagation.h"

/***************************************************************************
 * sBP implemetation 
 **************************************************************************/
sBP::sBP(string path, int k, int t, double alpha, double beta):LdaBase(path, k, t, alpha, beta, true){
	
}

sBP::~sBP() {

}

void sBP::init() {
	int topic, w;
	double x;

	/* allocate space to sufficient statistics */
	phi = new double[W*K];
	memset(phi, 0, sizeof(double)*W*K);
	phitot = new double[K];
	memset(phitot, 0, sizeof(double)*K);
	theta = new double[D*K];
	memset(theta, 0, sizeof(double)*D*K);
	thetatot = new double[D];
	memset(thetatot, 0, sizeof(double)*D);
	mu = new double[NNZ*K];
	memset(mu, 0, sizeof(double)*NNZ*K);

	srand(time(0));

	/* initialize suffcient statistics */
	for(int d = 0; d < D; d++) {
		for(int i = jc[d]; i < jc[d+1]; i++) {
			w = ir[i];
			x = pr[i];
			
			/* randomly pick a topic */
			topic = rand() % K;
			mu[i*K+topic] = 1.0;				// assign this word the topic
			phi[w*K+topic] += x;
			phitot[topic] += x;
			theta[d*K+topic] += x;
			thetatot[d] += x;
		}
	}

}

void sBP::LearnTopics() {
	int iter, w;
	double x, mu_tot, perplexity;
	myTimer *tm = new myTimer();
	
	tm->start();

	printf("sBP training begin!\n");
	init();

	for(int iter = 1; iter <= T; iter++) {

		/* passing message */
		for(int d = 0; d < D; d++) {
			for(int i = jc[d]; i < jc[d+1]; i++) {
				w = ir[i];
				x = pr[i];
				mu_tot = 0.0;
				/* update mu through phi and theta except current word */
				for(int k = 0; k < K; k++) {
					mu[i*K+k] = (phi[w*K+k] - x*mu[i*K+k] +BETA) / (phitot[k] - x*mu[i*K+k] + WBETA) * (theta[d*K+k] - x*mu[i*K+k] + ALPHA);
					mu_tot += mu[i*K+k];
				}
				for(int k = 0; k < K; k++) mu[i*K+k] /= mu_tot;
			}
		}

		/* update theta and phi after passing message */
		memset(phi, 0, sizeof(double)*W*K);
		memset(phitot, 0, sizeof(double)*K);
		memset(theta, 0, sizeof(double)*D*K);
		for(int d = 0; d < D; d++) {
			for(int i = jc[d]; i < jc[d+1]; i++) {
				w = ir[i];
				x = pr[i];
				for(int k = 0; k < K; k++) {
					phi[w*K+k] += x*mu[i*K+k];
					phitot[k] += x*mu[i*K+k];
					theta[d*K+k] += x*mu[i*K+k];
				}
			}
		}


		/* compute perplexity */
		if(iter % 10 == 0) {
			perplexity = 0.0;
			for(int d = 0; d < D; d++) {
				for(int i = jc[d]; i < jc[d+1]; i++) {
					w = ir[i];
					x = pr[i];
					mu_tot = 0.0;
					for(int k = 0; k < K; k++) {
						mu_tot += (phi[w*K+k] + BETA) / (phitot[k] + WBETA) * (theta[d*K+k] + ALPHA) / (thetatot[d] + KALPHA);
					}
					perplexity -= (x * log(mu_tot));
				}
			}
			printf("\tIteration %d, Perplexity is: %.5lf\n", iter, exp(perplexity / tokens));
		}
	}

	tm->end();
	printf("\nTotal training time: %lf seconds.\n", tm->getTime());

	ParameterSet = true;
}


/***************************************************************************
 * aBP implementation
 **************************************************************************/
aBP::aBP(string path, int k, int t, double alpha, double beta):LdaBase(path, k, t, alpha, beta, true) {

}

aBP::~aBP() {

}

void aBP::init() {
	int topic, w;
	double x;

	/* allocate space to sufficient statistics */
	phi = new double[W*K];
	memset(phi, 0, sizeof(double)*W*K);
	phitot = new double[K];
	memset(phitot, 0, sizeof(double)*K);
	theta = new double[D*K];
	memset(theta, 0, sizeof(double)*D*K);
	thetatot = new double[D];
	memset(thetatot, 0, sizeof(double)*D);
	mu = new double[NNZ*K];
	memset(mu, 0, sizeof(double)*NNZ*K);

	srand(time(0));

	/* initialize suffcient statistics */
	for(int d = 0; d < D; d++) {
		for(int i = jc[d]; i < jc[d+1]; i++) {
			w = ir[i];
			x = pr[i];
			
			/* randomly pick a topic */
			topic = rand() % K;
			mu[i*K+topic] = 1.0;				// assign this word the topic
			phi[w*K+topic] += x;
			phitot[topic] += x;
			theta[d*K+topic] += x;
			thetatot[d] += x;
		}
	}
	
}

void aBP::LearnTopics() {
	int w;
	double x, mu_tot = 0.0, perplexity = 0.0;
	myTimer *tm = new myTimer();

	tm->start();

	printf("aBP training begin!\n");
	init();
	
	for(int iter = 1; iter <= T; iter++) {
		/* passing message */
		for(int d = 0; d < D; d++) {
			for(int i = jc[d]; i < jc[d+1]; i++) {
				w = ir[i];
				x = pr[i];
				mu_tot = 0.0;
				for(int k = 0; k < K; k++) {
					/* update mu through phi and theta except current word */
					phi[w*K+k] -= x*mu[i*K+k];
					phitot[k] -= x*mu[i*K+k];
					theta[d*K+k] -= x*mu[i*K+k];
					mu[i*K+k] = (phi[w*K+k] + BETA) / (phitot[k] + WBETA) * (theta[d*K+k] + ALPHA);
					mu_tot += mu[i*K+k];
				}
				/* update phi and theta while passing message */
				for(int k = 0; k < K; k++) {
					mu[i*K+k] /= mu_tot;
					phi[w*K+k] += x*mu[i*K+k];
					phitot[k] += x*mu[i*K+k];
					theta[d*K+k] += x*mu[i*K+k];
				}
			}
		}

		/* compute perplexity */
		if(iter % 10 == 0) {
			perplexity = 0.0;
			for(int d = 0; d < D; d++) {
				for(int i = jc[d]; i < jc[d+1]; i++) {
					w = ir[i];
					x = pr[i];
					mu_tot = 0.0;
					for(int k = 0; k < K; k++) {
						mu_tot += (phi[w*K+k] + BETA) / (phitot[k] + WBETA) * (theta[d*K+k] + ALPHA) / (thetatot[d] + KALPHA);
					}
					perplexity -= (x * log(mu_tot));
				}
			}
			printf("\tIteration %d, Perplexity is: %.5lf\n", iter, exp(perplexity / tokens));
		}
	}

	tm->end();
	printf("\nTotal training time: %lf seconds.\n", tm->getTime());

	ParameterSet = true;
}


/***************************************************************************
 * RBP_doc implementation ---- residual accumulate by documents
 **************************************************************************/
RBP_doc::RBP_doc(string path, int k, int t, double alpha, double beta):LdaBase(path, k, t, alpha, beta, true) {
	residual = NULL;
	mu_new = NULL;
	seq = NULL;
}

RBP_doc::~RBP_doc() {
	if(residual != NULL) {
		delete residual;
		residual = NULL;
	}
	if(mu_new != NULL) {
		delete mu_new;
		mu_new = NULL;
	}
	if(seq != NULL) {
		delete seq;
		seq = NULL;
	}
}

void RBP_doc::init() {
	int topic, w;
	double x;

	/* allocate space to sufficient statistics */
	phi = new double[W*K];
	memset(phi, 0, sizeof(double)*W*K);
	phitot = new double[K];
	memset(phitot, 0, sizeof(double)*K);
	theta = new double[D*K];
	memset(theta, 0, sizeof(double)*D*K);
	thetatot = new double[D];
	memset(thetatot, 0, sizeof(double)*D);
	mu = new double[NNZ*K];
	memset(mu, 0, sizeof(double)*NNZ*K);
	mu_new = new double[K];
	memset(mu_new, 0, sizeof(double)*K);
	residual = new double[D];
	memset(residual, 0, sizeof(double)*D);
	seq = new int[D];
	for(int d = 0; d < D; d++) seq[d] = d;

	srand(time(0));

	/* initialize suffcient statistics */
	for(int d = 0; d < D; d++) {
		for(int i = jc[d]; i < jc[d+1]; i++) {
			w = ir[i];
			x = pr[i];
			
			/* randomly pick a topic */
			topic = rand() % K;
			mu[i*K+topic] = 1.0;				// assign this word the topic
			phi[w*K+topic] += x;
			phitot[topic] += x;
			theta[d*K+topic] += x;
			thetatot[d] += x;
		}
	}
}

void RBP_doc::LearnTopics() {
	int d, w;
	double x, mu_tot, perplexity;
	myTimer *tm = new myTimer();

	tm->start();

	init();

	for(int iter = 1; iter <= T; iter++) {

		/* clear residual */
		memset(residual, 0, sizeof(double)*D);

		/* passing message*/
		for(int j = 0; j < D; j++) {
			d = seq[j];
			for(int i = jc[d]; i < jc[d+1]; i++) {
				w = ir[i];
				x = pr[i];
				mu_tot = 0.0;
				for(int k = 0; k < K; k++) {
					/* exclude current word */
					phi[w*K+k] -= x*mu[i*K+k];
					phitot[k] -= x*mu[i*K+k];
					theta[d*K+k] -= x*mu[i*K+k];
					/* update mu stored in mu_new */
					mu_new[k] = (phi[w*K+k] + BETA) / (phitot[k] + WBETA) * (theta[d*K+k] + ALPHA);
					mu_tot += mu_new[k];
				}
				for(int k = 0; k < K; k++) {
					mu_new[k] /= mu_tot;
					/* compute residual */
					residual[d] += x*fabs(mu_new[k] - mu[i*K+k]);
					mu[i*K+k] = mu_new[k];
					/* update sufficient statistics ----> phi and theta */
					phi[w*K+k] += x*mu[i*K+k];
					phitot[k] += x*mu[i*K+k];
					theta[d*K+k] += x*mu[i*K+k];
				}
			}
		}

		/* sort residual and update seq to get the scan schedule in next iteration */
		dsort(D, residual, -1, seq);

		/* compute perplexity */
		if(iter % 10 == 0) {
			perplexity = 0.0;
			for(int d = 0; d < D; d++) {
				for(int i = jc[d]; i < jc[d+1]; i++) {
					w = ir[i];
					x = pr[i];
					mu_tot = 0.0;
					for(int k = 0; k < K; k++) {
						mu_tot += (phi[w*K+k] + BETA) / (phitot[k] + WBETA) * (theta[d*K+k] + ALPHA) / (thetatot[d] + KALPHA);
					}
					perplexity -= (x*log(mu_tot));
				}
			}
			printf("\tIteration %d, Perplexity is: %.5lf.\n", iter, exp(perplexity / tokens));
		}
	}

	tm->end();
	printf("\nTotal training time: %lf seconds.\n", tm->getTime());

	ParameterSet = true;

}


/***************************************************************************
 * RBP_voc implementation ---- residual accumulate by documents
 **************************************************************************/
RBP_voc::RBP_voc(string path, int k, int t, double alpha, double beta):LdaBase(path, k, t, alpha, beta, false) {
	residual = NULL;
	mu_new = NULL;
	seq = NULL;
}

RBP_voc::~RBP_voc() {
	if(residual != NULL) {
		delete residual;
		residual = NULL;
	}
	if(mu_new != NULL) {
		delete mu_new;
		mu_new = NULL;
	}
	if(seq != NULL) {
		delete seq;
		seq = NULL;
	}
}

void RBP_voc::init() {
	int topic, d;
	double x;

	/* allocate space to sufficient statistics */
	phi = new double[W*K];
	memset(phi, 0, sizeof(double)*W*K);
	phitot = new double[K];
	memset(phitot, 0, sizeof(double)*K);
	theta = new double[D*K];
	memset(theta, 0, sizeof(double)*D*K);
	thetatot = new double[D];
	memset(thetatot, 0, sizeof(double)*D);
	mu = new double[NNZ*K];
	memset(mu, 0, sizeof(double)*NNZ*K);
	mu_new = new double[K];
	memset(mu_new, 0, sizeof(double)*K);
	residual = new double[W];
	memset(residual, 0, sizeof(double)*W);
	seq = new int[W];
	for(int w = 0; w < D; w++) seq[w] = w;

	srand(time(0));

	/* initialize suffcient statistics */
	for(int w = 0; w < W; w++) {
		for(int i = jc[w]; i < jc[w+1]; i++) {
			d = ir[i];
			x = pr[i];
			
			/* randomly pick a topic */
			topic = rand() % K;
			mu[i*K+topic] = 1.0;				// assign this word the topic
			phi[w*K+topic] += x;
			phitot[topic] += x;
			theta[d*K+topic] += x;
			thetatot[d] += x;
		}
	}
}

void RBP_voc::LearnTopics() {
	int d, w;
	double x, mu_tot, perplexity;
	myTimer *tm = new myTimer();

	tm->start();

	init();

	for(int iter = 1; iter <= T; iter++) {

		/* clear residual */
		memset(residual, 0, sizeof(double)*D);

		/* passing message*/
		for(int j = 0; j < W; j++) {
			w = seq[j];
			for(int i = jc[w]; i < jc[w+1]; i++) {
				d = ir[i];
				x = pr[i];
				mu_tot = 0.0;
				for(int k = 0; k < K; k++) {
					/* exclude current word */
					phi[w*K+k] -= x*mu[i*K+k];
					phitot[k] -= x*mu[i*K+k];
					theta[d*K+k] -= x*mu[i*K+k];
					/* update mu stored in mu_new */
					mu_new[k] = (phi[w*K+k] + BETA) / (phitot[k] + WBETA) * (theta[d*K+k] + ALPHA);
					mu_tot += mu_new[k];
				}
				for(int k = 0; k < K; k++) {
					mu_new[k] /= mu_tot;
					/* compute residual */
					residual[w] += x*fabs(mu_new[k] - mu[i*K+k]);
					mu[i*K+k] = mu_new[k];
					/* update sufficient statistics ----> phi and theta */
					phi[w*K+k] += x*mu[i*K+k];
					phitot[k] += x*mu[i*K+k];
					theta[d*K+k] += x*mu[i*K+k];
				}
			}
		}

		/* sort residual and update seq to get the scan schedule in next iteration */
		dsort(W, residual, -1, seq);

		/* compute perplexity */
		if(iter % 10 == 0) {
			perplexity = 0.0;
			for(int w = 0; w < W; w++) {
				for(int i = jc[w]; i < jc[w+1]; i++) {
					d = ir[i];
					x = pr[i];
					mu_tot = 0.0;
					for(int k = 0; k < K; k++) {
						mu_tot += (phi[w*K+k] + BETA) / (phitot[k] + WBETA) * (theta[d*K+k] + ALPHA) / (thetatot[d] + KALPHA);
					}
					perplexity -= (x*log(mu_tot));
				}
			}
			printf("\tIteration %d, Perplexity is: %.5lf.\n", iter, exp(perplexity / tokens));
		}
	}

	tm->end();
	printf("\nTotal training time: %lf seconds.\n", tm->getTime());

	ParameterSet = true;

}
