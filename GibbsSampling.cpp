#include "GibbsSampling.h"
#include "util.h"

GibbsSampling::GibbsSampling(string path, int k, int t, int burn_in/*, int sample_lag*/, double alpha, double beta):LdaBase(path, k, t, alpha, beta, true) {
	BURN_IN = burn_in;
//	SAMPLE_LAG = sample_lag;
	tokens = 0;
	z = NULL;
	wd = NULL;
	doc = NULL;
	p = NULL;
}

GibbsSampling::~GibbsSampling() {
	if(z != NULL) {
		delete z;
		z = NULL;
	}
	if(wd != NULL) {
		delete wd;
		wd = NULL;
	}
	if(doc != NULL) {
		delete doc;
		doc = NULL;
	}
	if(p != NULL) {
		delete p;
		p = NULL;
	}
}

void GibbsSampling::init() {
	int d, w, x, k, topic;

	// allocate free space to parameters
	phitot = new double[K];
	phi = new double[K*W];
	memset(phitot, 0, sizeof(double)*K);
	memset(phi, 0, sizeof(double)*K*W);

	thetatot = new double[D];
	theta = new double[D*K];
	memset(thetatot, 0, sizeof(double)*D);
	memset(theta, 0, sizeof(double)*D*K);

	srand((unsigned)time(0));

	// count how many tokens in the these documents
	for(int i = 0; i < NNZ; i++) tokens += (int)pr[i];

	wd = new int[tokens];
	doc = new int[tokens];
	z = new int[tokens];

	p = new double[K];

	k = 0;
	for(d = 0; d < D; d++) {
		for(int i = jc[d]; i < jc[d+1]; i++) {
			w = ir[i];
			x = (int)pr[i];
			for(int j = 0; j < x; j++) {
				wd[k] = w;
				doc[k] = d;
				k++;
			}
		}
	}

	srand((unsigned)time(0));
	
	for(int i = 0; i < tokens; i++) {
		// random pick a topic from 0-K
		topic = rand() % K;
		w = wd[i];
		d = doc[i];

		z[i] = topic;
		phi[w*K + topic]++;
		theta[d*K + topic]++;
		phitot[topic]++;
		thetatot[d]++;
	}

}

/* token: index of the sampling word*/
int GibbsSampling::sampleTopic(int token) {
	int topic, w, d;
	topic = z[token];
	w = wd[token];
	d = doc[token];

	// decrement counts and sums of this token for the current topic
	phi[w*K+topic]--;
	phitot[topic]--;
	theta[d*K+topic]--;
	thetatot[d]--;
	
	for(int k = 0; k < K; k++) {
		// full conditional
		p[k] = (phi[w*K+k]+BETA) / (phitot[k]+WBETA) * (theta[d*K+k]+ALPHA);
	}
	

	// cumulate the prob
	for(int k = 1; k < K; k++) p[k] += p[k-1];

	// generate a random number from 0-sum(p)
	double r = (rand() / (double)(RAND_MAX)) * p[K-1];

	// index the topic in the cumulted prob array
	for(topic = 0; topic < K; topic++) {
		if(r < p[topic]) break;
	}

	// increment counts and sums of this token for the new sampled topic
	phi[w*K+topic]++;
	phitot[topic]++;
	theta[d*K+topic]++;
	thetatot[d]++;
	
	return topic;
}

void GibbsSampling::LearnTopics() {
	int topic, w, d;
	double perplexity, totprob;
	myTimer *tm = new myTimer();

	tm->start();
	/* initialize Markov chain */
	init();

	/* begin sample */
	for(int t = 0; t < T; t++) {
		
		// calculate perplexity
		if((t % 10 == 0) && (t != 0)) {
			perplexity = 0.0;
			for(int i = 0; i < tokens; i++) {
				w = wd[i];
				d = doc[i];
				totprob = 0.0;
				for(int k = 0; k < K; k++) {
					totprob += (phi[w*K+k]+BETA) / (phitot[k]+WBETA) * (theta[d*K+k]+ALPHA) / (thetatot[d]+KALPHA);
				}
				perplexity -= log(totprob);
			}
			printf("Iteration %d of %d:\t%.5lf\n", t, T, exp(perplexity/tokens));
		}

		for(int i = 0; i < tokens; i++) {
			topic = sampleTopic(i);
			z[i] = topic;
		}
	}

	tm->end();
	printf("Learning finished. Using time %.3lf seconds.\n", tm->getTime());
	ParameterSet = true;
}

double* GibbsSampling::getTheta() {
	double *p_theta = new double[W*K];

	if(!ParameterSet) {
		printf("Please call LearnTopic() first.\n");
		return NULL;
	}

	for(int w = 0; w < W; w++) {
		for(int k = 0; k < K; k++) {
			p_theta[w*K+k] = (phi[w*K+k]+BETA) / (phitot[k]+WBETA);
		}
	}
	return p_theta;
}

double* GibbsSampling::getPhi() {
	double *p_phi = new double[D*K];
	
	if(!ParameterSet) {
		printf("Please call LearnTopic() first.\n");
		return NULL;
	}

	for(int d = 0; d < D; d++) {
		for(int k = 0; k < K; k++) {
			p_phi[d*K+k] = (theta[d*K+k]+ALPHA) / (thetatot[d]+KALPHA);
		}
	}
	return p_phi;
}
