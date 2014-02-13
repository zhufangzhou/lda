#include "GibbsSampling.h"

GibbsSampling::GibbsSampling() {

}

GibbsSampling::GibbsSampling(string path, int k, int t, int burn_in, int sample_lag, double alpha, double beta):LdaBase(path, k, t, alpha, beta) {
	this->BURN_IN = burn_in;
	this->SAMPLE_LAG = sample_lag;
	this->tokens = 0;
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
}

void GibbsSampling::init() {
	int d, w, x, k, topic;


	// count how many tokens in the these documents
	for(int i = 0; i < this->NNZ; i++) this->tokens += (int)pr[i];

	this->wd = new int[this->tokens];
	this->doc = new int[this->tokens];
	this->z = new int[this->tokens];

	k = 0;
	for(d = 0; d < this->D; d++) {
		for(int i = jc[d]; i < jc[d+1]; i++) {
			w = ir[i];
			x = (int)pr[i];
			for(int j = 0; j < x; j++) {
				this->wd[k] = w;
				this->doc[k] = d;
				k++;
			}
		}
	}

	srand((unsigned)time(0));
	
	for(int i = 0; i < this->tokens; i++) {
		// random pick a topic from 0-K
		topic = rand() % this->K;
		w = wd[i];
		d = doc[i];

		z[i] = topic;
		phi[w*this->K + topic]++;
		theta[d*this->K + topic]++;
		phitot[topic]++;
		thetatot[d]++;
	}

}

/* token: index of the sampling word*/
int GibbsSampling::sampleTopic(int token) {
	int topic, w, d;
	double *p;
	topic = this->z[token];
	w = this->wd[token];
	d = this->doc[token];

	// decrement counts and sums of this token for the current topic
	this->phi[w*this->K+topic]--;
	this->phitot[topic]--;
	this->theta[d*this->K+topic]--;
	this->thetatot[d]--;
	
	p = new double[this->K];
	for(int k = 0; k < this->K; k++) {
		// full conditional
		p[k] = (this->phi[w*this->K+k]+this->BETA) / (this->phitot[k]+this->WBETA)
			 * (this->theta[d*this->K+k]+this->ALPHA);
	}
	
	// cumulate the probabilities
	for(int k = 1; k < this->K; k++) p[k] += p[k-1];

	// generate a random number from 0-sum(p)
	srand((unsigned)time(0));
	double r = (rand() / (double)(RAND_MAX)) * p[this->K-1];

	// index the topic in the cumulted prob array
	for(topic = 0; topic < this->K; topic++) {
		if(r < p[topic]) break;
	}

	// increment counts and sums of this token for the new sampled topic
	this->phi[w*this->K]+topic++;
	this->phitot[topic]++;
	this->theta[d*this->K+topic]++;
	this->thetatot[d]++;
	
	// free space
	if(p != NULL) {
		delete p;
		p = NULL;
	}
	return topic;
}

void GibbsSampling::LearnTopics() {
	int topic, w, d;
	double perplexity, totprob;
	/* initialize Markov chain */
	init();

	/* begin sample */
	for(int t = 0; t < this->T; t++) {
		
		// calculate perplexity
		if((t % 10 == 0) && (t != 0)) {
			perplexity = 0.0;
			for(int i = 0; i < this->tokens; i++) {
				w = this->wd[i];
				d = this->doc[i];
				totprob = 0.0;
				for(int k = 0; k < this->K; k++) {
					totprob += (this->phi[w*this->K + k]+this->BETA) / (this->phitot[k]+this->WBETA)
							 * (this->theta[d*this->K + k]+this->ALPHA) / (this->thetatot[d]+this->KALPHA);
				}
				perplexity -= log(totprob);
			}
			printf("Iteration %d of %d:\t%.5lf\n", t, this->T, exp(perplexity/this->tokens));
		}

		for(int i = 0; i < this->tokens; i++) {
			topic = sampleTopic(i);
			z[i] = topic;
		}
	}
}
