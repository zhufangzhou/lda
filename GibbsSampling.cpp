#include "GibbsSampling.h"
#include "util.h"

GibbsSampling::GibbsSampling(string path, int k, int t, int burn_in, double alpha, double beta):LdaBase(path, k, t, alpha, beta, true) {
	BURN_IN = burn_in;
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
	int d, w, x, topic;
	int k;

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
	
	for(int k = 0; k < K; k++) {
		// full conditional
		p[k] = (phi[w*K+k]+BETA) / (phitot[k]+WBETA) * (theta[d*K+k]+ALPHA);
	}
	

	// cumulate the probability
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
	for(int t = 1; t <= T; t++) {
		
		for(int i = 0; i < tokens; i++) {
			topic = sampleTopic(i);
			z[i] = topic;
		}

		// calculate perplexity
		if(t % 10 == 0) {
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
	}

	tm->end();
	printf("Learning finished. Using time %.3lf seconds.\n", tm->getTime());
	ParameterSet = true;
}


/*
 * fast Gibbs Sampling implementation
 * */
fastGibbsSampling::fastGibbsSampling(string path, int k, int t, int burn_in, double alpha, double beta):LdaBase(path, k, t, alpha, beta, true) {
	BURN_IN = burn_in;

	wd = NULL;
	doc = NULL;
	z = NULL;
	sum_p = NULL;
	phi_norm = NULL;
	theta_norm = NULL;
	phitot_idx = NULL;
	phitot_ridx = NULL;
}

fastGibbsSampling::~fastGibbsSampling() {
	if(wd != NULL) {
		delete wd;
		wd = NULL;
	}
	if(doc != NULL) {
		delete doc;
		doc = NULL;
	}
	if(z != NULL) {
		delete z;
		z = NULL;
	}
	if(sum_p != NULL) {
		delete sum_p;
		sum_p = NULL;
	}
	if(phi_norm != NULL) {
		delete phi_norm;
		phi_norm = NULL;
	}
	if(theta_norm != NULL) {
		delete theta_norm;
		theta_norm = NULL;
	}
	if(phitot_idx != NULL) {
		delete phitot_idx;
		phitot_idx = NULL;
	}
	if(phitot_ridx != NULL) {
		delete phitot_ridx;
		phitot_ridx = NULL;
	}
}

void fastGibbsSampling::init() {
	int w, d, topic;
	int k;
	double x;

	/* allocate space to suffcient statistics */
	phi = new double[W*K];
	memset(phi, 0, sizeof(double)*W*K);
	phitot = new double[K];
	memset(phitot, 0, sizeof(double)*K);
	theta = new double[D*K];
	memset(theta, 0, sizeof(double)*D*K);
	thetatot = new double[D];
	memset(thetatot, 0, sizeof(double)*D);
	sum_p = new double[K];
	memset(sum_p, 0, sizeof(double)*K);

	wd = new int[tokens];
	doc = new int[tokens];
	z = new int[tokens];

	k = 0;
	for(int d = 0; d < D; d++) {
		for(int i = jc[d]; i < jc[d+1]; i++) {
			w = ir[i];
			x = pr[i];
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
		phitot[topic]++;
		theta[d*K + topic]++;
		thetatot[d]++;
	}

	/* set up norm */
	phi_norm = new double[W];
	memset(phi_norm, 0, sizeof(double)*W);
	theta_norm = new double[D];
	memset(theta_norm, 0, sizeof(double)*D);
	phitot_idx = new int[K];
	phitot_ridx = new int[K];
	for(int k = 0; k < K; k++) phitot_idx[k] = phitot_ridx[k] = k;

	for(int k = 0; k < K; k++) {
		for(int w = 0; w < W; w++) {
			phi_norm[w] += SQUARE(phi[w*K+k] + BETA);
		}
		for(int d = 0; d < D; d++) {
			theta_norm[d] += SQUARE(theta[d*K+k] + ALPHA);
		}
	}

	quick_sort_des(phitot, phitot_idx, 0, K);
	quick_sort_asc(phitot_idx, phitot_ridx, 0, K);


}

void fastGibbsSampling::update_sort(int n, double *value, int *idx, int *ridx, int id, bool des) {
	int tmp;
	
	if(des) {
		// update after decrease
		id = ridx[id];
		while((id < n-1) && (value[idx[id]] < value[idx[id+1]])) {
			// update idx
			tmp = idx[id];
			idx[id] = idx[id+1];
			idx[id+1] = tmp;

			// update ridx
			tmp = ridx[idx[id]];
			ridx[idx[id]] = ridx[idx[id+1]];
			ridx[idx[id+1]] = tmp;

			id++;
		}
	} else {
		// update after increase
		id = ridx[id];
		while((id > 0) && (value[idx[id]] > value[idx[id-1]])) {
			// update idx
			tmp = idx[id];
			idx[id] = idx[id-1];
			idx[id-1] = tmp;

			// update ridx
			tmp = ridx[idx[id]];
			ridx[idx[id]] = ridx[idx[id-1]];
			ridx[idx[id-1]] = tmp;
		}
	}
}

int fastGibbsSampling::sampleTopic(int token) {
	int w, d, topic, topic_new;
	double u, tp, v_phi_norm, v_theta_norm, v_phitot, zk = 0.0, zk_old, r;

	srand((unsigned)time(0));

	w = wd[token];
	d = doc[token];

	topic = z[token];

	/* compute phi_norm and theta_norm */
	phi_norm[w] -= 2*(phi[w*K + topic] + BETA) - 1;
	theta_norm[d] -= 2*(theta[d*K + topic] + ALPHA) - 1;

	// exclude current token
	phi[w*K + topic]--;
	phitot[topic]--;
	theta[d*K + topic]--;

	update_sort(K, phitot, phitot_idx, phitot_ridx, topic, true);

	v_phi_norm = phi_norm[w];
	v_theta_norm = theta_norm[d];
	v_phitot = phitot[phitot_idx[K-1]]+WBETA;

	memset(sum_p, 0, sizeof(double)*K);	
	u = rand() / (double)(RAND_MAX);

	for(int k = 0; k < K; k++) {
		tp = (phi[w*K+k]+BETA)/(phitot[k]+WBETA)*(theta[d*K+k]+ALPHA);	
		if(k == 0) sum_p[k] = tp;
		else sum_p[k] = sum_p[k-1] + tp;

		/* phi_norm(l+1:K) and theta_norm(l+1:K) */
		v_phi_norm -= SQUARE(phi[w*K+k] + BETA);
		v_theta_norm -= SQUARE(theta[d*K+k] + ALPHA);

		zk_old = zk;
		zk = sum_p[k] + sqrt(v_phi_norm*v_theta_norm) / v_phitot;

		r = u*zk;
		if(r > sum_p[k]) continue;
		else {
			if((k == 0) || (r > sum_p[k-1])) {
				topic_new = k;
				break;
			} else {
				u = (u*zk_old - sum_p[k-1]) * zk / (zk_old - zk);
				for(int kk = 0; kk < k; kk++) {
					if(sum_p[kk] >= u) {
						topic_new = kk;
						break;
					}
				}
				break;
			}
		}
		
	}

	phi_norm[w] += 2*(phi[w*K+topic_new] + BETA) + 1;
	theta_norm[d] += 2*(phi[d*K+topic_new] + ALPHA) + 1;

	// increase current token 
	phi[w*K + topic_new]++;
	phitot[topic_new]++;
	theta[d*K + topic_new]++;

	update_sort(K, phitot, phitot_idx, phitot_ridx, topic_new, false);
	return topic_new;
}

void fastGibbsSampling::LearnTopics() {
	int w, d, topic;
	double x, perplexity, mu_tot;

	init();

	for(int iter = 1; iter <= T; iter++) {

		for(int i = 0; i < tokens; i++) {
			topic = sampleTopic(i);
			z[i] = topic;
		}

		/* compute perplexity */
		if(iter % 10 == 0) {
			perplexity = 0.0;
			for(int i = 0; i < tokens; i++) {
				w = wd[i];
				d = doc[i];
				mu_tot = 0.0;
				for(int k = 0; k < K; k++) {
					mu_tot += (phi[w*K+k] + BETA) / (phitot[k] + WBETA) * (theta[d*K+k] + ALPHA) / (thetatot[d] + KALPHA);
				}
				perplexity -= log(mu_tot);
			}
			printf("\tIteration %d, Perplexity is: %.5lf\n", iter, exp(perplexity / tokens));
		}
	}

	init();

}
