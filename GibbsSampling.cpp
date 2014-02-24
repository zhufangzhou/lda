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

	min_phitot = 0.0;
	wd = NULL;
	doc = NULL;
	z = NULL;
	sum_p = NULL;
	phi_norm = NULL;
	theta_norm = NULL;
	p = NULL;
	theta_idx = NULL;
	theta_ridx = NULL;
	d_last = NULL;
	w_last = NULL;
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
	if(theta_idx != NULL) {
		delete theta_idx;
		theta_idx = NULL;
	}
	if(theta_ridx != NULL) {
		delete theta_ridx;
		theta_ridx = NULL;
	}
	if(p != NULL) {
		delete p;
		p = NULL;
	}
	if(d_last != NULL) {
		delete d_last;
		d_last = NULL;
	}
	if(w_last != NULL) {
		delete w_last;
		w_last = NULL;
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
	theta_idx = new int[D*K];
	theta_ridx = new int[D*K];
	thetatot = new double[D];
	memset(thetatot, 0, sizeof(double)*D);
	sum_p = new double[K];
	memset(sum_p, 0, sizeof(double)*K);

	p = new double[K];
	memset(p, 0, sizeof(double)*K);

	wd = new int[tokens];
	doc = new int[tokens];
	z = new int[tokens];
	memset(z, 0, sizeof(int)*tokens);

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

	w_last = new int[W];
	memset(w_last, 0xff, sizeof(int)*W);
	d_last = new int[D];
	memset(d_last, 0xff, sizeof(int)*D);


	/* initial min_phitot */
	min_phitot = phitot[0];
	for(int k = 0; k < K; k++) min_phitot = min(min_phitot, phitot[k]);
}

void fastGibbsSampling::update_sort(int n, double *value, int *idx, int *ridx, int ip, int im) {
	int tmp;
		
	/*
	 * Increment
	 * did ++ get bigger than prev
	 */
	ip = ridx[ip];
	while((ip > 0) && (value[idx[ip]] > value[idx[ip-1]])) {
		/* swap idx */
		tmp = idx[ip];
		idx[ip] = idx[ip-1];
		idx[ip-1] = tmp;

		/* swap ridx */
		tmp = ridx[idx[ip]];
		ridx[idx[ip]] = ridx[idx[ip-1]];
		ridx[idx[ip-1]] = tmp;
	}

	/*
	 * Decrement
	 * did -- get smaller than next
	 */
	im = ridx[im];
	while((im < n) && (value[idx[im]] < value[idx[im+1]])) {
		/* swap idx */
		tmp = idx[im];
		idx[im] = idx[im+1];
		idx[im+1] = tmp;

		/* swap ridx */
		tmp = ridx[idx[im]];
		ridx[idx[im]] = ridx[idx[im+1]];
		ridx[idx[im+1]] = tmp;
	}
}

int fastGibbsSampling::sampleTopic(int token, int iter) {
	int w, d, topic;
	double p_tot, r, u, v_phi_norm, v_theta_norm, prob, zk = 0.0, zk_old;

	/* set up random seed */
	srand((unsigned)time(0));

	/* get current token infomation */
	w = wd[token];
	d = doc[token];
	topic = z[token];

	/* substract current token from counts before sample */
	phi[w*K+ topic]--;
	theta[d*K + topic]--;
	phitot[topic]--;

	if(iter == 1) { /* standard sample */
		/* set up phitot norm */
		min_phitot = min(min_phitot, phitot[topic]);

		/* sort theta (each K-theta array sort at the end of scanning its doc)*/
		if((token == tokens - 1) || (doc[token] != doc[token+1])) {
			dsort(K, theta + d*K, -1, theta_idx + d*K); /* descending */
			isort(K, theta_idx + d*K, 1, theta_ridx + d*K); /* ascending */
		}

		/* set up norms */
		phi_norm[w] = 0;
		theta_norm[d] = 0;
		for(int k = 0; k < K; k++) {
			phi_norm[w] += SQUARE(phi[w*K+k] + BETA);
			theta_norm[d] += SQUARE(theta[d*K+k] + ALPHA);
		}

		p_tot = 0.0;
		for(int k = 0; k < K; k++) {
			p[k] = (phi[w*K+k] + BETA) / (phitot[k] + WBETA) * (theta[d*K+k] + ALPHA);	
			p_tot += p[k];
		}

		// sample a topic 
		r = (rand() / (double)RAND_MAX) * p_tot;
		for(int k = 0; k < K; k++) {
			if(r < p[k]) {
				topic_new = k;
				break;
			} else {
				r -= p[k];
			}
		}

	} else { /* fast sample */
		/* prepare phitot_norm */
		if(topic_new != topic) {
			min_phitot = min(min_phitot, phitot[topic]);
			min_phitot = min(min_phitot, phitot[topic_new]);
		}
		
		/* prepare phi_norm */
		if(w_last[w] != topic) {
			phi_norm[w] += 2*(phi[w*K+ w_last[w]] - phi[w*K + topic] - 1);
		}
		v_phi_norm = phi_norm[w];

		/* prepare theta_norm */
		if(d_last[d] != topic) {
			theta_norm[d] += 2*(theta[d*K + d_last[d]] - theta[d*K + topic] - 1);
			update_sort(K, theta+d*K, theta_idx+d*K, theta_ridx+d*K, d_last[d], topic);
		}
		v_theta_norm = theta_norm[d];

		u = rand() / (double)RAND_MAX;
		for(int j = 0, k; j < K; j++) {
			/* scan topic in theta order descending */
			k = theta_idx[d*K + j];
			prob = (phi[w*K+k] + BETA) / (phitot[k] + WBETA) * (theta[d*K+k] + ALPHA);	

			if(j == 0) sum_p[j] = prob;
			else sum_p[j] = sum_p[j-1] + prob;

			/* update current norm-based bounds */
			v_phi_norm -= SQUARE(phi[w*K+k] + BETA);
			v_theta_norm -= SQUARE(theta[d*K+k] + ALPHA);

			/* make norm bigger or eaqul to zero */
			v_phi_norm = v_phi_norm < 0 ? 0: v_phi_norm;
			v_theta_norm = v_theta_norm < 0 ? 0 : v_theta_norm;

			zk_old = zk;
			zk = sum_p[j] + sqrt(v_phi_norm*v_theta_norm) / (min_phitot+ WBETA);

			if(u*zk > sum_p[j]) continue; /* conitnue to next k */
			else {
				if((j == 0) || (u*zk > sum_p[j-1])) {
					topic_new = theta_idx[d*K + j];
					break;
				} else {
					u = (u*zk_old - sum_p[j-1]) * zk / (zk_old - zk);
					for(int jj = 0; jj < j; jj++) {
						if(sum_p[jj] >= u) {
							topic_new = theta_idx[d*K + jj];
							break;
						}
					}
					break;
				}
			}
		}
	}

	phi[w*K + topic_new]++;
	theta[d*K + topic_new]++;
	phitot[topic_new]++;

	w_last[w] = topic_new;
	d_last[d] = topic_new;
	return topic_new;
}

void fastGibbsSampling::LearnTopics() {
	int w, d, topic;
	double x, perplexity, mu_tot;
	myTimer *tm = new myTimer();

	tm->start();

	init();

	for(int iter = 1; iter <= T; iter++) {

		for(int i = 0; i < tokens; i++) {
			topic = sampleTopic(i, iter);
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

	tm->end();
	printf("Learning finished. Using time %.3lf seconds.\n", tm->getTime());
	ParameterSet = true;
}
