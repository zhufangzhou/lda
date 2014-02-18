#include "util.h"
#include "VariationalBayes.h"

VariationalBayes::VariationalBayes(string path, int k, int t, double alpha, double beta, double em_converged, double var_converged):LdaBase(path, k, t, alpha, beta){
	p_log_phi = NULL;
	p_theta = NULL;
	var_gamma = NULL;
	var_digamma = NULL;
	var_phi = NULL;
	doc_len = NULL;
	DOC_MAX_LEN = 0;	

	EM_CONVERGED = em_converged;
	VAR_CONVERGED = var_converged;

	VAR_MAX_ITER = 10;
}

VariationalBayes::~VariationalBayes() {
	if(p_log_phi != NULL) {
		delete p_log_phi;
		p_log_phi = NULL;
	}
	if(p_theta != NULL) {
		delete p_theta;
		p_theta = NULL;
	}
	if(var_gamma != NULL) {
		delete var_gamma;
		var_gamma = NULL;
	}
	if(var_digamma != NULL) {
		delete var_digamma;
		var_digamma = NULL;
	}
	if(var_phi != NULL) {
		delete var_phi;
		var_phi = NULL;
	}
	if(doc_len != NULL) {
		delete doc_len;
		doc_len = NULL;
	}
}

void VariationalBayes::init() {
	int topic, w;
	double x;

	// allocate free space and initialize to zero (LDA model) 
	phi = new double[K*W];
	memset(phi, 0, sizeof(double)*K*W);
	phitot = new double[K];
	memset(phitot, 0, sizeof(double)*K);
	p_log_phi = new double[K*W];
	memset(p_log_phi, 0, sizeof(double)*K*W);
	p_theta = new double[K*D];
	memset(p_theta, 0, sizeof(double)*K*D);
	doc_len = new int[D];
	memset(doc_len, 0, sizeof(int)*D);

	srand(time(0));
	for(int d = 0; d < D; d++) {
		DOC_MAX_LEN = max(DOC_MAX_LEN, jc[d+1]-jc[d]);
		for(int i = jc[d]; i < jc[d+1]; i++) {
			w = ir[i];
			x = pr[i];
			tokens += x;
			doc_len[d] += x;

			// random pick a topic
			topic = rand() % K;
			phi[w*K+topic] += x;
			phitot[topic] += x;
		}
	}

	// allocate free space and initialize to zero (Variational distribution)
	var_gamma = new double[D*K];
	memset(var_gamma, 0, sizeof(double)*D*K);
	var_digamma = new double[K];
	memset(var_digamma, 0, sizeof(double)*K);
	var_phi = new double[K*DOC_MAX_LEN];
	memset(var_phi, 0, sizeof(double)*K*DOC_MAX_LEN);

	/* update p_phi from sufficient statistics */
	mle();
}

/*
 * Update parameter using suffcient statistics
 */
void VariationalBayes::mle() {
	for(int w = 0; w < W; w++) {
		for(int k = 0; k < K; k++) {
			if(phi[w*K+k] > 0)
				p_log_phi[w*K+k] = log(phi[w*K+k]+BETA) - log(phitot[k]+WBETA);
			else
				p_log_phi[w*K+k] = -100;
		}
	}
}

/*
 * Formula refer to paper 'Laten Dirichlet Allocation' (15)
 */
double VariationalBayes::compute_likelihood(int d) {
	double likelihood = 0.0, var_gamma_sum = 0.0, var_digamma_sum = 0.0, x;
	int w, word_count;

	/* number of words in document d*/
	word_count = jc[d+1] - jc[d];

	for(int k = 0; k < K; k++) var_gamma_sum += var_gamma[d*K+k];
	var_digamma_sum = digamma(var_gamma_sum);
	
	/* formula line 1 and line 4 [ lgamma is log(gamma(x)) ]*/
	likelihood = lgamma(KALPHA) - K*lgamma(ALPHA) - lgamma(var_gamma_sum); 
	for(int k = 0; k < K; k++) {
		likelihood += (ALPHA-1)*(var_digamma[k]-var_digamma_sum) + lgamma(var_gamma[d*K+k]) - (var_gamma[d*K+k]-1)*(var_digamma[k]-var_digamma_sum);
	}

	/* formula line 2, 3 and 5*/
	for(int n = 0; n < word_count; n++) {
		w = ir[jc[d]+n];				// word index
		x = pr[jc[d]+n];				// word count
		for(int k = 0; k < K; k++) {
			likelihood += x*var_phi[n*K+k]*(var_digamma[k]-var_digamma_sum + p_log_phi[w*K+k] - log(var_phi[n*K+k]));
		}
	}

	return likelihood;
}

/*
 * Get the tightest lower-bound
 */
double VariationalBayes::variational_inference(int d) {
	double converged = 1, likelihood = 0, old_likelihood = 0, var_phi_sum, x, var_digamma_sum;
	int w, iter, word_count;

	/* number of words in document d*/
	word_count = jc[d+1] - jc[d];

	/* initialize var_gamma and var_phi */
	for(int k = 0; k < K; k++) {
		var_gamma[d*K+k] = ((double)doc_len[d] / (double)K) + ALPHA;
		var_digamma[k] = digamma(var_gamma[d*K+k]);	
		for(int n = 0; n < word_count; n++) {
			var_phi[n*K+k] = 1.0/K;
		}
	}
	
	iter = 0;
	while((iter++ < VAR_MAX_ITER) && ((converged < 0) || (converged > VAR_CONVERGED))) {
		for(int n = 0; n < word_count; n++) {
			w = ir[jc[d]+n]; // get word index
			x = pr[jc[d]+n]; // get word count
			var_phi_sum = 0; // var_phi_sum is used to normalize var_phi
			for(int k = 0; k < K; k++) {
				/* phi = beta * exp(digamma(gamma)), here put it in log space, get back when normalize */
				var_phi[n*K+k] = p_log_phi[w*K+k] + var_digamma[k];
				if(k == 0) var_phi_sum = var_phi[n*K+k];
				else var_phi_sum = log_sum(var_phi_sum, var_phi[n*K+k]);			
			}
			
			for(int k = 0; k < K; k++) {
				/* normalize var_phi */
				var_phi[n*K+k] = exp(var_phi[n*K+k] - var_phi_sum);
			}
		}

		/* gamma = ALPHA + sigma(phi)*/
		for(int k = 0; k < K; k++)  {
			var_gamma[d*K+k] = ALPHA;
		}
		for(int n = 0; n < word_count; n++) {
			x = pr[jc[d]+n];
			for(int k = 0; k < K; k++) {
				var_gamma[d*K+k] += x*var_phi[n*K+k];	
			}
		}
		/* update digamma(gamma) and compute converged*/
		for(int k = 0; k < K; k++) {
			var_digamma[k] = digamma(var_gamma[d*K+k]);
		}
		likelihood = compute_likelihood(d);
		converged = (old_likelihood - likelihood) / old_likelihood;
		old_likelihood = likelihood;
	}


	/* recover theta using gamma */
	var_digamma_sum = 0.0;
	for(int k = 0; k < K; k++) var_digamma_sum += var_gamma[d*K+k];
	var_digamma_sum = digamma(var_digamma_sum);
	for(int k = 0; k < K; k++) {
		p_theta[d*K+k] = exp(var_digamma[k] - var_digamma_sum);
	}
	
	return likelihood;
}

/*
 * Run EM algorithm
 */
void VariationalBayes::LearnTopics() {
	double converged = 1.0, likelihood = 0, x, mutot, perplexity;
	int iter = 1, w;
	myTimer *tm = new myTimer();	

	/* malloc space for variables, initialize variables*/
	init();

	tm->start();
	while((converged > EM_CONVERGED) || (converged < 0) || (iter <= 3)) {
		/* clear suffcient statistics */
		memset(phi, 0, sizeof(double)*K*W);
		memset(phitot, 0, sizeof(double)*K);

		likelihood = 0;

		/* E-step */
		for(int d = 0; d < D; d++) {
			likelihood += variational_inference(d);

			/* update sufficient statistics */
			for(int i = jc[d], n = 0; i < jc[d+1]; i++, n++) {
				w = ir[i];
				x = pr[i];
				for(int k = 0; k < K; k++) {
					phi[w*K+k] += x*var_phi[n*K+k];
					phitot[k] += x*var_phi[n*K+k];
				}
			}
		}

		/* M-step */
		mle();


		/* compute perplexity */
		if(iter % 10 == 0) {
			printf("Likelihood: %.10lf\n", likelihood);
			perplexity = 0.0;
			for(int d = 0; d < D; d++) {
				for(int i = jc[d]; i < jc[d+1]; i++) {
					w = ir[i];
					x = pr[i];
					mutot = 0.0;
					for(int k = 0; k < K; k++) {
						mutot += exp(p_log_phi[w*K+k])*p_theta[d*K+k];
					}
					perplexity -= (x*log(mutot));	
				}
			}
			printf("\tIteration %d:\t%.5lf\n", iter, exp(perplexity / tokens));
		}
		iter++;
		if(iter >= 100) break;
	}
	tm->end();
	printf("\tIteration: %d last %lf seconds.\n", iter, tm->getTime());
}

double* VariationalBayes::getPhi() {
	double* p_phi = new double[W*K];
	for(int w = 0; w < W; w++) {
		for(int k = 0; k < K; k++) {
			p_phi[w*K+k] = exp(p_log_phi[w*K+k]);
		}
	}
	return p_phi;
}

double* VariationalBayes::getTheta() {
	return p_theta;
}
