#include "util.h"
#include "VariationalBayes.h"

VariationalBayes::VariationalBayes(string path, int k, int t, double alpha, double beta, double em_converged, double var_converged):LdaBase(path, k, t, alpha, beta){
	p_log_phi = NULL;
	var_gamma = NULL;
	var_old_gamma = NULL;
	var_digamma = NULL;
	var_phi = NULL;
	doc_len = NULL;
	DOC_MAX_LEN = 0;	

	EM_CONVERGED = em_converged;
	VAR_CONVERGED = var_converged;
}

VariationalBayes::~VariationalBayes() {
	if(p_log_phi != NULL) {
		delete p_log_phi;
		p_log_phi = NULL;
	}
	if(var_gamma != NULL) {
		delete var_gamma;
		var_gamma = NULL;
	}
	if(var_old_gamma != NULL) {
		delete var_old_gamma;
		var_old_gamma = NULL;
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
	var_phi = new double[K*DOC_MAX_LEN];
	memset(var_phi, 0, sizeof(double)*K*DOC_MAX_LEN);

	/* update p_phi from sufficient statistics */
	mle();
}

void VariationalBayes::mle() {
	for(int w = 0; w < W; w++) {
		for(int k = 0; k < K; k++) {
			p_log_phi[w*K+k] = log(phi[w*K+k]+BETA) - log(phitot[k]+WBETA);
		}
	}
}

double VariationalBayes::variational_inference(int d) {
	double converged = 1, likelihood = 0, old_likelihood = 0, var_phi_sum, x;
	int w;

	/* initialize var_gamma and var_phi*/
	for(int k = 0; k < K; k++) {
		var_gamma[d*K+k] = (doc_len[d] / (double)K) + ALPHA;
		var_digamma[k] = digamma(var_gamma[d*K+k]);	
		for(int n = 0; n < doc_len[d]; n++) {
			var_phi[n*K+k] = 1.0/K;
		}
	}
	
	while((converged < 0) || (converged > VAR_CONVERGED)) {
		for(int n = 0; n < doc_len[d]; n++) {
			w = ir[jc[d]+n]; // get word index
			x = pr[jc[d]+n]; // get word count
			var_phi_sum = 0; // var_phi_sum is used to normalize var_phi
			for(int k = 0; k < K; k++) {
				/* phi = beta * exp(digamma(gamma)), here put it in log space, get back when normalize */
				var_phi[n*K+k] = p_log_phi[w*K+k] + var_digamma[k];

				var_phi_sum = log_sum(var_phi_sum, var_phi[n*K+k]);			
			}
			
			for(int k = 0; k < K; k++) {
				/* normalize var_phi */
				var_phi[n*K+k] = exp(var_phi[n*K+k] - var_phi_sum);
			}
		}

		/* gamma = ALPHA + sigma(phi)*/
		for(int k = 0; k < K; k++)  {
			var_old_gamma[k] = var_gamma[d*K+k];
			var_gamma[d*K+k] = ALPHA;
		}
		for(int n = 0; n < doc_len[d]; n++) {
			for(int k = 0; k < K; k++) {
				var_gamma[d*K+k] += var_phi[n*K+k];	
			}
		}
		converged = 0.0;
		/* update digamma(gamma) and compute converged*/
		for(int k = 0; k < K; k++) {
			var_digamma[k] = digamma(var_gamma[d*K+k]);
			converged += abs(var_old_gamma[k] - var_gamma[d*K+k]);
		}
		// converged = abs(change in gamma) / K
		converged /= K;
	}
}

void VariationalBayes::LearnTopics() {
	double converged = 1.0, likelihood = 0, x;
	int iter = 0, w;

	/* malloc space for variables, initialize variables*/
	init();

	while((converged > EM_CONVERGED) || (converged < 0) || (iter++ <= 3)) {
		
		/* clear suffcient statistics */
		memset(phi, 0, sizeof(double)*K*W);
		memset(phitot, 0, sizeof(double)*K);

		/* E-step */
		for(int d = 0; d < D; d++) {
			likelihood += variational_inference(d);
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

		/* perplexity */
	}
}

double* VariationalBayes::getPhi() {

}

double* VariationalBayes::getTheta() {

}
