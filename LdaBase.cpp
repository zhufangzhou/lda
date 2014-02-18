#include "LdaBase.h"
#include "util.h"

/*
 * Constructors
 */
LdaBase::LdaBase(string path, int k, int t, double alpha, double beta, bool document_major) {
	this->path = path;
	K = k;
	T = t;
	ALPHA = alpha;
	BETA = beta;

	ParameterSet = false;
	tokens = 0;

	phi = NULL;
	theta = NULL;
	phitot = NULL;
	thetatot = NULL;
	mu = NULL;
	ir = NULL;
	jc = NULL;
	pr = NULL;
	ReadData(path, document_major);
}

/*
 * Destructor
 */
LdaBase::~LdaBase() {
	if(ir != NULL) {
		delete ir;
		ir = NULL;
	}
	if(jc != NULL) {
		delete jc;
		jc = NULL;
	}
	if(pr != NULL) {
		delete pr;
		pr = NULL;
	}
	if(mu != NULL) {
		delete mu;
		mu = NULL;
	}
	if(theta != NULL) {
		delete theta;
		theta = NULL;
	}
	if(thetatot != NULL) {
		delete thetatot;
		thetatot = NULL;
	}
	if(phi != NULL) {
		delete phi;
		phi = NULL;
	}
	if(phitot != NULL) {
		delete phitot;
		phitot = NULL;
	}
}

/*
 * Read data from dataset
 */
void LdaBase::ReadData(string path, bool document_major) {
	FILE *fp = fopen(path.c_str(), "r"); 	
	char line[MAXLEN], *pch;
	int nnz, line_count;

	if(fp == NULL) {
		fprintf(stderr, "**********Text file cannot be found.**********\n");
		exit(1);
	}

	// get three integers --->   D  W  NNZ
	fgets(line, sizeof(line), fp);
	
	pch = strtok(line, " ");
	if((pch != NULL) && (atoi(pch) != 0)) {
		D = atoi(pch);
		printf("Documents #: %d\n", D);
	} else {
		fprintf(stderr, "**********File header error.**********\n");
		exit(1);
	}

	pch = strtok(NULL, " ");
	if((pch != NULL) && (atoi(pch) != 0)) {
		W = atoi(pch);
		printf("Vocabulary #: %d\n", W);
	} else {
		fprintf(stderr, "**********File header error.**********\n");
	}

	pch = strtok(NULL, " ");
	if((pch != NULL) && (atoi(pch) != 0)) {
		NNZ = atoi(pch);
		printf("NNZ #: %d\n", NNZ);
	} else {
		fprintf(stderr, "**********File header error.**********\n");
	}

	if(this->NNZ == 0) fprintf(stderr, "Empty file.\n");	

	WBETA = W * BETA;
	KALPHA = K * ALPHA;

	line_count = document_major ? D : W;

	jc = new int[line_count+1];
	ir = new int[NNZ];
	pr = new double[NNZ];

	nnz = 0;
	for(int i = 0; i < line_count; i++) {
		jc[i] = nnz;
		if(fgets(line, sizeof(line), fp) == NULL) break;
		pch = strtok(line, " ");
		if((pch != NULL) && (atoi(pch) != 0)) {
			ir[nnz] = atoi(pch) - 1;
			pch = strtok(NULL, " ");
			pr[nnz] = atof(pch);
			nnz++;
		}

		while((pch != NULL) && (atoi(pch) != 0)) {
			pch = strtok(NULL, " ");
			if((pch != NULL) && (atoi(pch) != 0)) {
				ir[nnz] = atoi(pch) - 1;
				pch = strtok(NULL, " ");
				pr[nnz] = atof(pch);
				nnz++;
			}
		}
	}
	jc[line_count] = nnz;
}

double* LdaBase::getPhi() {
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

double* LdaBase::getTheta() {
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
