#include "LdaBase.h"
#include "util.h"

/*
 * Constructors
 */
LdaBase::LdaBase(string path, int k, int t, double alpha, double beta) {
	this->path = path;
	K = k;
	T = t;
	ALPHA = alpha;
	BETA = beta;

	ParameterSet = false;
	ReadData(path);
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
void LdaBase::ReadData(string path) {
	FILE *fp = fopen(path.c_str(), "r"); 	
	char line[MAXLEN], *pch;
	int nnz;

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
	ALPHA = K * ALPHA;

	// allocate free space to parameters
	phitot = new double[K];
	phi = new double[K*W];
	memset(phitot, 0, sizeof(double)*K);
	memset(phi, 0, sizeof(double)*K*W);

	thetatot = new double[D];
	theta = new double[D*K];
	memset(thetatot, 0, sizeof(double)*D);
	memset(theta, 0, sizeof(double)*D*K);

	jc = new int[D+1];
	ir = new int[NNZ];
	pr = new double[NNZ];

	nnz = 0;
	for(int d = 0; d < D; d++) {
		jc[d] = nnz;
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
	jc[D] = nnz;
}
