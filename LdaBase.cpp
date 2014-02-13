#include "LdaBase.h"
#include "util.h"

/*
 * Constructors
 */
LdaBase::LdaBase() {
}

LdaBase::LdaBase(string path, int k, int t, double alpha, double beta) {
	this->path = path;
	this->K = k;
	this->T = t;
	this->ALPHA = alpha;
	this->BETA = beta;

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
		this->D = atoi(pch);
		printf("Documents #: %d\n", this->D);
	} else {
		fprintf(stderr, "**********File header error.**********\n");
		exit(1);
	}

	pch = strtok(NULL, " ");
	if((pch != NULL) && (atoi(pch) != 0)) {
		this->W = atoi(pch);
		printf("Vocabulary #: %d\n", this->W);
	} else {
		fprintf(stderr, "**********File header error.**********\n");
	}

	pch = strtok(NULL, " ");
	if((pch != NULL) && (atoi(pch) != 0)) {
		this->NNZ = atoi(pch);
		printf("NNZ #: %d\n", this->NNZ);
	} else {
		fprintf(stderr, "**********File header error.**********\n");
	}

	if(this->NNZ == 0) fprintf(stderr, "Empty file.\n");	

	this->WBETA = this->W * this->BETA;
	this->ALPHA = this->K * this->ALPHA;

	// allocate free space to parameters
	this->phitot = new double[this->K];
	this->phi = new double[this->K * this->W];
	memset(this->phitot, 0, sizeof(double)*this->K);
	memset(this->phi, 0, sizeof(double)*this->K*this->W);

	this->thetatot = new double[this->D];
	this->theta = new double[this->D * this->K];
	memset(this->thetatot, 0, sizeof(double)*this->D);
	memset(this->theta, 0, sizeof(double)*this->D*this->K);

	this->jc = new int[this->D + 1];
	this->ir = new int[this->NNZ];
	this->pr = new double[this->NNZ];

	nnz = 0;
	for(int d = 0; d < this->D; d++) {
		this->jc[d] = nnz;
		if(fgets(line, sizeof(line), fp) == NULL) break;
		pch = strtok(line, " ");
		if((pch != NULL) && (atoi(pch) != 0)) {
			this->ir[nnz] = atoi(pch) - 1;
			pch = strtok(NULL, " ");
			this->pr[nnz] = atof(pch);
			nnz++;
		}

		while((pch != NULL) && (atoi(pch) != 0)) {
			pch = strtok(NULL, " ");
			if((pch != NULL) && (atoi(pch) != 0)) {
				this->ir[nnz] = atoi(pch) - 1;
				pch = strtok(NULL, " ");
				this->pr[nnz] = atof(pch);
				nnz++;
			}
		}
	}
	this->jc[this->D] = nnz;
}
