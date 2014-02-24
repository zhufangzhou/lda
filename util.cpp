#include "util.h"

/*
 * Class myTimer
 */
myTimer::myTimer() {
	counts = 0;
}

void myTimer::start() {
	if(counts != 0)  {
		printf("Opeartion failed.\n");
		return;
	}
	counts++;
	v_start = clock();
}

void myTimer::end() {
	if(counts != 1) {
		printf("Operation failed.\n");
		return;
	}
	counts++;
	v_end = clock();
}

double myTimer::getTime() {
	counts = 0;
	return ((double)(v_end - v_start) / CLOCKS_PER_SEC);
}

//==========================================================================

/* data used for isort*/
static int *icomp_vec;
static double *dcomp_vec;

/*
 * From NIST Handbook of Mathematical Functions, using formula 5.11.2
 */
/*
double digamma(double x) {
	double y = 1.0 / (x*x);
	return log(x) - 0.5/x - y*(0.08333333333333333333-y*(-0.00833333333333333333-y*(0.00396825396825396825+y*0.00416666666666666667)));
}*/
double digamma(double x) {
	double y;
	x += 6;
	y = 1.0 / (x*x);
	return log(x) - 0.5/x - y*(0.08333333333333333333-y*(-0.00833333333333333333-y*(0.00396825396825396825+y*0.00416666666666666667))) - 1/(x-1) - 1/(x-2) - 1/(x-3) - 1/(x-4) - 1/(x-5) - 1/(x-6);
}


/*
 * return log(a+b)
 */
double log_sum(double log_a, double log_b) {
	if(log_a < log_b) {
		return log_b + log(1 + exp(log_a - log_b));
	} else {
		return log_a + log(1 + exp(log_b - log_a));
	}
}

static int icomp(const void *a, const void *b) {
	int i = *(int*)a, j = *(int*)b;
	return icomp_vec[i] - icomp_vec[j];
}

/* dir = +1 ascending order, dir = -1 descending order */
void isort(int n, int *value, int dir, int *idx) {
	icomp_vec = new int[n];
	for(int i = 0; i < n; i++) {
		icomp_vec[i] = dir * value[i];
		idx[i] = i;
	}
	qsort(idx, n, sizeof(int), icomp);

	delete icomp_vec;
	icomp_vec = NULL;
}

static int dcomp(const void *a, const void *b) {
	int i = *(int*)a, j = *(int*)b;
	return dcomp_vec[i] > dcomp_vec[j] ? 1 : -1;
}

/* dir = +1 ascending order, dir = -1 descending order */
void dsort(int n, double *value, int dir, int *idx) {
	dcomp_vec = new double[n];
	for(int i = 0; i < n; i++) {
		dcomp_vec[i] = dir * value[i];
		idx[i] = i;
	}
	qsort(idx, n, sizeof(int), dcomp);

	delete dcomp_vec;
	dcomp_vec = NULL;
}
