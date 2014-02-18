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

/*
 * quick_sort (descend order)
 */
void quick_sort(double *value, int *seq, int left, int right) {
	int pivot, l = left+1, r = right-1, tmp, mid;
	if(right-left <= 1) return ;
	else {
		/* swap middle to left and choose middle as pivot */
		mid = (left + right) / 2;
		tmp = seq[mid];
		seq[mid] = seq[left];
		seq[left] = tmp;

		pivot = seq[left];
		while(l <= r) {
			while(l<=r && value[seq[l]] > value[pivot]) l++;
			while(l<=r && value[seq[r]] < value[pivot]) r--;
			if(l<=r) {
				tmp = seq[l];
				seq[l] = seq[r];
				seq[r] = tmp;
			}
		}
		tmp = seq[r];
		seq[r] = seq[left];
		seq[left] = tmp;
		quick_sort(value, seq, left, r);
		quick_sort(value, seq, r+1, right);
	}
}
