#include "util.h"

int main() {
	int idx[5] = {0,1,2,3,4}, ridx[5] = {0,1,2,3,4};
	double val[5] = {3,1,4,2,5};
	quick_sort_des(val, idx, 0, 5);
	quick_sort_asc(idx, ridx, 0, 5);
	for(int i = 0; i < 5; i++) printf("%lf ", val[i]);
	printf("\n");
	for(int i = 0; i < 5; i++) printf("%d ", idx[i]);
	printf("\n");
	for(int i = 0; i < 5; i++) printf("%d ", ridx[i]);
	printf("\n");
	return 0;
}
