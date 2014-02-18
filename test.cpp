#include "util.h"

int main() {
	double value[]={2,3,1,5,4};
	int seq[]={3,0,2,4,1};
	quick_sort(value, seq, 0, 5);
	for(int i = 0; i < 5; i++) printf("%d ", seq[i]);
	printf("\n");
	return 0;
}
