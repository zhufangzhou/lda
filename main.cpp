#include "GibbsSampling.h"
#include "util.h"


int main() {
	LdaBase *lda = new GibbsSampling("nips.txt", 10, 100, 20, 10, 0.01, 0.01);
	lda->LearnTopics();
	return 0;
}
