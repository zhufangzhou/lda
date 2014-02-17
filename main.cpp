#include "GibbsSampling.h"
#include "VariationalBayes.h"
#include "util.h"


int main() {
	//LdaBase *lda = new GibbsSampling("nips.txt", 10, 100, 20, 0.01, 0.01);
	LdaBase *lda = new VariationalBayes("nips.txt", 10, 100, 0.01, 0.01, 1e-3, 1e-5);
	lda->LearnTopics();
	return 0;
}
