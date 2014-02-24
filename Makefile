prog = util.cpp LdaBase.cpp GibbsSampling.cpp VariationalBayes.cpp BeliefPropagation.cpp main.cpp

lda: $(prog)
	g++ -O2 -g $(prog) -o out -lm
