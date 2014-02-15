prog = util.cpp LdaBase.cpp GibbsSampling.cpp VariationalBayes.cpp main.cpp

lda: $(prog)
	g++ -g $(prog) -o out -lm
