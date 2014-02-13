prog = util.cpp LdaBase.cpp GibbsSampling.cpp main.cpp

lda: $(prog)
	g++ -g $(prog) -o out
