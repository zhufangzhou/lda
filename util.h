#ifndef __UTIL
#define __UTIL

#include <iostream>
#include <cmath>
#include <string>
#include <cstdlib>
#include <string.h>
#include <vector>
#include <stdio.h>
#include <ctime>


using namespace std;

#define MAXLEN 88888				// max characters in one line
#define max(a,b) ((a)>(b)?(a):(b))
#define min(a,b) ((a)<(b)?(a):(b))

class myTimer {
	private:
		clock_t v_start;
		clock_t v_end;
		int counts;
	public:
		myTimer();
		void start();
		void end();
		double getTime();
}; 

double log_sum(double log_a, double log_b);
double digamma(double x);



#endif
