#include<iostream>
#include"my_data.h"
#include"generate_e.h"
#include<map>
#include<fstream>
#include<time.h>
using namespace std;
typedef vector<double> vd;


vd false_light(double t) {
	return vd{t,0};
}

//统计检验
void static_test() {
	double delta=5;
	int Ndata = 100000;
	ofstream os;
	os.open("temp.txt");
	map<double, double> result;
	for (int i = 0; i < Ndata; i++) {
		double x = generate_time(false_light);
		//double x = generate_V();
		int xk = (int) x/delta;
		result[xk*delta] += 1;
	}
	for (auto& x : result) {
		x.second = x.second / Ndata;
		os << x.first << "\t" << x.second << endl;
	}
}




int main() {
	srand((unsigned)time(NULL));
	static_test();
	system("pause");
}