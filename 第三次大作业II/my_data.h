#pragma once
#include<vector>
constexpr double E = 2.71828182845904523536028747;
constexpr double PI = 3.14159265358979323846264338;
const double omega = 0.057;
const double A0 = 1.325;
const double T = 110.23;
using namespace std;
typedef vector<double> vd;
//包络
double f(double t) {
	return pow(cos(PI*t / 4 / T), 2);
}


//线偏振光
vd line_light(double t) {
	vd result{ 0,0 };
	result[0] = -omega*f(t)*A0*cos(omega*t);
	return result;
}


//椭圆偏振光
vd elli_light(double t) {
	vd result{ 0,0 };
	result[0] = -omega*f(t)*A0*cos(omega*t)/sqrt(1.25);
	result[1] = omega*f(t)*A0*sin(omega*t)*0.5 / sqrt(1.25);
	return result;
}

//求向量的模
double abs(vd vec) {
	double result = 0;
	for (auto x : vec) {
		result += x * x;
	}
	return sqrt(result);
}
