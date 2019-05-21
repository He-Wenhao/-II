#pragma once
#include<vector>
#include<array>
constexpr double Exp = 2.71828182845904523536028747;
constexpr double PI = 3.14159265358979323846264338;
const double omega = 0.057;
const double A0 = 1.325;
const double T = 110.23;
using namespace std;

typedef array<double,4> vd4;
typedef array<double, 2> vd2;

//线偏振光
vd2 line_light(double t) {
	vd2 result{ 0,0 };
	result[0] = -omega*pow(cos(omega*t/8.),2)*A0*cos(omega*t)
		+omega/8.*sin(omega*t/4.)*A0*sin(omega*t);
	return result;
}


//椭圆偏振光
vd2 elli_light(double t) {
	vd2 result{ 0,0 };
	result[0] = (-omega * pow(cos(omega*t / 8.), 2)*A0*cos(omega*t)
		+ omega / 8.*sin(omega*t / 4.)*A0*sin(omega*t))/sqrt(1.25);
	result[1] =( omega*pow(cos(omega*t/8.),2)*A0*sin(omega*t)
		+ omega / 8.*sin(omega*t / 4.)*A0*cos(omega*t))*0.5 / sqrt(1.25);
	return result;
}

//求向量的模
template<int n>
double abs(array<double,n> vec) {
	double result = 0.;
	for (auto x : vec) {
		result += x * x;
	}
	return sqrt(result);
}
//求向量模平方
template<int n>
double norm(array<double,n> vec) {
	double result = 0.;
	for (auto x : vec) {
		result += x * x;
	}
	return result;
}