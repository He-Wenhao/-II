#pragma once
#include"my_data.h"
#include<array>
using namespace std;
//const double max_E = 1 / pow(E, 3. / 2.);
const double max_g = sqrt(2 * Exp / PI);


typedef array<double, 4> vd4;
typedef array<double, 2> vd2;

struct state {
	vd4 r_v;//坐标和动量
	double t;//时间
};


//用第二类舍选法
//对分布标准正态分布取样
double generate_V() {
	while(1) {
		//对exp^(-x)抽样
		double eta = -log((double)rand_u(rand_e));
		//决定是否取这个值
		double zeta2 = max_g*(double)rand_u(rand_e);
		if (pow(eta-1,2)<=-2*log(zeta2)) {
			if ((double)rand_u(rand_e) > 0.5) { 
				return eta; 
			}
			else {
				return -eta;
			}
		}
	}
}

//为第二步生成样品vd={x,y,vx,vy}
//选取v的方向使得角动量沿z轴
state generate_sample2D(vd2(*light)(double),double t) {
	//t为电离时刻

	//计算电场大小
	double Et = abs(light(t));
	//对标准正态分布取样
	double V = generate_V();
	//计算v垂直
	double v_ver = sqrt(Et/2. )*V;//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	//单位径向量
	vd2 er = { -light(t)[0] / Et,-light(t)[1] / Et };
	double r0 = 0.5 / Et;
	vd4 r_v{ r0*er[0],r0*er[1],-v_ver * er[1],v_ver*er[0] };
	return state{ r_v,t };
}