#pragma once
#include"my_data.h"
using namespace std;
const double max_E = 1 / pow(E, 3. / 2.);
const double max_g = sqrt(2 * E / PI);

struct state {
	vd r_v;//坐标和动量
	double t;//时间
};


//用第一类舍选法
//根据分布Exp^(-2/3/E)/E^(-3/2)生成时间t
//参数是光强函数
double generate_time(vd (*light)(double)) {
	while(1) {
		//生成+-2T之间的随机时间
		double t = -2 * T + 4 * T*(double)rand() / RAND_MAX;
		double Et = abs(light(t));
		//决定是否取这个值
		double zeta = max_E * (double)rand() / RAND_MAX;
		if (zeta <= pow(E, -2. / 3. / Et) / pow(Et, 3. / 2.)) {
			return t;
		}
	}
}

//用第二类舍选法
//对分布标准正态分布取样
double generate_V() {
	while(1) {
		//对exp^(-x)抽样
		double eta = -log((double)rand() / RAND_MAX);
		//决定是否取这个值
		double zeta2 = max_g*(double)rand() / RAND_MAX;
		if (pow(eta-1,2)<=-2*log(zeta2)) {
			return eta;
		}
	}
}

//(这里利用了电场仅有x,y分量的特性)
//生成样品 vd={x,y,z,vx,vy,vz}
state generate_sample3D(vd(*light)(double)) {
	//电离时刻
	double t = generate_time(light);
	//计算电场大小
	double Et = abs(light(t));
	//对标准正态分布取样
	double V = generate_V();
	//计算v垂直
	double v_ver = sqrt(Et / 2)*V;
	//随机确定v垂直的方向
	double vz = -1 + 2 * (double)rand() / RAND_MAX;
	double vxy = sqrt(1 - vz * vz);
	//v的z分量和xy分量
	vz = vz * v_ver;
	vxy = vxy * v_ver;

	//单位径向量
	vd er = { -light(t)[0] / Et,-light(t)[1] / Et };
	double r0 = 0.5 / Et;
	vd r_v{ r0*er[0],r0*er[1],0.,-vxy * er[1],vxy*er[0],vz };
	return state{ r_v,t };
}



//为第二步生成样品vd={x,y,vx,vy}
//选取v的方向使得角动量沿z轴
state generate_sample2D(vd(*light)(double)) {
	//电离时刻
	double t = generate_time(light);
	//计算电场大小
	double Et = abs(light(t));
	//对标准正态分布取样
	double V = generate_V();
	//计算v垂直
	double v_ver = sqrt(Et / 2)*V;
	//单位径向量
	vd er = { -light(t)[0] / Et,-light(t)[1] / Et };
	double r0 = 0.5 / Et;
	vd r_v{ r0*er[0],r0*er[1],-v_ver * er[1],v_ver*er[0] };
	return state{ r_v,t };
}