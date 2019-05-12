#pragma once
#include<vector>

using namespace std;
typedef vector<double> vd;


//求解二维粒子运动
//坐标储存方式为{x,y,vx,vy}
//受力F 是x,y,t的函数
//init是初始值,t,tf是初,末时间
//N是总迭代次数
vd solve_2D_Newton(double t0, double tf, vd init, vd(*F)(double, double, double), int N) {
	//初始化,坐标,动量,时间
	vd r = { init[0],init[1] };
	vd v = { init[2],init[3] };
	double t = t0;
	//步长
	double delta_t = (tf - t0) / N;
	//迭代,应用四阶龙格库塔法
	vd k1(4, 0.), k2(4, 0.), k3(4, 0.), k4(4, 0.);
	vd F1, F2, F3, F4;
	for (int n = 0; n < N; n++) {
		//计算龙格库塔参数k1,k2,k3,k4
		F1 = F(r[0], r[1], t);
		k1 = { v[0],v[1],F1[0],F1[1] };

		F2 = F(r[0] + k1[0] * delta_t / 2, r[1] + k1[1] * delta_t / 2, t + delta_t / 2);
		k2 = { v[0]+k1[2]*delta_t/2,v[1]+ k1[3] * delta_t / 2,F2[0],F2[1] };

		F3 = F(r[0] + k2[0] * delta_t / 2, r[1] + k2[1] * delta_t / 2, t + delta_t / 2);
		k3 = { v[0] + k2[2] * delta_t / 2,v[1] + k2[3] * delta_t / 2,F3[0],F3[1] };

		F4 = F(r[0] + k3[0] * delta_t, r[1] + k3[1] * delta_t, t + delta_t);
		k4 = { v[0] + k3[2] * delta_t,v[1] + k3[3] * delta_t,F4[0],F4[1] };
		//完成迭代
		r[0] = r[0] + delta_t / 6.*(k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]);
		r[1] = r[1] + delta_t / 6.*(k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]);
		v[0] = v[0] + delta_t / 6.*(k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2]);
		v[1] = v[1] + delta_t / 6.*(k1[3] + 2 * k2[3] + 2 * k3[3] + k4[3]);
		t += delta_t;
	}
	return vd{ r[0],r[1],v[0],v[1] };
}

