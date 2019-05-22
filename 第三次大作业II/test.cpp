#include<random>
using namespace std;
default_random_engine rand_e;
uniform_real_distribution<double> rand_u(0, 1);
#include<math.h>
#include<iostream>
#include"my_data.h"
#include"generate_e.h"
#include"solve_eq.h"
#include<map>
#include<fstream>
#include<time.h>
#include<string>
#include<thread>
typedef array<double, 4> vd4;
typedef array<double, 2> vd2;

//线偏振光的总力场
vd2 line_F(double x, double y, double t) {
	//库仑力
	vd2 Fk{ -x / pow(x*x + y * y+0.04,3. / 2.),-y / pow(x*x + y * y+0.04,3. / 2.) };
	return vd2{ Fk[0] - line_light(t)[0],Fk[1] - line_light(t)[1] };
}

//椭圆偏振光的总力场
vd2 elli_F(double x, double y, double t) {
	//库仑力
	vd2 Fk{ -x / pow(x*x + y * y + 0.04,3. / 2.),-y / pow(x*x + y * y + 0.04,3 / 2.) };
	return vd2{ Fk[0] - elli_light(t)[0],Fk[1] - elli_light(t)[1] };
}

//第二步函数:给定初态求末态
//N:总迭代步数

//线偏振 解牛顿方程
vd4 line_rf_vf(state init,int N) {
	return solve_2D_Newton(init.t, 2 * T, init.r_v, line_F, N);
}

//椭圆偏振 解牛顿方程
vd4  elli_rf_vf(state init, int N) {
	return solve_2D_Newton(init.t, 2 * T, init.r_v, elli_F, N);
}


//第三步函数:求无穷远动量值
vd2 infp(vd4 rf_vf) {
	vd2 rf = { rf_vf[0],rf_vf[1] };
	vd2 pf = { rf_vf[2],rf_vf[3] };
	//无穷远初动量
	double p_inf = sqrt(norm(pf) - 2 / abs(rf));
	//角动量及龙哥楞次矢量
	double L = rf[0] * pf[1]-rf[1] * pf[0] ;
	vd2 a = { pf[1] * L - rf[0] / abs(rf),-pf[0] * L - rf[1] / abs(rf) };
	//计算结果
	vd2 result = { -p_inf*L* a[1]-a[0],p_inf*L* a[0]-a[1] };
	result = { result[0] * p_inf / (1 + pow(p_inf*L,2)),result[1] * p_inf / (1 + pow(p_inf*L,2)) };
	return result;
}



//最后一步统计电子
//线偏振
void line_test(int N_data0,double delta_tsample){
	double timenow = clock();//记录运行时间
	map<double, int> px, py;//储存x,y方向数目
	map<vd2, int> pvec;//储存二维分布的数目
	//初始化pvec
	for (double x = -1.51; x <= 1.51; x += 0.02) {
		for (double y = -1.51; y < 1.51; y += 0.02) {
			pvec[vd2{ x,y }] = 0;
		}
	}
	//遍历时间t
	for (double tsample = -2 * T; tsample <= 2 * T; tsample += delta_tsample) {
		//计算权重
		double Et = abs(line_light(tsample));
		double N_data_t = N_data0 * pow(Exp, -2. / 3. / Et) / pow(Et, 3. / 2.);
		for (double i = 0; i < N_data_t; i++) {
			//生成电子
			state e_sample = generate_sample2D(line_light,tsample);
			//激光结束时的状态
			vd4 fstate = line_rf_vf(e_sample, int(2*T-tsample)*10);
			//无穷远初动量
			vd2 inf_p = infp(fstate);
			if (abs(inf_p[0]) < 1.5&&abs(inf_p[1]) < 1.5) {
				//离散化
				double px_ = 0.02*(floor(inf_p[0] * 50) + 0.5);
				double py_ = 0.02*(floor(inf_p[1] * 50) + 0.5);
				//统计个数
				px[px_]++;
				py[py_]++;
				pvec[vd2{ px_,py_ }]++;
			}

		}
		if (abs(tsample - floor(tsample) )< 0.11) {
			cout << tsample << endl;//进程可视化
		}
	}
	//将数据写入txt
	ofstream osx;
	osx.open("line_temp_x.txt");
	for (auto x : px) {
		osx << x.first << "\t" << x.second << endl;
	}
	ofstream osy;
	osy.open("line_temp_y.txt");
	for (auto y : py) {
		osy << y.first << "\t" << y.second << endl;
	}
	ofstream osvec;
	osvec.open("line_temp_vec.txt");
	for (auto vec : pvec) {
		osvec << vec.first[0]<<"\t"<<vec.first[1] << "\t" << vec.second << endl;
	}
	cout << clock() - timenow;//输出运行时间
}


//最后一步统计电子
//椭圆偏振
void elli_test(int N_data0, double delta_tsample) {
	double timenow = clock();//记录运行时间
	map<double, int> px, py;//储存x,y方向数目
	map<vd2, int> pvec;//储存二维分布的数目
	//初始化
	for (double x = -1.51; x <= 1.51; x += 0.02) {
		for (double y = -1.51; y <= 1.51; y += 0.02) {
			pvec[vd2{ x,y }] = 0;
		}
	}
	//遍历时间t,
	double Et, N_data_t;
	for (double tsample = -2 * T; tsample <= 2 * T; tsample += delta_tsample) {
		//计算权重
		Et = abs(elli_light(tsample));
		N_data_t = N_data0 * pow(Exp, -2. / 3. / Et) / pow(Et, 3. / 2.);
		for (int i = 0; i < N_data_t; i++) {
			//生成电子
			state e_sample = generate_sample2D(elli_light, tsample);
			//激光结束时的状态
			vd4 fstate = elli_rf_vf(e_sample, int(2 * T - tsample) * 10);
			//无穷远初动量
			vd2 inf_p = infp(fstate);
			if (abs(inf_p[0]) < 1.5&&abs(inf_p[1]) < 1.5) {
				//离散化
				double px_ = 0.02*(floor(inf_p[0] * 50) + 0.5);
				double py_ = 0.02*(floor(inf_p[1] * 50) + 0.5);
				//统计个数
				px[px_]++;
				py[py_]++;
				pvec[vd2{ px_,py_ }]++;
			}

		}
		if (tsample - floor(tsample) < 0.11) {
			cout << tsample << endl;//进程可视化
		}
	}
	//将数据写入txt
	ofstream osx;
	osx.open("elli_temp_x.txt");
	for (auto x : px) {
		osx << x.first << "\t" << x.second << endl;
	}
	ofstream osy;
	osy.open("elli_temp_y.txt");
	for (auto y : py) {
		osy << y.first << "\t" << y.second << endl;
	}
	ofstream osvec;
	osvec.open("elli_temp_vec.txt");
	for (auto vec : pvec) {
		osvec << vec.first[0] << "\t" << vec.first[1] << "\t" << vec.second << endl;
	}
	cout << clock() - timenow;//输出运行时间
}






int main() {
	elli_test(5e5, 0.04);
	line_test(5e5, 0.04);
	system("pause");
}  