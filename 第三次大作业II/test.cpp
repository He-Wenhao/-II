#include<iostream>
#include"my_data.h"
#include"generate_e.h"
#include"solve_eq.h"
#include<map>
#include<fstream>
#include<time.h>
#include<string>
#include<thread>
using namespace std;
typedef vector<double> vd;

/*
vd false_light(double t) {
	return vd{t,0};
}

//统计检验
void static_test() {
	double tnow = clock();
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
	cout << clock() - tnow;
}
*/
/*

vd solve_2D_Newton_fake(double t0, double tf, vd init, vd(*F)(double, double, double), int N) {
	ofstream os;
	os.open("temp.txt");
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
		os << r[0] << "\t" << r[1] << endl;
		//计算龙格库塔参数k1,k2,k3,k4
		F1 = F(r[0], r[1], t);
		k1 = { v[0],v[1],F1[0],F1[1] };

		F2 = F(r[0] + k1[0] * delta_t / 2, r[1] + k1[1] * delta_t / 2, t + delta_t / 2);
		k2 = { v[0] + k1[2] * delta_t / 2,v[1] + k1[3] * delta_t / 2,F2[0],F2[1] };

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

vd F_fake(double x, double y, double t) {
	return vd{ -x / pow(x*x + y * y,3 / 2.),-y / pow(x*x + y * y,3 / 2.) };
}*/

//线偏振光的总力场
vd line_F(double x, double y, double t) {
	//库仑力
	vd Fk{ -x / pow(x*x + y * y+0.04,3 / 2.),-y / pow(x*x + y * y+0.04,3 / 2.) };
	return vd{ Fk[0] + line_light(t)[0],Fk[1] + line_light(t)[1] };
}

//椭圆偏振光的总力场
vd elli_F(double x, double y, double t) {
	//库仑力
	vd Fk{ -x / pow(x*x + y * y + 0.04,3 / 2.),-y / pow(x*x + y * y + 0.04,3 / 2.) };
	return vd{ Fk[0] + elli_light(t)[0],Fk[1] + elli_light(t)[1] };
}

//第二步函数:给定初态求末态
//N:总迭代步数

//线偏振
vd line_rf_vf(state init,int N) {
	return solve_2D_Newton(init.t, 2 * T, init.r_v, line_F, N);
}

//椭圆偏振
vd  elli_rf_vf(state init, int N) {
	return solve_2D_Newton(init.t, 2 * T, init.r_v, elli_F, N);
}


//第三步函数:求无穷远动量值
vd infp(vd rf_vf) {
	vd rf = { rf_vf[0],rf_vf[1] };
	vd pf = { rf_vf[2],rf_vf[3] };
	//无穷远初动量
	double p_inf = sqrt(norm(pf) - 2 / abs(rf));
	//角动量及龙哥楞次矢量
	vd L = { rf[0] * pf[1],-rf[1] * pf[0] };
	vd a = { pf[0] * L[1] - rf[0] / abs(rf),-pf[1] * L[0] - rf[1] / abs(rf) };
	//计算结果
	vd result = { p_inf*L[0] * a[1]-a[0],-p_inf*L[1] * a[0]-a[1] };
	result = { result[0] * p_inf / (1 + pow(p_inf,2)*norm(L)),result[1] * p_inf / (1 + pow(p_inf,2)*norm(L)) };
	return result;
}



//最后一步统计电子
//线偏振
void line_test(int N_t,int N_data){
	double timenow = clock();//!!!!!!!!!
	map<double, int> px, py;
	map<vd, int> pvec;

	for (int i = 0; i < N_data; i++) {
			//生成电子
			state e_sample = generate_sample2D(line_light);
			//激光结束时的状态
			vd fstate = line_rf_vf(e_sample, N_t);
			//无穷远初动量
			vd inf_p = infp(fstate);
			if (abs(inf_p[0]) < 1.5&&abs(inf_p[1]) < 1.5) {
				//离散化
				double px_ = 0.02*((int)(inf_p[0] * 50) + 0.5);
				double py_ = 0.02*((int)(inf_p[1] * 50) + 0.5);
				//统计个数
				px[px_]++;
				py[py_]++;
				pvec[vd{ px_,py_ }]++;
				if (i % 1000 == 0) {
					cout << i << endl;//!!!!!!!!
				}
			}
	}

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
	cout << clock() - timenow;//!!!!!!!!!!!!!!!!!!!
}
//椭圆偏振
void elli_test(int N_t, int N_data) {
	map<double, int> px, py;
	map<vd, int> pvec;

	for (int i = 0; i < N_data; i++) {
		//生成电子
		state e_sample = generate_sample2D(elli_light);
		//激光结束时的状态
		vd fstate = elli_rf_vf(e_sample, N_t);
		//无穷远初动量
		vd inf_p = infp(fstate);
		if (abs(inf_p[0]) < 1.5&&abs(inf_p[1]) < 1.5) {
			//离散化
			double px_ = 0.02*((int)(inf_p[0] * 50) + 0.5);
			double py_ = 0.02*((int)(inf_p[1] * 50) + 0.5);
			//统计个数
			px[px_]++;
			py[py_]++;
			pvec[vd{ px_,py_ }]++;
			if (i % 1000 == 0) {
				cout << i << endl;//!!!!!!!!
			}
		}
	}

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
}


//线偏振 多线程版本
void line_test_thread(int N_t, int N_data) {
	double timenow = clock();//!!!!!!!!!
	map<double, int> px, py;
	map<vd, int> pvec;

	for (int i = -1.54 / 0.02; i < 1.54 / 0.02; i++) {
		px[0.02*i+0.01] = 0;
		py[0.02*i+0.01] = 0;
	}

	for (int i = -1.54 / 0.02; i < 1.54 / 0.02; i++) {
		for (int j = -1.54 / 0.02; j < 1.54 / 0.02; j++) {
			pvec[vd{ 0.02*i+0.01,0.02*j+0.01 }] = 0;
		}
	}

	for (int i = 0; i < N_data; i++) {
		thread t([&px,&py, &pvec, N_t, i] {
			srand(i);
			//生成电子
			state e_sample = generate_sample2D(line_light);
			//激光结束时的状态
			vd fstate = line_rf_vf(e_sample, N_t);
			//无穷远初动量
			vd inf_p = infp(fstate);
			if (abs(inf_p[0]) < 1.5&&abs(inf_p[1]) < 1.5) {
				//离散化
				double px_ = 0.02*((int)(inf_p[0] * 50) + 0.5);
				double py_ = 0.02*((int)(inf_p[1] * 50) + 0.5);
				//统计个数
				px[px_]++;
				py[py_]++;
				pvec[vd{ px_,py_ }]++;
				if (i % 1000 == 0) {
					cout << i << endl;//!!!!!!!!
				}
			}
		});
		t.detach();
	}
	cout << clock() - timenow;//!!!!!!!!!!!!!!!!!!!

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
		osvec << vec.first[0] << "\t" << vec.first[1] << "\t" << vec.second << endl;
	}
}


//椭圆偏振 多线程版本
void elli_test_thread(int N_t, int N_data) {
	double timenow = clock();//!!!!!!!!!
	map<double, int> px, py;
	map<vd, int> pvec;

	for (int i = -1.54 / 0.02; i < 1.54 / 0.02; i++) {
		px[0.02*i+0.01] = 0;
		py[0.02*i+0.01] = 0;
	}

	for (int i = -1.54 / 0.02; i < 1.54 / 0.02; i++) {
		for (int j = -1.54 / 0.02; j < 1.54 / 0.02; j++) {
			pvec[vd{ 0.02*i+0.01,0.02*j +0.01}] = 0;
		}
	}

	for (int i = 0; i < N_data; i++) {
		thread t([&px, &py, &pvec, N_t, i] {
			srand(i);
			//生成电子
			state e_sample = generate_sample2D(elli_light);
			//激光结束时的状态
			vd fstate = elli_rf_vf(e_sample, N_t);
			//无穷远初动量
			vd inf_p = infp(fstate);
			if (abs(inf_p[0]) < 1.5&&abs(inf_p[1]) < 1.5) {
				//离散化
				double px_ = 0.02*((int)(inf_p[0] * 50) + 0.5);
				double py_ = 0.02*((int)(inf_p[1] * 50) + 0.5);
				//统计个数
				px[px_]++;
				py[py_]++;
				pvec[vd{ px_,py_ }]++;
				if (i % 1000 == 0) {
					cout << i << endl;//!!!!!!!!
				}
			}
		});
		t.detach();
	}

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
	cout << clock() - timenow;//!!!!!!!!!!!!!!!!!!!
}
int main() {
	line_test_thread(400, 1e5);
	//line_test(400, 5e3);
	//elli_test(400, 1e5);
	//static_test();
	system("pause");
}