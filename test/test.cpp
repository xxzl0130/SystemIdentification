#include "LeastSquareEstimation.h"
#include "Model.h"
#include <iostream>
#include <fstream>
#include <vector>
using namespace std;

int main()
{
    fstream fin("data.txt", ios_base::in);
    double u, y, Ts;
    int na, nb, nd;
    // 输入阶数、采样周期
    cin >> na >> nb >> nd >> Ts;
    Eigen::VectorXd a, b;
    a.resize(na + 1);
    b.resize(nb + 1);
    for (auto i = 0; i <= na; ++i)
    {
        cin >> a[i];
    }
    for (auto i = 0; i <= nb; ++i)
    {
        cin >> b[i];
    }
    // 创建z模型
    ArxModel zModel(na, nb, nd);
    zModel.setCoefA(a);
    zModel.setCoefB(b);

    /*读取数据并辨识
    vector<double> input,output;
    while((bool)(fin >> u >> y))
    {
        input.push_back(u);
        input.push_back(y);
    }
    zModel = leastSquare(input,output,na,nb,nd);*/
    // 转换
    SModel sModel = d2c(zModel, Ts,Tustin);
    cout << sModel.num.transpose() << endl << sModel.den.transpose() << endl;
    system("pause");
}