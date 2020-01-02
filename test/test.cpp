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
    // �����������������
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
    // ����zģ��
    ArxModel zModel(na, nb, nd);
    zModel.setCoefA(a);
    zModel.setCoefB(b);

    /*��ȡ���ݲ���ʶ
    vector<double> input,output;
    while((bool)(fin >> u >> y))
    {
        input.push_back(u);
        input.push_back(y);
    }
    zModel = leastSquare(input,output,na,nb,nd);*/
    // ת��
    SModel sModel = d2c(zModel, Ts,Tustin);
    cout << sModel.num.transpose() << endl << sModel.den.transpose() << endl;
    system("pause");
}