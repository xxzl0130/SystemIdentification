#include "LeastSquareEstimation.h"
#include "Model.h"
#include <iostream>
#include <fstream>
using namespace std;

int main()
{
    ArxModel dModel(2,1);
    auto B = Eigen::Vector2d(1, 100);
    auto A = Eigen::Vector3d(1, 110,10);
    dModel.setCoefB(B);
    dModel.setCoefA(A);
    auto cModel = d2c(dModel);
    cout << cModel.num.transpose() << endl << cModel.den.transpose() << endl;
    system("pause");
}