#include "LeastSquareEstimation.h"
#include "Model.h"
#include <iostream>
#include <fstream>
using namespace std;

int main()
{
    SModel sModel;
    sModel.num = Eigen::Vector2d(1, 100);
    sModel.den = Eigen::Vector3d(1, 110,10);
    auto zModel = c2d(sModel);
    cout << zModel.getCoefB().transpose() << endl << zModel.getCoefA().transpose() << endl;
    system("pause");
}