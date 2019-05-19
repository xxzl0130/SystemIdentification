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
    SModel sModel;
    sModel.num = B;
    sModel.den = A;
    auto cModel = c2d(sModel,0.01,Tustin);
    cout << cModel.getCoefB().transpose() << endl << cModel.getCoefA().transpose() << endl;
    system("pause");
}