#include "LeastSquareEstimation.h"
#include "Model.h"
#include <iostream>
#include <fstream>
using namespace std;

int main()
{
    fstream fin("data.txt", ios_base::in);
    fstream fout("out.txt", ios::out);
    ArxModel model(2, 1);
    model.setCoefA(Eigen::Vector3d(1, -1.97791, 0.97824));
    model.setCoefB(Eigen::Vector2d(0.00999945, -0.00978186));
    double u, y;
    std::vector<double> input, output;
    while((bool)(fin >> u >> y))
    {
        input.push_back(u);
        output.push_back(y);
        fout << model.update(u) << endl;
    }
    auto m = leastSquare(input, output, 2, 1);
    cout << m.getCoefA().transpose() << endl << m.getCoefB().transpose() << endl;
    system("pause");
}