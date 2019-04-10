#include "LeastSquareEstimation.h"
#include <iostream>
#include <fstream>
using namespace std;

int main()
{
    fstream fin("data.txt", ios_base::in);
    double u, y;
    std::vector<double> input, output;
    while((bool)(fin >> u >> y))
    {
        input.push_back(u);
        output.push_back(y);
    }
    auto m = leastSquare(input, output, 2, 1);
    cout << m.coefA.transpose() << endl << m.coefB.transpose() << endl;
    system("pause");
}