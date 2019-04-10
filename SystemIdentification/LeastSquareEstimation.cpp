#include "LeastSquareEstimation.h"
using namespace Eigen;

ArxModel leastSquare(const std::vector<double>& inputs, const std::vector<double>& outputs, unsigned na, unsigned nb,
    unsigned nd)
{
    return leastSquare(inputs.data(), outputs.data(), std::min(inputs.size(), outputs.size()),na,nb,nd);
}

ArxModel leastSquare(const double* inputs, const double* outputs, const unsigned size, unsigned na, unsigned nb,
    unsigned nd)
{
    auto nCoef = na + nb + 1;
    VectorXd y;
    MatrixXd phi(size, nCoef);
    y.resize(size);
    memcpy_s(y.data(), size * sizeof(double), outputs, size * sizeof(double));

    if(size < 2 * nCoef)
    {//too less input data to solve
        return ArxModel();
    }

    // deal with index which may less than 0
    for(auto i = 0;i < nCoef;++i)
    {
        for(auto j = 1;j <= na;++j)
        {
            if(i - j < 0)
            {
                phi(i, j - 1) = 0.0;
            }
            else
            {
                phi(i, j - 1) = -outputs[i - j];
            }
        }
        for(auto j = 0;j <= nb;++j)
        {
            if(int(i - nd - j) < 0)
            {
                phi(i, na + j) = 0;
            }
            else
            {
                phi(i, na + j) = inputs[i - nd - j];
            }
        }
    }

    for(auto i = nCoef;i < size;++i)
    {
        for (auto j = 1; j <= na; ++j)
        {
            phi(i, j - 1) = -outputs[i - j];
        }
        for (auto j = 0; j <= nb; ++j)
        {
            phi(i, na + j) = inputs[i - nd - j];
        }
    }

    //svd分解求最小二乘解
    VectorXd x = phi.jacobiSvd(ComputeFullU | ComputeFullV).solve(y);
    //VectorXd x = (phi.transpose() * phi).inverse() * phi.transpose() * y;
    ArxModel model(na, nb, nd);
    model.coefA << 1,x.head(na);
    model.coefB << x.tail(nb + 1);

    return model;
}
