#include "LeastSquareEstimation.h"
using namespace Eigen;

ArxModel leastSquare(const std::vector<double>& inputs, const std::vector<double>& outputs, unsigned na, unsigned nb,
    unsigned nd)
{
    return leastSquare(inputs.data(), outputs.data(), std::min(inputs.size(), outputs.size()), na, nb, nd);
}

ArxModel leastSquare(const double* inputs, const double* outputs, const unsigned size, unsigned na, unsigned nb,
    unsigned nd)
{
    auto nCoef = na + nb + 1;
    VectorXd y;
    MatrixXd phi(size, nCoef);
    y.resize(size);
    memcpy_s(y.data(), size * sizeof(double), outputs, size * sizeof(double));

    if (size < 2 * nCoef)
    {//too less input data to solve
        return ArxModel();
    }

    // deal with index which may less than 0
    for (auto i = 0; i < nCoef; ++i)
    {
        for (auto j = 1; j <= na; ++j)
        {
            if (i - j < 0)
            {
                phi(i, j - 1) = 0.0;
            }
            else
            {
                phi(i, j - 1) = -outputs[i - j];
            }
        }
        for (auto j = 0; j <= nb; ++j)
        {
            if (int(i - nd - j) < 0)
            {
                phi(i, na + j) = 0;
            }
            else
            {
                phi(i, na + j) = inputs[i - nd - j];
            }
        }
    }

    for (auto i = nCoef; i < size; ++i)
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

    //svd for least square
    VectorXd x = phi.jacobiSvd(ComputeFullU | ComputeFullV).solve(y);
    //VectorXd x = (phi.transpose() * phi).inverse() * phi.transpose() * y;
    ArxModel model(na, nb, nd);
    Eigen::VectorXd A, B;
    A.resize(na + 1);
    B.resize(nb + 1);
    A << 1, x.head(na);
    B << x.tail(nb + 1);
    model.setCoefA(A);
    model.setCoefB(B);

    return model;
}

RLS::RLS(unsigned na, unsigned nb, unsigned nd) :
    nCoefs(na + nb + 1),
    nStep(0),
    na(na),
    nb(nb),
    nd(nd), lambda(1.0), miu(0),
    model(na, nb, nd),
    eps(0.0)
{
    phi = VectorXd::Zero(nCoefs, 1);
    theta0 = VectorXd::Zero(nCoefs, 1);
    theta1 = VectorXd::Ones(nCoefs, 1);
    theta2 = VectorXd::Zero(nCoefs, 1);
    err = VectorXd::Zero(nCoefs, 1);
    y = VectorXd::Zero(na, 1);
    u = VectorXd::Zero(nb + nd + 1, 1);
    coefA = VectorXd::Zero(na + 1, 1);
    coefB = VectorXd::Zero(nb + 1, 1);
    p1 = MatrixXd::Identity(nCoefs, nCoefs) * 1e3;
    p0 = MatrixXd::Zero(nCoefs, nCoefs);
}

bool RLS::update(double input, double output)
{
    double ny;
    ++nStep;
    if (nStep <= nb + nd)
    {
        goto exit;
    }
    for (auto i = 0u; i < na; ++i)
    {
        phi(i) = -y(i);
    }
    for (auto i = 0u; i <= nb; ++i)
    {
        phi(na + i) = u(nd + i - 1);
    }
    p0 = (miu * (1 - lambda) * MatrixXd::Identity(nCoefs, nCoefs) + lambda * p1.inverse() + phi * phi.transpose()).inverse();
    theta0 = theta1 + lambda * miu * p0 * (theta1 - theta2) +
        p0 * phi * (output - phi.transpose() * theta1);
    err = theta1 - theta0;
    for (auto i = 0; i < err.rows(); ++i)
    {
        for (auto j = 0; j < err.cols(); ++j)
        {
            err(i, j) /= theta0(i, j);
        }
    }
    theta2 << theta1;
    theta1 << theta0;
    p1 << p0;
    coefA << 1, theta0.head(na);
    coefB << theta0.tail(nb + 1);
    model.setCoefA(coefA);
    model.setCoefB(coefB);
exit:
    VectorXd uh(u.head(u.size() - 1));
    u << input, uh;
    VectorXd yh(y.head(y.size() - 1));
    y << output, yh;
    return err.norm() < eps;
}

void RLS::setForgetFactor(double l)
{

    lambda = l;
}

void RLS::setDampingFactor(double m)
{
    miu = m;
}

void RLS::setStopError(double err)
{
    eps = err;
}

void RLS::setInitCoefs(const Eigen::VectorXd& coefs)
{
    if (!nStep && coefs.size() == theta1.size())
    {
        theta1 = coefs;
    }
}
