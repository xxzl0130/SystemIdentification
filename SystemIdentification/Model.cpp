#include "Model.h"
#include <iostream>
#include <Polynomial/Polynomial.hpp>
using namespace std;
typedef polynomial::Polynomial<Eigen::Dynamic> Polynomial;

ArxModel::ArxModel():na(0),nb(0),nd(1),y(0)
{
}

ArxModel::ArxModel(unsigned na, unsigned nb, unsigned nd) : na(na),nb(nb),nd(nd), y(0)
{
    coefA.resize(na + 1);
    coefB.resize(nb + 1);
    uk = Eigen::VectorXd::Zero(nd + nb, 1);
    yk = Eigen::VectorXd::Zero(na, 1);
    theta.resize(na + nb + 1);
    theta << coefA.tail(coefA.size() - 1), coefB;
}

void ArxModel::setNa(unsigned a)
{
    coefA.resize(a);
    na = a;
    yk = Eigen::VectorXd::Zero(na, 1);
    theta.resize(na + nb + 1);
    theta << coefA.tail(coefA.size() - 1), coefB;
}

void ArxModel::setNb(unsigned b)
{
    coefB.resize(b);
    nb = b;
    uk = Eigen::VectorXd::Zero(nd + nb, 1);
    theta.resize(na + nb + 1);
    theta << coefA.tail(coefA.size() - 1), coefB;
}

void ArxModel::setNd(unsigned d)
{
    nd = d;
    uk = Eigen::VectorXd::Zero(nd + nb, 1);
}

void ArxModel::setCoefA(const Eigen::VectorXd& A)
{
    if (A.size() - 1 != na)
        return;
    coefA << A;
    theta.resize(na + nb + 1);
    theta << coefA.tail(coefA.size() - 1), coefB;
}

void ArxModel::setCoefB(const Eigen::VectorXd& B)
{
    if (B.size() - 1 != nb)
        return;
    coefB << B;
    theta.resize(na + nb + 1);
    theta << coefA.tail(coefA.size() - 1), coefB;
}

double ArxModel::update(double input)
{
    inputs.push_back(input);
    Eigen::VectorXd phi;

    phi.resize(yk.size() + nb + 1,1);
    phi << -yk, uk.middleRows(nd - 1, nb + 1);
    y = phi.dot(theta);
    outputs.push_back(y);
    
    //sequence shift
    Eigen::VectorXd uh(uk.head(uk.size() - 1));
    uk << input, uh;
    Eigen::VectorXd yh(yk.head(yk.size() - 1));
    yk << y, yh;

    return y;
}

SModel d2c(const ArxModel& dModel, double Ts, DiscretizationMethod method)
{
    return SModel();
}

ArxModel c2d(const SModel& sModel, double Ts, DiscretizationMethod method)
{
    switch (method)
    {
    case Zoh:
        {
            // s = (-z^-1 + 1) / Ts
            Polynomial zoh(Eigen::Vector2d(-1, 1) / Ts);
            //only preset constant
            Polynomial num(Eigen::Vector2d(0, sModel.num(sModel.num.size() - 1)));
            auto t = zoh;
            for(auto i = sModel.num.size() - 2;i >= 0;--i)
            {
                num = num + t * sModel.num(i);
                t = t * zoh;//plus one order
            }
            //only preset constant
            Polynomial den(Eigen::Vector2d(0, sModel.den(sModel.den.size() - 1)));
            t = zoh;
            for (auto i = sModel.den.size() - 2; i >= 0; --i)
            {
                den = den + t * sModel.den(i);
                t = t * zoh;//plus one order
            }
            ArxModel zModel(den.coefficients().size() - 1,  den.coefficients().size() - 1, 1);
            
            auto A = den.coefficients();
            //reverse
            for(auto i = 0;i < A.size() / 2;++i)
            {
                std::swap(A(i), A(A.size() - 1 - i));
            }
            // set coefficients and normalize
            zModel.setCoefA(A / A(0));

            Eigen::VectorXd B = Eigen::VectorXd::Zero(A.size());
            //reverse and *z^n
            for(auto i = 0;i < num.coefficients().size();++i)
            {
                B(i) = num.coefficients()(num.coefficients().size() - 1 - i);
            }
            // set coefficients and normalize
            zModel.setCoefB(B / A(0));
            return zModel;
        }
    default:
        break;
    }

    return ArxModel();
}

