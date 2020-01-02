#include "Model.h"
#include <iostream>
#include <vector>
#include <Polynomial/Polynomial.hpp>
using namespace std;
typedef polynomial::Polynomial<Eigen::Dynamic> Polynomial;

ArxModel::ArxModel() :na(0), nb(0), nd(1), y(0)
{
}

ArxModel::ArxModel(unsigned na, unsigned nb, unsigned nd) : na(na), nb(nb), nd(nd), y(0)
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

    phi.resize(yk.size() + nb + 1, 1);
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
    switch (method)
    {
    case Zoh:
    {
        //z = 1/(1 - s*Ts)
        Polynomial zoh(Eigen::Vector2d(-Ts, 1));
        int order = std::max(dModel.getCoefA().size(), dModel.getCoefB().size());
        auto t = Polynomial(Eigen::VectorXd::Identity(1, 1));    //constant 1
        auto den = Polynomial(Eigen::Vector2d(0, 0));
        for (auto i = 0; i < order - dModel.getCoefA().size(); ++i)
        {
            t = t * zoh;
        }
        for (auto i = 0; i < dModel.getCoefA().size(); ++i)
        {
            den = den + t * dModel.getCoefA()(i);
            t = t * zoh;
        }
        t = Polynomial(Eigen::VectorXd::Identity(1, 1));
        auto num = Polynomial(Eigen::Vector2d(0, 0));
        for (auto i = 0; i < order - dModel.getCoefB().size(); ++i)
        {
            t = t * zoh;
        }
        for (auto i = 0; i < dModel.getCoefB().size(); ++i)
        {
            num = num + t * dModel.getCoefB()(i);
            t = t * zoh;
        }
        SModel sModel;
        sModel.num = num.coefficients() / den.coefficients()(0);
        sModel.den = den.coefficients() / den.coefficients()(0);
        return sModel;
    }
    case Tustin:
        //z = (1+s*Ts/2)/(1-s*Ts/2)
    {
        Polynomial denP0(Eigen::Vector2d(-Ts / 2, 1));     //(1-s*Ts/2)
        Polynomial numP0(Eigen::Vector2d(Ts / 2, 1));     //(1+s*Ts/2)
        Polynomial den(Eigen::VectorXd::Zero(1, 1));
        Polynomial num(Eigen::VectorXd::Zero(1, 1));
        int order = std::max(dModel.getCoefA().size(), dModel.getCoefB().size());
        auto tmpDen = Polynomial(Eigen::VectorXd::Identity(1, 1));    //constant 1
        auto tmpNum = tmpDen;    //constant 1
        std::vector<Polynomial> dens, nums;
        dens.resize(dModel.getCoefA().size(), Polynomial(Eigen::Vector2d(0, 0)));
        nums.resize(dModel.getCoefA().size(), Polynomial(Eigen::Vector2d(0, 0)));
        for (auto i = 0; i < order - dModel.getCoefA().size(); ++i)
        {
            tmpDen = tmpDen * denP0;
        }
        for (auto i = 0; i < dModel.getCoefA().size(); ++i)
        {
            dens[i] = tmpDen * dModel.getCoefA()(i);
            tmpDen = tmpDen * denP0;
        }
        for (auto i = dModel.getCoefA().size() - 1; i >= 0; --i)
        {
            den = den + dens[i] * tmpNum;
            tmpNum = tmpNum * numP0;
        }

        tmpDen = Polynomial(Eigen::VectorXd::Identity(1, 1));    //constant 1
        tmpNum = tmpDen;    //constant 1
        for (auto i = 0; i < order - dModel.getCoefB().size(); ++i)
        {
            tmpDen = tmpDen * denP0;
        }
        for (auto i = 0; i < dModel.getCoefB().size(); ++i)
        {
            nums[i] = tmpDen * dModel.getCoefB()(i);
            tmpDen = tmpDen * denP0;
        }
        for (auto i = dModel.getCoefB().size() - 1; i >= 0; --i)
        {
            num = num + nums[i] * tmpNum;
            tmpNum = tmpNum * numP0;
        }
        SModel sModel;
        sModel.num = num.coefficients() / den.coefficients()(0);
        sModel.den = den.coefficients() / den.coefficients()(0);
        return sModel;
    }
    default:
        break;
    }
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
        for (auto i = sModel.num.size() - 2; i >= 0; --i)
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
        ArxModel zModel(den.coefficients().size() - 1, den.coefficients().size() - 1, 1);

        auto A = den.coefficients();
        //reverse
        for (auto i = 0; i < A.size() / 2; ++i)
        {
            std::swap(A(i), A(A.size() - 1 - i));
        }
        // set coefficients and normalize
        zModel.setCoefA(A / A(0));

        Eigen::VectorXd B = Eigen::VectorXd::Zero(A.size());
        //reverse and *z^n
        for (auto i = 0; i < num.coefficients().size(); ++i)
        {
            B(i) = num.coefficients()(num.coefficients().size() - 1 - i);
        }
        // set coefficients and normalize
        zModel.setCoefB(B / A(0));
        return zModel;
    }
    case Tustin:
        //s = (2/Ts) * (z - 1) / (z + 1)
    {
        Polynomial denP0(Eigen::Vector2d(1, 1));        //(z + 1)
        Polynomial numP0 = Polynomial(Eigen::Vector2d(1, -1)) * (2 / Ts);       //(2/Ts) * (z - 1)
        Polynomial den(Eigen::VectorXd::Zero(1, 1));
        Polynomial num(Eigen::VectorXd::Zero(1, 1));
        int order = std::max(sModel.num.size(), sModel.den.size());
        auto tmpDen = Polynomial(Eigen::VectorXd::Identity(1, 1));    //constant 1
        auto tmpNum = tmpDen;    //constant 1
        std::vector<Polynomial> dens, nums;
        dens.resize(sModel.den.size(), Polynomial(Eigen::Vector2d(0, 0)));
        nums.resize(sModel.num.size(), Polynomial(Eigen::Vector2d(0, 0)));
        for (auto i = 0; i < order - sModel.den.size(); ++i)
        {
            tmpDen = tmpDen * denP0;
        }
        for (auto i = 0; i < sModel.den.size(); ++i)
        {
            dens[i] = tmpDen * sModel.den(i);
            tmpDen = tmpDen * denP0;
        }
        for (auto i = sModel.den.size() - 1; i >= 0; --i)
        {
            den = den + dens[i] * tmpNum;
            tmpNum = tmpNum * numP0;
        }

        tmpDen = Polynomial(Eigen::VectorXd::Identity(1, 1));    //constant 1
        tmpNum = tmpDen;    //constant 1
        for (auto i = 0; i < order - sModel.num.size(); ++i)
        {
            tmpDen = tmpDen * denP0;
        }
        for (auto i = 0; i < sModel.num.size(); ++i)
        {
            nums[i] = tmpDen * sModel.num(i);
            tmpDen = tmpDen * denP0;
        }
        for (auto i = sModel.num.size() - 1; i >= 0; --i)
        {
            num = num + nums[i] * tmpNum;
            tmpNum = tmpNum * numP0;
        }
        ArxModel dModel(den.coefficients().size() - 1, den.coefficients().size() - 1, 1);
        dModel.setCoefB(num.coefficients() / den.coefficients()(0));
        dModel.setCoefA(den.coefficients() / den.coefficients()(0));
        return dModel;
    }
    default:
        break;
    }

    return ArxModel();
}
