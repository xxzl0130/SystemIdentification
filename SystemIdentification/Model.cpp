#include "Model.h"
#include <iostream>
using namespace std;

ArxModel::ArxModel():nd(1),y(0)
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
    cout << theta << endl;
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

