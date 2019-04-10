#include "Model.h"

ArxModel::ArxModel():nd(1),y(0)
{
}

ArxModel::ArxModel(unsigned na, unsigned nb, unsigned nd) : nd(nd), y(0)
{
    coefA.resize(na + 1);
    coefB.resize(nb + 1);
}

//TODO:test
double ArxModel::update(double input, double noise)
{
    inputs.push_back(input);
    noises.push_back(noise);
    auto a = 0.0, b = 0.0;
    unsigned i;
    auto it = inputs.rbegin();
    for(i = 1;i < coefA.size() && it != inputs.rend();++i,++it)
    {
        a += coefA(i) * *it;
    }
    for(i = 0,it = outputs.rbegin() + nd;i < coefB.size() && it != outputs.rend();++i,++it)
    {
        b += coefB(i) * *it;
    }
    y = b - a;
    outputs.push_back(y);
    return y;
}

