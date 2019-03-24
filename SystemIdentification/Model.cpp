#include "Model.h"

Model::Model():nd(0),y(0)
{
}

Model::Model(unsigned na, unsigned nb, unsigned nc) : nd(0), y(0)
{
    coefA.resize(na);
    coefB.resize(nb);
    coefC.resize(nc);
}

double Model::update(double input, double noise)
{
    inputs.push_back(input);
    noises.push_back(noise);
    auto a = 0.0, b = 0.0, c = 0.0;
    unsigned i;
    auto it = inputs.rbegin();
    for(i = 1;i < coefA.size() && it != inputs.rend();++i,++it)
    {
        a += coefA(i) * *it;
    }
    for(i = 0,it = outputs.rbegin();i < coefB.size() && it != outputs.rend();++i,++it)
    {
        b += coefB(i) * *it;
    }
    for (i = 0, it = noises.rbegin(); i < coefC.size() && it != noises.rend(); ++i, ++it)
    {
        c += coefC(i) * *it;
    }
    y = b + c - a;
    outputs.push_back(y);
    return y;
}

