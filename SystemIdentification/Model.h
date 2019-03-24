#pragma once
#include <Eigen\Dense>
#include <vector>
class Model
{
public:
    Model();
    Model(unsigned na, unsigned nb = 0, unsigned nc = 0);
    virtual ~Model() = default;

    // coefficients
    Eigen::VectorXd coefA, coefB, coefC;
    // pure delay
    unsigned nd;
    double y;   //output
    // update the model with input and noise.
    virtual double update(double input, double noise);

protected:
    // record inputs, outputs and noises
    std::vector<double> inputs, outputs, noises;
    // records max size
    static constexpr unsigned maxLogSize = 1024;
};

