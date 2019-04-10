#pragma once
#include <Eigen\Dense>
#include <vector>
class ArxModel
{
public:
    ArxModel();
    ArxModel(unsigned na, unsigned nb = 0, unsigned nd = 1);
    virtual ~ArxModel() = default;

    // coefficients
    Eigen::VectorXd coefA, coefB;
    // pure delay
    unsigned nd;
    double y;   //output
    // update the model with input and noise.
    virtual double update(double input, double noise = 0);

protected:
    // record inputs, outputs
    std::vector<double> inputs, outputs, noises;
};
