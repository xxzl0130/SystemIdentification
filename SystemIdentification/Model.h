#pragma once
#include <Eigen\Dense>
#include <vector>
class ArxModel
{
public:
    ArxModel();
    ArxModel(unsigned na, unsigned nb = 0, unsigned nd = 1);
    virtual ~ArxModel() = default;

    void setNa(unsigned a);
    void setNb(unsigned b);
    void setNd(unsigned d);
    void setCoefA(const Eigen::VectorXd &A);
    Eigen::VectorXd getCoefA() const
    {
        return coefA;
    }
    void setCoefB(const Eigen::VectorXd &B);
    Eigen::VectorXd getCoefB() const
    {
        return coefB;
    }

    double y;   //output
    // update the model
    virtual double update(double input);

protected:
    // record inputs, outputs
    std::vector<double> inputs, outputs, noises;
    Eigen::VectorXd uk, yk;
    unsigned na, nb, nd;
    // coefficients
    Eigen::VectorXd coefA, coefB;
    Eigen::VectorXd theta;
};
