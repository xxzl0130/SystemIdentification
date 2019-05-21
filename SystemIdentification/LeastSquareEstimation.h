#pragma once
#include "Model.h"

/**
 * \brief estimate arx model (Ay = z^-d*Bu) with least square algorithm.
 * \param inputs target system input data record, or output of test system.
 * \param outputs target system output data record, or input of test system.
 * \param na dimension of denominator, as coefficient B.
 * \param nb dimension of numerator, as coefficient A.
 * \param nd time of pure delay.
 * \return ArxModel
 */
ArxModel leastSquare(const std::vector<double>& inputs, const std::vector<double>& outputs, unsigned na, unsigned nb,
    unsigned nd = 1);
/**
 * \brief estimate arx model (Ay = z^-d*Bu) with least square algorithm. C-style interface
 * \param inputs target system input data record, or output of test system.
 * \param outputs target system output data record, or input of test system.
 * \param size data size.
 * \param na dimension of denominator, as coefficient B.
 * \param nb dimension of numerator, as coefficient A.
 * \param nd time of pure delay.
 * \return ArxModel
 */
ArxModel leastSquare(const double* inputs, const double* outputs, const unsigned size, unsigned na, unsigned nb,
    unsigned nd = 1);

/**
 * \brief Recursion Least Square with forgetting factor and damping factor.
 * RLS has to record inputs and outputs, so make it class instead of function.
 */
class RLS
{
public:
    RLS(unsigned na, unsigned nb, unsigned nd = 1);
    virtual ~RLS() = default;

    /**
     * \brief update RLS with system input and output value.
     * \param input target system input, or output of test system.
     * \param output target system output, or input of test system
     * \return whether the disturbance of coefficients is less than error threshold.
     */
    bool update(double input, double output);
    const ArxModel& getModel() const
    {
        return model;
    }
    void setForgetFactor(double l);
    void setDampingFactor(double m);
    /**
     * \brief set error threshold.
     * \param err error threshold.
     */
    void setStopError(double err);
    // set initial coefficients of model, only available before first update.
    void setInitCoefs(const Eigen::VectorXd& coefs);
private:
    //temp vars in update
    Eigen::VectorXd phi, theta1, theta0, theta2, err, coefA, coefB;
    Eigen::MatrixXd p0, p1;
    unsigned nCoefs, nStep;
    unsigned na, nb, nd;

    //factors
    double lambda, miu;

    //record input and output
    Eigen::VectorXd u, y;

    ArxModel model;
    double eps;
};