#pragma once
#include <Eigen\Dense>
#include <vector>

// ARX Model in discrete domain (z domain)
// Ay = Bu*z^-d
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

// Transfer function model in continuous domain (s domain).
class SModel
{
public:
    Eigen::VectorXd num; // Numerator
    Eigen::VectorXd den; // Denominator
};

enum DiscretizationMethod
{
    Zoh,
    Foh,
    Tustin,
    ZeroPole,
};

/**
 * \brief convert ArxModel to SModel
 * \param dModel discrete domain model
 * \param Ts sample time
 * \param method conversion method
 * \return continuous model
 */
SModel d2c(const ArxModel& dModel, double Ts = 0.01, DiscretizationMethod method = Zoh);
/**
 * \brief convert SModel to ArxModel
 * \param sModel continuous model
 * \param Ts sample time
 * \param method conversion method
 * \return discrete model
 */
ArxModel c2d(const SModel& sModel, double Ts = 0.01, DiscretizationMethod method = Zoh);