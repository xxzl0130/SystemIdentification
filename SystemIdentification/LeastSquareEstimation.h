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
