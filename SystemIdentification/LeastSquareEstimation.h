#pragma once
#include "Model.h"

ArxModel leastSquare(const std::vector<double>& inputs, const std::vector<double>& outputs, unsigned na, unsigned nb, unsigned nd = 1);
ArxModel leastSquare(const double* inputs, const double* outputs, const unsigned size, unsigned na, unsigned nb, unsigned nd = 1);