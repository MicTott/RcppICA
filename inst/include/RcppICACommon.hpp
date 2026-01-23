// RcppICACommon.hpp - Common header for RcppICA package
//
// Copyright (C) 2026 Michael Totty
//
// This file is part of RcppICA
//
// RcppICA is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// RcppICA is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

#ifndef RCPPICA_COMMON_HPP
#define RCPPICA_COMMON_HPP

// OpenMP support (conditional compilation)
//[[Rcpp::plugins(openmp)]]
#ifdef _OPENMP
#include <omp.h>
#endif

// Eigen configuration - must be defined before including Eigen
#ifndef EIGEN_NO_DEBUG
#define EIGEN_NO_DEBUG
#endif

// Core includes
#include <RcppEigen.h>

// Standard library
#include <cmath>
#include <algorithm>
#include <random>
#include <memory>
#include <vector>
#include <string>

// RcppICA algorithm components
#include "RcppICA/Nonlinearities.hpp"
#include "RcppICA/Whitening.hpp"
#include "RcppICA/ICADeflation.hpp"
#include "RcppICA/ICAParallel.hpp"

#endif // RCPPICA_COMMON_HPP
