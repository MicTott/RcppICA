#ifndef RCPPICA_H
#define RCPPICA_H

// Eigen optimizations - must be defined before including Eigen
#ifndef EIGEN_NO_DEBUG
#define EIGEN_NO_DEBUG
#endif

#ifndef EIGEN_INITIALIZE_MATRICES_BY_ZERO
#define EIGEN_INITIALIZE_MATRICES_BY_ZERO
#endif

// Core includes
#include <RcppEigen.h>

// OpenMP support (optional)
#ifdef _OPENMP
#include <omp.h>
#endif

// Standard library
#include <cmath>
#include <algorithm>
#include <random>
#include <memory>

// RcppICA components
#include "RcppICA/Nonlinearities.hpp"
#include "RcppICA/Whitening.hpp"
#include "RcppICA/ICADeflation.hpp"
#include "RcppICA/ICAParallel.hpp"

#endif // RCPPICA_H
