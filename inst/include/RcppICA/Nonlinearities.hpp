#ifndef RCPPICA_NONLINEARITIES_HPP
#define RCPPICA_NONLINEARITIES_HPP

#include <Eigen/Dense>
#include <cmath>

namespace RcppICA {

// Nonlinearity type enumeration
enum class NonlinearityType { LOGCOSH = 0, EXP = 1, CUBE = 2 };

// LogCosh nonlinearity: g(u) = tanh(alpha * u), g'(u) = alpha * (1 - tanh^2(alpha * u))
// Most robust choice, good for super- and sub-Gaussian sources
template<typename Scalar = double>
struct LogCosh {
    Scalar alpha;

    explicit LogCosh(Scalar a = 1.0) : alpha(a) {}

    // g(u) = tanh(alpha * u)
    template<typename Derived>
    Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic>
    g(const Eigen::MatrixBase<Derived>& x) const {
        return (alpha * x.array()).tanh();
    }

    // g'(u) = alpha * (1 - tanh^2(alpha * u))
    template<typename Derived>
    Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic>
    gPrime(const Eigen::MatrixBase<Derived>& x) const {
        auto t = (alpha * x.array()).tanh();
        return alpha * (1.0 - t.square());
    }

    // Combined computation for efficiency (avoids computing tanh twice)
    template<typename Derived>
    void compute(const Eigen::MatrixBase<Derived>& x,
                 Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>& gx,
                 Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>& gpx) const {
        auto t = (alpha * x.array()).tanh();
        gx = t.matrix();
        gpx = (alpha * (1.0 - t.square())).matrix();
    }
};

// Exponential nonlinearity: g(u) = u * exp(-u^2/2), g'(u) = (1 - u^2) * exp(-u^2/2)
// Good for super-Gaussian sources with heavy tails
template<typename Scalar = double>
struct Exp {
    // g(u) = u * exp(-u^2 / 2)
    template<typename Derived>
    Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic>
    g(const Eigen::MatrixBase<Derived>& x) const {
        return x.array() * (-0.5 * x.array().square()).exp();
    }

    // g'(u) = (1 - u^2) * exp(-u^2 / 2)
    template<typename Derived>
    Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic>
    gPrime(const Eigen::MatrixBase<Derived>& x) const {
        auto u2 = x.array().square();
        return (1.0 - u2) * (-0.5 * u2).exp();
    }

    // Combined computation
    template<typename Derived>
    void compute(const Eigen::MatrixBase<Derived>& x,
                 Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>& gx,
                 Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>& gpx) const {
        auto u2 = x.array().square();
        auto expTerm = (-0.5 * u2).exp();
        gx = (x.array() * expTerm).matrix();
        gpx = ((1.0 - u2) * expTerm).matrix();
    }
};

// Cubic nonlinearity: g(u) = u^3, g'(u) = 3u^2
// Simplest, equivalent to kurtosis, but less robust to outliers
template<typename Scalar = double>
struct Cube {
    // g(u) = u^3
    template<typename Derived>
    Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic>
    g(const Eigen::MatrixBase<Derived>& x) const {
        return x.array().cube();
    }

    // g'(u) = 3 * u^2
    template<typename Derived>
    Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic>
    gPrime(const Eigen::MatrixBase<Derived>& x) const {
        return 3.0 * x.array().square();
    }

    // Combined computation
    template<typename Derived>
    void compute(const Eigen::MatrixBase<Derived>& x,
                 Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>& gx,
                 Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>& gpx) const {
        auto xarr = x.array();
        gx = xarr.cube().matrix();
        gpx = (3.0 * xarr.square()).matrix();
    }
};

} // namespace RcppICA

#endif // RCPPICA_NONLINEARITIES_HPP
