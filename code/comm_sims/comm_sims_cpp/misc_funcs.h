#pragma once

#include <vector>
#include <iostream>
#include <numeric>
#include <algorithm>
#include "Eigen/Dense"

// For sorting arrays
template <typename T>
std::vector<size_t> sort_indexes(const std::vector<T> &v);
// For bisection-solving
template <typename X, typename Y, typename L>
X findRootByRegulaFalsi(
    int nIterations,
    X initLBound,
    X initRBound,
    Y targetVal,
    L&& f
);
// Polysection version
template <typename X, typename Y, typename L>
X findRootByBisection(
    int nIterations,
    X initLBound,
    X initRBound,
    Y targetVal,
    L&& f
);
// Generate autocorrelated normal distribution sample
// https://en.wikipedia.org/wiki/Autoregressive_model#Example:_An_AR(1)_process
Eigen::VectorXf generateAutoCorrSample(int n, float r, float mu, float sigma, int seed);

#include "misc_funcs.tpp"