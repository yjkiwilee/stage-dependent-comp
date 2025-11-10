#pragma once

#include <iostream>
#include <vector>
#include <numeric>
#include <algorithm>

template <typename T>

// https://stackoverflow.com/questions/targetVal577475/c-sorting-and-keeping-track-of-indexes
std::vector<size_t> sort_indexes(const std::vector<T> &v) {

  // initialize original index locations
  std::vector<size_t> idx(v.size());
  std::iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  // using std::stable_sort instead of std::sort
  // to avoid unnecessary index re-orderings
  // when v contains elements of equal values 
  std::stable_sort(idx.begin(), idx.end(),
       [&v](size_t itargetVal, size_t i2) {return v[itargetVal] < v[i2];});

  return idx;
}

// https://stackoverflow.com/questions/1903954/is-there-a-standard-sign-function-signum-sgn-in-c-c
template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

// For bisection-solving
template <typename X, typename Y, typename L>
X findRootByRegulaFalsi(
    int nIterations,
    X initLBound,
    X initRBound,
    Y targetVal,
    float tolerance,
    L&& f
) {
    // Bisection method, assuming monotonic increase in Val with fertility scaling
    X lBound = initLBound, rBound = initRBound, mBound;
    Y lVal = f(lBound), rVal = f(rBound), mVal;

    // Iterate
    int i;
    for(i = 0; i < nIterations; i++) {
        // Swap vals as needed
        if(lVal > rVal) {
            std::swap(lVal, rVal);
            std::swap(lBound, rBound);
        }

        // Calculate mBound
        mBound = lBound + (targetVal - lVal) / (rVal - lVal) * (rBound - lBound);
        mVal = f(mBound);

        // If value within tolerance range, break
        if(abs(mVal - targetVal) / abs(targetVal) < tolerance) {
            break;
        }
        
        if(lVal < targetVal && targetVal < rVal) {
            // If Val = targetVal lies between bounds
            if(mVal < targetVal) {
                // If Val = targetVal is in the right section
                // 'zoom in' on right section
                lBound = mBound;
                lVal = mVal;
            } else if (mVal > targetVal) {
                // If Val = targetVal is in the left section
                // 'zoom in' on left section
                rBound = mBound;
                rVal = mVal;
            } else { // Val = targetVal at mBound; unlikely edge case
                return mBound;
            }
        } else if(lVal > targetVal) {
            // If Val = targetVal is left of bound
            // Double range towards left
            lBound -= rBound - lBound;
            lVal = f(lBound);
        } else if(rVal < targetVal) {
            // If Val = targetVal is right of bound
            // Double range towards right
            rBound += rBound - lBound;
            rVal = f(rBound);
        }
        else if(lVal == targetVal) { return lBound; } // Unlikely edge cases
        else if(rVal == targetVal) { return rBound; } // "
    }

    // std::cout << i << std::endl;

    // Return mBound
    return mBound;
};

// For ITP-solving (https://en.wikipedia.org/wiki/ITP_method)
// https://math.stackexchange.com/questions/1464795/what-are-the-difference-between-some-basic-numerical-root-finding-methods
template <typename X, typename Y, typename L>
X findRootByITP(
    Y ytarg, // target
    X initA,
    X initB,
    X epsilon,
    X k1,
    X k2,
    int n0,
    L&& f
) {
    X a = initA, b = initB;
    Y ya = f(a) - ytarg, yb = f(b) - ytarg;

    // Swap vals as needed
    if(ya > yb) {
        std::swap(ya, yb);
        std::swap(a, b);
    }

    // Initially, find straddling range
    while (ya > 0 || yb < 0) {
        if(ya > 0) {
            // If Val = targetVal is left of bound
            // Double range towards left
            a -= b - a;
            ya = f(a) - ytarg;
        } else if(yb < 0) {
            // If Val = targetVal is right of bound
            // Double range towards right
            b += b - a;
            yb = f(b) - ytarg;
        }
    }

    // Preprocessing
    X range = b - a;
    int nhalf = ceil(log2(range / (2 * epsilon)));
    int nmax = nhalf + n0;
    int j = 0;

    X xhalf, r, delta, xf, xt, xitp;
    Y yitp;
    int sigma;

    while(range > 2 * epsilon) {
        // Calculating parameters
        xhalf = (a + b) / 2;
        r = epsilon * pow(2, nmax - j) - range / 2;
        delta = k1 * pow(range, k2);

        // Interpolation
        xf = (yb * a - ya * b) / (yb - ya);

        // Truncation
        sigma = sgn(xhalf - xf);

        if(delta <= abs(xhalf - xf)) {
            xt = xf + sigma * delta;
        } else {
            xt = xhalf;
        }

        // Projection
        if(abs(xt - xhalf) <= r) {
            xitp = xt;
        } else {
            xitp = xhalf - sigma * r;
        }

        // Updating interval
        yitp = f(xitp) - ytarg;
        if(yitp > 0) {
            b = xitp;
            yb = yitp;
        } else if(yitp < 0) {
            a = xitp;
            ya = yitp;
        } else {
            a = xitp;
            b = xitp;
        }

        j++;

        range = b - a;
    }

    return (a + b) / 2;
};

// Generate autocorrelated normal distribution sample
// https://en.wikipedia.org/wiki/Autoregressive_model#Example:_An_AR(1)_process
Eigen::VectorXf generateAutoCorrSample(int n, float r, float mu = 0, float sigma = 1, int seed = 1) {
    // Calculate standard deviation of epsilon
    float sigmaE = sqrt(sigma * sigma * (1 - r * r));
    
    // Seed randomengine
    std::default_random_engine randEngine(seed);
    std::normal_distribution<float> xDist(0.0, sigma);
    std::normal_distribution<float> epsilonDist(0.0, sigmaE);

    // Choose initial value
    Eigen::VectorXf xVector = Eigen::VectorXf::Zero(n);
    xVector.head(1) << xDist(randEngine);

    // Populate epsilon
    auto randGen = [&]() {
        return epsilonDist(randEngine);
    };
    Eigen::VectorXf epsilonVector(n - 1);
    std::generate(epsilonVector.begin(), epsilonVector.end(), randGen);

    // Populate X
    for(int i = 0; i < n - 1; i++) {
        xVector[i + 1] = r * xVector[i] + epsilonVector[i];
    }

    // Return X
    return xVector;
}