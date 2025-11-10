#pragma once

#include <iostream>
#include <vector>
#include <numeric>
#include <iterator>
#include <algorithm>
#include <optional>
#include "Eigen/Dense"

#include "misc_funcs.h"

namespace Demo {
    // ===== VitalRates =====

    // Update MPM based on vital rates
    // NB: progression & retrogression shortly before census,
    // reproduction shortly after census

    void VitalRates::updateMPM() {
        Eigen::MatrixXf zeroMatrix = Eigen::MatrixXf::Zero(numStages, numStages);
        Eigen::VectorXf oneCol = Eigen::VectorXf::Ones(numStages);

        // Survival matrix
        Eigen::MatrixXf survMatrix = oneCol * s.transpose();

        // Stationary matrix
        Eigen::MatrixXf nonProgMatrix = (1 - g.array()).matrix().reshaped().asDiagonal();
        Eigen::MatrixXf nonRetroMatrix = (1 - r.array()).matrix().reshaped().asDiagonal();
        Eigen::MatrixXf statMatrix = (nonProgMatrix.array() * nonRetroMatrix.array()).matrix();

        // Progression & retrogression matrix
        Eigen::MatrixXf nonStatMatrix = zeroMatrix;
        // Progression in first subdiagonal
        nonStatMatrix.diagonal(-1) = g(Eigen::seq(0, Eigen::last - 1));
        // Retrogression in first superdiagonal
        nonStatMatrix.diagonal(1) = r(Eigen::seq(0, Eigen::last - 1));
        
        // Reproductive matrix
        Eigen::MatrixXf fertMatrix(numStages, numStages);
        fertMatrix << f.transpose(), Eigen::MatrixXf(numStages - 1, numStages);

        /* std::cout << survMatrix << std::endl <<
            nonProgMatrix << std::endl <<
            nonRetroMatrix << std::endl <<
            nonStatMatrix << std::endl; */

        // Merge all matrices together & save
        mpm = (survMatrix.array() * (statMatrix + nonStatMatrix).array() +
            fertMatrix.array()).matrix();
    };

    void VitalRates::updateMPMAfterCheck() {
        // Update MPM if needed
        if(doesMatrixNeedUpdate) {
            updateMPM();
            doesMatrixNeedUpdate = false;
        }
    }

    // Calculate lambda & equilibrium pop structure
    void VitalRates::updateEigen() {
        // Use a getter to invoke update if needed
        Eigen::MatrixXf mpm_ = getMPM();

        Eigen::EigenSolver<Eigen::MatrixXf> eigenSolver(mpm_);

        if(eigenSolver.info() != Eigen::Success) {
            throw std::runtime_error("EigenSolver failed to converge");
        }

        // Store
        eigenvalues = eigenSolver.eigenvalues();
        eigenvectors = eigenSolver.eigenvectors();

        // Calculate size order of eigens
        Eigen::MatrixXf eigenValueMags = abs(eigenvalues.array()).matrix();
        std::vector<float> eigenvalueMagsCast(
            eigenValueMags.data(),
            eigenValueMags.data() + eigenValueMags.rows() * eigenValueMags.cols()
        );

        eigenOrder = sort_indexes(eigenvalueMagsCast);
        std::reverse(std::begin(eigenOrder), std::end(eigenOrder));

        // Set demo variables as indet
        lambda = std::nullopt;
        eqPopStruct = std::nullopt;
        dampingRatio = std::nullopt;
    }

    void VitalRates::updateEigenAfterCheck() {
        // Update eigen if needed
        if(doesEigenNeedUpdate) {
            updateEigen();
            doesEigenNeedUpdate = false;
        }
    }

    float VitalRates::calcLambda() {
        updateEigenAfterCheck();

        if(lambda) {
            return lambda.value();
        }

        return eigenvalues[eigenOrder[0]].real();
    }

    Eigen::VectorXf VitalRates::calcEqPopStruct() {
        updateEigenAfterCheck();

        if(eqPopStruct) {
            return eqPopStruct.value();
        }

        Eigen::VectorXf eqPopStructNonNorm = eigenvectors.col(eigenOrder[0]).real();

        return eqPopStructNonNorm / eqPopStructNonNorm.sum();
    }
    
    float VitalRates::calcDampingRatio() {
        updateEigenAfterCheck();

        if(dampingRatio) {
            return dampingRatio.value();
        }

        return eigenvalues[eigenOrder[0]].real() / abs(eigenvalues[eigenOrder[1]]);
    }

    // Function for calculating the scaling factor for some vital rate
    // such that lambda is a given value
    // NB: This function assumes that lambda changes monotonically with the vital rate in question.
    float VitalRates::findVRScaleAtLambda(
        float targetLambda,
        VRType adjustedVRType,
        std::optional<int> nIterations,
        std::optional<float> initLScale,
        std::optional<float> initRScale,
        std::optional<float> tolerance
    ) {
        // Run bisection root finder
        float targetScale = findRootByRegulaFalsi(
            nIterations.value_or(20),
            initLScale.value_or(0.0f),
            initRScale.value_or(2.0f),
            targetLambda,
            tolerance.value_or(0.0001f),
            [&] (float scale) {
                VitalRates copyVitalRates = VitalRates(*this);
                Eigen::VectorXf adjustedVRs = getVR(adjustedVRType) * scale;
                if(adjustedVRType != F) {
                    adjustedVRs = adjustedVRs.array().min(Eigen::VectorXf::Ones(numStages).array()).matrix();
                }
                copyVitalRates.setVR(adjustedVRType, adjustedVRs);
                return copyVitalRates.calcLambda();
            }
        );

        /* float targetScale = findRootByITP(
            targetLambda,
            initLScale.value_or(0.0f),
            initRScale.value_or(2.0f),
            0.0001f,
            0.001f * abs(initRScale.value_or(2.0f) - initLScale.value_or(0.0f)),
            2.0f,
            2,
            [&] (float scale) {
                VitalRates copyVitalRates = VitalRates(*this);
                copyVitalRates.setF(f * scale);
                return copyVitalRates.calcLambda();
            }
        ); */

        // Return target scale
        return std::max(0.0f, targetScale);
    }

    void VitalRates::adjustVRToLambda(
        float targetLambda,
        VRType adjustedVRType,
        std::optional<int> nIterations,
        std::optional<float> initLScale,
        std::optional<float> initRScale,
        std::optional<float> tolerance
    ) {
        // Find target scale
        float targetScale = findVRScaleAtLambda(
            targetLambda, adjustedVRType, nIterations, initLScale, initRScale, tolerance
        );

        // Rescale VR
        Eigen::VectorXf adjustedVRs = getVR(adjustedVRType) * targetScale;
        if(adjustedVRType != F) {
            adjustedVRs = adjustedVRs.array().min(Eigen::VectorXf::Ones(numStages).array()).matrix();
        }
        setVR(adjustedVRType, adjustedVRs.eval());

        // std::cout << getVR(adjustedVRType).transpose() << std::endl;
    }

    // The above, but with constraints
    // The constraint vector defines the multiplicative bias to be applied to each life stage
    // in terms of dynamic adjustment. e.g. (1.0, 2.0) means that the adult vital rate is scaled by
    // a factor 2 times greater in magnitude than juveniles.
    float VitalRates::findVRScaleAtLambdaWithConstraints(
        float targetLambda,
        VRType adjustedVRType,
        Eigen::VectorXf constraints,
        std::optional<int> nIterations,
        std::optional<float> initLScale,
        std::optional<float> initRScale,
        std::optional<float> tolerance
    ) {
        // Run bisection root finder
        float targetScale = findRootByRegulaFalsi(
            nIterations.value_or(20),
            initLScale.value_or(0.0f),
            initRScale.value_or(2.0f),
            targetLambda,
            tolerance.value_or(0.0001f),
            [&] (float scale) {
                VitalRates copyVitalRates = VitalRates(*this);
                Eigen::VectorXf adjustedVRs = (
                    (1 - (1 - scale) * constraints.array())
                    * getVR(adjustedVRType).array()
                ).matrix();
                if(adjustedVRType != F) { // Vital rates except fertility has an upper bound of 1
                    adjustedVRs = adjustedVRs.array().min(Eigen::VectorXf::Ones(numStages).array()).matrix();
                }
                copyVitalRates.setVR(adjustedVRType, adjustedVRs);
                return copyVitalRates.calcLambda();
            }
        );

        /* float targetScale = findRootByITP(
            targetLambda,
            initLScale.value_or(0.0f),
            initRScale.value_or(2.0f),
            0.0001f,
            0.001f * abs(initRScale.value_or(2.0f) - initLScale.value_or(0.0f)),
            2.0f,
            2,
            [&] (float scale) {
                VitalRates copyVitalRates = VitalRates(*this);
                copyVitalRates.setF(f * scale);
                return copyVitalRates.calcLambda();
            }
        ); */

        // Return target scale, with a lower bound of 0
        return std::max(0.0f, targetScale);
    }

    void VitalRates::adjustVRToLambdaWithConstraints(
        float targetLambda,
        VRType adjustedVRType,
        Eigen::VectorXf constraints,
        std::optional<int> nIterations,
        std::optional<float> initLScale,
        std::optional<float> initRScale,
        std::optional<float> tolerance
    ) {
        // Find target scale
        float targetScale = findVRScaleAtLambdaWithConstraints(
            targetLambda, adjustedVRType, constraints, nIterations, initLScale, initRScale, tolerance
        );

        Eigen::VectorXf adjustedVRs = (
            (1 - (1 - targetScale) * constraints.array())
            * getVR(adjustedVRType).array()
        ).matrix();
        if(adjustedVRType != F) { // Vital rates except fertility has an upper bound of 1
            adjustedVRs = adjustedVRs.array().min(Eigen::VectorXf::Ones(numStages).array()).matrix();
        }

        // Rescale VR
        setVR(adjustedVRType, adjustedVRs.eval());

        // std::cout << getVR(adjustedVRType).transpose() << std::endl;
    }
}


