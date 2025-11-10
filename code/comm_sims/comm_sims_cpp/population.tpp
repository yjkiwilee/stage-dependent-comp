#pragma once

#include <iostream>
#include <vector>
#include <cmath>
#include "Eigen/Dense"

#include "vital_rates.h"

namespace Demo {
    // ===== Population =====

    // Function for calculating population projection
    Eigen::VectorXf Population::calcPopProjection() {
        return vitalRates.getMPM() * popVector;
    }
    Eigen::VectorXf Population::calcPopProjection(Eigen::VectorXf currPopVector) {
        return vitalRates.getMPM() * currPopVector;
    }
    Eigen::VectorXf Population::calcPopProjectionN(int nTimeSteps) {
        Eigen::VectorXf resPopVector = popVector;
        for(int i = 0; i < nTimeSteps; i++) {
            resPopVector = vitalRates.getMPM() * popVector;
        }
        return resPopVector;
    }
    Eigen::VectorXf Population::calcPopProjectionN(Eigen::VectorXf currPopVector, int nTimeSteps) {
        Eigen::VectorXf resPopVector = currPopVector;
        for(int i = 0; i < nTimeSteps; i++) {
            resPopVector = vitalRates.getMPM() * currPopVector;
        }
        return resPopVector;
    }

    // Project population by 1 time step
    void Population::projectPopulation() {
        popVector = calcPopProjection(popVector);
    }
    // Project population by n time steps
    void Population::projectPopulationN(int nTimeSteps) {
        popVector = calcPopProjectionN(popVector, nTimeSteps);
    }

    // Adjust vital rates to given lambda
    void Population::adjustVRToLambda(
        float targetLambda,
        VRType adjustedVRType,
        std::optional<int> nIterations,
        std::optional<float> initLScale,
        std::optional<float> initRScale,
        std::optional<float> tolerance
    ) {
        vitalRates.adjustVRToLambda(
            targetLambda, adjustedVRType, nIterations,
            initLScale, initRScale, tolerance
        );
    }

    // ===== 'Dynamic' Population =====

    // Constructor
    DynamicPopulation::DynamicPopulation(
        int numStages_,
        VitalRates vitalRates_,
        Eigen::VectorXf popVector_,
        VRType dynamicVRType_,
        std::optional<Eigen::VectorXf> constraints_
    ) : Population(numStages_, vitalRates_, popVector_),
        dynamicVRType(dynamicVRType_),
        constraints(constraints_.value_or(Eigen::VectorXf::Ones(numStages_)))
    {
        dynamicVRVector = Eigen::VectorXf(numStages_);
    };

    // Override projection functions
    Eigen::VectorXf DynamicPopulation::calcPopProjection(
        float dynamicFactor
    ) {
        return calcDynamicVitalRates(dynamicFactor).getMPM() * popVector;
    }
    Eigen::VectorXf DynamicPopulation::calcPopProjection(
        Eigen::VectorXf currPopVector,
        float dynamicFactor
    ) {
        return calcDynamicVitalRates(dynamicFactor).getMPM() * currPopVector;
    }
    Eigen::VectorXf DynamicPopulation::calcPopProjectionN(
        int nTimeSteps,
        Eigen::VectorXf dynamicFactorVector // nrow = nTimeSteps
    ) {
        Eigen::VectorXf resPopVector = popVector;
        for(int i = 0; i < nTimeSteps; i++) {
            resPopVector = calcDynamicVitalRates(dynamicFactorVector.row(i).value()).getMPM()
                * popVector;
        }
        return resPopVector;
    }
    Eigen::VectorXf DynamicPopulation::calcPopProjectionN(
        Eigen::VectorXf currPopVector, int nTimeSteps,
        Eigen::VectorXf dynamicFactorVector // nrow = nTimeSteps
    ) {
        Eigen::VectorXf resPopVector = currPopVector;
        for(int i = 0; i < nTimeSteps; i++) {
            resPopVector = calcDynamicVitalRates(dynamicFactorVector.row(i).value()).getMPM()
                * currPopVector;
        }
        return resPopVector;
    }
    // Project population by 1 time step
    void DynamicPopulation::projectPopulation(float dynamicFactor) {
        popVector = calcPopProjection(popVector, dynamicFactor);
    }
    // Project population by n time steps
    void DynamicPopulation::projectPopulationN(
        int nTimeSteps, Eigen::VectorXf dynamicFactorVector
    ) {
        popVector = calcPopProjectionN(popVector, nTimeSteps, dynamicFactorVector);
    }

    // ====== Dynamic population where lambda is adjusted =====

    VitalRates LambdaDynPopulation::calcDynamicVitalRates(
        float dynamicFactor
    ) {
        // Calculate new lambda
        float targetLambda = vitalRates.calcLambda()
            * std::pow(lambdaBase, dynamicFactor);
        // Adjust vital rates to fit lambda
        VitalRates newVitalRates = vitalRates;
        newVitalRates.adjustVRToLambdaWithConstraints(targetLambda, dynamicVRType, constraints);

        dynamicVRVector = newVitalRates.getVR(dynamicVRType);

        // Return new vital rates
        return newVitalRates;
    }
}