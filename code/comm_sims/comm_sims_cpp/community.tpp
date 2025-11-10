#pragma once

#include <iostream>
#include <vector>
#include <cmath>
#include "Eigen/Dense"

#include "vital_rates.h"
#include "population.h"

namespace Demo {
    void Community::updateStageEqDensities() {
        stageEqDensities = Eigen::VectorXf(numStages * numPopulations);
        
        for(int i = 0; i < numPopulations; i++) {
            Eigen::VectorXf popEqStruct = getPopulation(i).getVitalRates().calcEqPopStruct();
            stageEqDensities.middleRows(i * numStages, numStages) = popEqStruct * eqDensities.row(i).value();
        }
    }
    void Community::updateInteractionAtEq() {
        interactionAtEq = Eigen::MatrixXf(numPopulations, numPopulations);

        Eigen::VectorXf tempCommVector(numStages * numPopulations);
        for(int i = 0; i < numPopulations; i++) {
            tempCommVector = Eigen::VectorXf::Zero(numStages * numPopulations);
            tempCommVector.segment(i * numStages, numStages) = getStageEqDensities().segment(i * numStages, numStages);
            tempCommVector /= getEqDensities().segment(i, 1).value();
            interactionAtEq.col(i) = interactionCoeffMatrix * tempCommVector;
        }
    }

    Eigen::VectorXf Community::calcCommProjection() {
        Community newComm = Community(*this);
        Eigen::VectorXf commVector = newComm.getCommVector();
        Eigen::VectorXf factorVector = interactionCoeffMatrix * (commVector - newComm.getStageEqDensities());
        
        for(int i = 0; i < numPopulations; i++) {
            newComm.getPopulation(i).projectPopulation(factorVector.row(i).value());
        }

        return newComm.getCommVector();
    }
    Eigen::VectorXf Community::calcCommProjectionN(int nTimeSteps) {
        Community newComm = Community(*this);
        Eigen::VectorXf commVector, factorVector;
        
        for(int i = 0; i < nTimeSteps; i++) {
            commVector = newComm.getCommVector();
            factorVector = interactionCoeffMatrix * (commVector - newComm.getStageEqDensities());

            for(int j = 0; j < numPopulations; j++) {
                newComm.getPopulation(j).projectPopulation(factorVector.row(i).value());
            }
        }
        
        return newComm.getCommVector();
    }
    void Community::projectCommunity() {
        setCommVector(calcCommProjection());
    }
    void Community::projectCommunityN(int nTimeSteps) {
        setCommVector(calcCommProjectionN(nTimeSteps));
    }
    // NB: THIS ONLY WORKS IF THERE ARE TWO STAGES
    void Community::adjustStageDependence(float stageDependence) { // SD: -1 ~ +1
        if(numStages != 2) {
            throw std::logic_error("Function is not yet implemented");
        }

        Eigen::RowVectorXf onesVector = Eigen::RowVectorXf::Ones(numPopulations);

        Eigen::MatrixXf newInteractionCoeffMatrix(numPopulations, numStages * numPopulations);

        const bool juvsAreMoreComp = stageDependence < 0;

        Eigen::MatrixXf scaledPartMatrix = getInteractionAtEq();
        // Scale part matrix according to stage dependence factor
        scaledPartMatrix *= 1 - abs(stageDependence);
        Eigen::VectorXf scaledPartEqVector = (getStageEqDensities()(Eigen::seqN(
            juvsAreMoreComp ? 1 : 0,
            numPopulations,
            2
        )).array() / getEqDensities().array()).matrix();

        Eigen::VectorXf adjustedPartEqVector = (getStageEqDensities()(Eigen::seqN(
            juvsAreMoreComp ? 0 : 1,
            numPopulations,
            2
        )).array() / getEqDensities().array()).matrix();
        
        Eigen::MatrixXf adjustedPartMatrix = ((getInteractionAtEq().array()
            - (scaledPartMatrix.array() * (scaledPartEqVector * onesVector).array()))
            / (adjustedPartEqVector * onesVector).array()).matrix();
        
        newInteractionCoeffMatrix(Eigen::all, Eigen::seqN(
            juvsAreMoreComp ? 1 : 0,
            numPopulations,
            2
        )) = scaledPartMatrix;

        newInteractionCoeffMatrix(Eigen::all, Eigen::seqN(
            juvsAreMoreComp ? 0 : 1,
            numPopulations,
            2
        )) = adjustedPartMatrix;

        interactionCoeffMatrix = newInteractionCoeffMatrix;

        // updateInteractionAtEq();
    }

    // StochasticCommunity
    Eigen::VectorXf StochasticCommunity::calcCommProjection(Eigen::VectorXf stochVector) {
        Community newComm = Community(*this);
        Eigen::VectorXf commVector = newComm.getCommVector();
        Eigen::VectorXf factorVector = interactionCoeffMatrix *
            (commVector - newComm.getStageEqDensities()) +
            stochVector;
        
        for(int i = 0; i < numPopulations; i++) {
            newComm.getPopulation(i).projectPopulation(factorVector.row(i).value());
        }

        return newComm.getCommVector();
    }
    
    Eigen::VectorXf StochasticCommunity::calcCommProjectionN(int nTimeSteps, Eigen::MatrixXf stochMatrix) {
        Community newComm = Community(*this);
        Eigen::VectorXf commVector, factorVector;
        
        for(int i = 0; i < nTimeSteps; i++) {
            commVector = newComm.getCommVector();
            factorVector = interactionCoeffMatrix *
                (commVector - newComm.getStageEqDensities()) +
                stochMatrix.col(i);

            for(int i = 0; i < numPopulations; i++) {
                newComm.getPopulation(i).projectPopulation(factorVector.row(i).value());
            }
        }
        
        return newComm.getCommVector();
    }

    void StochasticCommunity::projectCommunity(Eigen::VectorXf stochVector) {
        Eigen::VectorXf commVector = getCommVector();
        Eigen::VectorXf factorVector = interactionCoeffMatrix *
            (commVector - getStageEqDensities()) +
            stochVector;
        
        for(int i = 0; i < numPopulations; i++) {
            getPopulation(i).projectPopulation(factorVector.row(i).value());
        }
    }

    void StochasticCommunity::projectCommunityN(int nTimeSteps, Eigen::MatrixXf stochMatrix) {
        Eigen::VectorXf commVector, factorVector;
        
        for(int i = 0; i < nTimeSteps; i++) {
            commVector = getCommVector();
            factorVector = interactionCoeffMatrix *
                (commVector - getStageEqDensities()) +
                stochMatrix.col(i);

            for(int i = 0; i < numPopulations; i++) {
                getPopulation(i).projectPopulation(factorVector.row(i).value());
            }
        }
    }
}