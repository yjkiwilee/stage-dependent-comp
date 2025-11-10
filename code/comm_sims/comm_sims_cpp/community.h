#pragma once

#include <vector>
#include <random>
#include "Eigen/Dense"

#include "vital_rates.h"
#include "population.h"

namespace Demo {
    class Community {
        protected:
            int numStages; // Number of life cycle stages
            int numPopulations; // Number of populations
            std::vector<LambdaDynPopulation> populations; // Populations
            Eigen::MatrixXf interactionCoeffMatrix; // Interaction coefficients (row: n. pop., col: n. stages * n. pop)
            Eigen::VectorXf eqDensities; // Equilibrium densities
            Eigen::VectorXf stageEqDensities;
            Eigen::MatrixXf interactionAtEq; // Species-level per-capita interaction strengths at equilibrium population structures

            // Internal
            void updateInternal() {
                updateStageEqDensities();
                updateInteractionAtEq();
            };
            void updateStageEqDensities();
            void updateInteractionAtEq();
            
        public:
            // Constructor
            Community(
                int numStages_,
                int numPopulations_,
                std::vector<LambdaDynPopulation> populations_,
                Eigen::MatrixXf interactionCoeffMatrix_,
                Eigen::VectorXf eqDensities_
            ) : numStages(numStages_), numPopulations(numPopulations_),
            populations(populations_), interactionCoeffMatrix(interactionCoeffMatrix_),
            eqDensities(eqDensities_)
            {
                if(
                    populations_.size() != numPopulations_ ||
                    eqDensities_.size() != numPopulations_ ||
                    interactionCoeffMatrix_.cols() != numPopulations_ * numStages_
                ) {
                    throw std::invalid_argument("Invalid number of populations");
                } else if (interactionCoeffMatrix_.rows() != numPopulations_) {
                    throw std::invalid_argument("Invalid number of life stages");
                }
                updateInternal();
            };

            // Getters
            int getNumStages() const { return numStages; };
            int getNumPopulations() const { return numPopulations; };
            std::vector<LambdaDynPopulation>& getPopulations() { return populations; };
            LambdaDynPopulation& getPopulation(int i) { return populations.at(i); };
            Eigen::MatrixXf& getInteractionCoeffs() { return interactionCoeffMatrix; };
            Eigen::VectorXf getCommVector() {
                Eigen::VectorXf commVector(numStages * numPopulations);

                for(int i = 0; i < numPopulations; i++) {
                    commVector.middleRows(i * numStages, numStages) = getPopulation(i).getPopVector();
                }

                return commVector;
            };
            Eigen::VectorXf getVRVector() {
                Eigen::VectorXf vitalRateVector(numStages * numPopulations);

                for(int i = 0; i < numPopulations; i++) {
                    vitalRateVector.middleRows(i * numStages, numStages) = getPopulation(i).getDynamicVRVector();
                }

                return vitalRateVector;
            };
            Eigen::VectorXf getConstraints() {
                Eigen::VectorXf constraintsVector(numStages * numPopulations);

                for(int i = 0; i < numPopulations; i++) {
                    constraintsVector.middleRows(i * numStages, numStages) = getPopulation(i).getConstraints();
                }

                return constraintsVector;
            }
            Eigen::VectorXf getEqDensities() const { return eqDensities; };
            Eigen::VectorXf getStageEqDensities() const { return stageEqDensities; };
            Eigen::MatrixXf getInteractionAtEq() const { return interactionAtEq; };

            // Setters
            void setPopulations(std::vector<LambdaDynPopulation> populations_) {
                if(populations_.size() != numPopulations) {
                    throw std::invalid_argument("Invalid number of populations");
                } // TODO: CHECK FOR N STAGES IN POPULATIONS
                populations = populations_;
                updateInternal();
            }
            void setPopulation(int i, LambdaDynPopulation population) {
                if(population.getNumStages() != numStages) {
                    throw std::invalid_argument("Invalid number of life stages");
                }
                populations.at(i) = population;
                updateInternal();
            }
            void setInteractionCoeffs(Eigen::MatrixXf interactionCoeffMatrix_) {
                if(interactionCoeffMatrix_.cols() != numPopulations) {
                    throw std::invalid_argument("Invalid number of populations");
                } else if (interactionCoeffMatrix_.rows() != numStages * numPopulations) {
                    throw std::invalid_argument("Invalid number of life stages");
                }
                interactionCoeffMatrix = interactionCoeffMatrix_;
                updateInteractionAtEq();
            }
            void setCommVector(Eigen::VectorXf commVector) {
                if(commVector.size() != numPopulations * numStages) {
                    throw std::invalid_argument("Invalid dimensions");
                }
                for(int i = 0; i < numPopulations; i++) {
                    getPopulation(i).setPopVector(commVector.middleRows(i * numStages, numStages));
                }
            }
            void setConstraints(Eigen::VectorXf constraints) {
                if(constraints.size() != numPopulations * numStages) {
                    throw std::invalid_argument("Invalid dimensions");
                }
                for(int i = 0; i < numPopulations; i++) {
                    getPopulation(i).setConstraints(constraints.middleRows(i * numStages, numStages));
                }
            }
            void setEqDensities(Eigen::VectorXf eqDensities_) {
                eqDensities = eqDensities_;
                updateStageEqDensities();
            }

            // Projection functions
            Eigen::VectorXf calcCommProjection();
            Eigen::VectorXf calcCommProjectionN(int nTimeSteps);
            void projectCommunity();
            void projectCommunityN(int nTimeSteps);

            // Calculation function for stage-dependence
            void adjustStageDependence(float stageDependence); // +1: Only adults exert effect, -1: only juveniles exert effect
    };

    class StochasticCommunity: public Community {
        public:
            // Constructor
            using Community::Community;

            // Calculation functions
            Eigen::VectorXf calcCommProjection(Eigen::VectorXf stochVector);
            Eigen::VectorXf calcCommProjectionN(int nTimeSteps, Eigen::MatrixXf stochMatrix);
            void projectCommunity(Eigen::VectorXf stochVector);
            void projectCommunityN(int nTimeSteps, Eigen::MatrixXf stochMatrix);
    };
}

#include "community.tpp"