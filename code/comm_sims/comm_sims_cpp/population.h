#pragma once

#include <vector>
#include "Eigen/Dense"

#include "vital_rates.h"

namespace Demo {
    class Population {
        protected:
            int numStages; // Number of life cycle stages
            VitalRates vitalRates; // Vital rates of the population
            Eigen::VectorXf popVector; // Population vector

        public:
            // Constructor
            Population(
                int numStages_,
                VitalRates vitalRates_,
                Eigen::VectorXf popVector_
            )
                : numStages(numStages_)
                , vitalRates(vitalRates_)
                , popVector(popVector_)
            {
                if (
                    vitalRates_.getNumStages() != numStages ||
                    popVector_.size() != numStages
                ) {
                    throw std::invalid_argument("Invalid number of life stages");
                }
            };

            // Getters
            int getNumStages() const { return numStages; };
            VitalRates getVitalRates() const { return vitalRates; };
            Eigen::VectorXf getPopVector() const {return popVector; };

            // Setters
            void setVitalRates(VitalRates vitalRates_) {
                if(vitalRates_.getNumStages() != numStages) {
                    throw std::invalid_argument("Invalid number of life stages");
                }
                vitalRates = vitalRates_;
            };
            void setPopVector(Eigen::VectorXf popVector_) {
                if(popVector_.size() != numStages) {
                    throw std::invalid_argument("Invalid number of life stages");
                }
                popVector = popVector_;
            };

            // Projection functions
            // Internal function for calculating projection
            Eigen::VectorXf calcPopProjection();
            Eigen::VectorXf calcPopProjection(Eigen::VectorXf currPopVector);
            Eigen::VectorXf calcPopProjectionN(int nTimeSteps);
            Eigen::VectorXf calcPopProjectionN(Eigen::VectorXf currPopVector, int nTimeSteps);
            void projectPopulation(); // Project population by 1 step
            void projectPopulationN(int nTimeSteps); // Project population by a given number of steps

            // Calculation functions
            void adjustVRToLambda(
                float targetLambda,
                VRType adjustedVRType,
                std::optional<int> nIterations = std::nullopt,
                std::optional<float> initLScale = std::nullopt,
                std::optional<float> initRScale = std::nullopt,
                std::optional<float> tolerance = std::nullopt
            ); // Adjust VR so that lambda is a given value
    };

    class DynamicPopulation: public Population {
        protected:
            VRType dynamicVRType; // Type of dynamic vital rate
            Eigen::VectorXf dynamicVRVector; // Snapshot of dynamic VR vector
            Eigen::VectorXf constraints;

        public:
            // Constructor
            DynamicPopulation(
                int numStages_,
                VitalRates vitalRates_,
                Eigen::VectorXf popVector_,
                VRType dynamicVRType_,
                std::optional<Eigen::VectorXf> constraints_
            );

            // Getter
            Eigen::VectorXf getDynamicVRVector() { return dynamicVRVector; };
            Eigen::VectorXf getConstraints() { return constraints; };

            // Setter
            void setConstraints(Eigen::VectorXf constraints_) {
                // TODO: ADD SANITY CHECK
                constraints = constraints_;
            }
            
            // Internal function for calculating projection
            Eigen::VectorXf calcPopProjection(float dynamicFactor);
            Eigen::VectorXf calcPopProjection(
                Eigen::VectorXf currPopVector,
                float dynamicFactor
            );
            Eigen::VectorXf calcPopProjectionN(
                int nTimeSteps,
                Eigen::VectorXf dynamicFactorVector
            );
            Eigen::VectorXf calcPopProjectionN(
                Eigen::VectorXf currPopVector, int nTimeSteps,
                Eigen::VectorXf dynamicFactorVector // nrow = nTimeSteps
            );
            void projectPopulation(float dynamicFactor);
            void projectPopulationN(int nTimeSteps, Eigen::VectorXf dynamicFactorVector);

            // Dynamic calculation function
            virtual VitalRates calcDynamicVitalRates(
                float dynamicFactor
            ) = 0;
    };

    class LambdaDynPopulation: public DynamicPopulation {
        protected:
            float lambdaBase; // Base of adjusting function

        public:
            LambdaDynPopulation(
                int numStages_,
                VitalRates vitalRates_,
                Eigen::VectorXf popVector_,
                VRType dynamicVRType_,
                float lambdaBase_,
                std::optional<Eigen::VectorXf> constraints_ = std::nullopt
            ) : DynamicPopulation(
                numStages_, vitalRates_, popVector_,
                dynamicVRType_, constraints_
            ), lambdaBase(lambdaBase_)
            {};

            VitalRates calcDynamicVitalRates(
                float dynamicFactor
            ) override;
    };
};

#include "population.tpp"