#pragma once

#include <vector>
#include <stdexcept>
#include <optional>
#include "Eigen/Dense"
#include "misc_funcs.h"

namespace Demo {
    // four types of vital rates
    enum VRType { S, G, R, F };

    class VitalRates {
        private:
            int numStages; // Number of life cycle stages

            Eigen::VectorXf s; // Survival (size = mNumStages)
            Eigen::VectorXf g; // Progression (")
            Eigen::VectorXf r; // Retrogression (")
            Eigen::VectorXf f; // Fertility (")

            bool doesMatrixNeedUpdate = true; // Does MPM require an update?
            Eigen::MatrixXf mpm; // MPM
            void updateMPM(); // Update mMPM from vital rates
            void updateMPMAfterCheck();

            bool doesEigenNeedUpdate = true; // Do eigenvalue & eigenvector need an update?
            Eigen::VectorXcf eigenvalues;
            Eigen::MatrixXcf eigenvectors;
            std::vector<size_t> eigenOrder;
            void updateEigen(); // Update eigenvalues & eigenvectors
            void updateEigenAfterCheck();
            std::optional<float> lambda;
            std::optional<Eigen::VectorXf> eqPopStruct;
            std::optional<float> dampingRatio;
        
        public:
            // Default constructor
            VitalRates(
                int numStages_,
                Eigen::VectorXf s_,
                Eigen::VectorXf g_,
                Eigen::VectorXf r_,
                Eigen::VectorXf f_
            )
                : numStages(numStages_)
                , s(s_), g(g_), r(r_), f(f_)
            {
                if (
                    s_.size() != numStages ||
                    g_.size() != numStages ||
                    r_.size() != numStages ||
                    f_.size() != numStages
                ) {
                    throw std::invalid_argument("Invalid number of life stages");
                }
            };

            // Constructor for deep copies
            VitalRates(const VitalRates& vrs) {
                numStages = vrs.getNumStages();
                s = Eigen::VectorXf(vrs.getS());
                g = Eigen::VectorXf(vrs.getG());
                r = Eigen::VectorXf(vrs.getR());
                f = Eigen::VectorXf(vrs.getF());
            }

            // Getters
            int getNumStages() const { return numStages; };
            Eigen::VectorXf getS() const { return s; };
            Eigen::VectorXf getG() const { return g; };
            Eigen::VectorXf getR() const { return r; };
            Eigen::VectorXf getF() const { return f; };
            Eigen::VectorXf getVR(const VRType vrType) const {
                switch(vrType) {
                    case S:
                        return s;
                    case G:
                        return g;
                    case R:
                        return r;
                    case F:
                        return f;
                }
            };
            Eigen::MatrixXf getMPM() {
                // Update MPM if needed
                updateMPMAfterCheck();
                
                return mpm;
            };
            Eigen::VectorXcf getEigenvalues() {
                // Update eigen if needed
                updateEigenAfterCheck();

                return eigenvalues;
            }
            Eigen::MatrixXcf getEigenvectors() {
                // Update eigen if needed
                updateEigenAfterCheck();

                return eigenvectors;
            }

            // Setters
            void setS(Eigen::VectorXf s_) {
                if(s_.size() != numStages) {
                    throw std::invalid_argument("Invalid number of life stages");
                }
                s = s_; doesMatrixNeedUpdate = true; doesEigenNeedUpdate = true;
            };
            void setG(Eigen::VectorXf g_) {
                if(g_.size() != numStages) {
                    throw std::invalid_argument("Invalid number of life stages");
                }
                g = g_; doesMatrixNeedUpdate = true; doesEigenNeedUpdate = true;
            };
            void setR(Eigen::VectorXf r_) {
                if(r_.size() != numStages) {
                    throw std::invalid_argument("Invalid number of life stages");
                }
                r = r_; doesMatrixNeedUpdate = true; doesEigenNeedUpdate = true;
            };
            void setF(Eigen::VectorXf f_) {
                if(f_.size() != numStages) {
                    throw std::invalid_argument("Invalid number of life stages");
                }
                f = f_; doesMatrixNeedUpdate = true; doesEigenNeedUpdate = true;
            };
            void setVR(const VRType vrType, Eigen::VectorXf vr_) {
                if(vr_.size() != numStages) {
                    throw std::invalid_argument("Invalid number of life stages");
                }
                switch(vrType) {
                    case S:
                        s = vr_;
                        break;
                    case G:
                        g = vr_;
                        break;
                    case R:
                        r = vr_;
                        break;
                    case F:
                        f = vr_;
                        break;
                }

                doesMatrixNeedUpdate = true;
                doesEigenNeedUpdate = true;
            };

            // Calculating functions
            float calcLambda(); // Calculate lambda
            Eigen::VectorXf calcEqPopStruct(); // Calculate equilibrium pop struct
            float calcDampingRatio(); // Calculate damping ratio

            // Internal function for calculating the scaling factor of fertility for a given value of lambda
            float findVRScaleAtLambda(
                float targetLambda,
                VRType adjustedVRType,
                std::optional<int> nIterations = std::nullopt,
                std::optional<float> initLScale = std::nullopt,
                std::optional<float> initRScale = std::nullopt,
                std::optional<float> tolerance = std::nullopt
            );
            void adjustVRToLambda(
                float targetLambda,
                VRType adjustedVRType,
                std::optional<int> nIterations = std::nullopt,
                std::optional<float> initLScale = std::nullopt,
                std::optional<float> initRScale = std::nullopt,
                std::optional<float> tolerance = std::nullopt
            );

            // The above with stage-specific constraints
            float findVRScaleAtLambdaWithConstraints(
                float targetLambda,
                VRType adjustedVRType,
                Eigen::VectorXf constraints,
                std::optional<int> nIterations = std::nullopt,
                std::optional<float> initLScale = std::nullopt,
                std::optional<float> initRScale = std::nullopt,
                std::optional<float> tolerance = std::nullopt
            );
            void adjustVRToLambdaWithConstraints(
                float targetLambda,
                VRType adjustedVRType,
                Eigen::VectorXf constraints,
                std::optional<int> nIterations = std::nullopt,
                std::optional<float> initLScale = std::nullopt,
                std::optional<float> initRScale = std::nullopt,
                std::optional<float> tolerance = std::nullopt
            );
    };
};

#include "vital_rates.tpp"