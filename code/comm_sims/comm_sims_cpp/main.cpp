#include <iostream>
#include <optional>
#include <random>
#include <map>

#include "libs.h"

int main(int argc, char **argv) {
    argparse::ArgumentParser program("comm_sims");

    program.add_argument("tsteps").scan<'i', int>();
    program.add_argument("--pres").default_value(-1).scan<'i', int>();
    program.add_argument("--eqdens").default_value(200.0f).scan<'g', float>();
    program.add_argument("--stochsd").default_value(0.1f).scan<'g', float>();
    program.add_argument("--intra").default_value(-0.01f).scan<'g', float>();
    program.add_argument("--inter").default_value(-0.001f).scan<'g', float>();
    program.add_argument("--recipstagedep").default_value(1.0f).scan<'g', float>();
    program.add_argument("--lbase").default_value(1.01f).scan<'g', float>();
    program.add_argument("--stagedep").default_value(0.0f).scan<'g', float>();
    program.add_argument("--randseed").default_value(1).scan<'i', int>();
    program.add_argument("--initjprop").default_value(0.5f).scan<'g', float>();
    program.add_argument("-sj").default_value(0.2f).scan<'g', float>();
    program.add_argument("-sa").default_value(0.8f).scan<'g', float>();
    program.add_argument("-g").default_value(0.1f).scan<'g', float>();
    program.add_argument("-r").default_value(0.0f).scan<'g', float>();
    program.add_argument("--ddvr").default_value("f");
    program.add_argument("--autocorr").default_value(0.0f).scan<'g', float>();
    program.add_argument("--dispdiff").default_value(""); // a: absolute, r: relative
    
    try {
        program.parse_args(argc, argv);
    }
    catch (const std::exception& err) {
        std::cerr << err.what() << std::endl;
        std::cerr << program;
        return 1;
    }
    
    int NUM_STAGES = 2;
    int NUM_POPS = 2;
    auto EQ_DENSITY = program.get<float>("--eqdens");

    auto TIME_STEPS = program.get<int>("tsteps");
    auto PRINT_RES = program.get<int>("--pres");
    if(PRINT_RES == -1) { PRINT_RES = std::max(TIME_STEPS / 10, 1); }
    auto DISP_DIFF = program.get<std::string>("--dispdiff");

    auto STOCH_SD = program.get<float>("--stochsd");
    auto INTRA = program.get<float>("--intra");
    auto INTER = program.get<float>("--inter");
    auto RECIP_SD_LEVEL = program.get<float>("--recipstagedep");
    Eigen::VectorXf CONSTRAINTS {{1.0f, RECIP_SD_LEVEL}};

    auto LAMBDA_BASE = program.get<float>("--lbase");
    auto SD_LEVEL = program.get<float>("--stagedep");
    auto RAND_SEED = program.get<int>("--randseed");
    auto AUTOCORR_R = program.get<float>("--autocorr");
    auto INIT_JUV_PROP = program.get<float>("--initjprop");
    auto SJ = program.get<float>("-sj");
    auto SA = program.get<float>("-sa");
    auto G = program.get<float>("-g");
    auto R = program.get<float>("-r");
    Demo::VRType DDVR_TYPE;
    switch(program.get<std::string>("--ddvr")[0]) {
        case 's':
            DDVR_TYPE = Demo::S;
            break;
        case 'g':
            DDVR_TYPE = Demo::G;
            break;
        case 'r':
            DDVR_TYPE = Demo::R;
            break;
        case 'f':
            DDVR_TYPE = Demo::F;
            break;
    }

    Eigen::VectorXf s {{SJ, SA}};
    Eigen::VectorXf g {{G, 0}};
    Eigen::VectorXf r {{0, R}};
    Eigen::VectorXf f {{0, 1}};

    Demo::VitalRates vr(NUM_STAGES, s, g, r, f);

    float scaleAtEq = vr.findVRScaleAtLambda(
        1, Demo::F, 50, 0.0, 2.0, 0.000001
    );

    vr.setF(f * scaleAtEq);

    Eigen::VectorXf popStruct = vr.calcEqPopStruct() * EQ_DENSITY;
    Eigen::VectorXf initCommStruct {{
        EQ_DENSITY * INIT_JUV_PROP,
        EQ_DENSITY * (1 - INIT_JUV_PROP),
        EQ_DENSITY * INIT_JUV_PROP,
        EQ_DENSITY * (1 - INIT_JUV_PROP)
    }};

    Demo::LambdaDynPopulation pop1(
        NUM_STAGES, vr, initCommStruct.segment(0, 2), DDVR_TYPE, LAMBDA_BASE
    );
    Demo::LambdaDynPopulation pop2(
        NUM_STAGES, vr, initCommStruct.segment(2, 2), DDVR_TYPE, LAMBDA_BASE
    );

    std::vector<Demo::LambdaDynPopulation> pops = {
        pop1, pop2
    };

    Eigen::MatrixXf interactionMatrix {
        {INTRA, INTRA, INTER, INTER},
        {INTER, INTER, INTRA, INTRA}
    };

    Eigen::VectorXf eqDensities {{EQ_DENSITY, EQ_DENSITY}};

    Eigen::MatrixXf randFactorMatrix(NUM_POPS, TIME_STEPS);

    for(int i = 0; i < NUM_POPS; i++) {
        Eigen::VectorXf randFactorVector = generateAutoCorrSample(TIME_STEPS, AUTOCORR_R, 0.0f, STOCH_SD, i + RAND_SEED);
        randFactorMatrix.row(i) = randFactorVector;
    }

    Demo::StochasticCommunity comm(
        NUM_STAGES, NUM_POPS, pops, interactionMatrix, eqDensities
    );

    struct CommRecordRow {
        int timeStep;
        Eigen::VectorXf densities;
    };

    struct VRRecordRow {
        int timeStep;
        Eigen::VectorXf vitalRates;
    };

    std::vector<CommRecordRow> commRecords(TIME_STEPS / PRINT_RES);
    std::vector<VRRecordRow> vitalRateRecords(TIME_STEPS / PRINT_RES);
    std::cout << std::endl;

    for(int i = 0; i < TIME_STEPS; i++) {
        std::cout << "\33[2K\r" << "original " << i << std::flush;
        // std::cout << "original " << i << std::flush;
        if(i % PRINT_RES == 0) {
            commRecords[i / PRINT_RES] = {i, Eigen::VectorXf::Zero(NUM_STAGES * NUM_POPS * 2)};
            commRecords[i / PRINT_RES].densities(Eigen::seqN(0, NUM_POPS * NUM_STAGES, 2)) = comm.getCommVector();
        }

        comm.projectCommunity(randFactorMatrix.col(i));

        if(i % PRINT_RES == 0) {
            vitalRateRecords[i / PRINT_RES] = {i, Eigen::VectorXf::Zero(NUM_STAGES * NUM_POPS * 2)};
            vitalRateRecords[i / PRINT_RES].vitalRates(Eigen::seqN(0, NUM_POPS * NUM_STAGES, 2)) = comm.getVRVector();
        }
    }
    std::cout << std::endl;

    comm.adjustStageDependence(SD_LEVEL);
    comm.setCommVector(initCommStruct);
    Eigen::VectorXf constraintsVector(NUM_STAGES * NUM_POPS);
    constraintsVector << CONSTRAINTS, CONSTRAINTS;
    comm.setConstraints(constraintsVector);
    
    for(int i = 0; i < TIME_STEPS; i++) {
        std::cout << "\33[2K\r" << "SD " << i << std::flush;
        // std::cout << "SD " << i << std::flush;
        if(i % PRINT_RES == 0) {
            commRecords[i / PRINT_RES].densities(Eigen::seqN(1, NUM_POPS * NUM_STAGES, 2)) = comm.getCommVector();
        }
        comm.projectCommunity(randFactorMatrix.col(i));
        if(i % PRINT_RES == 0) {
            vitalRateRecords[i / PRINT_RES].vitalRates(Eigen::seqN(1, NUM_POPS * NUM_STAGES, 2)) = comm.getVRVector();
        }
    }

    std::cout << std::endl << "DENSITIES" << std::endl;

    for(int i = 0; i < TIME_STEPS / PRINT_RES; i++) {
        if(DISP_DIFF == "") {
            std::cout << commRecords[i].timeStep << "\t" << commRecords[i].densities.transpose() << std::endl;
        } else {
            Eigen::VectorXf diffVector;

            if(DISP_DIFF == "a") { // absolute diff
                diffVector =
                    (commRecords[i].densities(Eigen::seq(0, Eigen::last, 2)).array()
                    - commRecords[i].densities(Eigen::seq(1, Eigen::last, 2)).array()).matrix();
            } else if(DISP_DIFF == "r") { // relative diff
                diffVector =
                    (((commRecords[i].densities(Eigen::seq(0, Eigen::last, 2)).array()
                    - commRecords[i].densities(Eigen::seq(1, Eigen::last, 2)).array()))
                    / commRecords[i].densities(Eigen::seq(0, Eigen::last, 2)).array()).matrix();
            }

            std::cout << commRecords[i].timeStep << "\t" << diffVector.transpose() << std::endl;
        }
    }

    std::cout << std::endl << "VRS" << std::endl;

    for(int i = 0; i < TIME_STEPS / PRINT_RES; i++) {
        if(DISP_DIFF == "") {
            std::cout << vitalRateRecords[i].timeStep << "\t" << vitalRateRecords[i].vitalRates.transpose() << std::endl;
        } else {
            Eigen::VectorXf diffVector;

            if(DISP_DIFF == "a") { // absolute diff
                diffVector =
                    (vitalRateRecords[i].vitalRates(Eigen::seq(0, Eigen::last, 2)).array()
                    - vitalRateRecords[i].vitalRates(Eigen::seq(1, Eigen::last, 2)).array()).matrix();
            } else if(DISP_DIFF == "r") { // relative diff
                diffVector =
                    (((vitalRateRecords[i].vitalRates(Eigen::seq(0, Eigen::last, 2)).array()
                    - vitalRateRecords[i].vitalRates(Eigen::seq(1, Eigen::last, 2)).array()))
                    / vitalRateRecords[i].vitalRates(Eigen::seq(0, Eigen::last, 2)).array()).matrix();
            }

            std::cout << vitalRateRecords[i].timeStep << "\t" << diffVector.transpose() << std::endl;
        }
    }

    return 0;
}