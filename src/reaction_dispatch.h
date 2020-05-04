#ifndef CTF_REACTION_DISPATCH_H_
#define CTF_REACTION_DISPATCH_H_
#include <thread>
#include <mutex>
#include <functional>
#include <vector>
#include <Eigen/Dense>
#include "ChemThermo.h"

class ReactionDispatch
{
private:
    ChemThermo& primary_gas_;
    ChemThermo& gas_helper_;
    int nsp_;
    bool init_;
    bool completed_;
    std::thread thread_;
    std::mutex mutex1_;
    std::mutex mutex2_;


    // solution data
    double tc_helper;
    double deltaT_helper;
    int half_;
    int len_;
    Eigen::VectorXd hs_helper;
    std::vector<Eigen::VectorXd> Y_helper;
    std::vector<Eigen::VectorXd> wdot_helper;

    void thread_solve();

public:

    ReactionDispatch(ChemThermo& primary, ChemThermo& helper);
    ~ReactionDispatch();
    double solve(const double& deltaT, const Eigen::VectorXd& hs, const std::vector<Eigen::VectorXd>& Y,
        std::vector<Eigen::VectorXd>& wdot, Eigen::VectorXd& qdot);
    void complete();  // call from main thread
};

void ReactionDispatch::thread_solve()
{
    mutex2_.lock();
    while (true) {
        mutex1_.lock();
        // log_debug("[thread_solve] mutex1 locked");
        mutex2_.unlock();
        // log_debug("[thread_solve] mutex2 unlocked");

        if (completed_) {
            log_debug("[thread_solve] completed");
            break;
        }

        for (int j = half_; j < len_; j++) {
            double y[nsp_];
            gas_helper_.massFractions(Y_helper, y, j);
            gas_helper_.setY(y);
            gas_helper_.thermo_.p() = dimensionedScalar("p", dimPressure, gas_helper_.p0_);
            gas_helper_.thermo_.he() = dimensionedScalar("h", dimEnergy/dimMass, hs_helper(j));
            gas_helper_.thermo_.correct();
            tc_helper = min(tc_helper, gas_helper_.chemistry_.solve(deltaT_helper));
            for (int k = 0; k < nsp_; k++) {
                wdot_helper[k](j) = gas_helper_.chemistry_.RR(k)[0];
            }
        }

        mutex1_.unlock();
        // log_debug("[thread_solve] mutex1 unlocked");
        mutex2_.lock();
        // log_debug("[thread_solve] mutex2 locked");
    }
    log_debug("[thread_solve] exiting");
}

ReactionDispatch::ReactionDispatch(ChemThermo& primary, ChemThermo& helper)
    : primary_gas_(primary)
    , gas_helper_(helper)
    , nsp_(primary.nsp())
    , init_(false)
    , completed_(false)
    , Y_helper(primary.nsp())
    , wdot_helper(primary.nsp())
{
    mutex1_.lock();
    thread_ = std::thread(std::bind(&ReactionDispatch::thread_solve, this));
}

ReactionDispatch::~ReactionDispatch()
{

}

double ReactionDispatch::solve(const double& deltaT, const Eigen::VectorXd& hs, const std::vector<Eigen::VectorXd>& Y,
    std::vector<Eigen::VectorXd>& wdot, Eigen::VectorXd& qdot)
{
    tc_helper = 1.0;
    deltaT_helper = deltaT;
    double tc = 1.0;
    int nsp = Y.size();
    int len = hs.size();
    int half = len / 2;

    if (!init_) {
        half_ = half;
        len_ = len;
        hs_helper.resize(len);
        for (int k = 0; k < nsp_; k++) {
            Y_helper[k].resize(len);
            wdot_helper[k].resize(len);
        }
        init_ = true;
    }

    // copy into helpers
    for (int j = half; j < len; j++) {
        hs_helper(j) = hs(j);
    }
    for (int k = 0; k < nsp_; k++) {
        for (int j = half; j < len; j++) {
            Y_helper[k](j) = Y[k](j);
            wdot_helper[k](j) = 0.0;
        }
    }
    mutex1_.unlock();
    // log_debug("[solve] mutex1 unlocked");
    mutex2_.lock();
    // log_debug("[solve] mutex2 locked");
    for (int j = 0; j < half; j++) {
        double y[nsp];
        primary_gas_.massFractions(Y, y, j);
        primary_gas_.setY(y);
        primary_gas_.thermo_.p() = dimensionedScalar("p", dimPressure, primary_gas_.p0_);
        primary_gas_.thermo_.he() = dimensionedScalar("h", dimEnergy/dimMass, hs(j));
        primary_gas_.thermo_.correct();
        tc = min(tc, primary_gas_.chemistry_.solve(deltaT));
        for (int k = 0; k < nsp; k++) {
            wdot[k](j) = primary_gas_.chemistry_.RR(k)[0];
        }
    }
    mutex1_.lock();
    // log_debug("[solve] mutex1 locked");
    mutex2_.unlock();
    // log_debug("[solve] mutex2 unlocked");
    tc = min(tc, tc_helper);
    for (int k = 0; k < nsp; k++) {
        for (int j = half; j < len; j++) {
            wdot[k](j) = wdot_helper[k](j);
        }
    }
    // Compute qdot
    for (int j = 0; j < len; j++) {
        qdot(j) = 0.0;
        for (int k = 0; k < nsp_; k++) {
            qdot(j) -= primary_gas_.composition_.Hc(k)*wdot[k](j);
        }
    }
    return max(tc, small);
}

void ReactionDispatch::complete()
{
    completed_ = true;
    mutex1_.unlock();
    if (thread_.joinable()) {
        log_debug("[complete] join thread");
        thread_.join();
        log_debug("[complete] successfully joined");
    }
}

#endif