// ChemFlow
// -- ChemFlow
// -- A Segregated Solution to the Quasi-1D Counterflow Flame
// -- xu-zhang@sjtu.edu.cn
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <ctime>
#include <stdexcept>
#include <vector>
#include <string>
#include <map>
#include <Eigen/Dense>

#include "lininterp.h"
#include "tdma.h"
#include "ChemThermo.h"

const double numlimSmall        =       1e-08;
const double numlimSmallSmall   =       1e-14;
const double numlimGreat        =        1e12;
const double tprecision         =       1e-08;
const size_t WIDTH              =          18;
const double Le                 =         1.0;

struct gas_phase {
    // constructor
    gas_phase(int nx, int nsp) :
        u(nx), V(nx), T(nx), hs(nx), Y(nsp), wdot(nsp), qdot(nx),
        rho(nx), rhoPrev(nx), mu(nx), kappa(nx), alpha(nx), D(nx)
    {
        for (int k = 0; k < nsp; k++) {
            Y[k].resize(nx);
            wdot[k].resize(nx);
            for (int j = 0; j < nx; j++) {
                Y[k](j) = 0.0;
                wdot[k](j) = 0.0;
            }
        }
    }

    // public member data
    // primary variables
    // Both std::vector and Eigen::VectorXd are used for the purpose of
    // distinguishing between containers for species and spatial grid points
    Eigen::VectorXd u;  // x-direction velocity [m/s]
    Eigen::VectorXd V;  // v/y = dv/dy [1/s]
    Eigen::VectorXd T;  // temperature [K]
    Eigen::VectorXd hs;  // sensible enthalpy [J/kg]
    std::vector<Eigen::VectorXd> Y;  // species mass fractions [-]
    std::vector<Eigen::VectorXd> wdot;  // reaction rates [kg/m3 s]
    Eigen::VectorXd qdot;  // heat source [J/m3 s]
    // properties
    Eigen::VectorXd rho;
    Eigen::VectorXd rhoPrev;
    Eigen::VectorXd mu;
    Eigen::VectorXd kappa;
    Eigen::VectorXd alpha;
    Eigen::VectorXd D;
};

struct {
    long nx;
    double XBEG;
    double XEND;
    double tBEG;
    double tEND;
    double dtMax;
    double ainit;
    // BC
    double VL;
    double VR;
    double TI;
    double TL;
    double TR;
    double YO2Air;
    double YN2Air;
    double YARAir;
    std::vector<double> YFUELS;  // multi-component fuel
    std::vector<double> YL;
    std::vector<double> YR;
    std::vector<std::string> FUELNAMES; // multi-component fuel
    double hsL;
    double hsR;

    bool rstr;
    bool ign;
    bool strain;
    bool log;
    double ignBEGt;
    double ignENDt;
    double ignHs;
    double aBEG;
    double aEND;
    int writeFreq;
    double rhoInf;
    double p0;
    double lrRatio;
} input_data;


void fill_input(const std::string fname);
void restore(const ChemThermo& gas, const Eigen::VectorXd& x, gas_phase& gp);
void write(const double iter, const ChemThermo& gas,
    const Eigen::VectorXd& x, const gas_phase& gp);

int main(int argc, char *argv[])
{
    // Input
    fill_input("input.txt");
    const double dx = (input_data.XEND - input_data.XBEG)/(input_data.nx - 1);
    double a = input_data.ainit;
    Eigen::VectorXd x(input_data.nx);
    double dtChem = numlimSmall;
    double dt = dtChem;
    double time = input_data.tBEG;
    // Output
    std::ofstream fm("data/monitor.csv");
    fm << "time (s),temperature (K)" << std::endl;

    // Solution and initial conditions
    #include "createFields.H"
    ChemThermo gas(mesh, runTime, input_data.p0);
    ChemThermo* gas_helper = new ChemThermo(mesh_helper, runTime_helper, input_data.p0);
    const int nsp = gas.nsp();  // number of species
    // construct gas_phase
    gas_phase gp(input_data.nx, nsp);
    input_data.YL.resize(nsp, 0.0);
    input_data.YR.resize(nsp, 0.0);
    int nFuels = input_data.FUELNAMES.size();
    for (int k = 0; k < nFuels; k++) {
        string s(input_data.FUELNAMES[k]);
        input_data.YL[gas.speciesIndex(s)] = input_data.YFUELS[k];
    }
    input_data.YR[gas.speciesIndex("O2")] = input_data.YO2Air;
    input_data.YR[gas.speciesIndex("N2")] = input_data.YN2Air;
    input_data.YR[gas.speciesIndex("AR")] = input_data.YARAir;
    input_data.hsL = gas.calcHs(input_data.TL, input_data.YL.data());
    input_data.hsR = gas.calcHs(input_data.TR, input_data.YR.data());

    for (int j = 0; j < input_data.nx; j++) {
        x(j) = input_data.XBEG + dx*j;
        gp.u(j) = -a*(x(j) - 0.5*(input_data.XEND-input_data.XBEG));
        gp.V(j) = a;
        gp.T(j) = input_data.TI;
        gp.Y[gas.speciesIndex("O2")](j) = input_data.YO2Air;
        gp.Y[gas.speciesIndex("N2")](j) = input_data.YN2Air;
        gp.Y[gas.speciesIndex("AR")](j) = input_data.YARAir;
        double y[nsp];
        gas.massFractions(gp.Y, y, j);
        gp.hs(j) = gas.calcHs(gp.T(j), y);
        gp.qdot(j) = 0.0;
    }

    if (input_data.rstr) restore(gas, x, gp);
    gas.updateThermo(gp.hs, gp.Y, Le, gp.rho, gp.mu, gp.kappa, gp.alpha, gp.D);
    gp.rhoPrev = gp.rho;


    // Time marching
    clock_t startTime, endTime;
    startTime = std::clock();
    Eigen::MatrixXd A(input_data.nx,input_data.nx);
    Eigen::MatrixXd b(input_data.nx,1);
    // conservative form for continuity equation
    Eigen::VectorXd m(input_data.nx);
    Eigen::VectorXd::Index loc;
    int iter = 0;
    while (true) {
        if (input_data.strain) {
            a = input_data.aBEG +
                (input_data.aEND-input_data.aBEG)/
                (input_data.tEND-input_data.tBEG)*
                (time-input_data.tBEG);  // strain the flame
            if (iter++%input_data.writeFreq == 0) write(iter, gas, x, gp);
        } else if (iter++%input_data.writeFreq == 0) {
            write(time, gas, x, gp);
        }
        fm << std::setprecision(10) << time << ","
           << std::setprecision(10) << gp.T.maxCoeff(&loc) << std::endl;
        // V equation
        A.setZero();
        b.setZero();
        for (int j=1; j<input_data.nx-1; j++) {
            const double mul = 0.5*(gp.mu(j)+gp.mu(j-1));
            const double mur = 0.5*(gp.mu(j)+gp.mu(j+1));
            A(j,j-1) = -mul*dt/(gp.rho(j)*dx*dx) +
                       (gp.u(j) > 0.0 ? -dt*gp.u(j)/dx : 0.0);
            A(j,j) = 1.0 + dt*gp.V(j) + mul*dt/(gp.rho(j)*dx*dx) +
                     mur*dt/(gp.rho(j)*dx*dx) +
                     (gp.u(j) > 0.0 ? dt*gp.u(j)/dx : -dt*gp.u(j)/dx);
            A(j,j+1) = -mur*dt/(gp.rho(j)*dx*dx) +
                       (gp.u(j) > 0.0 ? 0.0 : dt*gp.u(j)/dx);
            b(j) = input_data.rhoInf*a*a*dt/gp.rho(j) + gp.V(j);
        }
        A(0,0) = 1.0;
        A(input_data.nx-1,input_data.nx-1) = 1.0;
        b(0) = input_data.VL;
        b(input_data.nx-1) = input_data.VR;
        gp.V = tdma(A,b);
        if (input_data.log) {
        std::cout << std::setw(WIDTH) << "V.max "
                  << std::setw(WIDTH) << gp.V.maxCoeff(&loc) << " @ position "
                  << loc << std::endl;
        }

        // Continuity equation
        // Propagate from left to right
        m.setZero();
        m(0) = gp.rho(0) * gp.u(0);
        for (int j=1; j<input_data.nx; j++) {
            double drhodt0 = (numlimGreat*(gp.rho(j-1) - gp.rhoPrev(j-1)))/
                             (numlimGreat*dt);
            double drhodt1 = (numlimGreat*(gp.rho(j) - gp.rhoPrev(j)))/
                             (numlimGreat*dt);
            drhodt0 = (dt > tprecision ? drhodt0 : 0.0);
            drhodt1 = (dt > tprecision ? drhodt1 : 0.0);
            // double drhodt0 = 0.0;
            // double drhodt1 = 0.0;
            m(j) = m(j-1) + dx*(-0.5*(drhodt0+drhodt1) -
                   0.5*(gp.rho(j-1)*gp.V(j-1)+gp.rho(j)*gp.V(j)));
        }
        const double rhouOffset = (-input_data.lrRatio*gp.rho(0)*
                                   m(input_data.nx-1) -
                                   gp.rho(input_data.nx-1)*m(0)) /
                                  (input_data.lrRatio*gp.rho(0) +
                                   gp.rho(input_data.nx-1));
        m = m.array() + rhouOffset;
        gp.u = m.cwiseQuotient(gp.rho);
        if (input_data.log) {
        std::cout << std::setw(WIDTH) << "u.max "
                  << std::setw(WIDTH) << gp.u.maxCoeff(&loc) << " @ position "
                  << loc << std::endl;
        }

        dtChem = gas.solve(dt, gp.hs, gp.Y, gp.wdot, gp.qdot, *gas_helper);
        // Y equations
        for (int k=0; k<nsp; k++) {
            A.setZero();
            b.setZero();
            for (int j=1; j<input_data.nx-1; j++) {
                const double rhoDl = 0.5*(gp.rho(j)*gp.D(j)+
                                          gp.rho(j-1)*gp.D(j-1));
                const double rhoDr = 0.5*(gp.rho(j)*gp.D(j)+
                                          gp.rho(j+1)*gp.D(j+1));
                A(j,j-1) = -rhoDl*dt/(gp.rho(j)*dx*dx) +
                           (gp.u(j) > 0.0 ? -dt*gp.u(j)/dx : 0.0);
                A(j,j) = 1.0 + rhoDl*dt/(gp.rho(j)*dx*dx) +
                         rhoDr*dt/(gp.rho(j)*dx*dx) +
                         (gp.u(j) > 0.0 ? dt*gp.u(j)/dx : -dt*gp.u(j)/dx);
                A(j,j+1) = -rhoDr*dt/(gp.rho(j)*dx*dx) +
                           (gp.u(j) > 0.0 ? 0.0 : dt*gp.u(j)/dx);
                b(j) = gp.Y[k](j) + dt*gp.wdot[k](j)/gp.rho(j);
            }
            A(0,0) = 1.0;
            A(input_data.nx-1,input_data.nx-1) = 1.0;
            b(0) = input_data.YL[k];
            b(input_data.nx-1) = input_data.YR[k];
            gp.Y[k] = tdma(A,b);
            // std::cout << std::setw(WIDTH)
            //           << "Y-" + gas.speciesName(k) + ".max "
            //           << std::setw(WIDTH)
            //           << gp.Y[k].maxCoeff(&loc)
            //           << " @ position "
            //           << loc << std::endl;
        }
        // Correct
        for (int j=0; j<input_data.nx; j++) {
            double sumY = 0.0;
            for (int k=0; k<nsp; k++) {
                gp.Y[k](j) = (gp.Y[k](j) > 0.0 ? gp.Y[k](j) : 0.0);
                gp.Y[k](j) = (gp.Y[k](j) < 1.0 ? gp.Y[k](j) : 1.0);
                sumY += gp.Y[k](j);
            }
            for (int k=0; k<nsp; k++) {
                gp.Y[k](j) /= sumY;
            }
        }

        // Energy eqaution
        A.setZero();
        b.setZero();
        for (int j=1; j<input_data.nx-1; j++) {
            const double rhoAlphal = 0.5*(gp.rho(j)*gp.alpha(j)+
                                          gp.rho(j-1)*gp.alpha(j-1));
            const double rhoAlphar = 0.5*(gp.rho(j)*gp.alpha(j)+
                                          gp.rho(j+1)*gp.alpha(j+1));
            A(j,j-1) = -rhoAlphal*dt/(gp.rho(j)*dx*dx) +
                       (gp.u(j) > 0.0 ? -dt*gp.u(j)/dx : 0.0);
            A(j,j) = 1.0 + rhoAlphal*dt/(gp.rho(j)*dx*dx) +
                     rhoAlphar*dt/(gp.rho(j)*dx*dx) +
                     (gp.u(j) > 0.0 ? dt*gp.u(j)/dx : -dt*gp.u(j)/dx);
            A(j,j+1) = -rhoAlphar*dt/(gp.rho(j)*dx*dx) +
                       (gp.u(j) > 0.0 ? 0.0 : dt*gp.u(j)/dx);
            b(j) = gp.hs(j) + dt*gp.qdot(j)/gp.rho(j);
        }
        A(0,0) = 1.0;
        A(input_data.nx-1,input_data.nx-1) = 1.0;
        b(0) = input_data.hsL;
        b(input_data.nx-1) = input_data.hsR;
        // Ignition
        if (input_data.ign && time > input_data.ignBEGt &&
            time < input_data.ignENDt) {
            A(input_data.nx/2,input_data.nx/2-1) = 0.0;
            A(input_data.nx/2,input_data.nx/2) = 1.0;
            A(input_data.nx/2,input_data.nx/2+1) = 0.0;
            b(input_data.nx/2) = input_data.ignHs;
        }
        gp.hs = tdma(A,b);
        gas.calcT(gp.T, gp.Y, gp.hs);
        if (input_data.log) {
        std::cout << std::setw(WIDTH) << "T.max "
                  << std::setw(WIDTH) << gp.T.maxCoeff(&loc) << " @ position "
                  << loc << std::endl;
        }

        gp.rhoPrev = gp.rho;
        gas.updateThermo(gp.hs, gp.Y, Le, gp.rho,
                         gp.mu, gp.kappa, gp.alpha, gp.D);

        time += dt;
        // Adjustable time step according to chemical time scale
        std::cout << "Time =  " << time << std::setprecision(10) << std::endl;
        dt = std::min(dtChem, input_data.dtMax);
        if (time+tprecision > input_data.tEND) break;
        if (time+dt > input_data.tEND) dt = input_data.tEND - time;
        std::cout << std::endl;
    }
    std::cout << "End" << std::endl;
    endTime = std::clock();
    std::cout << "Run time   " << double(endTime - startTime) / CLOCKS_PER_SEC
              << std::setprecision(6) << " s" << std::endl;
    write(time, gas, x, gp);

    return 0;
}

void fill_input(const std::string fname)
{
    std::ifstream finp(fname);
    if (!finp) throw std::runtime_error("input.txt NOT FOUND");
    std::map<std::string, std::string> dict;
    std::string name;
    std::string value;
    while (finp >> name) {
        finp >> value;
        dict[name] = value;
    }
    // Discretize space and time
    input_data.nx = std::stoi(dict["nPoints"]);
    input_data.XBEG = std::stod(dict["XBEG"]);
    input_data.XEND = std::stod(dict["XEND"]);
    input_data.tBEG = std::stod(dict["tBEG"]);
    input_data.tEND = std::stod(dict["tEND"]);
    input_data.dtMax = std::stod(dict["dtMax"]);
    // BC
    input_data.ainit = std::stod(dict["strainRate"]);  // prescribed strain rate
    input_data.VL = input_data.ainit*0;
    input_data.VR = input_data.ainit*0;
    input_data.TI = std::stod(dict["TI"]);
    input_data.TL = std::stod(dict["TL"]);
    input_data.TR = std::stod(dict["TR"]);
    input_data.YO2Air = 0.23197;
    input_data.YN2Air = 0.75425;
    input_data.YARAir = 0.01378;
    std::string str;
    std::string temp = dict["fuelNames"];
    std::stringstream buf1(temp);
    while (std::getline(buf1, str, ',')) {
        input_data.FUELNAMES.push_back(str);
    }
    temp = dict["YFUELS"];
    std::stringstream buf2(temp);
    while (std::getline(buf2, str, ',')) {
        input_data.YFUELS.push_back(std::stod(str));
    }
    for (string s : input_data.FUELNAMES) {
        std::cout << s << std::endl;
    }
    for (double yi : input_data.YFUELS) {
        std::cout << yi << std::endl;
    }
    // IC
    input_data.rstr = (dict["restore"] == "true" ? true : false);
    input_data.ign = (dict["ignition"] == "true" ? true : false);
    input_data.strain = (dict["strain"] == "true" ? true : false);
    input_data.log = (dict["log"] == "true" ? true : false);
    input_data.ignBEGt = std::stod(dict["ignBEGt"]);
    input_data.ignENDt = std::stod(dict["ignENDt"]);
    input_data.ignHs = std::stod(dict["ignHs"]);
    input_data.aBEG = std::stod(dict["aBEG"]);
    input_data.aEND = std::stod(dict["aEND"]);
    input_data.writeFreq = std::stoi(dict["writeFreq"]);
    input_data.rhoInf = std::stod(dict["rhoInf"]);
    input_data.p0 = std::stod(dict["pressure"]);
    input_data.lrRatio = std::stod(dict["lrRatio"]);
}

void restore(const ChemThermo& gas, const Eigen::VectorXd& x, gas_phase& gp)
{
    std::vector<double> x0;
    std::vector<double> u0;
    std::vector<double> V0;
    std::vector<double> T0;
    std::vector<std::vector<double> > Y0(gas.nsp());
    std::ifstream fin("initial_solution.csv");
    std::string line, str;
    std::getline(fin, line);
    if ((std::count(line.begin(),line.end(),',')+1) != gas.nsp()+6)
        throw std::runtime_error("Input x,u,V,rho,D,T,Y0,...,Ynsp-1");
    while (std::getline(fin, line)) {
        int n = 0;
        std::istringstream buffer(line);
        while(std::getline(buffer, str, ',')) {
            switch (n) {
            case 0:
                x0.push_back(std::stold(str));
                break;
            case 1:
                u0.push_back(std::stold(str));
                break;
            case 2:
                V0.push_back(std::stold(str));
                break;
            case 3:
                break;
            case 4:
                break;
            case 5:
                T0.push_back(std::stold(str));
                break;
            default:
                Y0[n-6].push_back(std::stold(str));
                break;
            }
            ++n;
        }
    }

    for (int j=0; j<x.size(); j++) {
        gp.u(j) = lininterp(x(j), x0, u0);
        gp.V(j) = lininterp(x(j), x0, V0);
        gp.T(j) = lininterp(x(j), x0, T0);
        for (int k=0; k<gas.nsp(); k++) {
            gp.Y[k](j) = lininterp(x(j), x0, Y0[k]);
        }
        double y[gas.nsp()];
        gas.massFractions(gp.Y, y, j);
        gp.hs(j) = gas.calcHs(gp.T(j), y);
    }
}

void write(const double iter, const ChemThermo& gas,
           const Eigen::VectorXd& x, const gas_phase& gp)
{
    std::stringstream ss;
    ss << iter;
    std::ofstream fout("data/output-"+ss.str()+".csv");
    std::ofstream rout("data/reaction-"+ss.str()+".csv");
    // Output
    fout << "x (m),u (m/s),V (1/s),rho (kg/m3),D (m2/s2),T (K)";
    for (int k=0; k<gas.nsp(); k++) {
        fout << "," << gas.speciesName(k);
    }
    fout << std::endl;
    for (int j=0; j<x.size(); j++) {
        fout << std::setprecision(8) << x(j) << ","
             << std::setprecision(8) << gp.u(j) << ","
             << std::setprecision(8) << gp.V(j) << ","
             << std::setprecision(8) << gp.rho(j) << ","
             << std::setprecision(8) << gp.D(j) << ","
             << std::setprecision(8) << gp.T(j);
        for (int k=0; k<gas.nsp(); k++) {
            fout << "," << std::setprecision(8) << gp.Y[k](j);
        }
        fout << std::endl;
    }
    // Output reactions related quantities (source terms)
    rout << "x (m),Qdot (J/m3 s)";
    for (int k = 0; k < gas.nsp(); k++) {
        rout << "," << gas.speciesName(k);
    }
    rout << std::endl;
    for (int j=0; j<x.size(); j++) {
        rout << std::setprecision(6) << x(j) << ","
             << std::setprecision(6) << gp.qdot(j);
        for (int k = 0; k < gas.nsp(); k++) {
            rout << "," << std::setprecision(6) << gp.wdot[k](j);
        }
        rout << std::endl;
    }
}
