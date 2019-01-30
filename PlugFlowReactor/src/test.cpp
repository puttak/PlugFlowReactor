//! @file AppPFR.cpp

// TODO add Cantera licence when contributing.
// TODO add Cantera licence when contributing.

#include "cantera/onedim/PlugFlowReactor.h"


void integrate(double L, Cantera::PlugFlowReactor &pfr)
{
    // Solver parameters.
    double rtol = 1.0e-02;
    double atol = 1.0e-15;
    double initstep = 1.0e-06;
    unsigned maxsteps = 100000;

    Cantera::IDA_Solver solver {pfr};
    solver.setTolerances(rtol, atol);
    solver.setMaxNumSteps(maxsteps);
    solver.setJacobianType(0);
    solver.setDenseLinearSolver();
    // solver.setInitialStepSize(initstep);
    solver.setStopTime(L);
    solver.init(0.0);

    double x = 0.0, dx = 0.010;

    while (x < L)
    {
        x = x + std::min(L-x, dx);
        solver.solve(x);
    }

    size_t space = 15;
    size_t id0 = pfr.getSpeciesIndex("C2H2");
    size_t id1 = pfr.getSpeciesIndex("H2");
    size_t neq = pfr.nEquations();

    for (int i = 0; i != 89; ++i) std::cout << "=";
    std::cout << std::scientific << "\n" << std::right
              << std::setw(space) << "C2H2 | "
              << std::setw(space) << "H2 | "
              << std::setw(space) << "u | "
              << std::setw(space) << "rho | "
              << std::setw(space) << "p | "
              << std::setw(space) << "T | "
              << "\n"
              << solver.solution(id0) << " | "
              << solver.solution(id1) << " | "
              << solver.solution(neq-4) << " | "
              << solver.solution(neq-3) << " | "
              << solver.solution(neq-2) << " | "
              << solver.solution(neq-1) << " | "
              << std::endl << std::endl;
}

// double get_mdot(std::string mech, std::string phase, double Q, std::string X)
// {
//     Cantera::IdealGasMix gas {mech, phase};
//     gas.setState_TPX(273.15, 101325.0, X);
//     return gas.density() * (Q / 6.0e+07);
// }

void test_gas()
{
    // Mechanism and phase.
    const std::string mech = "data/CT-hydrocarbon-dalmazsi-2017-mech.cti";
    const std::string phase = "gas";

    // Inlet gas state.
    const double T = 1173.0;
    const double P = 5000.0;
    const std::string X = "N2:0.64, C2H2:0.3528, CH3COCH3:6.48e-03, CH4:7.2e-04";

    // Reactor geometry.
    const double D = 0.028;
    const double L = 0.400;

    // Inlet flow rate.
    const double Q = 222.0;

    // Inlet mass flow rate.
    Cantera::IdealGasMix gas {mech, phase};
    gas.setState_TPX(273.15, 101325.0, X);
    const double mdot = gas.density() * (Q / 6.0e+07);
    //const double mdot = get_mdot(mech, phase, Q, X);

    try
    {
        std::cout << "--- TEST 0 ---" << std::endl;

        Cantera::PlugFlowReactor pfr;
        pfr.setMechanism(mech, phase);
        pfr.setState_TPX(T, P, X);
        pfr.setGeometry(D, L);
        pfr.setMassFlowRate(mdot);
        pfr.setTemperatureProfile(0.0, T);
        pfr.setViscosity(3.965315e-05);
        pfr.init();

        integrate(L, pfr);
    }
    catch (std::exception &err)
    {
        std::cerr << err.what() << std::endl;
    }

    try
    {
        std::cout << "--- TEST 1 ---" << std::endl;

        Cantera::PlugFlowReactor pfr;
        pfr.setMechanism(mech, phase);
        pfr.setState_TPX(T, P, X);
        pfr.setGeometry(D, L);
        pfr.setMassFlowRate(mdot);
        pfr.setTemperatureProfile(20.0, T);
        pfr.init();

        integrate(L, pfr);
    }
    catch (std::exception &err)
    {
        std::cerr << err.what() << std::endl;
    }

    try
    {
        std::cout << "--- TEST 2 ---" << std::endl;

        const double Ta = 300.0;
        const double Tc = 1173.0;
        const double Ts = 400.0;
        const double x1 = 0.02492942, x2 = 0.40810172;
        const double m1 = 0.78913918, m2 = 11.91548263;

        std::function<double(double)> Tw = [&](double x) {
            double term1 = 1 - std::exp(-std::pow(x / x1, m1));
            double term2 = 1 - std::exp(-std::pow(x / x2, m2));
            double wallT = Ta + (Tc - Ta) * term1 - (Tc - Ts) * term2;
            return 0.97 * wallT;
        };

        Cantera::PlugFlowReactor pfr;
        pfr.setMechanism(mech, phase);
        pfr.setState_TPX(Ta, P, X);
        pfr.setGeometry(D, L);
        pfr.setMassFlowRate(mdot);
        pfr.setTemperatureProfile(20.0, Tw);
        pfr.init();

        integrate(L, pfr);
    }
    catch (std::exception &err)
    {
        std::cerr << err.what() << std::endl;
    }
}

void test_sur()
{

    return;
}

int main()
{
    test_gas();
    test_sur();
    return 0;
}

//Cantera::Transport *trn = Cantera::newDefaultTransportMgr(&gas);
//std::cout << trn->viscosity() << std::endl;
