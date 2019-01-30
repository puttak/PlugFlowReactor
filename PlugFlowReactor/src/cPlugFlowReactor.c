//! @file cPlugFlowReactor.c

// TODO add Cantera licence when contributing.
// TODO add Cantera licence when contributing.

#include "cantera/onedim/PlugFlowReactor.h"


int CXX_PlugFlowReactor(
    const std::string & mech,
    const std::string & phase,
    const double T,
    const double P,
    const std::string & X,
    const double D,
    const double L,
    const double mdot,
    const double htc,
    const std::function<double(double)>& Tw,
    const std::string & saveas,
    const double dx,
    const double rtol,
    const double atol,
    const unsigned maxsteps
)
{
    try
    {
        std::cout << "--- CXX_PlugFlowReactor ---" << std::endl;

        Cantera::PlugFlowReactor pfr;
        pfr.setMechanism(mech, phase);
        pfr.setState_TPX(T, P, X);
        pfr.setGeometry(D, L);
        pfr.setMassFlowRate(mdot);
        pfr.setTemperatureProfile(htc, Tw);
        //pfr.setViscosity(mu);
        pfr.init();

        Cantera::IDA_Solver solver {pfr};
        solver.setTolerances(rtol, atol);
        solver.setMaxNumSteps(maxsteps);
        solver.setJacobianType(0);
        solver.setDenseLinearSolver();
        solver.setStopTime(L);
        solver.init(0.0);

        double x = 0.0;
        pfr.feedResults(x, &solver);

        while (x < L)
        {
            x = x + std::min(L - x, dx);
            // TODO Manage return codes from IDA_Solver.solve.
            int ret = solver.solve(x);
            pfr.feedResults(x, &solver);
        }

        pfr.writeResults(saveas);
    }
    catch (std::exception &err)
    {
        std::cerr << err.what() << std::endl;
        return -1;
    }

    return 0;
}


extern "C"
{
    typedef double (*HTCFunc_t)(double);

    int C_PlugFlowReactor(
        const char mech[],
        const char phase[],
        const double T,
        const double P,
        const char X[],
        const double D,
        const double L,
        const double mdot,
        const double htc,
        const HTCFunc_t Tw,
        const char saveas[],
        const double dx,
        const double rtol,
        const double atol,
        const unsigned maxsteps
    )
    {
        std::string cxx_mech {mech};
        std::string cxx_phase {phase};
        std::string cxx_X {X};
        std::string cxx_saveas {saveas};

        std::cout << "\nCXX_PlutFlowReactor C-interface"
                  << "\nMechanism .............. " << cxx_mech
                  << "\nPhase name ............. " << cxx_phase
                  << "\nInlet composition ...... " << cxx_X
                  << "\nInlet temperature ...... " << T
                  << "\nInlet pressure ......... " << P
                  << "\nInlet mass flow rate ... " << mdot
                  << "\nWall HTC ............... " << htc
                  << "\nInlet wall temperature . " << Tw(0.0)
                  << "\nReactor diameter ....... " << D
                  << "\nReactor length ......... " << L
                  << "\nSave step .............. " << dx
                  << "\nRelative tolerance ..... " << rtol
                  << "\nAbsolute tolerance ..... " << atol
                  << "\nMaximum no. of steps ... " << maxsteps
                  << "\nSaving results to ...... " << cxx_saveas
                  << std::endl;

        return CXX_PlugFlowReactor(cxx_mech, cxx_phase, T, P, cxx_X, D, L, mdot,
            htc, Tw, cxx_saveas, dx, rtol, atol, maxsteps);
    }
} // extern "C"
