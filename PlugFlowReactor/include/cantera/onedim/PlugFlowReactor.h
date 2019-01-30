//! @file PlugFlowReactor.h

// TODO add Cantera licence when contributing.
// TODO add Cantera licence when contributing.

#ifndef CT_PLUGFLOWREACTOR_H
#define CT_PLUGFLOWREACTOR_H

#include <chrono>
#include <functional>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <string>
#include <sstream>
#include <vector>
#include "cantera/IdealGasMix.h"
#include "cantera/numerics/IDA_Solver.h"
#include "cantera/numerics/eigen_dense.h"
#include "cantera/transport.h"

#ifndef SUPPRESS_WARNINGS
#define SUPPRESS_WARNINGS true
#endif

namespace Cantera
{

/**
 * Generig plug-flow reactor with surface chemistry.
 * TODO inherit from ReactorBase too.
 * TODO replace function by Func1.
 * TODO replace vector double by vectorfp.
 * TODO replace runtime_error with CanteraError.
 * @ingroup onedim
 */
class PlugFlowReactor : public ResidJacEval
{
public:
    PlugFlowReactor(doublereal atol=1.0e-13) : ResidJacEval{atol} {}

    ~PlugFlowReactor()
    {
        if (m_gas != nullptr) delete m_gas;
        if (m_trn != nullptr) delete m_trn;
        // TODO Should I call this?
        appdelete();
    }

    //! Set mechanism file and phase.
    void setMechanism(const std::string &mech, const std::string &phase) {
        m_mech = mech;
        m_phase = phase;
        has_mech = true;
    }

    //! Set the temperature (K), pressure (Pa), and mole fractions.
    /*!
     * TODO: understand the following from original class:
     * Note, the mole fractions are set first before the pressure is set.
     * Setting the pressure may involve the solution of a nonlinear equation.
     *
     * @param t    Temperature (K)
     * @param p    Pressure (Pa)
     * @param x    String containing a composition map of the mole fractions.
     *             Species not in the composition map are assumed to have zero
     *             mole fraction
     */
    void setState_TPX(doublereal t, doublereal p, const std::string& x) {
        m_t0 = t;
        m_p0 = p;
        m_x0 = x;
        has_state = true;
    }

    void setGeometry(const doublereal d, const doublereal l) {
        m_diam = d;
        m_length = l;
        has_geo = true;

        m_ac = Cantera::Pi * m_diam * m_diam / 4.0;
        m_poa = 4.0 / m_diam;
    }

    virtual void setMassFlowRate(const doublereal mdot) {
        m_mdot = mdot;
        has_mdot = true;
    }

    virtual void setTemperatureProfile(const doublereal htc,
        const doublereal tw)
    {
        m_htc = htc;
        m_tw = [&](double x){ return tw; };
        has_profile = true;
    }

    virtual void setTemperatureProfile(const doublereal htc,
        std::function<doublereal(doublereal)> tw)
    {
        m_htc = htc;
        m_tw = tw;
        has_profile = true;
    }

    void setViscosity(const doublereal mu) {
        m_mu = [&](){ return mu; };
        has_viscosity = true;
    }

    void setViscosity(std::function<doublereal()> const& mu) {
        m_mu = mu;
        has_viscosity = true;
    }

    unsigned getSpeciesIndex(std::string name) const {
        return m_gas->kineticsSpeciesIndex(name);
    }

    doublereal getIntEnergyMass() const {
        return m_gas->intEnergy_mass();
    }

    std::vector<std::string> variablesNames() const {
        return m_var;
    }

    void init() {
        // Avoid verbosity with this reactor.
        suppress_thermo_warnings(SUPPRESS_WARNINGS);

        std::cout << std::boolalpha
                  << "\nIntegrating PFR"
                  << "\nUsing Sundials : " << CT_SUNDIALS_VERSION
                  << "\nUsing LAPACK   : " << bool(CT_SUNDIALS_USE_LAPACK)
                  << std::endl;

        if (has_init) {
            throw std::runtime_error("Reactor already initialized.");
        }
        if (!has_mech) {
            throw std::runtime_error("First call `setMechanism`.");
        }
        if (!has_state) {
            throw std::runtime_error("First call `setState_TPX`.");
        }
        if (!has_geo) {
            throw std::runtime_error("First call `setGeometry`.");
        }
        if (!has_mdot) {
            throw std::runtime_error("First call `setMassFlowRate`.");
        }
        if (!has_profile) {
            throw std::runtime_error("First call `setTemperatureProfile`.");
        }

        // Create ideal gas phase object.
        m_gas = new IdealGasMix {m_mech, m_phase};
        m_gas->setState_TPX(m_t0, m_p0, m_x0);

        // Compute initial speed.
        m_u0 = m_mdot / (m_gas->density() * m_ac);

        if (!has_viscosity) {
            try
            {
                std::cout << "Getting transport manager..." << std::endl;
                m_trn = newDefaultTransportMgr(m_gas);
                m_mu = [&](){ return m_trn->viscosity(); };
                m_mu(); // Try to call it!
                std::cout << "Using mechanism viscosity..." << std::endl;
            }
            catch (CanteraError &err)
            {
                std::cerr << err.what() << std::endl;
                std::cout << "Assuming inviscid flow" << std::endl;
                m_mu = [&](){ return 0.0; };
            }
        }

        // Initialize indexers.
        nspec_gas_ = m_gas->nSpecies();
        nspec_sur_ = 0;
        neq_ = nspec_gas_ + nspec_sur_ + neqs_ext_;
        idx0 = nspec_gas_ + nspec_sur_ + 0;
        idx1 = nspec_gas_ + nspec_sur_ + 1;
        idx2 = nspec_gas_ + nspec_sur_ + 2;
        idx3 = nspec_gas_ + nspec_sur_ + 3;

        // TODO add nspec_sur_.
        m_mw.resize(nspec_gas_);
        m_wdot.resize(nspec_gas_);
        m_hbar.resize(nspec_gas_);
        m_gas->getMolecularWeights(m_mw);

        // Start variables names with species.
        m_var = m_gas->speciesNames();
        // TODO add surface species.
        m_var.push_back("u");
        m_var.push_back("rho");
        m_var.push_back("P");
        m_var.push_back("T");
        has_init = true;
    }

    void feedResults(const double x, Cantera::IDA_Solver *solver)
    {
        if (!has_ssbuffer)
        {
            for (auto var : m_var) { m_ssbuffer << var << ","; }
            m_ssbuffer << "x" << std::endl;
            has_ssbuffer = true;
        }

        for (unsigned i = 0; i != neq_; ++i)
        {
            m_ssbuffer << solver->solution(i) << ",";
        }
        m_ssbuffer << x << std::endl;
    }

    void writeResults(const std::string & saveas)
    {
        // TODO test if m_ssbuffer is not empty.
        std::ofstream ofs(saveas, std::ios::out);
        if (!ofs) { throw std::runtime_error("Cannot write to " + saveas); }
        ofs << m_ssbuffer.str();
        ofs.close();
    }

    int getInitialConditions(const doublereal t0,
                             doublereal *const y,
                             doublereal *const ydot);

    int evalResidNJ(const doublereal t,
                    const doublereal delta_t,
                    const doublereal* const y,
                    const doublereal* const ydot,
                    doublereal* const resid,
                    const ResidEval_Type_Enum evalType = Base_ResidEval,
                    const int id_x = -1,
                    const doublereal delta_x = 0.0);

protected:

    // ---------------------

    //! Number of extra equations.
    static constexpr unsigned neqs_ext_ = 4;

    //! Number of gas phase species.
    unsigned nspec_gas_;

    //! Number of surface species.
    unsigned nspec_sur_;

    //! Index of extra equations.
    unsigned idx0 = 0, idx1 = 0, idx2 = 0, idx3 = 0;

    //! Global initialization check.
    bool has_init = false;

    // ---------------------

    //! Mechanism path.
    std::string m_mech;

    //! Mechanism phase.
    std::string m_phase;

    //! Mechanism check.
    bool has_mech = false;

    // ---------------------

    //! Inlet gas temperature [K].
    doublereal m_t0;

    //! Inlet gas pressure [Pa].
    doublereal m_p0;

    //! Inlet gas composition [mole fractions].
    std::string m_x0;

    //! State check.
    bool has_state = false;

    // ---------------------

    //! Reactor diameter [m].
    doublereal m_diam;

    //! Reactor length [m].
    doublereal m_length;

    //! Reactor cross-sectional area.
    doublereal m_ac;

    //! Ratio perimeter per cross section.
    doublereal m_poa;

    //! Geometry check.
    bool has_geo = false;

    // ---------------------

    //! Inlet mass flow rate.
    doublereal m_mdot;

    //! Inlet velocity.
    doublereal m_u0;

    //! Flow rate check.
    bool has_mdot = false;

    // ---------------------

    //! Global heat transfer coefficient.
    doublereal m_htc = 0.0;

    //! Wall temperature in terms of position.
    std::function<doublereal(doublereal)> m_tw;

    //! Wall profile check.
    bool has_profile = false;

    // ---------------------

    //! Viscosity function interface.
    std::function<doublereal()> m_mu;

    //! Viscosity check.
    bool has_viscosity = false;

    // ---------------------

    //! Pointer to the gas phase object.
    IdealGasMix *m_gas = nullptr;

    //! Pointer to the transport manager.
    Transport *m_trn = nullptr;

    // ---------------------

    //! Species molar weights.
    std::vector<doublereal> m_mw;

    //! Species net production rates.
    std::vector<doublereal> m_wdot;

    //! Species molar enthalpies.
    std::vector<doublereal> m_hbar;

    //! Species names and variables.
    std::vector<std::string> m_var;

    // ---------------------

    //! Buffer for results file.
    std::stringstream m_ssbuffer;

    //! Results buffer check.
    bool has_ssbuffer = false;

    // ---------------------

    //! Pressure drop model of viscous loss.
    virtual const doublereal viscousLoss(doublereal u) const {
        // TODO Re number test.
        return 8 * m_mu() * u * Cantera::Pi / m_ac;
    }

    //! Wall heat exchange term.
    virtual const doublereal wallExchange(doublereal x, doublereal T) const {
        return m_htc * m_poa * (m_tw(x) - T);
    }

}; // (class CanteraPFR)

}

int Cantera::PlugFlowReactor::getInitialConditions(const doublereal t0,
    doublereal *const y, doublereal *const ydot)
{
    std::cout << "Finding initial conditions..." << std::endl;

    const doublereal T0 = m_gas->temperature();
    const doublereal P0 = m_gas->pressure();
    const doublereal rho0 = m_gas->density();
    const doublereal Wavg = m_gas->meanMolecularWeight();
    const doublereal RT = T0 * Cantera::GasConstant;
    const doublereal rho0R = rho0 * Cantera::GasConstant;
    const doublereal rhoUCp = rho0 * m_u0 * m_gas->cp_mass();
    doublereal hdot = 0;
    doublereal wall = wallExchange(t0, T0);

    m_gas->getMassFractions(y);
    m_gas->getNetProductionRates(&m_wdot[0]);
    m_gas->getPartialMolarEnthalpies(&m_hbar[0]);

    Eigen::MatrixXd A(neq_, neq_);
    Eigen::VectorXd b(neq_);

    y[idx0] = m_u0;
    y[idx1] = rho0;
    y[idx2] = P0;
    y[idx3] = T0;

    for (unsigned k = 0; k != idx0; ++k)
    {
        // Energy contribution.
        hdot += m_wdot[k] * m_hbar[k];

        // For species equations.
        A(k, k) = rho0 * m_u0;
        b(k) = m_wdot[k] * m_mw[k];

        // Yk' for other equations, exceptionally here!
        A(idx2, k) = P0 * Wavg * Wavg / m_mw[k];
    }

    // Continuity equation elements.
    A(idx0, idx0) = rho0;           // u'
    A(idx0, idx1) = m_u0;           // rho'
    A(idx0, idx2) = 0;              // p'
    A(idx0, idx3) = 0;              // T'

    // Momentum equation elements.
    A(idx1, idx0) = rho0 * m_u0;    // u'
    A(idx1, idx1) = 0;              // rho'
    A(idx1, idx2) = 1;              // p'
    A(idx1, idx3) = 0;              // T'

    // State equation elements.
    A(idx2, idx0) = 0;              // u'
    A(idx2, idx1) = RT;             // rho'
    A(idx2, idx2) = -Wavg;          // p'
    A(idx2, idx3) = rho0R;          // T'

    // Energy equation elements.
    A(idx3, idx0) = 0;              // u'
    A(idx3, idx1) = 0;              // rho'
    A(idx3, idx2) = 0;              // p'
    A(idx3, idx3) = rhoUCp;         // T'

    b(idx0) = 0;                    // RHS continuity
    b(idx1) = -viscousLoss(m_u0);   // RHS momentum
    b(idx2) = 0;                    // RHS state
    b(idx3) = -hdot + wall;         // RHS energy

    Eigen::VectorXd x = A.fullPivLu().solve(b);
    Eigen::VectorXd::Map(ydot, x.rows()) = x;

    //for (int i = 0; i < neq_; ++i) std::cout << ydot[i] << "\n";

    std::cout << std::scientific << std::showpos
              << "Initial conditions found..."
              << "\n u = " << y[idx0] << " u' = " << ydot[idx0]
              << "\n r = " << y[idx1] << " r' = " << ydot[idx1]
              << "\n p = " << y[idx2] << " p' = " << ydot[idx2]
              << "\n T = " << y[idx3] << " T' = " << ydot[idx3]
              << std::noshowpos << std::endl;

    return 0;
}

int Cantera::PlugFlowReactor::evalResidNJ(const doublereal t,
    const doublereal delta_t, const doublereal* const y,
    const doublereal* const ydot, doublereal* const resid,
    const Cantera::ResidEval_Type_Enum evalType, const int id_x,
    const doublereal delta_x)
{
    const doublereal u = y[idx0], dudz = ydot[idx0];
    const doublereal r = y[idx1], drdz = ydot[idx1];
    const doublereal p = y[idx2], dpdz = ydot[idx2];
    const doublereal T = y[idx3], dTdz = ydot[idx3];

    const doublereal Cp = m_gas->cp_mass();
    doublereal hdot = 0.0;

    m_gas->setMassFractions_NoNorm(y);
    m_gas->setState_TP(T, p);
    m_gas->getNetProductionRates(&m_wdot[0]);
    m_gas->getPartialMolarEnthalpies(&m_hbar[0]);

    // TODO change index when adding surface processes.
    for (unsigned k = 0; k != idx0; ++k)
    {
        resid[k] = u * r * ydot[k] - m_wdot[k] * m_mw[k];
        hdot += m_wdot[k] * m_hbar[k];
    }

    resid[idx0] = r * dudz + u * drdz;
    resid[idx1] = u * r * dudz + dpdz + viscousLoss(u);
    resid[idx2] = m_gas->density() - r;
    resid[idx3] = r * u * Cp * dTdz + hdot - wallExchange(t, T);

    return 0;
}

#endif
//
// //! Set the temperature (K), pressure (Pa), and mole fractions.
// /*!
//  * Note, the mole fractions are set first before the pressure is set.
//  * Setting the pressure may involve the solution of a nonlinear equation.
//  *
//  * @param t    Temperature (K)
//  * @param p    Pressure (Pa)
//  * @param x    Vector of mole fractions.
//  *             Length is equal to m_kk.
//  */
// virtual void setState_TPX(doublereal t, doublereal p, const doublereal* x);
//
// //! Set the temperature (K), pressure (Pa), and mole fractions.
// /*!
//  * Note, the mole fractions are set first before the pressure is set.
//  * Setting the pressure may involve the solution of a nonlinear equation.
//  *
//  * @param t    Temperature (K)
//  * @param p    Pressure (Pa)
//  * @param x    Composition map of mole fractions. Species not in
//  *             the composition map are assumed to have zero mole fraction
//  */
// virtual void setState_TPX(doublereal t, doublereal p, const compositionMap& x);
