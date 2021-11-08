/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP_
#define SIMULATIONPARAMETERS_HPP_

// General includes
#include "GRParmParse.hpp"
#include "SimulationParametersBase.hpp"

// Problem specific includes:
#include "InitialScalarData.hpp"
#include "KerrBH.hpp"
#include "Potential.hpp"
#include "Oscilloton.hpp"

class SimulationParameters : public SimulationParametersBase
{
  public:
    SimulationParameters(GRParmParse &pp) : SimulationParametersBase(pp)
    {
        // read the problem specific params
        read_params(pp);
        check_params();
    }

    void read_params(GRParmParse &pp)
    {
        // Initial scalar field data

        pp.load("G_Newton", G_Newton,
                0.0); // for now the example neglects backreaction

        pp.load("scalar_mass", potential_params.scalar_mass, 0.1);
        pp.load("threshold_rho", threshold_rho, 0.1); // this param for regridding
        pp.load("threshold_K", threshold_K, 0.1);  // this param for regridding

        // Initial Kerr data

        pp.load("amplitudeSF", initial_params.amplitudeSF);
        pp.load("widthSF", initial_params.widthSF);
        pp.load("r_zero", initial_params.r_zero);

        // set the velocity
        pp.get("vx",vx);
        pp.get("vx2",vx2);

        //Boost Parameters
        initial_params.vx = vx;
        initial_params.vx2 = vx2;



        for(int i=0; i < CH_SPACEDIM; i++)
        {
            pp.get("centerSF", centerSF[i], i);
        }
        for(int i=0; i < CH_SPACEDIM; i++)
        {
            pp.get("centerSF2", centerSF2[i], i);
        }
        
        //set the center
        initial_params.centerSF = centerSF;
        initial_params.centerSF2 = centerSF2;
    }

    void check_params()
    {
        warn_parameter("scalar_mass", potential_params.scalar_mass,
                       potential_params.scalar_mass <
                           0.2 / coarsest_dx / dt_multiplier,
                       "oscillations of scalar field do not appear to be "
                       "resolved on coarsest level");
    }

   // this is the Oscilloton params
    Oscilloton::params_t initial_params;
    std::array<double, CH_SPACEDIM> centerSF;
    std::array<double, CH_SPACEDIM> centerSF2;


    // Initial data for matter and potential and BH
    double G_Newton;
    Potential::params_t potential_params;



    double threshold_rho;
    double threshold_K;

    double vx;
    double vx2;
};

#endif /* SIMULATIONPARAMETERS_HPP_ */
