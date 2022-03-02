/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

// General includes common to most GR problems
#include "ScalarFieldLevel.hpp"
#include "BoxLoops.hpp"
#include "NanCheck.hpp"
#include "PositiveChiAndAlpha.hpp"
#include "SixthOrderDerivatives.hpp"
#include "TraceARemoval.hpp"
#include "AMRReductions.hpp"
#include "ChiExtractionTaggingCriterion.hpp"
#include "ChiPunctureExtractionTaggingCriterion.hpp"
#include "SetValue.hpp"
#include "SmallDataIO.hpp"
#include "TraceARemoval.hpp"
#include "TwoPuncturesInitialData.hpp"
#include "Weyl4.hpp"
#include "WeylExtraction.hpp"
#include "Oscilloton.hpp"
#include "MatterEnergy.hpp"
#include "FluxExtraction.hpp"

// For RHS update
#include "MatterCCZ4RHS.hpp"

// For constraints calculation
#include "NewMatterConstraints.hpp"

// For tag cells
#include "RhoAndKExtractionTaggingCriterion.hpp"

// Problem specific includes
#include "ComputePack.hpp"
#include "GammaCalculator.hpp"
#include "InitialScalarData.hpp"
#include "KerrBH.hpp"
#include "Potential.hpp"
#include "ScalarField.hpp"
#include "SetValue.hpp"
#include "MatterEnergyDensity.hpp"
#include "MovingPunctureGauge.hpp"

// Things to do at each advance step, after the RK4 is calculated
void ScalarFieldLevel::specificAdvance()
{
    // Enforce trace free A_ij and positive chi and alpha
    BoxLoops::loop(
        make_compute_pack(TraceARemoval(),
                          PositiveChiAndAlpha(m_p.min_chi, m_p.min_lapse)),
        m_state_new, m_state_new, INCLUDE_GHOST_CELLS);

    // Check for nan's
    if (m_p.nan_check)
        BoxLoops::loop(NanCheck(), m_state_new, m_state_new,
                       EXCLUDE_GHOST_CELLS, disable_simd());
}

// Initial data for field and metric variables
void ScalarFieldLevel::initialData()
{
    CH_TIME("ScalarFieldLevel::initialData");
    if (m_verbosity)
        pout() << "ScalarFieldLevel::initialData " << m_level << endl;

    // First set everything to zero then initial conditions for scalar field -
    // here a Kerr BH and a scalar field profile
    
    double spacing = 0.01;
    BoxLoops::loop(make_compute_pack(SetValue(0.0), Oscilloton(m_p.initial_params, m_dx, spacing)),
                   m_state_new, m_state_new, FILL_GHOST_CELLS, disable_simd());

    fillAllGhosts();
    BoxLoops::loop(GammaCalculator(m_dx), m_state_new, m_state_new,
                   EXCLUDE_GHOST_CELLS);

    // compute rho
    Potential potential(m_p.potential_params);
    ScalarFieldWithPotential scalar_field(potential);
    fillAllGhosts();
    BoxLoops::loop(
            MatterEnergyDensity<ScalarFieldWithPotential>(scalar_field, m_dx),
            m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);

}

#ifdef CH_USE_HDF5
// Things to do before outputting a checkpoint file
void ScalarFieldLevel::prePlotLevel()
{
    fillAllGhosts();
    Potential potential(m_p.potential_params);
    ScalarFieldWithPotential scalar_field(potential);
    BoxLoops::loop(
        MatterConstraints<ScalarFieldWithPotential>(
            scalar_field, m_dx, m_p.G_Newton, c_Ham, Interval(c_Mom, c_Mom)),
        m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);
}
#endif

// Things to do in RHS update, at each RK4 step
void ScalarFieldLevel::specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                       const double a_time)
{
    // Enforce trace free A_ij and positive chi and alpha
    BoxLoops::loop(
        make_compute_pack(TraceARemoval(),
                          PositiveChiAndAlpha(m_p.min_chi, m_p.min_lapse)),
        a_soln, a_soln, INCLUDE_GHOST_CELLS);

    // Calculate MatterCCZ4 right hand side with matter_t = ScalarField
    Potential potential(m_p.potential_params);
    ScalarFieldWithPotential scalar_field(potential);
    if (m_p.max_spatial_derivative_order == 4)
    {
        MatterCCZ4RHS<ScalarFieldWithPotential, MovingPunctureGauge,
                      FourthOrderDerivatives>
            my_ccz4_matter(scalar_field, m_p.ccz4_params, m_dx, m_p.sigma,
                           m_p.formulation, m_p.G_Newton);
        BoxLoops::loop(my_ccz4_matter, a_soln, a_rhs, EXCLUDE_GHOST_CELLS);
    }
    else if (m_p.max_spatial_derivative_order == 6)
    {
        MatterCCZ4RHS<ScalarFieldWithPotential, MovingPunctureGauge,
                      SixthOrderDerivatives>
            my_ccz4_matter(scalar_field, m_p.ccz4_params, m_dx, m_p.sigma,
                           m_p.formulation, m_p.G_Newton);
        BoxLoops::loop(my_ccz4_matter, a_soln, a_rhs, EXCLUDE_GHOST_CELLS);
    }
}

// Things to do at ODE update, after soln + rhs
void ScalarFieldLevel::specificUpdateODE(GRLevelData &a_soln,
                                         const GRLevelData &a_rhs, Real a_dt)
{  // compute rho
    Potential potential(m_p.potential_params);
    ScalarFieldWithPotential scalar_field(potential);
    fillAllGhosts();
    BoxLoops::loop(
        MatterEnergyDensity<ScalarFieldWithPotential>(scalar_field, m_dx),
        m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);

    // Enforce trace free A_ij
    BoxLoops::loop(TraceARemoval(), a_soln, a_soln, INCLUDE_GHOST_CELLS);
}

void ScalarFieldLevel::preTagCells()
{
    // we don't need any ghosts filled for the fixed grids tagging criterion
    // used here so don't fill any
}

void ScalarFieldLevel::computeTaggingCriterion(FArrayBox &tagging_criterion,
                                               const FArrayBox &current_state)
{
    BoxLoops::loop(
            RhoAndKExtractionTaggingCriterion(m_dx, m_level, m_p.threshold_rho, m_p.threshold_K, m_p.extraction_params, m_p.activate_extraction, m_p.max_matter_level),
        current_state, tagging_criterion);
}


void ScalarFieldLevel::specificPostTimeStep()
{

      double coarsest_dt = m_p.coarsest_dx * m_p.dt_multiplier;

    const double remainder = fmod(m_time, coarsest_dt);
    if (min(abs(remainder), abs(remainder - coarsest_dt)) < 1.0e-8)
    {
        // calculate the density of the PF, but excise the BH region completely
        fillAllGhosts();
        Potential potential(m_p.potential_params);
        ScalarFieldWithPotential scalar_field(potential);

        BoxLoops::loop(MatterEnergy<ScalarFieldWithPotential>(scalar_field,
                                                              m_dx, m_p.center),
                       m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);

    }

    CH_TIME("ScalarFieldLevel::specificPostTimeStep");

    bool first_step =
            (m_time == 0.); // this form is used when 'specificPostTimeStep' was
    // called during setup at t=0 from Main
    // bool first_step = (m_time == m_dt); // if not called in Main

    // write out the integral after each coarse timestep
    if (m_level == 0)
    {
        bool first_step = (m_time == m_dt);

        // integrate the densities and write to a file
        AMRReductions<VariableType::diagnostic> amr_reductions(m_gr_amr);

        double rho2_sum = amr_reductions.sum(c_rho2);
        double source2_sum = amr_reductions.sum(c_source2);

        SmallDataIO integral_file("VolumeIntegrals", m_dt, m_time,
                                  m_restart_time, SmallDataIO::APPEND,
                                  first_step);
        // remove any duplicate data if this is post restart
        integral_file.remove_duplicate_time_data();
        std::vector<double> data_for_writing = { rho2_sum, source2_sum};
        // write data
        if (first_step)
        {
            integral_file.write_header_line({ "rho2", "source2"});
        }
        integral_file.write_time_data_line(data_for_writing);

        // Now refresh the interpolator and do the interpolation
        bool fill_ghosts = false;
        m_gr_amr.m_interpolator->refresh(fill_ghosts);
        m_gr_amr.fill_multilevel_ghosts(VariableType::diagnostic,
                                        Interval (c_flux2,c_flux2));
        FluxExtraction my_extraction(m_p.extraction_params, m_dt, m_time,
                                     m_restart_time);
        my_extraction.execute_query(m_gr_amr.m_interpolator);
    }

    if (m_p.activate_extraction == 1) {
        int min_level = m_p.extraction_params.min_extraction_level();
        bool calculate_weyl = at_level_timestep_multiple(min_level);
        if (calculate_weyl) {
            // Populate the Weyl Scalar values on the grid
            fillAllGhosts();
            BoxLoops::loop(
                    Weyl4(m_p.extraction_params.center, m_dx, m_p.formulation),
                    m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);


            // Do the extraction on the min extraction level
            if (m_level == min_level) {
                CH_TIME("WeylExtraction");
                // Now refresh the interpolator and do the interpolation
                // fill ghosts manually to minimise communication
                bool fill_ghosts = false;



                m_gr_amr.m_interpolator->refresh(fill_ghosts);



                m_gr_amr.fill_multilevel_ghosts(
                        VariableType::diagnostic, Interval(c_Weyl4_Re, c_Weyl4_Im),
                        min_level);



                WeylExtraction my_extraction(m_p.extraction_params, m_dt,
                                             m_time, first_step,
                                             m_restart_time);


                my_extraction.execute_query(m_gr_amr.m_interpolator);


            }
        }
    }
}