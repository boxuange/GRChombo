/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef RHOANDKEXTRACTIONTAGGINGCRITERION_HPP_
#define RHOANDKEXTRACTIONTAGGINGCRITERION_HPP_

#include "Cell.hpp"
#include "DimensionDefinitions.hpp"
#include "FourthOrderDerivatives.hpp"
#include "Tensor.hpp"

class RhoAndKExtractionTaggingCriterion
{
  protected:
    const double m_dx;
    const int m_level;
    const FourthOrderDerivatives m_deriv;
    const double m_threshold_rho;
    const double m_threshold_K;
    const SphericalExtraction::params_t m_params;
    const bool m_activate_extraction;
    const int m_max_matter_level;

  public:
    RhoAndKExtractionTaggingCriterion(double dx, const int a_level, double threshold_rho,
                                      double threshold_K, const SphericalExtraction::params_t a_params, const bool activate_extraction = false, int a_max_matter_level =1000)

        : m_dx(dx), m_level(a_level), m_deriv(dx), m_threshold_rho(threshold_rho),
          m_threshold_K(threshold_K), m_params(a_params), m_activate_extraction(activate_extraction), m_max_matter_level(a_max_matter_level) {};


    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        Tensor<1, data_t> d1_rho;
        FOR(idir) m_deriv.diff1(d1_rho, current_cell, idir, c_rho);

        Tensor<1, data_t> d1_K;
        FOR(idir) m_deriv.diff1(d1_K, current_cell, idir, c_K);

        data_t mod_d1_rho = 0;
        data_t mod_d1_K = 0;
        FOR(idir)
        {

            if(m_level < m_max_matter_level){mod_d1_rho += d1_rho[idir] * d1_rho[idir];}
            mod_d1_K += d1_K[idir] * d1_K[idir];
        }

        data_t criterion = m_dx * (sqrt(mod_d1_rho) / m_threshold_rho +
                                   sqrt(mod_d1_K) / m_threshold_K);
        

        if (m_activate_extraction)
        {
            for (int iradius = 0; iradius < m_params.num_extraction_radii;
                 ++iradius)
            {
                // regrid if within extraction level and not at required
                // refinement
                if (m_level < m_params.extraction_levels[iradius])
                {
                    const Coordinates<data_t> coords(current_cell, m_dx,
                                                     m_params.center);
                    const data_t r = coords.get_radius();
                    // add a 20% buffer to extraction zone so not too near to
                    // boundary
                    auto regrid = simd_compare_lt(
                            r, 1.2 * m_params.extraction_radii[iradius]);
                    criterion = simd_conditional(regrid, 100.0, criterion);
                }
            }
        }
        // Write back into the flattened Chombo box
        current_cell.store_vars(criterion, 0);
    }
};

#endif /* RHOANDKEXTRACTIONTAGGINGCRITERION_HPP_ */
