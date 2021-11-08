/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(MATTERENERGYDENSITY_HPP_)
#error "This file should only be included through NewMatterConstraints.hpp"
#endif

#ifndef MATTERENERGYDENSITY_IMPL_HPP_
#define MATTERENERGYDENSITY_IMPL_HPP_
#include "DimensionDefinitions.hpp"

template <class matter_t>
MatterEnergyDensity<matter_t>::MatterEnergyDensity(
    const matter_t a_matter, double dx)
    :my_matter(a_matter), m_deriv(dx)
{
}

template <class matter_t>
template <class data_t>
void MatterEnergyDensity<matter_t>::compute(Cell<data_t> current_cell) const
{
    // Load local vars and calculate derivs
    const auto vars = current_cell.template load_vars<BSSNMatterVars>();
    const auto d1 = m_deriv.template diff1<BSSNMatterVars>(current_cell);


    // Inverse metric and Christoffel symbol
    const auto h_UU = TensorAlgebra::compute_inverse_sym(vars.h);
    const auto chris = TensorAlgebra::compute_christoffel(d1.h, h_UU);


    // Energy Momentum Tensor
    const auto emtensor = my_matter.compute_emtensor(vars, d1, h_UU, chris.ULL);

    data_t rho = emtensor.rho;

    // Write the constraints into the output FArrayBox
    current_cell.store_vars(rho, c_rho);
}

#endif /* MATTERENERGYDENSITY_IMPL_HPP_ */
