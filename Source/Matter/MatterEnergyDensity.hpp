/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef MATTERENERGYDENSITY_HPP_
#define MATTERENERGYDENSITY_HPP_

#include "CCZ4Geometry.hpp"
#include "Cell.hpp"
#include "FourthOrderDerivatives.hpp"
#include "GRInterval.hpp"
#include "Tensor.hpp"
#include "simd.hpp"
#include <array>

//!  Calculates the Energy density with matter fields
/*!
     The class calculates the Energy density  at each point in a box.
*/
template <class matter_t> class MatterEnergyDensity
{
  public:

    /// Matter variables
    template <class data_t>
    using MatterVars = typename matter_t::template Vars<data_t>;

    /// CCZ4 variables
    template <class data_t> using MetricVars = BSSNVars::VarsNoGauge<data_t>;

    /// Vars object for MatterEnergyDensity
    
    template <class data_t>
    struct BSSNMatterVars : public MetricVars<data_t>,
                            public MatterVars<data_t>
    {
        /// Defines the mapping between members of Vars and Chombo grid
        /// variables (enum in User_Variables)
        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t mapping_function)
        {
            MetricVars<data_t>::enum_mapping(mapping_function);
            MatterVars<data_t>::enum_mapping(mapping_function);
        }
    };

    //! Constructor of class MatterConstraints
    /*!
        Can specify the vars of the constraint vars instead of using the
        hardcoded ones.
    */
    MatterEnergyDensity(const matter_t a_matter, double dx);

    //! The compute member which calculates the constraints at each point in the
    //! box
    template <class data_t> void compute(Cell<data_t> current_cell) const;

  protected:
    matter_t my_matter; //!< The matter object, e.g. a scalar field
    const FourthOrderDerivatives m_deriv;
};

#include "MatterEnergyDensity.impl.hpp"

#endif /* MATTERENERGYDENSITY_HPP_ */
