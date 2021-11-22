/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef DIAGNOSTICVARIABLES_HPP
#define DIAGNOSTICVARIABLES_HPP

// assign an enum to each variable
enum
{
    c_rho,

    c_Ham,

    c_Mom,

    c_Weyl4_Re,
    c_Weyl4_Im,

    // the rho2 flux2 source2 come from the killing vector quantity arxiv:2104.13420(katy) equ (10-11 , 18)
    c_rho2,
    c_flux2,
    c_source2,

    NUM_DIAGNOSTIC_VARS
};

namespace DiagnosticVariables
{
static const std::array<std::string, NUM_DIAGNOSTIC_VARS> variable_names = {
    "rho",

    "Ham",

    "Mom",
    
    "Weyl4_Re", "Weyl4_Im",

    "rho2", "flux2", "source2"
};
}

#endif /* DIAGNOSTICVARIABLES_HPP */
