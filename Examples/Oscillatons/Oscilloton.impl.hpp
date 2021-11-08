// Last edited James Cook 3. August 2017
#if !defined(OSCILLOTON_HPP_)
#error "This file should only be included through Oscilloton.hpp"
#endif

#ifndef OSCILLOTON_IMPL_HPP_
#define OSCILLOTON_IMPL_HPP_

#include <fstream>

#include "ScalarField.hpp"
#include "simd.hpp"

inline
Oscilloton::Oscilloton(params_t a_params, double a_dx, double spacing)
    :m_dx (a_dx), m_spacing(spacing) ,m_params (a_params)
{

  // Define some variables for the constructor
  double input_grid_spacing; // Input file spacing
  double d_alpha2_dr; // Radial derivative of lapse^2 at last position in data file

  // Set values of protected class objects
  m_row_max = m_file_length; // Data file row max (CHANGE ME)
  m_column_max = 8; // Data file column max (CHANGE ME)
  m_gamma = 1.0 / sqrt( 1.0 - pow(m_params.vx,2)); // Gamma factor calculated from vx
  m_gamma2 = 1.0 / sqrt( 1.0 - pow(m_params.vx2,2)); // Gamma factor calculated from vx2
  m_time = 0.0; // Set the 'unboosted' time (Usually set to 0 to simplify)

  // Load the data files
  // Notice ifs1 does not exist
  ifstream ifs ("Oscilloton_data/general_a202.dat");
  ifstream ifs2 ("Oscilloton_data/general_c202.dat");
  ifstream ifs3 ("Oscilloton_data/general_phi202.dat");

  // Load the data into 2D arrays
  // Notice phi has one less column due to the way it is resummed
  for (int i = 0; i < m_row_max; ++i){
    for (int j=0; j < m_column_max; j++){
      if(j<m_column_max-1){
        ifs >> a202[i][j];
        ifs2 >> c202[i][j];
        ifs3 >> phi202[i][j];
      }
      else{
        ifs >> a202[i][j];
        ifs2 >> c202[i][j];
      }
    }
  }

  // Close the data files
  ifs.close();
  ifs2.close();
  ifs3.close();

  // Calculate the grid spacing of the data files
  input_grid_spacing = a202[1][0]-a202[0][0];

  // Split the 2D array into 1D arrays
  for (int i = 0; i < m_row_max; i++){

    // Radial Components
    a_rad[i] = a202[i][0];
    c_rad[i] = c202[i][0];
    phi_rad[i] = phi202[i][0];

    // Components of A
    a0[i] = a202[i][1];
    a2[i] = a202[i][2];
    a4[i] = a202[i][3];
    a6[i] = a202[i][4];
    a8[i] = a202[i][5];
    a10[i] = a202[i][6];
    a12[i] = a202[i][7];

    // Components of C
    c0[i] = c202[i][1];
    c2[i] = c202[i][2];
    c4[i] = c202[i][3];
    c6[i] = c202[i][4];
    c8[i] = c202[i][5];
    c10[i] = c202[i][6];
    c12[i] = c202[i][7];

    // Components of Phi
    phi1[i] = phi202[i][1];
    phi3[i] = phi202[i][2];
    phi5[i] = phi202[i][3];
    phi7[i] = phi202[i][4];
    phi9[i] = phi202[i][5];
    phi11[i] = phi202[i][6];
  }

  // Calculate radial derivatives of the components of A, C, and Phi
  // We use the 2D arrays to do this
  for (int i = 0; i < m_row_max; i++){
    for (int j=0; j < m_column_max; j++){
      // At the center we have even functions
      // Remembering that Phi has a different column length....
      if (i==0){
        if(j<m_column_max-1){
          d_a202_dr[i][j] = 0.0;
          d_c202_dr[i][j] = 0.0;
          d_phi202_dr[i][j] = 0.0;
        }
        else{
          d_a202_dr[i][j] = 0.0;
          d_c202_dr[i][j] = 0.0;
        }
      }
      else if (i==1){
        if(j<m_column_max-1){
          d_a202_dr[i][j] = (1.0/(12.0*input_grid_spacing)) * (-a202[i+2][j]+8.0*a202[i+1][j] - 8.0*a202[i-1][j]+a202[i][j]);
          d_c202_dr[i][j] = (1.0/(12.0*input_grid_spacing)) * (-c202[i+2][j]+8.0*c202[i+1][j] - 8.0*c202[i-1][j]+c202[i][j]);
          d_phi202_dr[i][j] = (1.0/(12.0*input_grid_spacing)) * (-phi202[i+2][j]+8.0*phi202[i+1][j] - 8.0*phi202[i-1][j]+phi202[i][j]);
        }
        else{
          d_a202_dr[i][j] = (1.0/(12.0*input_grid_spacing)) * (-a202[i+2][j]+8.0*a202[i+1][j] - 8.0*a202[i-1][j]+a202[i][j]);
          d_c202_dr[i][j] = (1.0/(12.0*input_grid_spacing)) * (-c202[i+2][j]+8.0*c202[i+1][j] - 8.0*c202[i-1][j]+c202[i][j]);
        }
      }
      // Approximate the last differential
      else if (i==m_row_max-1){
        if(j<m_column_max-1){
          d_a202_dr[i][j] = (1.0/input_grid_spacing) * (a202[i][j] - a202[i-1][j]);
          d_c202_dr[i][j] = (1.0/input_grid_spacing) * (c202[i][j] - c202[i-1][j]);
          d_phi202_dr[i][j] = (1.0/input_grid_spacing) * (phi202[i][j] - phi202[i-1][j]);
        }
        else{
          d_a202_dr[i][j] = (1.0/input_grid_spacing) * (a202[i][j] - a202[i-1][j]);
          d_c202_dr[i][j] = (1.0/input_grid_spacing) * (c202[i][j] - c202[i-1][j]);
        }
      }
      // For second to last differential drop order of the finite differencing
      else if (i==m_row_max-2){
        if(j<m_column_max-1){
          d_a202_dr[i][j] = (1.0/(2.0*input_grid_spacing)) * (a202[i+1][j] - a202[i-1][j]);
          d_c202_dr[i][j] = (1.0/(2.0*input_grid_spacing)) * (c202[i+1][j] - c202[i-1][j]);
          d_phi202_dr[i][j] = (1.0/(2.0*input_grid_spacing)) * (phi202[i+1][j] - phi202[i-1][j]);
        }
        else{
          d_a202_dr[i][j] = (1.0/(2.0*input_grid_spacing)) * (a202[i+1][j] - a202[i-1][j]);
          d_c202_dr[i][j] = (1.0/(2.0*input_grid_spacing)) * (c202[i+1][j] - c202[i-1][j]);
        }
      }
      // Standard 2nd Order Finite Differencing (higher orders not required due to grid spacing of input data)
      else{
        if(j<m_column_max-1){
          d_a202_dr[i][j] = (1.0/(12.0*input_grid_spacing)) * (-a202[i+2][j]+8.0*a202[i+1][j] - 8.0*a202[i-1][j]+a202[i-2][j]);
          d_c202_dr[i][j] = (1.0/(12.0*input_grid_spacing)) * (-c202[i+2][j]+8.0*c202[i+1][j] - 8.0*c202[i-1][j]+c202[i-2][j]);
          d_phi202_dr[i][j] = (1.0/(12.0*input_grid_spacing)) * (-phi202[i+2][j]+8.0*phi202[i+1][j] - 8.0*phi202[i-1][j]+phi202[i-2][j]);
        }
        else{
          d_a202_dr[i][j] = (1.0/(12.0*input_grid_spacing)) * (-a202[i+2][j]+8.0*a202[i+1][j] - 8.0*a202[i-1][j]+a202[i-2][j]);
          d_c202_dr[i][j] = (1.0/(12.0*input_grid_spacing)) * (-c202[i+2][j]+8.0*c202[i+1][j] - 8.0*c202[i-1][j]+c202[i-2][j]);
        }
      }
    }
  }

  // Split the 2D arrays into 1D arrays
  for (int i = 0; i < m_row_max; i++){
    // Radial components
    d_a_rad_dr[i] = d_a202_dr[i][0];
    d_c_rad_dr[i] = d_c202_dr[i][0];
    d_phi_rad_dr[i] = d_phi202_dr[i][0];

    // Radial derivatives of A
    d_a0_dr[i] = d_a202_dr[i][1];
    d_a2_dr[i] = d_a202_dr[i][2];
    d_a4_dr[i] = d_a202_dr[i][3];
    d_a6_dr[i] = d_a202_dr[i][4];
    d_a8_dr[i] = d_a202_dr[i][5];
    d_a10_dr[i] = d_a202_dr[i][6];
    d_a12_dr[i] = d_a202_dr[i][7];

    // Radial derivatives of C
    d_c0_dr[i] = d_c202_dr[i][1];
    d_c2_dr[i] = d_c202_dr[i][2];
    d_c4_dr[i] = d_c202_dr[i][3];
    d_c6_dr[i] = d_c202_dr[i][4];
    d_c8_dr[i] = d_c202_dr[i][5];
    d_c10_dr[i] = d_c202_dr[i][6];
    d_c12_dr[i] = d_c202_dr[i][7];

    // Radial derivatives of Phi
    d_phi1_dr[i] = d_phi202_dr[i][1];
    d_phi3_dr[i] = d_phi202_dr[i][2];
    d_phi5_dr[i] = d_phi202_dr[i][3];
    d_phi7_dr[i] = d_phi202_dr[i][4];
    d_phi9_dr[i] = d_phi202_dr[i][5];
    d_phi11_dr[i] = d_phi202_dr[i][6];
  }

  // Calculate m_omega (Notice how we need only the data files and not the
  // differentials)
  m_omega =sqrt(c0[m_row_max-1])/a0[m_row_max-1];

  // Calculate the parameters used to extend the initial data beyond the data files
  // lapse = (m_b - m_rs_alpha/rr)
  // grr = (m_c - m_rs_g_rr/rr)^(-1)
  d_alpha2_dr = m_omega*m_omega*(d_a202_dr[m_row_max-1][1]/c202[m_row_max-1][1] - (a202[m_row_max-1][1]*d_c202_dr[m_row_max-1][1])/pow(c202[m_row_max-1][1],2));
  m_b = m_omega*m_omega*a202[m_row_max-1][1]/c202[m_row_max-1][1] + a202[m_row_max-1][0] * d_alpha2_dr;
  m_rs_alpha = pow(a202[m_row_max-1][0],2) * d_alpha2_dr;
  m_rs_g_rr = - pow(a202[m_row_max-1][0],2) * d_a202_dr[m_row_max-1][1] * (1./pow(a202[m_row_max-1][1],2));
  m_c = m_rs_g_rr/a202[m_row_max-1][0] + 1./a202[m_row_max-1][1];
}

// Compute the value of the initial vars on the grid
template <class data_t>
void Oscilloton::compute(Cell<data_t> current_cell) const {


  MatterCCZ4<ScalarField<>>::Vars<data_t> vars;
  VarsTools::assign(vars, 0.); 

 // ScalarField<>::Vars<double> vars;
  // Define Coordinates
  Coordinates<double> coords(current_cell,m_dx);
  // metric in radial coordinates
  Tensor<2, double> h_radial;
  // metric in cartesian coordinates
  Tensor<2, double> h_cartesian;
  // inverse metric in cartesian coordinates
  Tensor<2, double> h_cartesian_UU;
  //h_cartesian
  double deth;
  // jacobian
  Tensor<2, double> jacobian ;
  //extrinsic curvature in radial coordinates
  Tensor<2, double> K_radial;
  //extrinsic curvature in cartesian coordinates
  Tensor<2, double> K_cartesian;
  // Trace of cartesian extrinsic curvature
  double Trace_K;
  // traceless extrinsic curvature in cartesian coordiantes
  Tensor<2, double> A_cartesian;
  // conformal decomposition of A_cartesian
  Tensor<2, double> A_cartesian_conformal;
  // radial component of extrinsic curcature in radial coordinates
  double Krr ;
  // rr compontent of metric in radial coordinates
  double grr;
  // Values of the scalar field
  double phi;
  // Values of the conjugate momentum of the scalar field
  double Pi;
  // Values of the lapse
  double lapse, lapse_unboosted;
  // Values of the decomposed metric
  double chi;
  Tensor<2,double> h_cartesian_conformal;

  //Vacuum Modification Variables (Require Dilithium Crystals)
  double A_crit;
  double lapse_crit;
  double seperation;
  double lam1, lam2, lam3;

  double beta_U[3]; //Beta is non-zero as we boost
  double beta[3];

  // We set these values to zero to initalise
  for(int i=0; i<3;i++){
    beta_U[i]=0.0;
    beta[i]=0.0;
  }

  // Set doubles that are used to calculated initial data
  // The value of phi, A and C (components) linearly interpolated at this coord
  double phi_components[6];
  double A_components[7];
  double C_components[7];

  // The value of dr phi, dr A and dr C (components) linearly interpolated at this coord
  double d_phi_dr_components[6];
  double d_A_dr_components[7];
  double d_C_dr_components[7];

  // Value of dx phi, A, C, lapse, grr, at this coord (fully resummed)
  double d_phi_dx[3];
  double d_A_dx[3];
  double d_C_dx[3];
  double d_lapse_dx[3];
  double d_grr_dx[3];

  // We set these values to zero to initalise
  for(int i=0; i<3;i++){
    d_phi_dx[i] = 0.0;
    d_A_dx[3]   = 0.0;
    d_C_dx[3]   = 0.0;
    d_lapse_dx[3]  = 0.0;
    d_grr_dx[3] = 0.0;
  }

  // dt and dr of Lapse, Phi (only dt), grr, A and C
  double d_phi_dt = 0.0;
  double d_lapse_dt = 0.0;
  double d_lapse_dr = 0.0;
  double d_grr_dt = 0.0;
  double d_grr_dr = 0.0;
  double d_A_dt = 0.0;
  double d_C_dt = 0.0;
  double d_A_dr = 0.0;
  double d_C_dr = 0.0;

  // Resummed values of phi, A, C, lapse and grr (non-boosted coords)
  double phi_resummed = 0.0;
  double A_resummed = 0.0;
  double C_resummed = 0.0;
  double lapse_resummed;
  double g_rr_resummed;

  // Setting up coordinates
  // For reference we are in the 'prime' coordinate system
  double t = m_time;
  double x = coords.x - m_params.centerSF[0];
  double y = coords.y - m_params.centerSF[1];
  double z = coords.z - m_params.centerSF[2];

  // We also define the coordinate transform between the prime coordinates
  // and the old coordinates
  double tboost = m_gamma*(t - m_params.vx * x);
  double xboost = m_gamma*(x - m_params.vx * t);

  // Here although we are in the 'prime' coordinate system we need
  // the 'tilde' version of the r coordinate
  double rr2 = pow(xboost,2)+ pow(y,2) + pow(z,2);

  // Regularising radius near 0
  double minimum_dis = 1e-12;
  auto r_is_too_small = simd_compare_lt(rr2, minimum_dis);
  rr2 = simd_conditional(r_is_too_small, minimum_dis, rr2);

  double rr = sqrt(rr2);

  // Radial and Time derivatives of the 'tilde' coordinates
  double dr_dx[3];
  double dr_dt;

  dr_dx[0] = (1.0/rr) * m_gamma * xboost;
  dr_dx[1] = (1.0/rr) * y;
  dr_dx[2] = (1.0/rr) * z;
  dr_dt = -(1.0/rr) * m_gamma * m_params.vx * xboost;

  // Calculate variables from the data file
  if(floor(rr/m_spacing)<m_row_max-1){

    // Interpolate phi, A, and C (and derivatives)
    phi_components[0] = linear_interpolation_new(coords,phi1);
    phi_components[1] = linear_interpolation_new(coords,phi3);
    phi_components[2] = linear_interpolation_new(coords,phi5);
    phi_components[3] = linear_interpolation_new(coords,phi7);
    phi_components[4] = linear_interpolation_new(coords,phi9);
    phi_components[5] = linear_interpolation_new(coords,phi11);

    A_components[0] = linear_interpolation_new(coords,a0);
    A_components[1] = linear_interpolation_new(coords,a2);
    A_components[2] = linear_interpolation_new(coords,a4);
    A_components[3] = linear_interpolation_new(coords,a6);
    A_components[4] = linear_interpolation_new(coords,a8);
    A_components[5] = linear_interpolation_new(coords,a10);
    A_components[6] = linear_interpolation_new(coords,a12);

    C_components[0] = linear_interpolation_new(coords,c0);
    C_components[1] = linear_interpolation_new(coords,c2);
    C_components[2] = linear_interpolation_new(coords,c4);
    C_components[3] = linear_interpolation_new(coords,c6);
    C_components[4] = linear_interpolation_new(coords,c8);
    C_components[5] = linear_interpolation_new(coords,c10);
    C_components[6] = linear_interpolation_new(coords,c12);

    d_phi_dr_components[0] = linear_interpolation_new(coords,d_phi1_dr);
    d_phi_dr_components[1] = linear_interpolation_new(coords,d_phi3_dr);
    d_phi_dr_components[2] = linear_interpolation_new(coords,d_phi5_dr);
    d_phi_dr_components[3] = linear_interpolation_new(coords,d_phi7_dr);
    d_phi_dr_components[4] = linear_interpolation_new(coords,d_phi9_dr);
    d_phi_dr_components[5] = linear_interpolation_new(coords,d_phi11_dr);

    d_A_dr_components[0] = linear_interpolation_new(coords,d_a0_dr);
    d_A_dr_components[1] = linear_interpolation_new(coords,d_a2_dr);
    d_A_dr_components[2] = linear_interpolation_new(coords,d_a4_dr);
    d_A_dr_components[3] = linear_interpolation_new(coords,d_a6_dr);
    d_A_dr_components[4] = linear_interpolation_new(coords,d_a8_dr);
    d_A_dr_components[5] = linear_interpolation_new(coords,d_a10_dr);
    d_A_dr_components[6] = linear_interpolation_new(coords,d_a12_dr);

    d_C_dr_components[0] = linear_interpolation_new(coords,d_c0_dr);
    d_C_dr_components[1] = linear_interpolation_new(coords,d_c2_dr);
    d_C_dr_components[2] = linear_interpolation_new(coords,d_c4_dr);
    d_C_dr_components[3] = linear_interpolation_new(coords,d_c6_dr);
    d_C_dr_components[4] = linear_interpolation_new(coords,d_c8_dr);
    d_C_dr_components[5] = linear_interpolation_new(coords,d_c10_dr);
    d_C_dr_components[6] = linear_interpolation_new(coords,d_c12_dr);

    // Resum phi, and the calculate d_phi_dx(dt)
    for(int i=0; i<6;i++){
      phi_resummed += phi_components[i]*cos((2*i+1)*m_omega* tboost);
      d_phi_dx[0]  += d_phi_dr_components[i]*dr_dx[0]*cos((2*i+1)*m_omega*tboost) + phi_components[i]*(2*i+1)*m_omega*m_gamma*m_params.vx*sin((2*i+1)*m_omega*tboost);
      d_phi_dx[1]  += d_phi_dr_components[i]*dr_dx[1]*cos((2*i+1)*m_omega*tboost);
      d_phi_dx[2]  += d_phi_dr_components[i]*dr_dx[2]*cos((2*i+1)*m_omega*tboost);
      d_phi_dt     += d_phi_dr_components[i]*dr_dt*cos((2*i+1)*m_omega*tboost) - phi_components[i]*(2*i+1)*m_omega*m_gamma*sin((2*i+1)*m_omega*tboost);
    }

    // Correct the value of d_phi_dx by sqrt 8 pi
    for(int i=0; i<3; i++){
      d_phi_dx[i] = d_phi_dx[i]/sqrt(8.0 * M_PI);
    }

    // Correct the value of phi and d_phi_dt by sqrt 8 pi
    d_phi_dt = d_phi_dt/sqrt(8.0 * M_PI);
    phi_resummed = phi_resummed/sqrt(8.0 * M_PI);

    // Resum A, C and calculate dr, dt, A, C
    for(int i=0; i<7;i++){
      A_resummed += A_components[i]*cos((2*i)*m_omega*tboost);
      d_A_dr     += d_A_dr_components[i]*cos((2*i)*m_omega*tboost);
      d_A_dt     += - A_components[i]*(2*i)*m_omega*sin((2*i)*m_omega*tboost);
      C_resummed += C_components[i]*cos((2*i)*m_omega*tboost);
      d_C_dr     += d_C_dr_components[i]*cos((2*i)*m_omega*tboost);
      d_C_dt     += - C_components[i]*(2*i)*m_omega*sin((2*i)*m_omega*tboost);
    }

    // Calculate dr, dt of grr and lapse, and then resum then resum grr and lapse
    d_grr_dr = d_A_dr;
    d_grr_dt = d_A_dt;
    d_grr_dx[0] = d_A_dx[0];
    d_grr_dx[1] = d_A_dx[1];
    d_grr_dx[2] = d_A_dx[2];
    d_lapse_dr = m_omega*0.5 * sqrt(C_resummed/A_resummed) * ( (d_A_dr/C_resummed) - (d_C_dr * A_resummed)/pow(C_resummed,2) );
    d_lapse_dt = m_omega*0.5 * sqrt(C_resummed/A_resummed) * ( (d_A_dt/C_resummed) - (d_C_dt * A_resummed)/pow(C_resummed,2) );
    d_lapse_dx[0] = m_omega*0.5 * sqrt(C_resummed/A_resummed) * ( (d_A_dx[0]/C_resummed) - (d_C_dx[0] * A_resummed)/pow(C_resummed,2) );
    d_lapse_dx[1] = m_omega*0.5 * sqrt(C_resummed/A_resummed) * ( (d_A_dx[1]/C_resummed) - (d_C_dx[1] * A_resummed)/pow(C_resummed,2) );
    d_lapse_dx[2] = m_omega*0.5 * sqrt(C_resummed/A_resummed) * ( (d_A_dx[2]/C_resummed) - (d_C_dx[2] * A_resummed)/pow(C_resummed,2) );
    lapse_resummed = sqrt(A_resummed/C_resummed);
    g_rr_resummed = A_resummed;

    // Set values of grr, phi, Pi, lapse_unboosted, and Krr
    grr = g_rr_resummed;
    phi = phi_resummed;
    Pi = 0.0;
    lapse_unboosted = m_omega*lapse_resummed;
    Krr = 0.0;
  }
  // Repeat as above but beyond the data files
  else{

    //Values of the derivatives of phi in the extention
    d_phi_dt = 0.0;
    d_phi_dx[0] = 0.0;
    d_phi_dx[1] = 0.0;
    d_phi_dx[2] = 0.0;

    // Calculate dr, dt of grr and lapse
    // TO DO FIX!
    grr = (1.0)/(m_c-m_rs_g_rr/rr);
    d_grr_dr = - m_rs_g_rr/pow(rr,2) * pow(grr,2);
    d_grr_dt = - pow(grr,2) * (m_rs_g_rr/(pow(rr,2))) * dr_dt; //((1.)/(m_c-m_rs_g_rr/rr))*((1.)/(m_c-m_rs_g_rr/rr)) * m_rs_g_rr/pow(rr,3)* m_gamma * m_params.vx * x;
    lapse_unboosted = sqrt(m_b - m_rs_alpha/rr);
    d_lapse_dr = m_rs_alpha/pow(rr,2) * (1./(2.*lapse_unboosted));
    d_lapse_dt = (1./(2.*lapse_unboosted)) * (m_rs_alpha/pow(rr,2)) * dr_dt; //m_rs_alpha/pow(rr,3) * m_gamma * m_params.vx * x * (1./(2.*lapse_unboosted));

    // Set values of grr, phi, Pi, lapse_unboosted, and Krr
    phi = 0.0;
    Pi = 0.0;
    Krr = 0.0;
  }

  // Variables for Mathematica CForm
  double lapseunboosted = lapse_unboosted;
  double paramgamma = m_gamma;
  double vel = m_params.vx;
  double aaa = grr;
  double daldt = d_lapse_dt;
  double daldr = d_lapse_dr;
  double daldx = d_lapse_dx[0];
  double daldy = d_lapse_dx[1];
  double daldz = d_lapse_dx[2];
  double daaadt = d_grr_dt;
  double daaadr = d_grr_dr;
  double daaadx = d_grr_dx[0];
  double daaady = d_grr_dx[1];
  double daaadz = d_grr_dx[2];

  // Calculate lapse (now boosted)
  lapse = pow(pow(paramgamma,2)*pow(-1 + pow(vel,2),2)*pow(pow(lapseunboosted,-2) - pow(aaa,-1)*pow(paramgamma,2)*pow(vel,2)*pow(-(t*vel) + x,2)*pow(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2),-1) -
      pow(vel,2)*(pow(y,2) + pow(z,2))*pow(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2),-1),-1),0.5);

  // Calculated beta_i (now boosted)
  beta[0] = vel*pow(paramgamma,2)*(pow(lapseunboosted,2) - (aaa*pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2))*pow(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2),-1));
  beta[1] = (-1 + aaa)*vel*(t*vel - x)*y*pow(paramgamma,2)*pow(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2),-1);
  beta[2] = (-1 + aaa)*vel*(t*vel - x)*z*pow(paramgamma,2)*pow(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2),-1);

  // Initalising
  FOR2(i,j){
	  h_radial[i][j]    = 0.0;
	  h_cartesian[i][j] = 0.0;
	  h_cartesian_conformal[i][j] = 0.0;
	  jacobian[i][j]    = 0.0;
    K_radial[i][j]    = 0.0;
    K_cartesian[i][j] = 0.0;
    A_cartesian[i][j] = 0.0;
 	  A_cartesian_conformal[i][j] = 0.0;
  }

  // Spatial Metric

  h_cartesian[0][0] = pow(paramgamma,2)*(-(pow(lapseunboosted,2)*pow(vel,2)) + (aaa*pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2))*pow(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2),-1));
  h_cartesian[0][1] = (-1 + aaa)*(-(t*vel) + x)*y*pow(paramgamma,2)*pow(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2),-1);
  h_cartesian[0][2] = (-1 + aaa)*(-(t*vel) + x)*z*pow(paramgamma,2)*pow(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2),-1);
  h_cartesian[1][1] = (pow(paramgamma,2)*pow(-(t*vel) + x,2) + aaa*pow(y,2) + pow(z,2))*pow(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2),-1);
  h_cartesian[2][2] = (pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + aaa*pow(z,2))*pow(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2),-1);
  h_cartesian[1][2] = (-1 + aaa)*y*z*pow(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2),-1);
  h_cartesian[2][0] = h_cartesian[0][2];
  h_cartesian[1][0] = h_cartesian[0][1];
  h_cartesian[2][1] = h_cartesian[1][2];

  // Extrinsic Curvature

  K_cartesian[0][0] = -(pow(paramgamma,3)*(-1 + pow(vel,2))*pow(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2),-1.5)*
     (-(paramgamma*vel*(t*vel - x)*pow(lapseunboosted,2)*(daaadr*pow(paramgamma,2)*pow(-(t*vel) + x,2)*(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2)) +
            2*(daaadt*vel*pow(paramgamma,3)*pow(t*vel - x,3) + (-1 + daaadt*paramgamma*vel*(t*vel - x))*pow(y,2) + (-1 + daaadt*paramgamma*vel*(t*vel - x))*pow(z,2))*
             pow(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2),0.5) + 2*daldr*lapseunboosted*pow(vel,2)*pow(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2),2))) +
       aaa*(4*daldr*lapseunboosted*paramgamma*vel*(t*vel - x)*pow(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2),2) +
          pow(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2),0.5)*(paramgamma*(t*vel - x)*
              (-3*daaadt*x*pow(paramgamma,3)*pow(t,2)*pow(vel,2) + daaadt*pow(paramgamma,3)*pow(t,3)*pow(vel,3) - daaadt*paramgamma*x*(pow(paramgamma,2)*pow(x,2) + pow(y,2) + pow(z,2)) +
                vel*(3*daaadt*t*pow(paramgamma,3)*pow(x,2) + (daaadt*paramgamma*t - 2*pow(lapseunboosted,2))*pow(y,2) + (daaadt*paramgamma*t - 2*pow(lapseunboosted,2))*pow(z,2))) +
             2*daldt*lapseunboosted*pow(vel,2)*pow(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2),2))))*
     pow(pow(lapseunboosted,2)*pow(paramgamma,2)*pow(vel,2)*pow(-(t*vel) + x,2) + aaa*(-(pow(paramgamma,2)*pow(-(t*vel) + x,2)) + (-1 + pow(lapseunboosted,2)*pow(vel,2))*pow(y,2) +
          (-1 + pow(lapseunboosted,2)*pow(vel,2))*pow(z,2)),-1));

  K_cartesian[0][1] = y*pow(paramgamma,2)*(-1 + pow(vel,2))*pow(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2),-1.5)*
   (-(paramgamma*vel*(t*vel - x)*pow(lapseunboosted,2)*(daaadr*paramgamma*(t*vel - x)*(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2)) +
          (2*paramgamma*(t*vel - x) + daaadt*vel*(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2)))*pow(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2),0.5))) +
     aaa*(paramgamma*(t*vel - x)*(2*paramgamma*vel*(t*vel - x)*pow(lapseunboosted,2) + daaadt*(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2)))*pow(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2),0.5) +
        2*daldr*lapseunboosted*vel*pow(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2),2)))*
   pow(pow(lapseunboosted,2)*pow(paramgamma,2)*pow(vel,2)*pow(-(t*vel) + x,2) + aaa*(-(pow(paramgamma,2)*pow(-(t*vel) + x,2)) + (-1 + pow(lapseunboosted,2)*pow(vel,2))*pow(y,2) + (-1 + pow(lapseunboosted,2)*pow(vel,2))*pow(z,2)),
    -1);

  K_cartesian[0][2] = z*pow(paramgamma,2)*(-1 + pow(vel,2))*pow(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2),-1.5)*
   (-(paramgamma*vel*(t*vel - x)*pow(lapseunboosted,2)*(daaadr*paramgamma*(t*vel - x)*(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2)) +
          (2*paramgamma*(t*vel - x) + daaadt*vel*(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2)))*pow(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2),0.5))) +
     aaa*(paramgamma*(t*vel - x)*(2*paramgamma*vel*(t*vel - x)*pow(lapseunboosted,2) + daaadt*(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2)))*pow(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2),0.5) +
        2*daldr*lapseunboosted*vel*pow(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2),2)))*
   pow(pow(lapseunboosted,2)*pow(paramgamma,2)*pow(vel,2)*pow(-(t*vel) + x,2) + aaa*(-(pow(paramgamma,2)*pow(-(t*vel) + x,2)) + (-1 + pow(lapseunboosted,2)*pow(vel,2))*pow(y,2) + (-1 + pow(lapseunboosted,2)*pow(vel,2))*pow(z,2)),
    -1);

  K_cartesian[1][1] = paramgamma*(-1 + pow(vel,2))*pow(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2),-1.5)*
   (paramgamma*vel*(t*vel - x)*pow(lapseunboosted,2)*(daaadr*pow(y,2)*(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2)) -
        2*(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(z,2))*pow(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2),0.5) +
        2*aaa*(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(z,2))*pow(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2),0.5)) - aaa*daaadt*pow(y,2)*pow(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2),1.5))*
   pow(pow(lapseunboosted,2)*pow(paramgamma,2)*pow(vel,2)*pow(-(t*vel) + x,2) + aaa*(-(pow(paramgamma,2)*pow(-(t*vel) + x,2)) + (-1 + pow(lapseunboosted,2)*pow(vel,2))*pow(y,2) + (-1 + pow(lapseunboosted,2)*pow(vel,2))*pow(z,2)),
    -1);

  K_cartesian[1][2] = -(paramgamma*y*z*(-1 + pow(vel,2))*pow(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2),-1.5)*
     (-(paramgamma*vel*(t*vel - x)*pow(lapseunboosted,2)*(daaadr*(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2)) + 2*pow(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2),0.5) -
            2*aaa*pow(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2),0.5))) + aaa*daaadt*pow(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2),1.5))*
     pow(pow(lapseunboosted,2)*pow(paramgamma,2)*pow(vel,2)*pow(-(t*vel) + x,2) + aaa*(-(pow(paramgamma,2)*pow(-(t*vel) + x,2)) + (-1 + pow(lapseunboosted,2)*pow(vel,2))*pow(y,2) +
          (-1 + pow(lapseunboosted,2)*pow(vel,2))*pow(z,2)),-1));

  K_cartesian[2][2] = paramgamma*(-1 + pow(vel,2))*pow(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2),-1.5)*
   (paramgamma*vel*(t*vel - x)*pow(lapseunboosted,2)*(daaadr*pow(z,2)*(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2)) -
        2*(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2))*pow(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2),0.5) +
        2*aaa*(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2))*pow(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2),0.5)) - aaa*daaadt*pow(z,2)*pow(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2),1.5))*
   pow(pow(lapseunboosted,2)*pow(paramgamma,2)*pow(vel,2)*pow(-(t*vel) + x,2) + aaa*(-(pow(paramgamma,2)*pow(-(t*vel) + x,2)) + (-1 + pow(lapseunboosted,2)*pow(vel,2))*pow(y,2) + (-1 + pow(lapseunboosted,2)*pow(vel,2))*pow(z,2)),
    -1);

  K_cartesian[1][0] = K_cartesian[0][1];
  K_cartesian[2][0] = K_cartesian[0][2];
  K_cartesian[2][1] = K_cartesian[1][2];

  FOR2(i,j){
    K_cartesian[i][j]=1.0/(2.0*lapse) * K_cartesian[i][j];
  }

  // Inverse of metric
  h_cartesian_UU = TensorAlgebra::compute_inverse(h_cartesian);

  FOR2(i,j){
    beta_U[i]+=beta[j]*h_cartesian_UU[i][j];
  }

  // We calculate Pi
  for(int i=0; i<3;i++){
    Pi += -(1.0/lapse)*beta_U[i]*d_phi_dx[i];
  }
  Pi += (1.0/lapse) * (d_phi_dt);

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// STAR 2
///////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

double phi2, Pi2, lapse2;

double beta2_U[3];
double beta2[3];

Tensor<2, double> h_cartesian2;
Tensor<2, double> K_cartesian2;
Tensor<2, double> h_cartesian2_UU;

// We set these values to zero to initalise
for(int i=0; i<3;i++){
  beta2_U[i]=0.0;
  beta2[i]=0.0;
}

t = 0.0;

x = coords.x - m_params.centerSF2[0];
y = coords.y - m_params.centerSF2[1];
z = coords.z - m_params.centerSF2[2];

tboost = m_gamma2*(t - m_params.vx2 * x);
xboost = m_gamma2*(x - m_params.vx2 * t);

rr2 = pow(xboost,2)+ pow(y,2) + pow(z,2);
r_is_too_small = simd_compare_lt(rr2, minimum_dis);
rr2 = simd_conditional(r_is_too_small, minimum_dis, rr2);
rr = sqrt(rr2);

dr_dx[0] = (1.0/rr) * m_gamma2 * xboost;
dr_dx[1] = (1.0/rr) * y;
dr_dx[2] = (1.0/rr) * z;
dr_dt = -(1.0/rr) * m_gamma2 * m_params.vx2 * xboost;

// Calculate variables from the data file
if(floor(rr/m_spacing)<m_row_max-1){

  // Interpolate phi, A, and C (and derivatives)
  phi_components[0] = linear_interpolation_new2(coords,phi1);
  phi_components[1] = linear_interpolation_new2(coords,phi3);
  phi_components[2] = linear_interpolation_new2(coords,phi5);
  phi_components[3] = linear_interpolation_new2(coords,phi7);
  phi_components[4] = linear_interpolation_new2(coords,phi9);
  phi_components[5] = linear_interpolation_new2(coords,phi11);

  A_components[0] = linear_interpolation_new2(coords,a0);
  A_components[1] = linear_interpolation_new2(coords,a2);
  A_components[2] = linear_interpolation_new2(coords,a4);
  A_components[3] = linear_interpolation_new2(coords,a6);
  A_components[4] = linear_interpolation_new2(coords,a8);
  A_components[5] = linear_interpolation_new2(coords,a10);
  A_components[6] = linear_interpolation_new2(coords,a12);

  C_components[0] = linear_interpolation_new2(coords,c0);
  C_components[1] = linear_interpolation_new2(coords,c2);
  C_components[2] = linear_interpolation_new2(coords,c4);
  C_components[3] = linear_interpolation_new2(coords,c6);
  C_components[4] = linear_interpolation_new2(coords,c8);
  C_components[5] = linear_interpolation_new2(coords,c10);
  C_components[6] = linear_interpolation_new2(coords,c12);

  d_phi_dr_components[0] = linear_interpolation_new2(coords,d_phi1_dr);
  d_phi_dr_components[1] = linear_interpolation_new2(coords,d_phi3_dr);
  d_phi_dr_components[2] = linear_interpolation_new2(coords,d_phi5_dr);
  d_phi_dr_components[3] = linear_interpolation_new2(coords,d_phi7_dr);
  d_phi_dr_components[4] = linear_interpolation_new2(coords,d_phi9_dr);
  d_phi_dr_components[5] = linear_interpolation_new2(coords,d_phi11_dr);

  d_A_dr_components[0] = linear_interpolation_new2(coords,d_a0_dr);
  d_A_dr_components[1] = linear_interpolation_new2(coords,d_a2_dr);
  d_A_dr_components[2] = linear_interpolation_new2(coords,d_a4_dr);
  d_A_dr_components[3] = linear_interpolation_new2(coords,d_a6_dr);
  d_A_dr_components[4] = linear_interpolation_new2(coords,d_a8_dr);
  d_A_dr_components[5] = linear_interpolation_new2(coords,d_a10_dr);
  d_A_dr_components[6] = linear_interpolation_new2(coords,d_a12_dr);

  d_C_dr_components[0] = linear_interpolation_new2(coords,d_c0_dr);
  d_C_dr_components[1] = linear_interpolation_new2(coords,d_c2_dr);
  d_C_dr_components[2] = linear_interpolation_new2(coords,d_c4_dr);
  d_C_dr_components[3] = linear_interpolation_new2(coords,d_c6_dr);
  d_C_dr_components[4] = linear_interpolation_new2(coords,d_c8_dr);
  d_C_dr_components[5] = linear_interpolation_new2(coords,d_c10_dr);
  d_C_dr_components[6] = linear_interpolation_new2(coords,d_c12_dr);

  phi_resummed = 0.0;
  d_phi_dx[0]  = 0.0;
  d_phi_dx[1]  = 0.0;
  d_phi_dx[2]  = 0.0;
  d_phi_dt     = 0.0;

  // Resum phi, and the calculate d_phi_dx(dt)
  for(int i=0; i<6;i++){
    phi_resummed += phi_components[i]*cos((2*i+1)*m_omega* tboost);
    d_phi_dx[0]  += d_phi_dr_components[i]*dr_dx[0]*cos((2*i+1)*m_omega*tboost) + phi_components[i]*(2*i+1)*m_omega*m_gamma2*m_params.vx2*sin((2*i+1)*m_omega*tboost);
    d_phi_dx[1]  += d_phi_dr_components[i]*dr_dx[1]*cos((2*i+1)*m_omega*tboost);
    d_phi_dx[2]  += d_phi_dr_components[i]*dr_dx[2]*cos((2*i+1)*m_omega*tboost);
    d_phi_dt     += d_phi_dr_components[i]*dr_dt*cos((2*i+1)*m_omega*tboost) - phi_components[i]*(2*i+1)*m_omega*m_gamma2*sin((2*i+1)*m_omega*tboost);
  }

  // Correct the value of d_phi_dx by sqrt 8 pi
  for(int i=0; i<3; i++){
    d_phi_dx[i] = d_phi_dx[i]/sqrt(8.0 * M_PI);
  }

  // Correct the value of phi and d_phi_dt by sqrt 8 pi
  d_phi_dt = d_phi_dt/sqrt(8.0 * M_PI);
  phi_resummed = phi_resummed/sqrt(8.0 * M_PI);

  // Resum A, C and calculate dr, dt, A, C

  A_resummed = 0.0;
  d_A_dr     = 0.0;
  d_A_dt     = 0.0;
  C_resummed = 0.0;
  d_C_dr     = 0.0;
  d_C_dt     = 0.0;

  for(int i=0; i<7;i++){
    A_resummed += A_components[i]*cos((2*i)*m_omega*tboost);
    d_A_dr     += d_A_dr_components[i]*cos((2*i)*m_omega*tboost);
    d_A_dt     += - A_components[i]*(2*i)*m_omega*sin((2*i)*m_omega*tboost);
    C_resummed += C_components[i]*cos((2*i)*m_omega*tboost);
    d_C_dr     += d_C_dr_components[i]*cos((2*i)*m_omega*tboost);
    d_C_dt     += - C_components[i]*(2*i)*m_omega*sin((2*i)*m_omega*tboost);
  }

  // Calculate dr, dt of grr and lapse, and then resum then resum grr and lapse
  d_grr_dr = d_A_dr;
  d_grr_dt = d_A_dt;
  d_grr_dx[0] = d_A_dx[0];
  d_grr_dx[1] = d_A_dx[1];
  d_grr_dx[2] = d_A_dx[2];
  d_lapse_dr = m_omega*0.5 * sqrt(C_resummed/A_resummed) * ( (d_A_dr/C_resummed) - (d_C_dr * A_resummed)/pow(C_resummed,2) );
  d_lapse_dt = m_omega*0.5 * sqrt(C_resummed/A_resummed) * ( (d_A_dt/C_resummed) - (d_C_dt * A_resummed)/pow(C_resummed,2) );
  d_lapse_dx[0] = m_omega*0.5 * sqrt(C_resummed/A_resummed) * ( (d_A_dx[0]/C_resummed) - (d_C_dx[0] * A_resummed)/pow(C_resummed,2) );
  d_lapse_dx[1] = m_omega*0.5 * sqrt(C_resummed/A_resummed) * ( (d_A_dx[1]/C_resummed) - (d_C_dx[1] * A_resummed)/pow(C_resummed,2) );
  d_lapse_dx[2] = m_omega*0.5 * sqrt(C_resummed/A_resummed) * ( (d_A_dx[2]/C_resummed) - (d_C_dx[2] * A_resummed)/pow(C_resummed,2) );
  lapse_resummed = sqrt(A_resummed/C_resummed);
  g_rr_resummed = A_resummed;

  // Set values of grr, phi, Pi, lapse_unboosted, and Krr
  grr = g_rr_resummed;
  phi2 = phi_resummed;
  Pi2 = 0.0;
  lapse_unboosted = m_omega*lapse_resummed;
  Krr = 0.0;
}
// Repeat as above but beyond the data files
else{

  //Values of the derivatives of phi in the extention
  d_phi_dt = 0.0;
  d_phi_dx[0] = 0.0;
  d_phi_dx[1] = 0.0;
  d_phi_dx[2] = 0.0;

  // Calculate dr, dt of grr and lapse
  // TO DO FIX!
  grr = (1.0)/(m_c-m_rs_g_rr/rr);
  d_grr_dr = - m_rs_g_rr/pow(rr,2) * pow(grr,2);
  d_grr_dt = - pow(grr,2) * (m_rs_g_rr/(pow(rr,2))) * dr_dt; //((1.)/(m_c-m_rs_g_rr/rr))*((1.)/(m_c-m_rs_g_rr/rr)) * m_rs_g_rr/pow(rr,3)* m_gamma2 * m_params.vx2 * x;
  lapse_unboosted = sqrt(m_b - m_rs_alpha/rr);
  d_lapse_dr = m_rs_alpha/pow(rr,2) * (1./(2.*lapse_unboosted));
  d_lapse_dt = (1./(2.*lapse_unboosted)) * (m_rs_alpha/pow(rr,2)) * dr_dt; //m_rs_alpha/pow(rr,3) * m_gamma2 * m_params.vx2 * x * (1./(2.*lapse_unboosted));

  // Set values of grr, phi, Pi, lapse_unboosted, and Krr
  phi2 = 0.0;
  Pi2 = 0.0;
  Krr = 0.0;
}

// Variables for Mathematica CForm
lapseunboosted = lapse_unboosted;
paramgamma = m_gamma2;
vel = m_params.vx2;
aaa = grr;
daldt = d_lapse_dt;
daldr = d_lapse_dr;
daldx = d_lapse_dx[0];
daldy = d_lapse_dx[1];
daldz = d_lapse_dx[2];
daaadt = d_grr_dt;
daaadr = d_grr_dr;
daaadx = d_grr_dx[0];
daaady = d_grr_dx[1];
daaadz = d_grr_dx[2];

// Calculate lapse (now boosted)
lapse2 = pow(pow(paramgamma,2)*pow(-1 + pow(vel,2),2)*pow(pow(lapseunboosted,-2) - pow(aaa,-1)*pow(paramgamma,2)*pow(vel,2)*pow(-(t*vel) + x,2)*pow(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2),-1) -
    pow(vel,2)*(pow(y,2) + pow(z,2))*pow(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2),-1),-1),0.5);

beta2[0] = vel*pow(paramgamma,2)*(pow(lapseunboosted,2) - (aaa*pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2))*pow(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2),-1));
beta2[1] = (-1 + aaa)*vel*(t*vel - x)*y*pow(paramgamma,2)*pow(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2),-1);
beta2[2] = (-1 + aaa)*vel*(t*vel - x)*z*pow(paramgamma,2)*pow(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2),-1);

FOR2(i,j){
  h_cartesian2[i][j] = 0.0;
  K_cartesian2[i][j] = 0.0;
}

h_cartesian2[0][0] = pow(paramgamma,2)*(-(pow(lapseunboosted,2)*pow(vel,2)) + (aaa*pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2))*pow(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2),-1));
h_cartesian2[0][1] = (-1 + aaa)*(-(t*vel) + x)*y*pow(paramgamma,2)*pow(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2),-1);
h_cartesian2[0][2] = (-1 + aaa)*(-(t*vel) + x)*z*pow(paramgamma,2)*pow(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2),-1);
h_cartesian2[1][1] = (pow(paramgamma,2)*pow(-(t*vel) + x,2) + aaa*pow(y,2) + pow(z,2))*pow(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2),-1);
h_cartesian2[2][2] = (pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + aaa*pow(z,2))*pow(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2),-1);
h_cartesian2[1][2] = (-1 + aaa)*y*z*pow(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2),-1);
h_cartesian2[2][0] = h_cartesian2[0][2];
h_cartesian2[1][0] = h_cartesian2[0][1];
h_cartesian2[2][1] = h_cartesian2[1][2];

K_cartesian2[0][0] = -(pow(paramgamma,3)*(-1 + pow(vel,2))*pow(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2),-1.5)*
   (-(paramgamma*vel*(t*vel - x)*pow(lapseunboosted,2)*(daaadr*pow(paramgamma,2)*pow(-(t*vel) + x,2)*(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2)) +
          2*(daaadt*vel*pow(paramgamma,3)*pow(t*vel - x,3) + (-1 + daaadt*paramgamma*vel*(t*vel - x))*pow(y,2) + (-1 + daaadt*paramgamma*vel*(t*vel - x))*pow(z,2))*
           pow(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2),0.5) + 2*daldr*lapseunboosted*pow(vel,2)*pow(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2),2))) +
     aaa*(4*daldr*lapseunboosted*paramgamma*vel*(t*vel - x)*pow(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2),2) +
        pow(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2),0.5)*(paramgamma*(t*vel - x)*
            (-3*daaadt*x*pow(paramgamma,3)*pow(t,2)*pow(vel,2) + daaadt*pow(paramgamma,3)*pow(t,3)*pow(vel,3) - daaadt*paramgamma*x*(pow(paramgamma,2)*pow(x,2) + pow(y,2) + pow(z,2)) +
              vel*(3*daaadt*t*pow(paramgamma,3)*pow(x,2) + (daaadt*paramgamma*t - 2*pow(lapseunboosted,2))*pow(y,2) + (daaadt*paramgamma*t - 2*pow(lapseunboosted,2))*pow(z,2))) +
           2*daldt*lapseunboosted*pow(vel,2)*pow(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2),2))))*
   pow(pow(lapseunboosted,2)*pow(paramgamma,2)*pow(vel,2)*pow(-(t*vel) + x,2) + aaa*(-(pow(paramgamma,2)*pow(-(t*vel) + x,2)) + (-1 + pow(lapseunboosted,2)*pow(vel,2))*pow(y,2) +
        (-1 + pow(lapseunboosted,2)*pow(vel,2))*pow(z,2)),-1));

K_cartesian2[0][1] = y*pow(paramgamma,2)*(-1 + pow(vel,2))*pow(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2),-1.5)*
 (-(paramgamma*vel*(t*vel - x)*pow(lapseunboosted,2)*(daaadr*paramgamma*(t*vel - x)*(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2)) +
        (2*paramgamma*(t*vel - x) + daaadt*vel*(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2)))*pow(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2),0.5))) +
   aaa*(paramgamma*(t*vel - x)*(2*paramgamma*vel*(t*vel - x)*pow(lapseunboosted,2) + daaadt*(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2)))*pow(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2),0.5) +
      2*daldr*lapseunboosted*vel*pow(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2),2)))*
 pow(pow(lapseunboosted,2)*pow(paramgamma,2)*pow(vel,2)*pow(-(t*vel) + x,2) + aaa*(-(pow(paramgamma,2)*pow(-(t*vel) + x,2)) + (-1 + pow(lapseunboosted,2)*pow(vel,2))*pow(y,2) + (-1 + pow(lapseunboosted,2)*pow(vel,2))*pow(z,2)),
  -1);

K_cartesian2[0][2] = z*pow(paramgamma,2)*(-1 + pow(vel,2))*pow(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2),-1.5)*
 (-(paramgamma*vel*(t*vel - x)*pow(lapseunboosted,2)*(daaadr*paramgamma*(t*vel - x)*(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2)) +
        (2*paramgamma*(t*vel - x) + daaadt*vel*(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2)))*pow(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2),0.5))) +
   aaa*(paramgamma*(t*vel - x)*(2*paramgamma*vel*(t*vel - x)*pow(lapseunboosted,2) + daaadt*(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2)))*pow(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2),0.5) +
      2*daldr*lapseunboosted*vel*pow(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2),2)))*
 pow(pow(lapseunboosted,2)*pow(paramgamma,2)*pow(vel,2)*pow(-(t*vel) + x,2) + aaa*(-(pow(paramgamma,2)*pow(-(t*vel) + x,2)) + (-1 + pow(lapseunboosted,2)*pow(vel,2))*pow(y,2) + (-1 + pow(lapseunboosted,2)*pow(vel,2))*pow(z,2)),
  -1);

K_cartesian2[1][1] = paramgamma*(-1 + pow(vel,2))*pow(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2),-1.5)*
 (paramgamma*vel*(t*vel - x)*pow(lapseunboosted,2)*(daaadr*pow(y,2)*(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2)) -
      2*(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(z,2))*pow(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2),0.5) +
      2*aaa*(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(z,2))*pow(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2),0.5)) - aaa*daaadt*pow(y,2)*pow(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2),1.5))*
 pow(pow(lapseunboosted,2)*pow(paramgamma,2)*pow(vel,2)*pow(-(t*vel) + x,2) + aaa*(-(pow(paramgamma,2)*pow(-(t*vel) + x,2)) + (-1 + pow(lapseunboosted,2)*pow(vel,2))*pow(y,2) + (-1 + pow(lapseunboosted,2)*pow(vel,2))*pow(z,2)),
  -1);

K_cartesian2[1][2] = -(paramgamma*y*z*(-1 + pow(vel,2))*pow(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2),-1.5)*
   (-(paramgamma*vel*(t*vel - x)*pow(lapseunboosted,2)*(daaadr*(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2)) + 2*pow(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2),0.5) -
          2*aaa*pow(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2),0.5))) + aaa*daaadt*pow(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2),1.5))*
   pow(pow(lapseunboosted,2)*pow(paramgamma,2)*pow(vel,2)*pow(-(t*vel) + x,2) + aaa*(-(pow(paramgamma,2)*pow(-(t*vel) + x,2)) + (-1 + pow(lapseunboosted,2)*pow(vel,2))*pow(y,2) +
        (-1 + pow(lapseunboosted,2)*pow(vel,2))*pow(z,2)),-1));

K_cartesian2[2][2] = paramgamma*(-1 + pow(vel,2))*pow(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2),-1.5)*
 (paramgamma*vel*(t*vel - x)*pow(lapseunboosted,2)*(daaadr*pow(z,2)*(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2)) -
      2*(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2))*pow(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2),0.5) +
      2*aaa*(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2))*pow(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2),0.5)) - aaa*daaadt*pow(z,2)*pow(pow(paramgamma,2)*pow(-(t*vel) + x,2) + pow(y,2) + pow(z,2),1.5))*
 pow(pow(lapseunboosted,2)*pow(paramgamma,2)*pow(vel,2)*pow(-(t*vel) + x,2) + aaa*(-(pow(paramgamma,2)*pow(-(t*vel) + x,2)) + (-1 + pow(lapseunboosted,2)*pow(vel,2))*pow(y,2) + (-1 + pow(lapseunboosted,2)*pow(vel,2))*pow(z,2)),
  -1);

K_cartesian2[1][0] = K_cartesian2[0][1];
K_cartesian2[2][0] = K_cartesian2[0][2];
K_cartesian2[2][1] = K_cartesian2[1][2];

FOR2(i,j){
  K_cartesian2[i][j]=1.0/(2.0*lapse2) * K_cartesian2[i][j];
}

// Inverse of metric
h_cartesian2_UU = TensorAlgebra::compute_inverse(h_cartesian2);

FOR2(i,j){
  beta2_U[i] += beta2[j]*h_cartesian2_UU[i][j];
}

// We calculate Pi
for(int i=0; i<3;i++){
  Pi2 += -(1.0/lapse2)*beta2_U[i]*d_phi_dx[i];
}
Pi2 += (1.0/lapse2) * (d_phi_dt);

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// COMBINING STAR 1 and STAR 2
///////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

phi = phi + phi2;
Pi = Pi + Pi2;
lapse = lapse + lapse2 - 1.0;

for(int i=0; i<3;i++){
  beta[i] = beta[i] + beta2[i];
}

FOR2(i,j){
  h_cartesian[i][j] = h_cartesian[i][j] + h_cartesian2[i][j];
  K_cartesian[i][j] = K_cartesian[i][j] + K_cartesian2[i][j];
}


double xpos = abs(m_params.centerSF2[0] - m_params.centerSF[0]);
double ypos = abs(m_params.centerSF2[1] - m_params.centerSF[1]);
double zpos = abs(m_params.centerSF2[2] - m_params.centerSF[2]);

seperation = sqrt(xpos*xpos+ypos*ypos+zpos*zpos);


A_crit = (1.0)/(m_c-m_rs_g_rr/seperation);
lapse_crit = sqrt(m_b - m_rs_alpha/seperation);

double r_normal, rr2_normal, x_normal, y_normal, z_normal;
double eul, shift_thomas;

x_normal = coords.x - (m_params.centerSF2[0] + m_params.centerSF[0])/2.0;
y_normal = coords.y - (m_params.centerSF2[1] + m_params.centerSF[1])/2.0;
z_normal = coords.z - (m_params.centerSF2[2] + m_params.centerSF[2])/2.0;


rr2_normal = pow(x_normal,2)+ pow(y_normal,2) + pow(z_normal,2);
r_normal = sqrt(rr2_normal);

shift_thomas = (seperation/2.+5.);

eul = 2.7182818284590452353602874713527;

CH_assert(m_params.centerSF[2]==m_params.centerSF2[2]);
CH_assert(m_params.vx==-m_params.vx2);
CH_assert(t==0);

lam1  = pow(m_gamma,2)*(-(pow(lapse_crit,2)*pow(m_params.vx,2)) + (A_crit*pow(m_gamma,2)*pow(xpos,2) + pow(ypos,2))*pow(pow(m_gamma,2)*pow(xpos,2) + pow(ypos,2),-1.));
lam2 = (pow(m_gamma,2)*pow(xpos,2) + A_crit*pow(ypos,2))*pow(pow(m_gamma,2)*pow(xpos,2) + pow(ypos,2),-1.);
lam3 = (-1. + A_crit)*xpos*ypos*pow(m_gamma,2)*pow(pow(m_gamma,2)*pow(xpos,2) + pow(ypos,2),-1);
if(r_normal>shift_thomas){
    lam1 = (lam1-1.)*pow(eul,-pow((r_normal-shift_thomas)/shift_thomas,2))+1.;
    lam2 = (lam2-1.)*pow(eul,-pow((r_normal-shift_thomas)/shift_thomas,2))+1.;
    lam3 = (lam3)*pow(eul,-pow((r_normal-shift_thomas)/shift_thomas,2));
}
h_cartesian[0][0] = h_cartesian[0][0] - lam1;
h_cartesian[1][1] = h_cartesian[1][1] - lam2;
h_cartesian[0][1] = h_cartesian[0][1] - lam3;
h_cartesian[1][0] = h_cartesian[1][0] - lam3;
h_cartesian[2][2] = h_cartesian[2][2] - 1.0;



h_cartesian_UU = TensorAlgebra::compute_inverse(h_cartesian);

FOR1(i){
  beta_U[i] = 0.0;
}

FOR2(i,j){
  beta_U[i]+=beta[j]*h_cartesian_UU[i][j];
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// CONFORMAL DECOMPOSITION
///////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

  // Calculating Trace of extrinsic curvature K
  Trace_K = 0 ;

  FOR2(i,j){
   	Trace_K += K_cartesian[i][j]*h_cartesian_UU[i][j];
  }
  FOR2(i,j){
   	A_cartesian[i][j] += K_cartesian[i][j]-1.0/3.0*h_cartesian[i][j]*Trace_K;
  }

//  conformal decomposition

  deth = h_cartesian[0][0]*(h_cartesian[1][1]*h_cartesian[2][2]-h_cartesian[1][2]*h_cartesian[2][1])-
         h_cartesian[0][1]*(h_cartesian[2][2]*h_cartesian[1][0]-h_cartesian[1][2]*h_cartesian[2][0])+
         h_cartesian[0][2]*(h_cartesian[1][0]*h_cartesian[2][1]-h_cartesian[1][1]*h_cartesian[2][0]);


  chi = pow(deth,-1.0/3.0);

  FOR2(i,j)
  {
    h_cartesian_conformal[i][j] = h_cartesian[i][j]*chi;
    A_cartesian_conformal[i][j] = A_cartesian[i][j]*chi;
  }

  // making A tracefree (conformal decomposition later)

//Send to Global Array
  vars.chi = chi;
  vars.lapse = lapse;
  vars.phi = phi;
  vars.Pi = Pi;
  vars.K = Trace_K;
  FOR1(i){
    vars.shift[i] = beta_U[i];
  }
  FOR2(i,j){
     vars.h[i][j] = h_cartesian_conformal[i][j];
     vars.A[i][j] = A_cartesian_conformal[i][j];
  }
  current_cell.store_vars(vars);
}

// Interpolate data
template <class data_t>
double Oscilloton::linear_interpolation_new(Coordinates<data_t> coords, const double (*data_jt)) const {

  double t = m_time;
  double x = coords.x - m_params.centerSF[0];
  double y = coords.y - m_params.centerSF[1];
  double z = coords.z - m_params.centerSF[2];

  double rr2 =   pow((x - m_params.vx * t)*m_gamma,2)+ pow(y,2) + pow(z,2);

  double minimum_rr2 = 1e-12;
  auto r_is_too_small = simd_compare_lt(rr2, minimum_rr2);
  rr2 = simd_conditional(r_is_too_small, minimum_rr2, rr2);

  double rr = sqrt(rr2);

  int indxL = static_cast<int>(floor(rr/m_spacing));
  int indxH = static_cast<int>(ceil(rr/m_spacing));

  double data_L = *(data_jt +indxL);
  double data_H = *(data_jt +indxH);

  double out =  data_L+ (rr/m_spacing - indxL) * (data_H - data_L);

  return out;
}

// Interpolate data
template <class data_t>
double Oscilloton::linear_interpolation_new2(Coordinates<data_t> coords, const double (*data_jt)) const  {

  double t = m_time;
  double x = coords.x - m_params.centerSF2[0];
  double y = coords.y - m_params.centerSF2[1];
  double z = coords.z - m_params.centerSF2[2];

  double rr2 =   pow((x - m_params.vx2 * t)*m_gamma2,2)+ pow(y,2) + pow(z,2);

  double minimum_rr2 = 1e-12;
  auto r_is_too_small = simd_compare_lt(rr2, minimum_rr2);
  rr2 = simd_conditional(r_is_too_small, minimum_rr2, rr2);

  double rr = sqrt(rr2);

  int indxL = static_cast<int>(floor(rr/m_spacing));
  int indxH = static_cast<int>(ceil(rr/m_spacing));

  double data_L = *(data_jt +indxL);
  double data_H = *(data_jt +indxH);

  double out =  data_L+ (rr/m_spacing - indxL) * (data_H - data_L);

  return out;
}

#endif /* OSCILLOTON_IMPL_HPP_ */
