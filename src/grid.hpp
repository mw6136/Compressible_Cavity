#include <cmath>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "H5Cpp.h"

#ifndef GRID_HPP_
#define GRID_HPP_

std::string
ZeroPadNumber (int num)
{
  std::ostringstream ss;
  ss << std::setw (4) << std::setfill ('0') << num;
  std::string result = ss.str ();
  if (result.length () > 7)
    {
      result.erase (0, result.length () - 7);
    }
  return result;
}

std::string
getFname (int iter, std::string prefix)
{

  std::string fname;
  std::stringstream ss;
  std::stringstream ssfname;
  std::string newnum = ZeroPadNumber (iter);

  ss << prefix << newnum << ".h5";

  ss >> fname;

  return fname;
}

class Array2D
{
  int x1s_, x2s_, nvar_;

public:
  std::vector<double> m_data;

  Array2D (size_t nvar, size_t x1, size_t x2, double init = 0.0)
      : nvar_ (nvar), x1s_ (x1), x2s_ (x2), m_data (nvar_ * x1s_ * x2s_, init)
  {
  }
  double &
  operator() (size_t n, size_t x, size_t y)
  {
    return m_data[y + x1s_ * (x + x2s_ * n)];
  }
  double
  operator() (size_t n, size_t x, size_t y) const
  {
    return m_data[y + x1s_ * (x + x2s_ * n)];
  }
};

class SimpleMesh
{

public:
  std::string m_wall_cond;

  // Number of grid points in x1, x2 direction and number of ghost zones
  int m_nx1, m_nx2, m_ng;
  double m_x1max, m_x1min, m_x2max, m_x2min, m_dx, m_dy;

  // Fluid quantities to track
  double gamma, R;

  // Vector for Un, Unp1, x direction flux (F), y direction flux (F)
  Array2D Un, Unp1, G_i_jpoh, G_i_jmoh, F_ipoh_j, F_imoh_j;

  // Arrays for other important intermediate quantities
  Array2D thermlam, mu, rho, u, v, et, e, p, T;

  // Arrays to store the derivative information for viscosity, temp diffusion
  // at the locations of where the fluxes are

  SimpleMesh (int nx1, int nx2, double x1max, double x1min, double x2max,
              double x2min, int ng, double i_gamma, double i_R,
              std::string wall_cond)
      : Un (4, nx2 + 2 * ng, nx1 + 2 * ng),
        G_i_jpoh (4, nx2 + (2 * ng), nx1 + (2 * ng)),
        G_i_jmoh (4, nx2 + (2 * ng), nx1 + (2 * ng)),
        F_ipoh_j (4, nx2 + (2 * ng), nx1 + (2 * ng)),
        F_imoh_j (4, nx2 + (2 * ng), nx1 + (2 * ng)),
        Unp1 (4, nx2 + 2 * ng, nx1 + 2 * ng),
        thermlam (1, nx2 + 2 * ng, nx1 + 2 * ng),
        mu (1, nx2 + 2 * ng, nx1 + 2 * ng), rho (1, nx2 + 2 * ng, nx1 + 2 * ng),
        u (1, nx2 + 2 * ng, nx1 + 2 * ng), v (1, nx2 + 2 * ng, nx1 + 2 * ng),
        et (1, nx2 + 2 * ng, nx1 + 2 * ng), e (1, nx2 + 2 * ng, nx1 + 2 * ng),
        p (1, nx2 + 2 * ng, nx1 + 2 * ng), T (1, nx2 + 2 * ng, nx1 + 2 * ng)
  {

    m_ng = ng;
    m_nx1 = nx1;
    m_nx2 = nx2;
    m_wall_cond = wall_cond;

    m_x1max = x1max;
    m_x1min = x1min;
    m_x2max = x2max;
    m_x2min = x2min;

    m_dx = (x1max - x1min) / nx1;
    m_dy = (x2max - x2min) / nx2;

    gamma = i_gamma;
    R = i_R;
  }

  void saveOutput (int iter, double t, std::string fpref);
  // void SingleTimestep();
  void swapArrays ();
  void fillPrimitive (double nu, std::string errcheck);
  // void enforceBCs(double t, double Uw, double omega);

  void new_SingleTimestep ();
  void new_enforceBCs (double t, double Uw, double omega, double nu);
  // void fillDerivativeArrays();
  // void tempFillFluxesInefficient();
};

void
SimpleMesh::new_enforceBCs (double t, double Uw, double omega, double nu)
{

  double T_inner, e_inner;
  int nx1t = m_nx1 + 2 * m_ng;
  int nx2t = m_nx2 + 2 * m_ng;

  // Left wall
  for (int j = 1; j <= m_nx2; ++j)
    {

      rho (0, 0, j) = rho (0, 1, j);
      u (0, 0, j) = -u (0, 1, j);
      v (0, 0, j) = -v (0, 1, j);

      if (m_wall_cond == "isothermal")
        {
          T (0, 0, j) = 600 - T (0, 1, j);
        }
      else
        {
          T (0, 0, j) = T (0, 1, j);
        }

      e (0, 0, j) = 1 / (gamma - 1) * R * T (0, 0, j);
      et (0, 0, j)
          = e (0, 0, j)
            + 1 / 2 * (u (0, 0, j) * u (0, 0, j) + v (0, 0, j) * v (0, 0, j));
      p (0, 0, j) = rho (0, 0, j) * R * T (0, 0, j);
      thermlam (0, 0, j) = rho (0, 0, j) * nu / 0.7;
      mu (0, 0, j) = rho (0, 0, j) * nu;
    }

  // Right wall
  for (int j = 1; j <= m_nx2; ++j)
    {

      rho (0, nx1t - 1, j) = rho (0, nx1t - 2, j);
      u (0, nx1t - 1, j) = -u (0, nx1t - 2, j);
      v (0, nx1t - 1, j) = -v (0, nx1t - 2, j);

      if (m_wall_cond == "isothermal")
        {
          T (0, nx1t - 1, j) = 600 - T (0, nx1t - 2, j);
        }
      else
        {
          T (0, nx1t - 1, j) = T (0, nx1t - 2, j);
        }

      e (0, nx1t - 1, j) = 1 / (gamma - 1) * R * T (0, nx1t - 1, j);
      et (0, nx1t - 1, j) = e (0, nx1t - 1, j)
                            + 1 / 2
                                  * (u (0, nx1t - 1, j) * u (0, nx1t - 1, j)
                                     + v (0, nx1t - 1, j) * v (0, nx1t - 1, j));
      p (0, nx1t - 1, j) = rho (0, nx1t - 1, j) * R * T (0, nx1t - 1, j);
      thermlam (0, nx1t - 1, j) = rho (0, nx1t - 1, j) * nu / 0.7;
      mu (0, nx1t - 1, j) = rho (0, nx1t - 1, j) * nu;
    }

  // Bottom wall
  for (int i = 1; i <= m_nx1; ++i)
    {

      rho (0, i, 0) = rho (0, i, 1);
      u (0, i, 0) = -u (0, i, 1);
      v (0, i, 0) = -v (0, i, 1);

      if (m_wall_cond == "isothermal")
        {
          T (0, i, 0) = 600 - T (0, i, 1);
        }
      else
        {
          T (0, i, 0) = T (0, i, 1);
        }

      e (0, i, 0) = 1 / (gamma - 1) * R * T (0, i, 0);
      et (0, i, 0)
          = e (0, i, 0)
            + 1 / 2 * (u (0, i, 0) * u (0, i, 0) + v (0, i, 0) * v (0, i, 0));
      p (0, i, 0) = rho (0, i, 0) * R * T (0, i, 0);
      thermlam (0, i, 0) = rho (0, i, 0) * nu / 0.7;
      mu (0, i, 0) = rho (0, i, 0) * nu;
    }

  // Top wall
  for (int i = 1; i <= m_nx1; ++i)
    {

      rho (0, i, nx2t - 1) = rho (0, i, nx2t - 2);
      u (0, i, nx2t - 1) = 2 * Uw * std::sin (omega * t) - u (0, i, nx2t - 2);
      v (0, i, nx2t - 1) = -v (0, i, nx2t - 2);

      if (m_wall_cond == "isothermal")
        {
          T (0, i, nx2t - 1) = 600 - T (0, i, nx2t - 2);
        }
      else
        {
          T (0, i, nx2t - 1) = T (0, i, nx2t - 2);
        }

      e (0, i, nx2t - 1) = 1 / (gamma - 1) * R * T (0, i, nx2t - 1);
      et (0, i, nx2t - 1) = e (0, i, nx2t - 1)
                            + 1 / 2
                                  * (u (0, i, nx2t - 1) * u (0, i, nx2t - 1)
                                     + v (0, i, nx2t - 1) * v (0, i, nx2t - 1));
      p (0, i, nx2t - 1) = rho (0, i, nx2t - 1) * R * T (0, i, nx2t - 1);
      thermlam (0, i, nx2t - 1) = rho (0, i, nx2t - 1) * nu / 0.7;
      mu (0, i, nx2t - 1) = rho (0, i, nx2t - 1) * nu;
    }
}

void
SimpleMesh::new_SingleTimestep ()
{

  // Each of the empty values
  double rho_ipoh_j, rho_imoh_j, rho_i_jpoh, rho_i_jmoh;
  double u_ipoh_j, u_imoh_j, u_i_jpoh, u_i_jmoh;
  double v_ipoh_j, v_imoh_j, v_i_jpoh, v_i_jmoh;
  double p_ipoh_j, p_imoh_j, p_i_jpoh, p_i_jmoh;
  double et_ipoh_j, et_imoh_j, et_i_jpoh, et_i_jmoh;

  double u_ipoh_jpoh, u_ipoh_jmoh, u_imoh_jpoh, u_imoh_jmoh;
  double v_ipoh_jpoh, v_ipoh_jmoh, v_imoh_jpoh, v_imoh_jmoh;

  double sigma_11_ipoh_j, sigma_11_imoh_j;
  double sigma_12_ipoh_j, sigma_12_imoh_j, sigma_12_i_jpoh, sigma_12_i_jmoh;
  double sigma_22_i_jpoh, sigma_22_i_jmoh;

  double lamb_ipoh_j, lamb_imoh_j, lamb_i_jpoh, lamb_i_jmoh;
  double mu_ipoh_j, mu_imoh_j, mu_i_jpoh, mu_i_jmoh;

  double dTdx_ipoh_j, dTdx_imoh_j, dTdy_i_jpoh, dTdy_i_jmoh;

  double dudx_ipoh_j, dudx_imoh_j, dudx_i_jpoh, dudx_i_jmoh;
  double dudy_ipoh_j, dudy_imoh_j, dudy_i_jpoh, dudy_i_jmoh;

  double dvdx_ipoh_j, dvdx_imoh_j, dvdx_i_jpoh, dvdx_i_jmoh;
  double dvdy_ipoh_j, dvdy_imoh_j, dvdy_i_jpoh, dvdy_i_jmoh;

  for (int i = 1; i <= m_nx1; ++i)
    {
      for (int j = 1; j <= m_nx2; ++j)
        {

          rho_ipoh_j = (rho (0, i + 1, j) + rho (0, i, j)) / 2;
          rho_imoh_j = (rho (0, i, j) + rho (0, i - 1, j)) / 2;
          rho_i_jpoh = (rho (0, i, j + 1) + rho (0, i, j)) / 2;
          rho_i_jmoh = (rho (0, i, j) + rho (0, i, j - 1)) / 2;

          u_ipoh_j = (u (0, i + 1, j) + u (0, i, j)) / 2;
          u_imoh_j = (u (0, i, j) + u (0, i - 1, j)) / 2;
          u_i_jpoh = (u (0, i, j + 1) + u (0, i, j)) / 2;
          u_i_jmoh = (u (0, i, j) + u (0, i, j - 1)) / 2;

          v_ipoh_j = (v (0, i + 1, j) + v (0, i, j)) / 2;
          v_imoh_j = (v (0, i, j) + v (0, i - 1, j)) / 2;
          v_i_jpoh = (v (0, i, j + 1) + v (0, i, j)) / 2;
          v_i_jmoh = (v (0, i, j) + v (0, i, j - 1)) / 2;

          p_ipoh_j = (p (0, i + 1, j) + p (0, i, j)) / 2;
          p_imoh_j = (p (0, i, j) + p (0, i - 1, j)) / 2;
          p_i_jpoh = (p (0, i, j + 1) + p (0, i, j)) / 2;
          p_i_jmoh = (p (0, i, j) + p (0, i, j - 1)) / 2;

          et_ipoh_j = (et (0, i + 1, j) + et (0, i, j)) / 2;
          et_imoh_j = (et (0, i, j) + et (0, i - 1, j)) / 2;
          et_i_jpoh = (et (0, i, j + 1) + et (0, i, j)) / 2;
          et_i_jmoh = (et (0, i, j) + et (0, i, j - 1)) / 2;

          mu_ipoh_j = (mu (0, i + 1, j) + mu (0, i, j)) / 2;
          mu_imoh_j = (mu (0, i, j) + mu (0, i - 1, j)) / 2;
          mu_i_jpoh = (mu (0, i, j + 1) + mu (0, i, j)) / 2;
          mu_i_jmoh = (mu (0, i, j) + mu (0, i, j - 1)) / 2;

          lamb_ipoh_j = (thermlam (0, i + 1, j) + thermlam (0, i, j)) / 2;
          lamb_imoh_j = (thermlam (0, i, j) + thermlam (0, i - 1, j)) / 2;
          lamb_i_jpoh = (thermlam (0, i, j + 1) + thermlam (0, i, j)) / 2;
          lamb_i_jmoh = (thermlam (0, i, j) + thermlam (0, i, j - 1)) / 2;

          u_ipoh_jpoh = 1 / 4
                        * (u (0, i, j) + u (0, i + 1, j) + u (0, i + 1, j + 1)
                           + u (0, i, j + 1));
          u_ipoh_jmoh = 1 / 4
                        * (u (0, i, j) + u (0, i + 1, j) + u (0, i + 1, j - 1)
                           + u (0, i, j - 1));
          u_imoh_jpoh = 1 / 4
                        * (u (0, i, j) + u (0, i, j + 1) + u (0, i - 1, j + 1)
                           + u (0, i - 1, j));
          u_imoh_jmoh = 1 / 4
                        * (u (0, i, j) + u (0, i - 1, j) + u (0, i - 1, j - 1)
                           + u (0, i, j - 1));

          v_ipoh_jpoh = 1 / 4
                        * (v (0, i, j) + v (0, i + 1, j) + v (0, i + 1, j + 1)
                           + v (0, i, j + 1));
          v_ipoh_jmoh = 1 / 4
                        * (v (0, i, j) + v (0, i + 1, j) + v (0, i + 1, j - 1)
                           + v (0, i, j - 1));
          v_imoh_jpoh = 1 / 4
                        * (v (0, i, j) + v (0, i, j + 1) + v (0, i - 1, j + 1)
                           + v (0, i - 1, j));
          v_imoh_jmoh = 1 / 4
                        * (v (0, i, j) + v (0, i - 1, j) + v (0, i - 1, j - 1)
                           + v (0, i, j - 1));

          // left wall
          if (i == 1)
            {

              dudy_imoh_j = 0.0;
              dvdy_imoh_j = 0.0;

              dudx_imoh_j = (u (0, i, j) - u (0, i - 1, j)) / m_dx;
              dvdx_imoh_j = (v (0, i, j) - v (0, i - 1, j)) / m_dx;
            }
          else
            {

              dudy_imoh_j = (u_imoh_jpoh - u_imoh_jmoh) / m_dy;
              dvdy_imoh_j = (v_imoh_jpoh - v_imoh_jmoh) / m_dy;

              dudx_imoh_j = (u (0, i, j) - u (0, i - 1, j)) / m_dx;
              dvdx_imoh_j = (v (0, i, j) - v (0, i - 1, j)) / m_dx;
            }

          // Right wall
          if (i == m_nx1)
            {

              dudx_ipoh_j = (u (0, i + 1, j) - u (0, i, j)) / m_dx;
              dvdx_ipoh_j = (v (0, i + 1, j) - v (0, i, j)) / m_dx;

              dudy_ipoh_j = 0.0;
              dvdy_ipoh_j = 0.0;
            }
          else
            {

              dudx_ipoh_j = (u (0, i + 1, j) - u (0, i, j)) / m_dx;
              dvdx_ipoh_j = (v (0, i + 1, j) - v (0, i, j)) / m_dx;

              dudy_ipoh_j = (u_ipoh_jpoh - u_ipoh_jmoh) / m_dy;
              dvdy_ipoh_j = (v_ipoh_jpoh - v_ipoh_jmoh) / m_dy;
            }

          // Top wall
          if (j == m_nx2)
            {

              dudx_i_jpoh = 0.0;
              dvdx_i_jpoh = 0.0;

              dudy_i_jpoh = (u (0, i, j + 1) - u (0, i, j)) / m_dy;
              dvdy_i_jpoh = (v (0, i, j + 1) - v (0, i, j)) / m_dy;
            }
          else
            {

              dudx_i_jpoh = (u_ipoh_jpoh - u_imoh_jpoh);
              dvdx_i_jpoh = (v_ipoh_jpoh - v_imoh_jpoh);

              dudy_i_jpoh = (u (0, i, j + 1) - u (0, i, j)) / m_dy;
              dvdy_i_jpoh = (v (0, i, j + 1) - v (0, i, j)) / m_dy;
            }

          // Bottom wall
          if (j == 1)
            {

              dudx_i_jmoh = 0.0;
              dvdx_i_jmoh = 0.0;

              dudy_i_jmoh = (u (0, i, j) - u (0, i, j - 1)) / m_dy;
              dvdy_i_jmoh = (v (0, i, j) - v (0, i, j - 1)) / m_dy;
            }
          else
            {

              dudx_i_jmoh = (u_ipoh_jmoh - u_imoh_jmoh) / m_dx;
              dvdx_i_jmoh = (v_ipoh_jmoh - v_imoh_jmoh) / m_dx;

              dudy_i_jmoh = (u (0, i, j) - u (0, i, j - 1)) / m_dy;
              dvdy_i_jmoh = (v (0, i, j) - v (0, i, j - 1)) / m_dy;
            }

          sigma_11_ipoh_j
              = -p_ipoh_j
                + mu_ipoh_j * (4 / 3 * dudx_ipoh_j - 2 / 3 * dvdy_ipoh_j);
          sigma_11_imoh_j
              = -p_imoh_j
                + mu_imoh_j * (4 / 3 * dudx_imoh_j - 2 / 3 * dvdy_imoh_j);

          sigma_22_i_jpoh
              = -p_i_jpoh
                + mu_i_jpoh * (4 / 3 * dvdy_i_jpoh - 2 / 3 * dudx_i_jpoh);
          sigma_22_i_jmoh
              = -p_i_jmoh
                + mu_i_jmoh * (4 / 3 * dvdy_i_jmoh - 2 / 3 * dudx_i_jmoh);

          sigma_12_ipoh_j = mu_ipoh_j * (dudy_ipoh_j + dvdx_ipoh_j);
          sigma_12_imoh_j = mu_imoh_j * (dudy_imoh_j + dvdx_imoh_j);
          sigma_12_i_jpoh = mu_i_jpoh * (dudy_i_jpoh + dvdx_i_jpoh);
          sigma_12_i_jmoh = mu_i_jmoh * (dudy_i_jmoh + dvdx_i_jmoh);

          dTdx_ipoh_j = (T (0, i + 1, j) - T (0, i, j)) / m_dx;
          dTdx_imoh_j = (T (0, i, j) - T (0, i - 1, j)) / m_dx;
          dTdy_i_jpoh = (T (0, i, j + 1) - T (0, i, j)) / m_dy;
          dTdy_i_jmoh = (T (0, i, j) - T (0, i, j - 1)) / m_dy;

          F_ipoh_j (0, i, j) = rho_ipoh_j * u_ipoh_j;
          F_ipoh_j (1, i, j)
              = rho_ipoh_j * u_ipoh_j * u_ipoh_j - sigma_11_ipoh_j;
          F_ipoh_j (2, i, j)
              = rho_ipoh_j * u_ipoh_j * v_ipoh_j - sigma_12_ipoh_j;
          F_ipoh_j (3, i, j)
              = rho_ipoh_j * u_ipoh_j * et_ipoh_j - sigma_11_ipoh_j * u_ipoh_j
                - sigma_12_ipoh_j * v_ipoh_j - lamb_ipoh_j * dTdx_ipoh_j;

          F_imoh_j (0, i, j) = rho_imoh_j * u_imoh_j;
          F_imoh_j (1, i, j)
              = rho_imoh_j * u_imoh_j * u_imoh_j - sigma_11_imoh_j;
          F_imoh_j (2, i, j)
              = rho_imoh_j * u_imoh_j * v_imoh_j - sigma_12_imoh_j;
          F_imoh_j (3, i, j)
              = rho_imoh_j * u_imoh_j * et_imoh_j - sigma_11_imoh_j * u_imoh_j
                - sigma_12_imoh_j * v_imoh_j - lamb_imoh_j * dTdx_imoh_j;

          G_i_jpoh (0, i, j) = rho_i_jpoh * v_i_jpoh;
          G_i_jpoh (1, i, j)
              = rho_i_jpoh * u_i_jpoh * v_i_jpoh - sigma_12_i_jpoh;
          G_i_jpoh (2, i, j)
              = rho_i_jpoh * v_i_jpoh * v_i_jpoh - sigma_22_i_jpoh;
          G_i_jpoh (3, i, j)
              = rho_i_jpoh * v_i_jpoh * et_i_jpoh - sigma_22_i_jpoh * v_i_jpoh
                - sigma_12_i_jpoh * u_i_jpoh - lamb_i_jpoh * dTdy_i_jpoh;

          G_i_jmoh (0, i, j) = rho_i_jmoh * v_i_jmoh;
          G_i_jmoh (1, i, j)
              = rho_i_jmoh * u_i_jmoh * v_i_jmoh - sigma_12_i_jmoh;
          G_i_jmoh (2, i, j)
              = rho_i_jmoh * v_i_jmoh * v_i_jmoh - sigma_22_i_jmoh;
          G_i_jmoh (3, i, j)
              = rho_i_jmoh * v_i_jmoh * et_i_jmoh - sigma_22_i_jmoh * v_i_jmoh
                - sigma_12_i_jmoh * u_i_jmoh - lamb_i_jmoh * dTdy_i_jmoh;
        }
    }
}

void
SimpleMesh::fillPrimitive (double nu, std::string errcheck)
{
  /*
  Fills arrays with their values to be used to calculate thermal
  diffusion and viscous forces
  */

  double vrho, vu, vv, vet, ve, vP, vT;

#pragma omp simd
  for (int i = 0; i < m_nx1 + 2 * m_ng; ++i)
    {
      for (int j = 0; j < m_nx2 + 2 * m_ng; ++j)
        {

          vrho = Un (0, i, j);
          vu = Un (1, i, j) / Un (0, i, j);
          vv = Un (2, i, j) / Un (0, i, j);
          vet = Un (3, i, j) / Un (0, i, j);

          ve = vet - 1 / 2 * (vu * vu + vv * vv);
          vT = (gamma - 1.) * ve / R;
          vP = vrho * R * vT;

          rho (0, i, j) = vrho;
          u (0, i, j) = vu;
          v (0, i, j) = vv;
          et (0, i, j) = vet;

          e (0, i, j) = ve;
          T (0, i, j) = vT;
          p (0, i, j) = vP;

          mu (0, i, j) = vrho * nu;
          thermlam (0, i, j) = vrho * nu * (gamma * R / (gamma - 1.0)) / 0.7;
        }
    }
}

void
SimpleMesh::saveOutput (int iter, double t, std::string fpref)
{

  std::string fname = getFname (iter, fpref);

  int nvar = 4;
  int nx1 = m_nx1 + 2 * m_ng;
  int nx2 = m_nx2 + 2 * m_ng;

  H5::H5File file (fname.c_str (), H5F_ACC_TRUNC);

  hsize_t dims[1];
  dims[0] = nx2 * nx1;

  H5::DataSpace dataspace (1, dims);

  H5::DataSet RHO
      = file.createDataSet ("Density", H5::PredType::NATIVE_DOUBLE, dataspace);
  RHO.write (rho.m_data.data (), H5::PredType::NATIVE_DOUBLE);

  H5::DataSet XVEL
      = file.createDataSet ("VelX1", H5::PredType::NATIVE_DOUBLE, dataspace);
  XVEL.write (u.m_data.data (), H5::PredType::NATIVE_DOUBLE);

  H5::DataSet YVEL
      = file.createDataSet ("VelX2", H5::PredType::NATIVE_DOUBLE, dataspace);
  YVEL.write (v.m_data.data (), H5::PredType::NATIVE_DOUBLE);

  H5::DataSet PRES
      = file.createDataSet ("Press", H5::PredType::NATIVE_DOUBLE, dataspace);
  PRES.write (p.m_data.data (), H5::PredType::NATIVE_DOUBLE);

  H5::DataSet TEMP
      = file.createDataSet ("Temp", H5::PredType::NATIVE_DOUBLE, dataspace);
  TEMP.write (T.m_data.data (), H5::PredType::NATIVE_DOUBLE);
}

void
SimpleMesh::swapArrays ()
{
#pragma omp simd
  for (int i = 0; i <= m_nx1 + 1; ++i)
    {
      for (int j = 0; j <= m_nx2 + 1; ++j)
        {
          for (int n = 0; n < 4; ++n)
            {

              Un (n, i, j) = Unp1 (n, i, j);
            }
        }
    }
}

#endif