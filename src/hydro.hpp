#include "grid.hpp"

#ifndef HYDRO_HPP_
#define HYDRO_HPP_

void
enforceBCs (SimpleMesh *pmesh, double t, double Uw, double omega)
{

  double T0;
  int nx1t = pmesh->m_nx1 + 2 * pmesh->m_ng;
  int nx2t = pmesh->m_nx2 + 2 * pmesh->m_ng;

  // Left wall cond
  for (int i = 0; i < pmesh->m_nx1 + 2 * pmesh->m_ng; ++i)
    {
      pmesh->Un (0, i, 0) = pmesh->Un (0, i, 1);
      pmesh->Un (1, i, 0) = 0.0;
      pmesh->Un (2, i, 0) = 0.0;

      T0 = (pmesh->gamma - 1) / pmesh->R
           * (pmesh->Un (4, i, 1) / pmesh->Un (0, i, 1)
              - 1 / 2 * pmesh->Un (0, i, 1)
                    * (pow (pmesh->Un (1, i, 1) / pmesh->Un (0, i, 1), 2.0)
                       + pow (pmesh->Un (2, i, 1) / pmesh->Un (0, i, 1), 2.0)));

      std::cout << "Left wall cond T0 = " << pmesh->Un (4, i, 1) << std::endl;

      if (pmesh->m_wall_cond == "isothermal")
        {
          pmesh->Un (3, i, 0)
              = pmesh->Un (0, i, 1) * (1 / (pmesh->gamma - 1) * pmesh->R * (600 - T0));
        }
      else if (pmesh->m_wall_cond == "adiabatic")
        {
          pmesh->Un (3, i, 0) = pmesh->Un (0, i, 1) * (1 / (pmesh->gamma - 1) * pmesh->R * (T0));
        }
    }

  // Bottom wall cond
  for (int i = 0; i < pmesh->m_nx2 + 2 * pmesh->m_ng; ++i)
    {
      pmesh->Un (0, nx1t - 1, i) = pmesh->Un (0, nx1t - 2, i);
      pmesh->Un (1, nx1t - 1, i) = 0.0;
      pmesh->Un (2, nx1t - 1, i) = 0.0;

      T0 = (pmesh->gamma - 1) / pmesh->R
           * (pmesh->Un (4, nx1t - 2, i) / pmesh->Un (0, nx1t - 2, i)
              - 1 / 2 * pmesh->Un (0, nx1t - 2, i)
                    * (pow (pmesh->Un (1, nx1t - 2, i) / pmesh->Un (0, nx1t - 2, i), 2.0)
                       + pow (pmesh->Un (2, nx1t - 2, i) / pmesh->Un (0, nx1t - 2, i), 2.0)));

      std::cout << "Bottom wall cond T0 = " << pmesh->Un (4, nx1t - 2, i) << std::endl;

      if (pmesh->m_wall_cond == "isothermal")
        {
          pmesh->Un (3, nx1t - 1, i)
              = pmesh->Un (0, nx1t - 1, i) * (1 / (pmesh->gamma - 1) * pmesh->R * (600 - T0));
        }
      else if (pmesh->m_wall_cond == "adiabatic")
        {
          pmesh->Un (3, nx1t - 1, i)
              = pmesh->Un (0, nx1t - 1, i) * (1 / (pmesh->gamma - 1) * pmesh->R * (T0));
        }
    }

  // Right wall cond
  for (int i = 0; i < pmesh->m_nx1 + 2 * pmesh->m_ng; ++i)
    {
      pmesh->Un (0, i, nx2t - 1) = pmesh->Un (0, i, nx2t - 2);
      pmesh->Un (1, i, nx2t - 1) = 0.0;
      pmesh->Un (2, i, nx2t - 1) = 0.0;

      T0 = (pmesh->gamma - 1) / pmesh->R
           * (pmesh->Un (4, i, nx2t - 2) / pmesh->Un (0, i, nx2t - 2)
              - 1 / 2 * pmesh->Un (0, i, nx2t - 2)
                    * (pow (pmesh->Un (1, i, nx2t - 2) / pmesh->Un (0, i, nx2t - 2), 2.0)
                       + pow (pmesh->Un (2, i, nx2t - 2) / pmesh->Un (0, i, nx2t - 2), 2.0)));

      std::cout << "Right wall cond T0 = " << pmesh->Un (4, i, nx2t - 2) << std::endl;

      if (pmesh->m_wall_cond == "isothermal")
        {
          pmesh->Un (3, i, nx2t - 1)
              = pmesh->Un (0, i, nx2t - 1) * (1 / (pmesh->gamma - 1) * pmesh->R * (600 - T0));
        }
      else if (pmesh->m_wall_cond == "adiabatic")
        {
          pmesh->Un (3, i, nx2t - 1)
              = pmesh->Un (0, i, nx2t - 1) * (1 / (pmesh->gamma - 1) * pmesh->R * (T0));
        }
    }

  // Top lid cond
  for (int i = 0; i < pmesh->m_nx2 + 2 * pmesh->m_ng; ++i)
    {
      pmesh->Un (0, 0, i) = pmesh->Un (0, 1, i);
      pmesh->Un (1, 0, i) = Uw * std::sin (omega * t);
      pmesh->Un (2, 0, i) = 0.0;

      T0 = (pmesh->gamma - 1) / pmesh->R
           * (pmesh->Un (4, 1, i) / pmesh->Un (0, 1, i)
              - 1 / 2 * pmesh->Un (0, 1, i)
                    * (pow (pmesh->Un (1, 1, i) / pmesh->Un (0, 1, i), 2.0)
                       + pow (pmesh->Un (2, 1, i) / pmesh->Un (0, 1, i), 2.0)));

      std::cout << "Top wall cond T0 = " << pmesh->Un (4, 1, i) << std::endl;

      if (pmesh->m_wall_cond == "isothermal")
        {
          pmesh->Un (3, 0, i) = pmesh->Un (0, 0, i)
                                * (1 / (pmesh->gamma - 1) * pmesh->R * (600 - T0)
                                   + 1 / 2 * pmesh->Un (0, 0, i)
                                         * pow (pmesh->Un (1, 0, i) / pmesh->Un (0, 0, i), 2.0));
        }
      else if (pmesh->m_wall_cond == "adiabatic")
        {
          pmesh->Un (3, 0, i) = pmesh->Un (0, 0, i)
                                * (1 / (pmesh->gamma - 1) * pmesh->R * (T0)
                                   + 1 / 2 * pmesh->Un (0, 0, i)
                                         * pow (pmesh->Un (1, 0, i) / pmesh->Un (0, 0, i), 2.0));
        }
    }
}

void
enforceBCs_T (SimpleMesh *pmesh, double t, double Uw, double omega)
{

  double T0;
  int nx1t = pmesh->m_nx1 + 2 * pmesh->m_ng;
  int nx2t = pmesh->m_nx2 + 2 * pmesh->m_ng;

  // Left wall cond
  for (int i = 0; i < pmesh->m_nx1 + 2 * pmesh->m_ng; ++i)
    {
      pmesh->Un (0, 0, i) = pmesh->Un (0, 1, i);
      pmesh->Un (1, 0, i) = 0.0;
      pmesh->Un (2, 0, i) = 0.0;

      T0 = (pmesh->gamma - 1) / pmesh->R
           * (pmesh->Un (4, 1, i)
              - 1 / 2 * pmesh->Un (0, 1, i)
                    * (pow (pmesh->Un (1, 1, i) / pmesh->Un (0, 1, i), 2.0)
                       + pow (pmesh->Un (2, 1, i) / pmesh->Un (0, 1, i), 2.0)));

      if (pmesh->m_wall_cond == "isothermal")
        {
          pmesh->Un (3, 0, i)
              = pmesh->Un (0, 0, i) * (1 / (pmesh->gamma - 1) * pmesh->R * (600 - T0));
        }
      else if (pmesh->m_wall_cond == "adiabatic")
        {
          pmesh->Un (3, 0, i) = pmesh->Un (0, 0, i) * (1 / (pmesh->gamma - 1) * pmesh->R * (T0));
        }
    }

  // Bottom wall cond
  for (int i = 0; i < pmesh->m_nx2 + 2 * pmesh->m_ng; ++i)
    {
      pmesh->Un (0, i, nx1t - 1) = pmesh->Un (0, i, nx1t - 2);
      pmesh->Un (1, i, nx1t - 1) = 0.0;
      pmesh->Un (2, i, nx1t - 1) = 0.0;

      T0 = (pmesh->gamma - 1) / pmesh->R
           * (pmesh->Un (4, i, nx1t - 2)
              - 1 / 2 * pmesh->Un (0, i, nx1t - 2)
                    * (pow (pmesh->Un (1, i, nx1t - 2) / pmesh->Un (0, i, nx1t - 2), 2.0)
                       + pow (pmesh->Un (2, i, nx1t - 2) / pmesh->Un (0, i, nx1t - 2), 2.0)));

      if (pmesh->m_wall_cond == "isothermal")
        {
          pmesh->Un (3, i, nx1t - 1)
              = pmesh->Un (0, i, nx1t - 1) * (1 / (pmesh->gamma - 1) * pmesh->R * (600 - T0));
        }
      else if (pmesh->m_wall_cond == "adiabatic")
        {
          pmesh->Un (3, i, nx1t - 1)
              = pmesh->Un (0, i, nx1t - 1) * (1 / (pmesh->gamma - 1) * pmesh->R * (T0));
        }
    }

  // Right wall cond
  for (int i = 0; i < pmesh->m_nx1 + 2 * pmesh->m_ng; ++i)
    {
      pmesh->Un (0, nx2t - 1, i) = pmesh->Un (0, nx2t - 2, i);
      pmesh->Un (1, nx2t - 1, i) = 0.0;
      pmesh->Un (2, nx2t - 1, i) = 0.0;

      T0 = (pmesh->gamma - 1) / pmesh->R
           * (pmesh->Un (4, nx2t - 2, i)
              - 1 / 2 * pmesh->Un (0, nx2t - 2, i)
                    * (pow (pmesh->Un (1, nx2t - 2, i) / pmesh->Un (0, nx2t - 2, i), 2.0)
                       + pow (pmesh->Un (2, nx2t - 2, i) / pmesh->Un (0, nx2t - 2, i), 2.0)));

      if (pmesh->m_wall_cond == "isothermal")
        {
          pmesh->Un (3, nx2t - 1, i)
              = pmesh->Un (0, nx2t - 1, i) * (1 / (pmesh->gamma - 1) * pmesh->R * (600 - T0));
        }
      else if (pmesh->m_wall_cond == "adiabatic")
        {
          pmesh->Un (3, nx2t - 1, i)
              = pmesh->Un (0, nx2t - 1, i) * (1 / (pmesh->gamma - 1) * pmesh->R * (T0));
        }
    }

  // Top lid cond
  for (int i = 0; i < pmesh->m_nx2 + 2 * pmesh->m_ng; ++i)
    {
      pmesh->Un (0, i, 0) = pmesh->Un (0, i, 1);
      pmesh->Un (1, i, 0) = pmesh->Un (0, i, 1) * Uw * std::sin (omega * t);
      pmesh->Un (2, i, 0) = 0.0;

      T0 = (pmesh->gamma - 1) / pmesh->R
           * (pmesh->Un (4, i, 1)
              - 1 / 2 * pmesh->Un (0, i, 1)
                    * (pow (pmesh->Un (1, i, 1) / pmesh->Un (0, i, 1), 2.0)
                       + pow (pmesh->Un (2, i, 1) / pmesh->Un (0, i, 1), 2.0)));

      if (pmesh->m_wall_cond == "isothermal")
        {
          pmesh->Un (3, i, 0) = pmesh->Un (0, i, 0)
                                * (1 / (pmesh->gamma - 1) * pmesh->R * (600 - T0)
                                   + 1 / 2 * pmesh->Un (0, i, 0)
                                         * pow (pmesh->Un (1, i, 0) / pmesh->Un (0, i, 0), 2.0));
        }
      else if (pmesh->m_wall_cond == "adiabatic")
        {
          pmesh->Un (3, i, 0) = pmesh->Un (0, i, 0)
                                * (1 / (pmesh->gamma - 1) * pmesh->R * (T0)
                                   + 1 / 2 * pmesh->Un (0, i, 0)
                                         * pow (pmesh->Un (1, i, 0) / pmesh->Un (0, i, 0), 2.0));
        }
    }
}

double
calculateTimestep (SimpleMesh *pmesh, double CFL)
{

  double maxvel = 0.0;

  for (int i = 0; i < pmesh->m_nx1 + 2 * pmesh->m_ng; ++i)
    {
      for (int j = 0; j < pmesh->m_nx2 + 2 * pmesh->m_ng; ++j)
        {

          double a = std::sqrt (pmesh->gamma * pmesh->gamma * pmesh->T (0, i, j));

          double cellvelx
              = std::max (std::abs (pmesh->u (0, i, j)), std::abs (pmesh->u (0, i, j)) + a);
          double cellvely
              = std::max (std::abs (pmesh->v (0, i, j)), std::abs (pmesh->v (0, i, j)) + a);

          if (cellvelx > maxvel)
            {
              maxvel = cellvelx;
            }

          if (cellvely > maxvel)
            {
              maxvel = cellvely;
            }
        }
    }

  return CFL * pmesh->m_dx / maxvel;
}

#endif