#include <iomanip>
#include <iostream>

#include "grid.hpp"
#include "hydro.hpp"

#include "omp.h"

void
ProblemGenerator (SimpleMesh *pmesh, double p0, double T0)
{

  double u0 = 0.0;
  double v0 = 0.0;

  double rho0 = p0 / (pmesh->R * T0);
  double e0 = 1 / (pmesh->gamma - 1) * pmesh->R * T0;

  for (int i = 0; i < pmesh->m_nx1 + 2 * pmesh->m_ng; ++i)
    {
      for (int j = 0; j < pmesh->m_nx2 + 2 * pmesh->m_ng; ++j)
        {

          pmesh->Un (0, i, j) = rho0;
          pmesh->Un (1, i, j) = rho0 * u0;
          pmesh->Un (2, i, j) = rho0 * v0;
          pmesh->Un (3, i, j) = rho0 * (e0 + 1 / 2 * (u0 * u0 + v0 * v0));
        }
    }
}

int
main (int argc, char *argv[])
{

  int nx1 = 64;
  int nx2 = 64;
  int ng = 1;

  double gamma = 1.4;
  double R = 287;

  double CFL = 0.01;
  int printfreq = 100;

  double Re = 100;
  double T0 = 300;   // K
  double Ma = 0.025; // Mach number

  double L = 1.0 / (Ma / 0.025);
  double x1max = L;
  double x1min = 0.0;
  double x2max = L;
  double x2min = 0.0;

  double a = std::sqrt (gamma * R * T0);
  double Uw = a * Ma;
  double nu = Uw * L / Re;
  double omega = 1 / (L * L) * 2 * nu;

  double omegat_max = 1.0;
  double savedt = 0.1; // in omega ts

  omp_set_num_threads (16);

  std::string time_int;
  time_int = "euler";

  std::string savepref = "test.output.";

  // --------------------------------------------------------------------------------------------//

  SimpleMesh *pmesh;
  pmesh = new SimpleMesh (nx1, nx2, x1max, x1min, x2max, x2min, ng, gamma, R, "adiabatic");

  ProblemGenerator (pmesh, 101330, 300);

  pmesh->fillPrimitive (nu, "true");

  pmesh->saveOutput (69, 0.0, savepref);

  double t = 0;
  int saveiter = 0;
  double dt = CFL * pmesh->m_dx / (Uw + a);
  int iter = 0;
  int savefreq = (int)(savedt / (omega * dt));

  while (omega * t <= omegat_max)
    {

      if (omega * t > omegat_max)
        {
          break;
        }

      if (iter % printfreq == 0)
        {
          std::cout << "Iteration : " << iter << ", Time (omega * t): " << omega * t
                    << ", dt = " << dt << std::endl;
        }

      pmesh->fillPrimitive (nu, "false");

      pmesh->new_enforceBCs (t, Uw, omega, nu);

      // Get fluxes for timestep
      pmesh->new_SingleTimestep ();

      if (time_int == "euler")
        {
          for (int i = 1; i <= pmesh->m_nx1; ++i)
            {
              for (int j = 1; j <= pmesh->m_nx2; ++j)
                {
                  for (int n = 0; n < 4; ++n)
                    {

                      pmesh->Unp1 (n, i, j)
                          = pmesh->Un (n, i, j)
                            - dt
                                  * ((pmesh->F_ipoh_j (n, i, j) - pmesh->F_imoh_j (n, i, j))
                                         / pmesh->m_dx
                                     + (pmesh->G_i_jpoh (n, i, j) - pmesh->G_i_jmoh (n, i, j))
                                           / pmesh->m_dy);
                    }
                }
            }
        }

      pmesh->swapArrays ();

      if (iter % savefreq == 0)
        {
          std::cout << "Iteration : " << iter << ", Time (omega * t): " << omega * t
                    << ", dt = " << dt << ", Data saved" << std::endl;
          pmesh->fillPrimitive (nu, "false");
          pmesh->saveOutput (saveiter, t, savepref);
          saveiter += 1;
        }

      t += dt;
      iter += 1;
    }

  return 0;
}