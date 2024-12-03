#include "core.hpp"
#include "poissonSolver.hpp"

void IcoNS::solve_time_step(Real time)
{
#ifdef PERIODIC
    const size_t start = 0;
#endif

#ifdef DIRICHELET
    const size_t start = 1;
#endif

#ifdef DIRICHELET
    const size_t end = -1;
#endif

#ifdef PERIODIC
    const size_t end = 0;
#endif

    // TODO Should be global 
    PoissonSolver poissonSolver(true,true,true);

    for (size_t i = start; i < NX + end; i++)
    {
        for (size_t j = start; j < NY + 1 + end; j++)
        {
            for (size_t k = start; k < NZ + 1 + end; k++)
            {
#ifdef PERIODIC
                Y2_x[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k] = grid.u[indexingPeriodicx(i, j, k)] + 64.0 / 120.0 * DT * functionF_u(grid.u, grid.v, grid.w, i, j, k, time) -
                                                                   64.0 / 120.0 * DT * (grid.p[indexingPeriodicp(i + 1, j, k)] - grid.p[indexingPeriodicp(i, j, k)]) / (DX);
#endif
#ifdef DIRICHELET
                Y2_x[indexingDiricheletx(i, j, k)] = grid.u[indexingDiricheletx(i, j, k)] + 64.0 / 120.0 * DT * functionF_u(grid.u, grid.v, grid.w, i, j, k, time);
#endif
            }
        }
    }

    for (size_t i = start; i < NX + 1 + end; i++)
    {
        for (size_t j = start; j < NY + end; j++)
        {
            for (size_t k = start; k < NZ + 1 + end; k++)
            {
#ifdef PERIODIC
                Y2_y[i * NY * (NZ + 1) + j * (NZ + 1) + k] = grid.v[indexingPeriodicy(i, j, k)] + 64.0 / 120.0 * DT * functionF_v(grid.u, grid.v, grid.w, i, j, k, time) -
                                                             64.0 / 120.0 * DT * (grid.p[indexingPeriodicp(i, j + 1, k)] - grid.p[indexingPeriodicp(i, j, k)]) / (DY);
#endif
#ifdef DIRICHELET
                Y2_y[indexingDirichelety(i, j, k)] = grid.v[indexingDirichelety(i, j, k)] + 64.0 / 120.0 * DT * functionF_v(grid.u, grid.v, grid.w, i, j, k, time);
#endif
            }
        }
    }

    for (size_t i = start; i < NX + 1 + end; i++)
    {
        for (size_t j = start; j < NY + 1 + end; j++)
        {
            for (size_t k = start; k < NZ + end; k++)
            {
#ifdef PERIODIC
                Y2_z[i * (NY + 1) * NZ + j * NZ + k] = grid.w[indexingPeriodicz(i, j, k)] + 64.0 / 120.0 * DT * functionF_w(grid.u, grid.v, grid.w, i, j, k, time) -
                                                       64.0 / 120.0 * DT * (grid.p[indexingPeriodicp(i, j, k + 1)] - grid.p[indexingPeriodicp(i, j, k)]) / (DZ);
#endif
#ifdef DIRICHELET
                Y2_z[indexingDiricheletz(i, j, k)] = grid.w[indexingDiricheletz(i, j, k)] + 64.0 / 120.0 * DT * functionF_w(grid.u, grid.v, grid.w, i, j, k, time);
#endif
            }
        }
    }

    for (size_t i = start; i < NX + end; i++)
    {
        for (size_t j = start; j < NY + end; j++)
        {
            for (size_t k = start; k < NZ + end; k++)
            {
#ifdef PERIODIC
                Y2_p[i * NY * NZ + j * NZ + k] = 120.0 / (64.0 * DT) * ((Y2_x[indexingPeriodicx(i, j, k)] - Y2_x[indexingPeriodicx(i - 1, j, k)]) / (DX) + (Y2_y[indexingPeriodicy(i, j, k)] - Y2_y[indexingPeriodicy(i, j - 1, k)]) / (DY) + (Y2_z[indexingPeriodicz(i, j, k)] - Y2_z[indexingPeriodicz(i, j, k - 1)]) / (DZ));
#endif
            }
        }
    }

    /////////////////////poisson_solver.solvePoisson(Y2_p);//////////////////////////////////// -> the solution is stored in Sol_p

    
    poissonSolver.solveDirichletPoisson(Y2_p,helper);

    for (size_t i = start; i < NX + end; i++)
    {
        for (size_t j = start; j < NY + 1 + end; j++)
        {
            for (size_t k = start; k < NZ + 1 + end; k++)
            {
#ifdef PERIODIC
                Y2_x[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k] = Y2_x[indexingPeriodicx(i, j, k)] - 64.0 * DT / (120.0) * (Y2_p[indexingPeriodicp(i + 1, j, k)] - Y2_p[indexingPeriodicp(i, j, k)]) / (DX);
#endif
            }
        }
    }

    for (size_t i = start; i < NX + 1 + end; i++)
    {
        for (size_t j = start; j < NY + end; j++)
        {
            for (size_t k = start; k < NZ + 1 + end; k++)
            {
#ifdef PERIODIC
                Y2_y[i * NY * (NZ + 1) + j * (NZ + 1) + k] = Y2_y[indexingPeriodicy(i, j, k)] - 64.0 * DT / (120.0) * (Y2_p[indexingPeriodicp(i, j + 1, k)] - Y2_p[indexingPeriodicp(i, j, k)]) / (DY);
#endif
            }
        }
    }

    for (size_t i = start; i < NX + 1 + end; i++)
    {
        for (size_t j = start; j < NY + 1 + end; j++)
        {
            for (size_t k = start; k < NZ + end; k++)
            {
#ifdef PERIODIC
                Y2_z[i * (NY + 1) * NZ + j * NZ + k] = Y2_z[indexingPeriodicz(i, j, k)] - 64.0 * DT / (120.0) * (Y2_p[indexingPeriodicp(i, j, k + 1)] - Y2_p[indexingPeriodicp(i, j, k)]) / (DZ);
#endif
            }
        }
    }

    for (size_t i = start; i < NX + end; i++)
    {
        for (size_t j = start; j < NY + end; j++)
        {
            for (size_t k = start; k < NZ + end; k++)
            {
#ifdef PERIODIC
                Phi_p[i * NY * NZ + j * NZ + k] = Y2_p[indexingPeriodicp(i, j, k)] + grid.p[indexingPeriodicp(i, j, k)]; // phi^2
#endif
            }
        }
    }

#ifdef DIRICHELET
    boundary.update_boundary(Y2_x, Y2_y, Y2_z, time + 64.0 / 120.0 * DT);
#endif

    for (size_t i = start; i < NX + end; i++)
    {
        for (size_t j = start; j < NY + 1 + end; j++)
        {
            for (size_t k = start; k < NZ + 1 + end; k++)
            {
#ifdef PERIODIC
                Y3_x[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k] = Y2_x[indexingPeriodicx(i, j, k)] +
                                                                   50.0 / 120.0 * DT * functionF_u(Y2_x, Y2_y, Y2_z, i, j, k, time + 64.0 / 120.0 * DT) -
                                                                   34.0 / 120.0 * DT * functionF_u(grid.u, grid.v, grid.w, i, j, k, time) -
                                                                   16.0 / 120.0 * DT * (Phi_p[indexingPeriodicp(i + 1, j, k)] - Phi_p[indexingPeriodicp(i, j, k)]) / (DX);
#endif
#ifdef DIRICHELET
                Y3_x[indexingDiricheletx(i, j, k)] = Y2_x[indexingDiricheletx(i, j, k)] +
                                                     50.0 / 120.0 * DT * functionF_u(Y2_x, Y2_y, Y2_z, i, j, k, time + 64.0 / 120.0 * DT) -
                                                     34.0 / 120.0 * DT * functionF_u(grid.u, grid.v, grid.w, i, j, k, time);
#endif
            }
        }
    }

    for (size_t i = start; i < NX + 1 + end; i++)
    {
        for (size_t j = start; j < NY + end; j++)
        {
            for (size_t k = start; k < NZ + 1 + end; k++)
            {
#ifdef PERIODIC
                Y3_y[i * NY * (NZ + 1) + j * (NZ + 1) + k] = Y2_y[indexingPeriodicy(i, j, k)] +
                                                             50.0 / 120.0 * DT * functionF_v(Y2_x, Y2_y, Y2_z, i, j, k, time + 64.0 / 120.0 * DT) -
                                                             34.0 / 120.0 * DT * functionF_v(grid.u, grid.v, grid.w, i, j, k, time) -
                                                             16.0 / 120.0 * DT * (Phi_p[indexingPeriodicp(i, j + 1, k)] - Phi_p[indexingPeriodicp(i, j, k)]) / (DY);
#endif
#ifdef DIRICHELET
                Y3_y[indexingDirichelety(i, j, k)] = Y2_y[indexingDirichelety(i, j, k)] +
                                                     50.0 / 120.0 * DT * functionF_v(Y2_x, Y2_y, Y2_z, i, j, k, time + 64.0 / 120.0 * DT) -
                                                     34.0 / 120.0 * DT * functionF_v(grid.u, grid.v, grid.w, i, j, k, time);
#endif
            }
        }
    }

    for (size_t i = start; i < NX + 1 + end; i++)
    {
        for (size_t j = start; j < NY + 1 + end; j++)
        {
            for (size_t k = start; k < NZ + end; k++)
            {
#ifdef PERIODIC
                Y3_z[i * (NY + 1) * NZ + j * NZ + k] = Y2_z[indexingPeriodicz(i, j, k)] +
                                                       50.0 / 120.0 * DT * functionF_w(Y2_x, Y2_y, Y2_z, i, j, k, time + 64.0 / 120.0 * DT) -
                                                       34.0 / 120.0 * DT * functionF_w(grid.u, grid.v, grid.w, i, j, k, time) -
                                                       16.0 / 120.0 * DT * (Phi_p[indexingPeriodicp(i, j, k + 1)] - Phi_p[indexingPeriodicp(i, j, k)]) / (DZ);
#endif
#ifdef DIRICHELET
                Y3_z[indexingDiricheletz(i, j, k)] = Y2_z[indexingDiricheletz(i, j, k)] +
                                                     50.0 / 120.0 * DT * functionF_w(Y2_x, Y2_y, Y2_z, i, j, k, time + 64.0 / 120.0 * DT) -
                                                     34.0 / 120.0 * DT * functionF_w(grid.u, grid.v, grid.w, i, j, k, time);
#endif
            }
        }
    }

    for (size_t i = start; i < NX + end; i++)
    {
        for (size_t j = start; j < NY + end; j++)
        {
            for (size_t k = start; k < NZ + end; k++)
            {
#ifdef PERIODIC
                Y3_p[i * NY * NZ + j * NZ + k] = 120.0 / (16.0 * DT) * ((Y3_x[indexingPeriodicx(i, j, k)] - Y3_x[indexingPeriodicx(i - 1, j, k)]) / (DX) + (Y3_y[indexingPeriodicy(i, j, k)] - Y3_y[indexingPeriodicy(i, j - 1, k)]) / (DY) + (Y3_z[indexingPeriodicz(i, j, k)] - Y3_z[indexingPeriodicz(i, j, k - 1)]) / (DZ));
#endif
            }
        }
    }

    /////////////////////poisson_solver.solvePoisson(Y3_p);//////////////////////////////////// -> the solution is stored in Sol_p=Y3_p
    // using Y3_p to store solution but we can use Y2_p
    poissonSolver.solveDirichletPoisson(Y3_p,helper);

    for (size_t i = start; i < NX + end; i++)
    {
        for (size_t j = start; j < NY + end + 1; j++)
        {
            for (size_t k = start; k < NZ + end + 1; k++)
            {
#ifdef PERIODIC
                Y3_x[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k] = Y3_x[indexingPeriodicx(i, j, k)] - 16.0 * DT / (120.0) * (Y3_p[indexingPeriodicp(i + 1, j, k)] - Y3_p[indexingPeriodicp(i, j, k)]) / (DX);
#endif
            }
        }
    }

    for (size_t i = start; i < NX + end + 1; i++)
    {
        for (size_t j = start; j < NY + end; j++)
        {
            for (size_t k = start; k < NZ + end + 1; k++)
            {
#ifdef PERIODIC
                Y3_y[i * NY * (NZ + 1) + j * (NZ + 1) + k] = Y3_y[indexingPeriodicy(i, j, k)] - 16.0 * DT / (120.0) * (Y3_p[indexingPeriodicp(i, j + 1, k)] - Y3_p[indexingPeriodicp(i, j, k)]) / (DY);
#endif
            }
        }
    }

    for (size_t i = start; i < NX + end + 1; i++)
    {
        for (size_t j = start; j < NY + end + 1; j++)
        {
            for (size_t k = start; k < NZ + end; k++)
            {
#ifdef PERIODIC
                Y3_z[i * (NY + 1) * NZ + j * NZ + k] = Y3_z[indexingPeriodicz(i, j, k)] - 16.0 * DT / (120.0) * (Y3_p[indexingPeriodicp(i, j, k + 1)] - Y3_p[indexingPeriodicp(i, j, k)]) / (DZ);
#endif
            }
        }
    }

    for (size_t i = start; i < NX + end; i++)
    {
        for (size_t j = start; j < NY + end; j++)
        {
            for (size_t k = start; k < NZ + end; k++)
            {
#ifdef PERIODIC
                Phi_p[i * NY * NZ + j * NZ + k] += Y3_p[indexingPeriodicp(i, j, k)]; // Phi_p=phi^3
#endif
            }
        }
    }

#ifdef DIRICHELET
    boundary.update_boundary(Y3_x, Y3_y, Y3_z, time + 80.0 / 120.0 * DT);
#endif

    for (size_t i = start; i < NX + end; i++)
    {
        for (size_t j = start; j < NY + 1 + end; j++)
        {
            for (size_t k = start; k < NZ + 1 + end; k++)
            {
#ifdef PERIODIC
                grid.u[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k] = Y3_x[indexingPeriodicx(i, j, k)] +
                                                                     90.0 / 120.0 * DT * functionF_u(Y3_x, Y3_y, Y3_z, i, j, k, time + 80.0 / 120.0 * DT) -
                                                                     50.0 / 120.0 * DT * functionF_u(Y2_x, Y2_y, Y2_z, i, j, k, time + 64.0 / 120.0 * DT) -
                                                                     40.0 / 120.0 * DT * (Phi_p[indexingPeriodicp(i + 1, j, k)] - Phi_p[indexingPeriodicp(i, j, k)]) / (DX);
#endif
#ifdef DIRICHELET
                grid.u[indexingDiricheletx(i, j, k)] = Y3_x[indexingDiricheletx(i, j, k)] +
                                                       90.0 / 120.0 * DT * functionF_u(Y3_x, Y3_y, Y3_z, i, j, k, time + 80.0 / 120.0 * DT) -
                                                       50.0 / 120.0 * DT * functionF_u(Y2_x, Y2_y, Y2_z, i, j, k, time + 64.0 / 120.0 * DT);
#endif
            }
        }
    }

    for (size_t i = start; i < NX + 1 + end; i++)
    {
        for (size_t j = start; j < NY + end; j++)
        {
            for (size_t k = start; k < NZ + 1 + end; k++)
            {
#ifdef PERIODIC
                grid.v[i * NY * (NZ + 1) + j * (NZ + 1) + k] = Y3_y[indexingPeriodicy(i, j, k)] +
                                                               90.0 / 120.0 * DT * functionF_v(Y3_x, Y3_y, Y3_z, i, j, k, time + 80.0 / 120.0 * DT) -
                                                               50.0 / 120.0 * DT * functionF_v(Y2_x, Y2_y, Y2_z, i, j, k, time + 64.0 / 120.0 * DT) -
                                                               40.0 / 120.0 * DT * (Phi_p[indexingPeriodicp(i, j + 1, k)] - Phi_p[indexingPeriodicp(i, j, k)]) / (DY);
#endif
#ifdef DIRICHELET
                grid.v[indexingDirichelety(i, j, k)] = Y3_y[indexingDirichelety(i, j, k)] +
                                                       90.0 / 120.0 * DT * functionF_v(Y3_x, Y3_y, Y3_z, i, j, k, time + 80.0 / 120.0 * DT) -
                                                       50.0 / 120.0 * DT * functionF_v(Y2_x, Y2_y, Y2_z, i, j, k, time + 64.0 / 120.0 * DT);
#endif
            }
        }
    }

    for (size_t i = start; i < NX + 1 + end; i++)
    {
        for (size_t j = start; j < NY + 1 + end; j++)
        {
            for (size_t k = start; k < NZ + end; k++)
            {
#ifdef PERIODIC
                grid.w[i * (NY + 1) * NZ + j * NZ + k] = Y3_z[indexingPeriodicz(i, j, k)] +
                                                         90.0 / 120.0 * DT * functionF_w(Y3_x, Y3_y, Y3_z, i, j, k, time + 80.0 / 120.0 * DT) -
                                                         50.0 / 120.0 * DT * functionF_w(Y2_x, Y2_y, Y2_z, i, j, k, time + 64.0 / 120.0 * DT) -
                                                         40.0 / 120.0 * DT * (Phi_p[indexingPeriodicp(i, j, k + 1)] - Phi_p[indexingPeriodicp(i, j, k)]) / (DZ);
#endif
#ifdef DIRICHELET
                grid.w[i * (NY + 1) * (NZ) + j * (NZ) + k] = Y3_z[indexingDiricheletz(i, j, k)] +
                                                             90.0 / 120.0 * DT * functionF_w(Y3_x, Y3_y, Y3_z, i, j, k, time + 80.0 / 120.0 * DT) -
                                                             50.0 / 120.0 * DT * functionF_w(Y2_x, Y2_y, Y2_z, i, j, k, time + 64.0 / 120.0 * DT);

#endif
            }
        }
    }

    for (size_t i = start; i < NX + end; i++)
    {
        for (size_t j = start; j < NY + end; j++)
        {
            for (size_t k = start; k < NZ + end; k++)
            {
                #ifdef PERIODIC
                Y2_p[i * NY * NZ + j * NZ + k] = 120.0 / (40.0 * DT) * ((grid.u[indexingPeriodicx(i, j, k)] - grid.u[indexingPeriodicx(i - 1, j, k)]) / (DX) + (grid.v[indexingPeriodicy(i, j, k)] - grid.v[indexingPeriodicy(i, j - 1, k)]) / (DY) + (grid.w[indexingPeriodicz(i, j, k)] - grid.w[indexingPeriodicz(i, j, k - 1)]) / (DZ));
                #endif
            }
        }
    }
    /////////////////////poisson_solver.solvePoisson(Y2_p);//////////////////////////////////// -> the solution is stored in Sol_p
    // using Y2_p to store solution
    poissonSolver.solveDirichletPoisson(Y2_p,helper);

    for(size_t i = start; i < NX + end; i++)
    {
        for(size_t j = start; j < NY + end + 1; j++)
        {
            for(size_t k = start; k < NZ + end + 1; k++)
            {
                #ifdef PERIODIC
                grid.u[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k] -= 40.0 * DT / (120.0) * (Y2_p[indexingPeriodicp(i + 1, j, k)] - Y2_p[indexingPeriodicp(i, j, k)]) / (DX); 
                #endif
            }
        }
    }
    for(size_t i = start; i < NX + end + 1; i++)
    {
        for(size_t j = start; j < NY + end; j++)
        {
            for(size_t k = start; k < NZ + end + 1; k++)
            {
                #ifdef PERIODIC
                grid.v[i * NY * (NZ + 1) + j * (NZ + 1) + k] -= 40.0 * DT / (120.0) * (Y2_p[indexingPeriodicp(i, j + 1, k)] - Y2_p[indexingPeriodicp(i, j, k)]) / (DY); 
                #endif
            }
        }
    }
    for(size_t i = start; i < NX + end + 1; i++)
    {
        for(size_t j = start; j < NY + end + 1; j++)
        {
            for(size_t k = start; k < NZ + end; k++)
            {
                #ifdef PERIODIC
                grid.w[i * (NY + 1) * NZ + j * NZ + k] -= 40.0 * DT / (120.0) * (Y2_p[indexingPeriodicp(i, j, k + 1)] - Y2_p[indexingPeriodicp(i, j, k)]) / (DZ); 
                #endif
            }
        }
    }
    for(size_t i = start; i < NX + end+1; i++)
    {
        for(size_t j = start; j < NY + end+1; j++)
        {
            for(size_t k = start; k < NZ + end+1; k++)
            {
                #ifdef PERIODIC
                grid.p[i * NY * NZ + j * NZ + k] += Y2_p[indexingPeriodicp(i, j, k)]; 
                #endif
            }
        }
    }
#ifdef DIRICHELET
    boundary.update_boundary(grid.u, grid.v, grid.w, time + DT);
#endif
}

Real IcoNS::functionF_u(const std::array<Real, NX *(NY + 1) * (NZ + 1)> &u, const std::array<Real, (NX + 1) * NY *(NZ + 1)> &v, const std::array<Real, (NX + 1) * (NY + 1) * NZ> &w, size_t i, size_t j, size_t k, Real t)
{
#ifdef PERIODIC
    return -(u[indexingPeriodicx(i, j, k)] * (u[indexingPeriodicx(i + 1, j, k)] - u[indexingPeriodicx(i - 1, j, k)]) / (2.0 * DX) +
             (v[indexingPeriodicy(i, j, k)] + v[indexingPeriodicy(i + 1, j, k)] + v[indexingPeriodicy(i, j - 1, k)] + v[indexingPeriodicy(i + 1, j - 1, k)]) / 4.0 * (u[indexingPeriodicx(i, j + 1, k)] - u[indexingPeriodicx(i, j - 1, k)]) / (2.0 * DY) +
             (w[indexingPeriodicz(i, j, k)] + w[indexingPeriodicz(i + 1, j, k)] + w[indexingPeriodicz(i, j, k - 1)] + w[indexingPeriodicz(i + 1, j, k - 1)]) / 4.0 * (u[indexingPeriodicx(i, j, k + 1)] - u[indexingPeriodicx(i, j, k - 1)]) / (2.0 * DZ)) +
           1 / RE * ((u[indexingPeriodicx(i + 1, j, k)] - 2 * u[indexingPeriodicx(i, j, k)] + u[indexingPeriodicx(i - 1, j, k)]) / (DX * DX) + (u[indexingPeriodicx(i, j + 1, k)] - 2 * u[indexingPeriodicx(i, j, k)] + u[indexingPeriodicx(i, j - 1, k)]) / (DY * DY) + (u[indexingPeriodicx(i, j, k + 1)] - 2 * u[indexingPeriodicx(i, j, k)] + u[indexingPeriodicx(i, j, k - 1)]) / (DZ * DZ)) +
           functionG_u(i, j, k, t);
#endif
#ifdef DIRICHELET
    return -(u[indexingDiricheletx(i, j, k)] * (u[indexingDiricheletx(i + 1, j, k)] - u[indexingDiricheletx(i - 1, j, k)]) / (2.0 * DX) +
             (v[indexingDirichelety(i, j, k)] + v[indexingDirichelety(i + 1, j, k)] + v[indexingDirichelety(i, j - 1, k)] + v[indexingDirichelety(i + 1, j - 1, k)]) / 4.0 * (u[indexingDiricheletx(i, j + 1, k)] - u[indexingDiricheletx(i, j - 1, k)]) / (2.0 * DY) +
             (w[indexingDiricheletz(i, j, k)] + w[indexingDiricheletz(i + 1, j, k)] + w[indexingDiricheletz(i, j, k - 1)] + w[indexingDiricheletz(i + 1, j, k - 1)]) / 4.0 * (u[indexingDiricheletx(i, j, k + 1)] - u[indexingDiricheletx(i, j, k - 1)]) / (2.0 * DZ)) +
           1 / RE * ((u[indexingDiricheletx(i + 1, j, k)] - 2 * u[indexingDiricheletx(i, j, k)] + u[indexingDiricheletx(i - 1, j, k)]) / (DX * DX) + (u[indexingDiricheletx(i, j + 1, k)] - 2 * u[indexingDiricheletx(i, j, k)] + u[indexingDiricheletx(i, j - 1, k)]) / (DY * DY) + (u[indexingDiricheletx(i, j, k + 1)] - 2 * u[indexingDiricheletx(i, j, k)] + u[indexingDiricheletx(i, j, k - 1)]) / (DZ * DZ)) +
           functionG_u(i, j, k, t);
#endif
}

Real IcoNS::functionF_v(const std::array<Real, NX *(NY + 1) * (NZ + 1)> &u, const std::array<Real, (NX + 1) * NY *(NZ + 1)> &v, const std::array<Real, (NX + 1) * (NY + 1) * NZ> &w, size_t i, size_t j, size_t k, Real t)
{
#ifdef PERIODIC
    return -((u[indexingPeriodicx(i, j, k)] + u[indexingPeriodicx(i, j + 1, k)] + u[indexingPeriodicx(i - 1, j, k)] + u[indexingPeriodicx(i - 1, j + 1, k)]) / 4.0 * (v[indexingPeriodicy(i + 1, j, k)] - v[indexingPeriodicy(i - 1, j, k)]) / (2.0 * DX) +
             v[indexingPeriodicy(i, j, k)] * (v[indexingPeriodicy(i, j + 1, k)] - v[indexingPeriodicy(i, j - 1, k)]) / (2.0 * DY) +
             (w[indexingPeriodicz(i, j, k)] + w[indexingPeriodicz(i, j + 1, k)] + w[indexingPeriodicz(i, j, k - 1)] + w[indexingPeriodicz(i, j + 1, k - 1)]) / 4.0 * (v[indexingPeriodicy(i, j, k + 1)] - v[indexingPeriodicy(i, j, k - 1)]) / (2.0 * DZ)) +
           1 / RE * ((v[indexingPeriodicy(i + 1, j, k)] - 2 * v[indexingPeriodicy(i, j, k)] + v[indexingPeriodicy(i - 1, j, k)]) / (DX * DX) + (v[indexingPeriodicy(i, j + 1, k)] - 2 * v[indexingPeriodicy(i, j, k)] + v[indexingPeriodicy(i, j - 1, k)]) / (DY * DY) + (v[indexingPeriodicy(i, j, k + 1)] - 2 * v[indexingPeriodicy(i, j, k)] + v[indexingPeriodicy(i, j, k - 1)]) / (DZ * DZ)) +
           functionG_v(i, j, k, t);
#endif
#ifdef DIRICHELET
    return -((u[indexingDiricheletx(i, j, k)] + u[indexingDiricheletx(i, j + 1, k)] + u[indexingDiricheletx(i - 1, j, k)] + u[indexingDiricheletx(i - 1, j + 1, k)]) / 4.0 * (v[indexingDirichelety(i + 1, j, k)] - v[indexingDirichelety(i - 1, j, k)]) / (2.0 * DX) +
             v[indexingDirichelety(i, j, k)] * (v[indexingDirichelety(i, j + 1, k)] - v[indexingDirichelety(i, j - 1, k)]) / (2.0 * DY) +
             (w[indexingDiricheletz(i, j, k)] + w[indexingDiricheletz(i, j + 1, k)] + w[indexingDiricheletz(i, j, k - 1)] + w[indexingDiricheletz(i, j + 1, k - 1)]) / 4.0 * (v[indexingDirichelety(i, j, k + 1)] - v[indexingDirichelety(i, j, k - 1)]) / (2.0 * DZ)) +
           1 / RE * ((v[indexingDirichelety(i + 1, j, k)] - 2 * v[indexingDirichelety(i, j, k)] + v[indexingDirichelety(i - 1, j, k)]) / (DX * DX) + (v[indexingDirichelety(i, j + 1, k)] - 2 * v[indexingDirichelety(i, j, k)] + v[indexingDirichelety(i, j - 1, k)]) / (DY * DY) + (v[indexingDirichelety(i, j, k + 1)] - 2 * v[indexingDirichelety(i, j, k)] + v[indexingDirichelety(i, j, k - 1)]) / (DZ * DZ)) +
           functionG_v(i, j, k, t);
#endif
}

Real IcoNS::functionF_w(const std::array<Real, NX *(NY + 1) * (NZ + 1)> &u, const std::array<Real, (NX + 1) * NY *(NZ + 1)> &v, const std::array<Real, (NX + 1) * (NY + 1) * NZ> &w, size_t i, size_t j, size_t k, Real t)

{
#ifdef PERIODIC
    return -((u[indexingPeriodicx(i, j, k)] + u[indexingPeriodicx(i, j, k + 1)] + u[indexingPeriodicx(i - 1, j, k)] + u[indexingPeriodicx(i - 1, j, k + 1)]) / 4.0 * (w[indexingPeriodicz(i + 1, j, k)] - w[indexingPeriodicz(i - 1, j, k)]) / (2.0 * DX) +
             (v[indexingPeriodicy(i, j, k)] + v[indexingPeriodicy(i, j - 1, k)] + v[indexingPeriodicy(i, j, k + 1)] + v[indexingPeriodicy(i, j - 1, k + 1)]) / 4.0 * (w[indexingPeriodicz(i, j + 1, k)] - w[indexingPeriodicz(i, j - 1, k)]) / (2.0 * DY) +
             w[indexingPeriodicz(i, j, k)] * (w[indexingPeriodicz(i, j, k + 1)] - w[indexingPeriodicz(i, j, k - 1)]) / (2.0 * DZ)) +
           1 / RE * ((w[indexingPeriodicz(i + 1, j, k)] - 2 * w[indexingPeriodicz(i, j, k)] + w[indexingPeriodicz(i - 1, j, k)]) / (DX * DX) + (w[indexingPeriodicz(i, j + 1, k)] - 2 * w[indexingPeriodicz(i, j, k)] + w[indexingPeriodicz(i, j - 1, k)]) / (DY * DY) + (w[indexingPeriodicz(i, j, k + 1)] - 2 * w[indexingPeriodicz(i, j, k)] + w[indexingPeriodicz(i, j, k - 1)]) / (DZ * DZ)) +
           functionG_w(i, j, k, t);
#endif
#ifdef DIRICHELET
    return -((u[indexingDiricheletx(i, j, k)] + u[indexingDiricheletx(i, j, k + 1)] + u[indexingDiricheletx(i - 1, j, k)] + u[indexingDiricheletx(i - 1, j, k + 1)]) / 4.0 * (w[indexingDiricheletz(i + 1, j, k)] - w[indexingDiricheletz(i - 1, j, k)]) / (2.0 * DX) +
             (v[indexingDirichelety(i, j, k)] + v[indexingDirichelety(i, j - 1, k)] + v[indexingDirichelety(i, j, k + 1)] + v[indexingDirichelety(i, j - 1, k + 1)]) / 4.0 * (w[indexingDiricheletz(i, j + 1, k)] - w[indexingDiricheletz(i, j - 1, k)]) / (2.0 * DY) +
             w[indexingDiricheletz(i, j, k)] * (w[indexingDiricheletz(i, j, k + 1)] - w[indexingDiricheletz(i, j, k - 1)]) / (2.0 * DZ)) +
           1 / RE * ((w[indexingDiricheletz(i + 1, j, k)] - 2 * w[indexingDiricheletz(i, j, k)] + w[indexingDiricheletz(i - 1, j, k)]) / (DX * DX) + (w[indexingDiricheletz(i, j + 1, k)] - 2 * w[indexingDiricheletz(i, j, k)] + w[indexingDiricheletz(i, j - 1, k)]) / (DY * DY) + (w[indexingDiricheletz(i, j, k + 1)] - 2 * w[indexingDiricheletz(i, j, k)] + w[indexingDiricheletz(i, j, k - 1)]) / (DZ * DZ)) +
           functionG_w(i, j, k, t);
#endif
}

// w[((i + 1) % (NX) + NX) * (NY + 1) * NZ + (j) % (NY)*NZ + (k) % (NZ)] - 2.0 * w[lw] + w[((i - 1) % (NX) + NX) * (NY + 1) * NZ + (j) % (NY)*NZ + (k) % (NZ)]) / (DX * DX);
// v[((i + 1) % (NX) + NX) * NY * (NZ + 1) + (j) % (NY) * (NZ + 1) + (k) % (NZ)] - 2.0 * v[lv] + v[((i - 1) % (NX) + NX) * NY * (NZ + 1) + (j) % (NY) * (NZ + 1) + (k) % (NZ)]) / (DX * DX);

/* Real IcoNS::functionF_u(const std::array<Real, NX *(NY + 1) * (NZ + 1)> &u, const std::array<Real, (NX + 1) * NY *(NZ + 1)> &v, const std::array<Real, (NX + 1) * (NY + 1) * NZ> &w, size_t i, size_t j, size_t k, Real t)
{
    size_t lu = i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k;
    size_t lv = i * NY * (NZ + 1) + j * (NZ + 1) + k;
    size_t lw = i * (NY + 1) * NZ + j * NZ + k;

    return -(u[lu] * (u[lu + (NY + 1) * (NZ + 1)] - u[lu - (NY + 1) * (NZ + 1)]) / (2.0 * DX) +
             (v[lv] + v[lv + NY * (NZ + 1)] + v[lv - (NZ + 1)] + v[lv + NY * (NZ + 1) - (NZ + 1)]) / 4.0 * (u[lu + (NZ + 1)] - u[lu - (NZ + 1)]) / (2.0 * DY) +
             (w[lw] + w[lw + (NY + 1) * NZ] + w[lw - 1] + w[lw + (NY + 1) * NZ - 1]) / 4.0 * (u[lu + 1] - u[lu - 1]) / (2.0 * DZ)) +
           1 / RE * ((u[lu + (NY + 1) * (NZ + 1)] - 2 * u[lu] + u[lu - (NY + 1) * (NZ + 1)]) / (DX * DX) + (u[lu + (NZ + 1)] - 2 * u[lu] + u[lu - (NZ + 1)]) / (DY * DY) + (u[lu + 1] - 2 * u[lu] + u[lu - 1]) / (DZ * DZ)) +
             functionG_u(i, j, k, t);
}

Real IcoNS::functionF_v(const std::array<Real, NX *(NY + 1) * (NZ + 1)> &u, const std::array<Real, (NX + 1) * NY *(NZ + 1)> &v, const std::array<Real, (NX + 1) * (NY + 1) * NZ> &w, size_t i, size_t j, size_t k, Real t)
{
    size_t lu = i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k;
    size_t lv = i * NY * (NZ + 1) + j * (NZ + 1) + k;
    size_t lw = i * (NY + 1) * NZ + j * NZ + k;

    return -((u[lu] + u[lu + (NZ + 1)] + u[lu - (NY + 1) * (NZ + 1)] + u[lu - (NY + 1) * (NZ + 1) + (NZ + 1)]) / 4.0 * ((v[lv + NY * (NZ + 1)] - v[lv - NY * (NZ + 1)]) / (2.0 * DX)) +
             v[lv] * (v[lv + (NZ + 1)] - v[lv - (NZ + 1)]) / (2.0 * DY) +
             (w[lw] + w[lw - 1] + w[lw + (NY + 1)] + w[lw + (NY + 1) - 1]) / 4.0 * (v[lv + 1] - v[lv - 1]) / (2.0 * DZ)) +
           (1.0 / RE) * ((v[lv + NY * (NZ + 1)] - 2.0 * v[lv] + v[lv - NY * (NZ + 1)]) / (DX * DX) +
                         (v[lv + (NZ + 1)] - 2.0 * v[lv] + v[lv - (NZ + 1)]) / (DY * DY) +
                         (v[lv + 1] - 2.0 * v[lv] + v[lv - 1]) / (DZ * DZ)) +
           functionG_v(i, j, k, t);
}

Real IcoNS::functionF_w(const std::array<Real, NX *(NY + 1) * (NZ + 1)> &u, const std::array<Real, (NX + 1) * NY *(NZ + 1)> &v, const std::array<Real, (NX + 1) * (NY + 1) * NZ> &w, size_t i, size_t j, size_t k, Real t)

{
    size_t lu = i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k;
    size_t lv = i * NY * (NZ + 1) + j * (NZ + 1) + k;
    size_t lw = i * (NY + 1) * NZ + j * NZ + k;

    return -((u[lu] + u[lu - (NY + 1) * (NZ + 1)] + u[lu + 1] + u[lu - (NZ + 1) * (NY + 1) + 1]) / 4.0 * (w[lw + (NY + 1) * NZ] - w[lw - (NY + 1) * NZ]) / (2.0 * DX) +
             (v[lv + 1] + v[lv - (NZ + 1) + 1] + v[lv] + v[lv - (NZ + 1)]) / 4.0 * (w[lw + NZ] - w[lw - NZ]) / (2.0 * DY) +
             w[lw] * (w[lw + 1] - w[lw - 1]) / (2.0 * DZ)) +
           (1.0 / RE) * ((w[lw + (NY + 1) * NZ] - 2.0 * w[lw] + w[lw - (NY + 1) * NZ]) / (DX * DX) +
                         (w[lw + NZ] - 2.0 * w[lw] + w[lw - NZ]) / (DY * DY) +
                         (w[lw + 1] - 2.0 * w[lw] + w[lw - 1]) / (DZ * DZ)) +
           functionG_w(i, j, k, t);
} */

Real IcoNS::functionG_u(size_t i, size_t j, size_t k, Real t)
{
    Real x = i * DX + DX / 2;
    Real y = j * DY;
    Real z = k * DZ;
    return std::sin(x) * std::cos(y) * std::sin(z) * std::cos(t) +
           std::sin(x) * std::cos(x) * std::cos(y) * std::cos(y) * std::sin(z) * std::sin(z) * std::sin(t) * std::sin(t) -
           std::sin(x) * std::cos(x) * std::sin(y) * std::sin(y) * std::sin(z) * std::sin(z) * std::sin(t) * std::sin(t) + 2.0 * std::sin(x) * std::cos(x) * std::cos(y) * std::cos(y) * std::cos(z) * std::cos(z) * std::sin(t) * std::sin(t) +
           3.0 / RE * std::sin(x) * std::cos(y) * std::sin(z) * std::sin(t) - std::sin(x) * std::cos(y) * std::sin(z) * std::sin(t);
    /*     return std::sin(z)*std::cos(t) + std::sin(y)*std::sin(t)*std::cos(z)*std::sin(t)+1.0 / RE *std::sin(z)*std::sin(t);
     */
    // return std::cos(t);
}

Real IcoNS::functionG_v(size_t i, size_t j, size_t k, Real t)
{
    Real x = i * DX;
    Real y = j * DY + DY / 2;
    Real z = k * DZ;
    return std::cos(x) * std::sin(y) * std::sin(z) * std::cos(t) -
           std::sin(x) * std::sin(x) * std::sin(y) * std::cos(y) * std::sin(z) * std::sin(z) * std::sin(t) * std::sin(t) +
           std::cos(x) * std::cos(x) * std::sin(y) * std::cos(y) * std::sin(z) * std::sin(z) * std::sin(t) * std::sin(t) +
           2.0 * std::cos(x) * std::cos(x) * std::sin(y) * std::cos(y) * std::cos(z) * std::cos(z) * std::sin(t) * std::sin(t) +
           3.0 / RE * std::cos(x) * std::sin(y) * std::sin(z) * std::sin(t) - std::cos(x) * std::sin(y) * std::sin(z) * std::sin(t);
    /*     return std::sin(x)*std::cos(t) + std::sin(z)*std::sin(t)*std::cos(x)*std::sin(t)+1.0 / RE *std::sin(x)*std::sin(t);
     */
    // return std::cos(t);
}

Real IcoNS::functionG_w(size_t i, size_t j, size_t k, Real t)
{
    Real x = i * DX;
    Real y = j * DY;
    Real z = k * DZ + DZ / 2;
    return 2.0 * std::cos(x) * std::cos(y) * std::cos(z) * std::cos(t) -
           2.0 * std::sin(x) * std::sin(x) * std::cos(y) * std::cos(y) * std::sin(z) * std::cos(z) * std::sin(t) * std::sin(t) -
           2.0 * std::cos(x) * std::cos(x) * std::sin(y) * std::sin(y) * std::sin(z) * std::cos(z) * std::sin(t) * std::sin(t) -
           4.0 * std::cos(x) * std::cos(x) * std::cos(y) * std::cos(y) * std::sin(z) * std::cos(z) * std::sin(t) * std::sin(t) +
           6.0 / RE * std::cos(x) * std::cos(y) * std::cos(z) * std::sin(t) + std::cos(x) * std::cos(y) * std::cos(z) * std::sin(t);
    /*     return std::sin(y)*std::cos(t) + std::sin(x)*std::sin(t)*std::cos(y)*std::sin(t)+1.0 / RE *std::sin(y)*std::sin(t);
     */
    // return std::cos(t);
}
