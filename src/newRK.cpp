#include "core.hpp"
#include "poissonSolver.hpp"

void IcoNS::solve_time_step(Real time)
{
// do we need to exlude the boundaries from the computation?
// we could compute, use for poisson solver and then overwrite the boundaries
#ifdef PERIODIC
    const int start = 0;
#endif

#ifdef DIRICHELET
    const int start = 0;
#endif

#ifdef DIRICHELET
    const int end = 0;
#endif

#ifdef PERIODIC
    const int end = 0;
#endif

    // TODO Should be global
    PoissonSolver poissonSolver(true, true, true);

    for (int i = start; i < NX + end; i++)
    {
        for (int j = start; j < NY + 1 + end; j++)
        {
            for (int k = start; k < NZ + 1 + end; k++)
            {
#ifdef PERIODIC
                Y2_x[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k] = grid.u[indexingPeriodicx(i, j, k)] + 64.0 / 120.0 * DT * functionF_u(grid.u, grid.v, grid.w, i, j, k, time) -
                                                                   64.0 / 120.0 * DT * /*d_Px(i,j,k,time);*/ (grid.p[(i + 1) * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k] - grid.p[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k]) / (DX);
#endif
#ifdef DIRICHELET
                Y2_x[indexingDiricheletx(i, j, k)] = grid.u[indexingDiricheletx(i, j, k)] + 64.0 / 120.0 * DT * functionF_u(grid.u, grid.v, grid.w, i, j, k, time) -
                                                     64.0 / 120.0 * DT * /*d_Px(i,j,k,time);*/ (grid.p[(i + 1) * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k] - grid.p[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k]) / (DX);
#endif
            }
        }
    }

    for (int i = start; i < NX + 1 + end; i++)
    {
        for (int j = start; j < NY + end; j++)
        {
            for (int k = start; k < NZ + 1 + end; k++)
            {
#ifdef PERIODIC
                Y2_y[i * NY * (NZ + 1) + j * (NZ + 1) + k] = grid.v[indexingPeriodicy(i, j, k)] + 64.0 / 120.0 * DT * functionF_v(grid.u, grid.v, grid.w, i, j, k, time) -
                                                             64.0 / 120.0 * DT * /*d_Py(i,j,k,time);*/ (grid.p[i * (NY + 1) * (NZ + 1) + (j + 1) * (NZ + 1) + k] - grid.p[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k]) / (DY);
#endif
#ifdef DIRICHELET
                Y2_y[indexingDirichelety(i, j, k)] = grid.v[indexingDirichelety(i, j, k)] + 64.0 / 120.0 * DT * functionF_v(grid.u, grid.v, grid.w, i, j, k, time) -
                                                     64.0 / 120.0 * DT * /*d_Py(i,j,k,time);*/ (grid.p[i * (NY + 1) * (NZ + 1) + (j + 1) * (NZ + 1) + k] - grid.p[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k]) / (DY);
#endif
            }
        }
    }

    for (int i = start; i < NX + 1 + end; i++)
    {
        for (int j = start; j < NY + 1 + end; j++)
        {
            for (int k = start; k < NZ + end; k++)
            {
#ifdef PERIODIC
                Y2_z[i * (NY + 1) * NZ + j * NZ + k] = grid.w[indexingPeriodicz(i, j, k)] + 64.0 / 120.0 * DT * functionF_w(grid.u, grid.v, grid.w, i, j, k, time) -
                                                       64.0 / 120.0 * DT * /*d_Pz(i,j,k,time);*/ (grid.p[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k + 1] - grid.p[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k]) / (DZ);
#endif
#ifdef DIRICHELET
                Y2_z[indexingDiricheletz(i, j, k)] = grid.w[indexingDiricheletz(i, j, k)] + 64.0 / 120.0 * DT * functionF_w(grid.u, grid.v, grid.w, i, j, k, time) -
                                                     64.0 / 120.0 * DT * /*d_Pz(i,j,k,time);*/ (grid.p[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k + 1] - grid.p[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k]) / (DZ);
#endif
            }
        }
    }

#ifdef PERIODIC
    for (int i = 0; i < NX; i++)
    {
        for (int j = 0; j < NY; j++)
        {
            for (int k = 0; k < NZ; k++)
            {
                Y2_p[i * NY * NZ + j * NZ + k] = 120.0 / (64.0 * DT) * ((Y2_x[indexingPeriodicx(i, j, k)] - Y2_x[indexingPeriodicx(i - 1, j, k)]) / (DX) + (Y2_y[indexingPeriodicy(i, j, k)] - Y2_y[indexingPeriodicy(i, j - 1, k)]) / (DY) + (Y2_z[indexingPeriodicz(i, j, k)] - Y2_z[indexingPeriodicz(i, j, k - 1)]) / (DZ));
            }
        }
    }
    /////////////////////poisson_solver.solvePoisson(step 1);//////////////////////////////////// -> the solution is stored in Sol_p=Y2_p
    poissonSolver.solveDirichletPoisson(Y2_p, helper);
#endif

#ifdef DIRICHELET
    // boundary.update_boundary(Y2_x, Y2_y, Y2_z, time + 64.0 / 120.0 * DT);
    //   boundary.divergence(Y2_x, Y2_y, Y2_z, Y2_p, time + 64.0 / 120.0 * DT, 64.0);
    for (int i = 0; i < NX + 1; i++)
    {
        for (int j = 0; j < NY + 1; j++)
        {
            for (int k = 0; k < NZ + 1; k++)
            {
                // Y2_p[i * (NY+1) * (NZ+1) + j * (NZ+1) + k] = 120.0 / (64.0 * DT) * ((Y2_x[indexingDiricheletx(i, j, k)] - Y2_x[indexingDiricheletx(i - 1, j, k)]) / (DX) + (Y2_y[indexingDirichelety(i, j, k)] - Y2_y[indexingDirichelety(i, j - 1, k)]) / (DY) + (Y2_z[indexingDiricheletz(i, j, k)] - Y2_z[indexingDiricheletz(i, j, k - 1)]) / (DZ));
                Real diffx = 0.0;
                Real diffy = 0.0;
                Real diffz = 0.0;

                if (i == 0)
                {
                    // continue;
                    diffx = (-8 * boundary.boundary_value_u[0]->value(i - 0.5, j, k, time + 64.0 / 120.0 * DT) + 9 * Y2_x[indexingDiricheletx(i, j, k)] - Y2_x[indexingDiricheletx(i + 1, j, k)]) / (3 * DX);
                    // diffx =  64.0 * DT / (120.0) * std::cos(i * DX) * std::cos(j * DY) * std::cos(k * DZ) * (sin(time) - sin(time + 64.0 * DT / (120.0)));
                }
                else if (i == NX)
                {
                    // continue;
                    diffx = (8 * boundary.boundary_value_u[0]->value(NX - 0.5, j, k, time + 64.0 / 120.0 * DT) - 9 * Y2_x[indexingDiricheletx(i - 1, j, k)] + Y2_x[indexingDiricheletx(i - 2, j, k)]) / (3 * DX);
                    // diffx =  64.0 * DT / (120.0) * std::cos(i * DX) * std::cos(j * DY) * std::cos(k * DZ) * (sin(time) - sin(time + 64.0 * DT / (120.0)));
                }
                else
                {
                    diffx = (Y2_x[indexingDiricheletx(i, j, k)] - Y2_x[indexingDiricheletx(i - 1, j, k)]) / (DX);
                    // diffx =  64.0 * DT / (120.0) * std::cos(i * DX) * std::cos(j * DY) * std::cos(k * DZ) * (sin(time) - sin(time + 64.0 * DT / (120.0)));
                }

                if (j == 0)
                {
                    // continue;
                    diffy = (-8 * boundary.boundary_value_v[0]->value(i, j - 0.5, k, time + 64.0 / 120.0 * DT) + 9 * Y2_y[indexingDirichelety(i, j, k)] - Y2_y[indexingDirichelety(i, j + 1, k)]) / (3 * DY);
                    // diffy =  64.0 * DT / (120.0) * std::cos(i * DX) * std::cos(j * DY) * std::cos(k * DZ) * (sin(time) - sin(time + 64.0 * DT / (120.0)));
                }
                else if (j == NY)
                {
                    // continue;
                    diffy = (8 * boundary.boundary_value_v[0]->value(i, j - 0.5, k, time + 64.0 / 120.0 * DT) - 9 * Y2_y[indexingDirichelety(i, j - 1, k)] + Y2_y[indexingDirichelety(i, j - 2, k)]) / (3 * DY);
                    // diffy =  64.0 * DT / (120.0) * std::cos(i * DX) * std::cos(j * DY) * std::cos(k * DZ) * (sin(time) - sin(time + 64.0 * DT / (120.0)));
                }
                else
                {
                    diffy = (Y2_y[indexingDirichelety(i, j, k)] - Y2_y[indexingDirichelety(i, j - 1, k)]) / (DY);
                    // diffy =  64.0 * DT / (120.0) * std::cos(i * DX) * std::cos(j * DY) * std::cos(k * DZ) * (sin(time) - sin(time + 64.0 * DT / (120.0)));
                }

                if (k == 0)
                {
                    // continue;
                    diffz = (-8 * boundary.boundary_value_w[0]->value(i, j, k - 0.5, time + 64.0 / 120.0 * DT) + 9 * Y2_z[indexingDiricheletz(i, j, k)] - Y2_z[indexingDiricheletz(i, j, k + 1)]) / (3 * DZ);
                    // diffz =  64.0 * DT / (120.0) * std::cos(i * DX) * std::cos(j * DY) * std::cos(k * DZ) * (sin(time) - sin(time + 64.0 * DT / (120.0)));
                }
                else if (k == NZ)
                {
                    // continue;
                    diffz = (8 * boundary.boundary_value_w[0]->value(i, j, k - 0.5, time + 64.0 / 120.0 * DT) - 9 * Y2_z[indexingDiricheletz(i, j, k - 1)] + Y2_z[indexingDiricheletz(i, j, k - 2)]) / (3 * DZ);
                    // diffz =  64.0 * DT / (120.0) * std::cos(i * DX) * std::cos(j * DY) * std::cos(k * DZ) * (sin(time) - sin(time + 64.0 * DT / (120.0)));
                }
                else
                {
                    diffz = (Y2_z[indexingDiricheletz(i, j, k)] - Y2_z[indexingDiricheletz(i, j, k - 1)]) / (DZ);
                    // diffz =  64.0 * DT / (120.0) * std::cos(i * DX) * std::cos(j * DY) * std::cos(k * DZ) * (sin(time) - sin(time + 64.0 * DT / (120.0)));
                }

                Y2_p[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k] = 120.0 / (64.0 * DT) * (diffx + diffy + diffz);
            }
        }
    }
    // boundary.divergence(Y2_x, Y2_y, Y2_z, Y2_p, time + 64.0 / 120.0 * DT, 64.0);

    poissonSolver.solveNeumannPoisson(Y2_p);
#endif

    for (int i = start; i < NX + end; i++)
    {
        for (int j = start; j < NY + 1 + end; j++)
        {
            for (int k = start; k < NZ + 1 + end; k++)
            {
#ifdef PERIODIC
                Y2_x[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k] = Y2_x[indexingPeriodicx(i, j, k)] - 64.0 * DT / (120.0) * /*(d_Px(i,j,k,time+64.0 * DT / (120.0))-d_Px(i,j,k,time));*/ (Y2_p[indexingPeriodicp(i + 1, j, k)] - Y2_p[indexingPeriodicp(i, j, k)]) / (DX);
#endif
#ifdef DIRICHELET
                Y2_x[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k] = Y2_x[indexingDiricheletx(i, j, k)] - 64.0 * DT / (120.0) * /*(d_Px(i,j,k,time+64.0 * DT / (120.0))-d_Px(i,j,k,time));*/ (Y2_p[indexingDiricheletp(i + 1, j, k)] - Y2_p[indexingDiricheletp(i, j, k)]) / (DX);
#endif
            }
        }
    }

    for (int i = start; i < NX + 1 + end; i++)
    {
        for (int j = start; j < NY + end; j++)
        {
            for (int k = start; k < NZ + 1 + end; k++)
            {
#ifdef PERIODIC
                Y2_y[i * NY * (NZ + 1) + j * (NZ + 1) + k] = Y2_y[indexingPeriodicy(i, j, k)] - 64.0 * DT / (120.0) * /*(d_Py(i,j,k,time+64.0 * DT / (120.0))-d_Py(i,j,k,time));*/ (Y2_p[indexingPeriodicp(i, j + 1, k)] - Y2_p[indexingPeriodicp(i, j, k)]) / (DY);
#endif
#ifdef DIRICHELET
                Y2_y[i * NY * (NZ + 1) + j * (NZ + 1) + k] = Y2_y[indexingDirichelety(i, j, k)] - 64.0 * DT / (120.0) * /*(d_Py(i,j,k,time+64.0 * DT / (120.0))-d_Py(i,j,k,time));*/ (Y2_p[indexingDiricheletp(i, j + 1, k)] - Y2_p[indexingDiricheletp(i, j, k)]) / (DY);
#endif
            }
        }
    }

    for (int i = start; i < NX + 1 + end; i++)
    {
        for (int j = start; j < NY + 1 + end; j++)
        {
            for (int k = start; k < NZ + end; k++)
            {
#ifdef PERIODIC
                Y2_z[i * (NY + 1) * NZ + j * NZ + k] = Y2_z[indexingPeriodicz(i, j, k)] - 64.0 * DT / (120.0) * /*(d_Pz(i,j,k,time+64.0 * DT / (120.0))-d_Pz(i,j,k,time));*/ (Y2_p[indexingPeriodicp(i, j, k + 1)] - Y2_p[indexingPeriodicp(i, j, k)]) / (DZ);
#endif
#ifdef DIRICHELET
                Y2_z[i * (NY + 1) * NZ + j * NZ + k] = Y2_z[indexingDiricheletz(i, j, k)] - 64.0 * DT / (120.0) * /*(d_Pz(i,j,k,time+64.0 * DT / (120.0))-d_Pz(i,j,k,time));*/ (Y2_p[indexingDiricheletp(i, j, k + 1)] - Y2_p[indexingDiricheletp(i, j, k)]) / (DZ);
#endif
            }
        }
    }

    for (int i = 0; i < NX + 1; i++)
    {
        for (int j = 0; j < NY + 1; j++)
        {
            for (int k = 0; k < NZ + 1; k++)
            {
#ifdef PERIODIC
                Phi_p[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k] = Y2_p[indexingPeriodicp(i, j, k)] + grid.p[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k]; // phi^2
#endif
#ifdef DIRICHELET
                Phi_p[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k] = Y2_p[indexingDiricheletp(i, j, k)] + grid.p[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k]; // phi^2
#endif
            }
        }
    }

#ifdef DIRICHELET
    pressionCorrection(Phi_p);
    boundary.update_boundary(Y2_x, Y2_y, Y2_z, time + 64.0 / 120.0 * DT);
#endif

    for (int i = start; i < NX + end; i++)
    {
        for (int j = start; j < NY + 1 + end; j++)
        {
            for (int k = start; k < NZ + 1 + end; k++)
            {
#ifdef PERIODIC
                Y3_x[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k] = Y2_x[indexingPeriodicx(i, j, k)] +
                                                                   50.0 / 120.0 * DT * functionF_u(Y2_x, Y2_y, Y2_z, i, j, k, time + 64.0 / 120.0 * DT) -
                                                                   34.0 / 120.0 * DT * functionF_u(grid.u, grid.v, grid.w, i, j, k, time) -
                                                                   16.0 / 120.0 * DT * /*d_Px(i,j,k,time + 64.0/120.0*DT);*/ (Phi_p[(i + 1) * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k] - Phi_p[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k]) / (DX);
#endif
#ifdef DIRICHELET
                Y3_x[indexingDiricheletx(i, j, k)] = Y2_x[indexingDiricheletx(i, j, k)] +
                                                     50.0 / 120.0 * DT * functionF_u(Y2_x, Y2_y, Y2_z, i, j, k, time + 64.0 / 120.0 * DT) -
                                                     34.0 / 120.0 * DT * functionF_u(grid.u, grid.v, grid.w, i, j, k, time) -
                                                     16.0 / 120.0 * DT * /*d_Px(i,j,k,time + 64.0/120.0*DT);*/ (Phi_p[(i + 1) * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k] - Phi_p[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k]) / (DX);
#endif
            }
        }
    }

    for (int i = start; i < NX + 1 + end; i++)
    {
        for (int j = start; j < NY + end; j++)
        {
            for (int k = start; k < NZ + 1 + end; k++)
            {
#ifdef PERIODIC
                Y3_y[i * NY * (NZ + 1) + j * (NZ + 1) + k] = Y2_y[indexingPeriodicy(i, j, k)] +
                                                             50.0 / 120.0 * DT * functionF_v(Y2_x, Y2_y, Y2_z, i, j, k, time + 64.0 / 120.0 * DT) -
                                                             34.0 / 120.0 * DT * functionF_v(grid.u, grid.v, grid.w, i, j, k, time) -
                                                             16.0 / 120.0 * DT * /*d_Py(i,j,k,time + 64.0/120.0*DT);*/ (Phi_p[i * (NY + 1) * (NZ + 1) + (j + 1) * (NZ + 1) + k] - Phi_p[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k]) / (DY);
#endif
#ifdef DIRICHELET
                Y3_y[indexingDirichelety(i, j, k)] = Y2_y[indexingDirichelety(i, j, k)] +
                                                     50.0 / 120.0 * DT * functionF_v(Y2_x, Y2_y, Y2_z, i, j, k, time + 64.0 / 120.0 * DT) -
                                                     34.0 / 120.0 * DT * functionF_v(grid.u, grid.v, grid.w, i, j, k, time) -
                                                     16.0 / 120.0 * DT * /*d_Py(i,j,k,time + 64.0/120.0*DT);*/ (Phi_p[i * (NY + 1) * (NZ + 1) + (j + 1) * (NZ + 1) + k] - Phi_p[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k]) / (DY);
#endif
            }
        }
    }

    for (int i = start; i < NX + 1 + end; i++)
    {
        for (int j = start; j < NY + 1 + end; j++)
        {
            for (int k = start; k < NZ + end; k++)
            {
#ifdef PERIODIC
                Y3_z[i * (NY + 1) * NZ + j * NZ + k] = Y2_z[indexingPeriodicz(i, j, k)] +
                                                       50.0 / 120.0 * DT * functionF_w(Y2_x, Y2_y, Y2_z, i, j, k, time + 64.0 / 120.0 * DT) -
                                                       34.0 / 120.0 * DT * functionF_w(grid.u, grid.v, grid.w, i, j, k, time) -
                                                       16.0 / 120.0 * DT * /*d_Pz(i,j,k,time + 64.0/120.0*DT);*/ (Phi_p[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k + 1] - Phi_p[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k]) / (DZ);
#endif
#ifdef DIRICHELET
                Y3_z[indexingDiricheletz(i, j, k)] = Y2_z[indexingDiricheletz(i, j, k)] +
                                                     50.0 / 120.0 * DT * functionF_w(Y2_x, Y2_y, Y2_z, i, j, k, time + 64.0 / 120.0 * DT) -
                                                     34.0 / 120.0 * DT * functionF_w(grid.u, grid.v, grid.w, i, j, k, time) -
                                                     16.0 / 120.0 * DT * /*d_Pz(i,j,k,time + 64.0/120.0*DT);*/ (Phi_p[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k + 1] - Phi_p[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k]) / (DZ);
#endif
            }
        }
    }

#ifdef PERIODIC
    for (int i = 0; i < NX; i++)
    {
        for (int j = 0; j < NY; j++)
        {
            for (int k = 0; k < NZ; k++)
            {
                Y2_p[i * NY * NZ + j * NZ + k] = 120.0 / (16.0 * DT) * ((Y3_x[indexingPeriodicx(i, j, k)] - Y3_x[indexingPeriodicx(i - 1, j, k)]) / (DX) + (Y3_y[indexingPeriodicy(i, j, k)] - Y3_y[indexingPeriodicy(i, j - 1, k)]) / (DY) + (Y3_z[indexingPeriodicz(i, j, k)] - Y3_z[indexingPeriodicz(i, j, k - 1)]) / (DZ));
            }
        }
    }
    /////////////////////poisson_solver.solvePoisson(step 2);//////////////////////////////////// -> the solution is stored in Sol_p=Y2_p
    poissonSolver.solveDirichletPoisson(Y2_p, helper);
#endif
#ifdef DIRICHELET
    // boundary.update_boundary(Y3_x, Y3_y, Y3_z, time + 80.0 / 120.0 * DT);
    //   boundary.divergence(Y3_x, Y3_y, Y3_z, Y2_p, time + 80.0 / 120.0 * DT, 64.0);
    for (int i = 0; i < NX + 1; i++)
    {
        for (int j = 0; j < NY + 1; j++)
        {
            for (int k = 0; k < NZ + 1; k++)
            {
                // Y2_p[i * (NY+1) * (NZ+1) + j * (NZ+1) + k] = 120.0 / (64.0 * DT) * ((Y3_x[indexingDiricheletx(i, j, k)] - Y3_x[indexingDiricheletx(i - 1, j, k)]) / (DX) + (Y3_y[indexingDirichelety(i, j, k)] - Y3_y[indexingDirichelety(i, j - 1, k)]) / (DY) + (Y3_z[indexingDiricheletz(i, j, k)] - Y3_z[indexingDiricheletz(i, j, k - 1)]) / (DZ));
                Real diffx = 0.0;
                Real diffy = 0.0;
                Real diffz = 0.0;

                if (i == 0)
                {
                    // continue;
                    diffx = (-8 * boundary.boundary_value_u[0]->value(i - 0.5, j, k, time + 80.0 / 120.0 * DT) + 9 * Y3_x[indexingDiricheletx(i, j, k)] - Y3_x[indexingDiricheletx(i + 1, j, k)]) / (3 * DX);
                }
                else if (i == NX)
                {
                    // continue;
                    diffx = (8 * boundary.boundary_value_u[0]->value(NX - 0.5, j, k, time + 80.0 / 120.0 * DT) - 9 * Y3_x[indexingDiricheletx(i - 1, j, k)] + Y3_x[indexingDiricheletx(i - 2, j, k)]) / (3 * DX);
                }
                else
                {
                    diffx = (Y3_x[indexingDiricheletx(i, j, k)] - Y3_x[indexingDiricheletx(i - 1, j, k)]) / (DX);
                }

                if (j == 0)
                {
                    // continue;
                    diffy = (-8 * boundary.boundary_value_v[0]->value(i, j - 0.5, k, time + 80.0 / 120.0 * DT) + 9 * Y3_y[indexingDirichelety(i, j, k)] - Y3_y[indexingDirichelety(i, j + 1, k)]) / (3 * DY);
                }
                else if (j == NY)
                {
                    // continue;
                    diffy = (8 * boundary.boundary_value_v[0]->value(i, j - 0.5, k, time + 80.0 / 120.0 * DT) - 9 * Y3_y[indexingDirichelety(i, j - 1, k)] + Y3_y[indexingDirichelety(i, j - 2, k)]) / (3 * DY);
                }
                else
                {
                    diffy = (Y3_y[indexingDirichelety(i, j, k)] - Y3_y[indexingDirichelety(i, j - 1, k)]) / (DY);
                }

                if (k == 0)
                {
                    // continue;
                    diffz = (-8 * boundary.boundary_value_w[0]->value(i, j, k - 0.5, time + 80.0 / 120.0 * DT) + 9 * Y3_z[indexingDiricheletz(i, j, k)] - Y3_z[indexingDiricheletz(i, j, k + 1)]) / (3 * DZ);
                }
                else if (k == NZ)
                {
                    // continue;
                    diffz = (8 * boundary.boundary_value_w[0]->value(i, j, k - 0.5, time + 80.0 / 120.0 * DT) - 9 * Y3_z[indexingDiricheletz(i, j, k - 1)] + Y3_z[indexingDiricheletz(i, j, k - 2)]) / (3 * DZ);
                }
                else
                {
                    diffz = (Y3_z[indexingDiricheletz(i, j, k)] - Y3_z[indexingDiricheletz(i, j, k - 1)]) / (DZ);
                }

                Y2_p[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k] = 120.0 / (16.0 * DT) * (diffx + diffy + diffz);
            }
        }
    }
    // boundary.divergence(Y3_x, Y3_y, Y3_z, Y2_p, time + 80.0 / 120.0 * DT, 64.0);
    poissonSolver.solveNeumannPoisson(Y2_p);

#endif

    for (int i = start; i < NX + end; i++)
    {
        for (int j = start; j < NY + 1 + end; j++)
        {
            for (int k = start; k < NZ + 1 + end; k++)
            {
#ifdef PERIODIC
                Y3_x[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k] = Y3_x[indexingPeriodicx(i, j, k)] - 16.0 * DT / (120.0) * /*(d_Px(i,j,k,time+80.0 * DT / (120.0))-d_Px(i,j,k,time+64.0 * DT / (120.0)));*/ (Y2_p[indexingPeriodicp(i + 1, j, k)] - Y2_p[indexingPeriodicp(i, j, k)]) / (DX);
#endif
#ifdef DIRICHELET
                Y3_x[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k] = Y3_x[indexingDiricheletx(i, j, k)] - 16.0 * DT / (120.0) * /*(d_Px(i,j,k,time+80.0 * DT / (120.0))-d_Px(i,j,k,time+64.0 * DT / (120.0)));*/ (Y2_p[indexingDiricheletp(i + 1, j, k)] - Y2_p[indexingDiricheletp(i, j, k)]) / (DX);
#endif
            }
        }
    }

    for (int i = start; i < NX + 1 + end; i++)
    {
        for (int j = start; j < NY + end; j++)
        {
            for (int k = start; k < NZ + 1 + end; k++)
            {
#ifdef PERIODIC
                Y3_y[i * NY * (NZ + 1) + j * (NZ + 1) + k] = Y3_y[indexingPeriodicy(i, j, k)] - 16.0 * DT / (120.0) * /*(d_Py(i,j,k,time+80.0 * DT / (120.0))-d_Py(i,j,k,time+64.0 * DT / (120.0)));*/ (Y2_p[indexingPeriodicp(i, j + 1, k)] - Y2_p[indexingPeriodicp(i, j, k)]) / (DY);
#endif
#ifdef DIRICHELET
                Y3_y[i * NY * (NZ + 1) + j * (NZ + 1) + k] = Y3_y[indexingDirichelety(i, j, k)] - 16.0 * DT / (120.0) * /*(d_Py(i,j,k,time+80.0 * DT / (120.0))-d_Py(i,j,k,time+64.0 * DT / (120.0)));*/ (Y2_p[indexingDiricheletp(i, j + 1, k)] - Y2_p[indexingDiricheletp(i, j, k)]) / (DY);
#endif
            }
        }
    }

    for (int i = start; i < NX + 1 + end; i++)
    {
        for (int j = start; j < NY + 1 + end; j++)
        {
            for (int k = start; k < NZ + end; k++)
            {
#ifdef PERIODIC
                Y3_z[i * (NY + 1) * NZ + j * NZ + k] = Y3_z[indexingPeriodicz(i, j, k)] - 16.0 * DT / (120.0) * /*(d_Pz(i,j,k,time+80.0 * DT / (120.0))-d_Pz(i,j,k,time+64.0 * DT / (120.0)));*/ (Y2_p[indexingPeriodicp(i, j, k + 1)] - Y2_p[indexingPeriodicp(i, j, k)]) / (DZ);
#endif
#ifdef DIRICHELET
                Y3_z[i * (NY + 1) * NZ + j * NZ + k] = Y3_z[indexingDiricheletz(i, j, k)] - 16.0 * DT / (120.0) * /*(d_Pz(i,j,k,time+80.0 * DT / (120.0))-d_Pz(i,j,k,time+64.0 * DT / (120.0)));*/ (Y2_p[indexingDiricheletp(i, j, k + 1)] - Y2_p[indexingDiricheletp(i, j, k)]) / (DZ);
#endif
            }
        }
    }

    for (int i = 0; i < NX + 1; i++)
    {
        for (int j = 0; j < NY + 1; j++)
        {
            for (int k = 0; k < NZ + 1; k++)
            {
#ifdef PERIODIC
                Phi_p[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k] = Y2_p[indexingPeriodicp(i, j, k)] + Phi_p[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k]; // Phi_p=phi^3
#endif
#ifdef DIRICHELET
                Phi_p[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k] = Y2_p[indexingDiricheletp(i, j, k)] + Phi_p[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k]; // Phi_p=phi^3
#endif
            }
        }
    }

#ifdef DIRICHELET
    pressionCorrection(Phi_p);
    boundary.update_boundary(Y3_x, Y3_y, Y3_z, time + 80.0 / 120.0 * DT);
#endif

    for (int i = start; i < NX + end; i++)
    {
        for (int j = start; j < NY + 1 + end; j++)
        {
            for (int k = start; k < NZ + 1 + end; k++)
            {
#ifdef PERIODIC
                grid.u[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k] = Y3_x[indexingPeriodicx(i, j, k)] +
                                                                     90.0 / 120.0 * DT * functionF_u(Y3_x, Y3_y, Y3_z, i, j, k, time + 80.0 / 120.0 * DT) -
                                                                     50.0 / 120.0 * DT * functionF_u(Y2_x, Y2_y, Y2_z, i, j, k, time + 64.0 / 120.0 * DT) -
                                                                     40.0 / 120.0 * DT * /*d_Px(i,j,k,time + 80.0/120.0*DT);*/ (Phi_p[(i + 1) * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k] - Phi_p[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k]) / (DX);
#endif
#ifdef DIRICHELET
                grid.u[indexingDiricheletx(i, j, k)] = Y3_x[indexingDiricheletx(i, j, k)] +
                                                       90.0 / 120.0 * DT * functionF_u(Y3_x, Y3_y, Y3_z, i, j, k, time + 80.0 / 120.0 * DT) -
                                                       50.0 / 120.0 * DT * functionF_u(Y2_x, Y2_y, Y2_z, i, j, k, time + 64.0 / 120.0 * DT) -
                                                       40.0 / 120.0 * DT * /*d_Px(i,j,k,time + 80.0/120.0*DT);*/ (Phi_p[(i + 1) * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k] - Phi_p[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k]) / (DX);
#endif
            }
        }
    }

    for (int i = start; i < NX + 1 + end; i++)
    {
        for (int j = start; j < NY + end; j++)
        {
            for (int k = start; k < NZ + 1 + end; k++)
            {
#ifdef PERIODIC
                grid.v[i * NY * (NZ + 1) + j * (NZ + 1) + k] = Y3_y[indexingPeriodicy(i, j, k)] +
                                                               90.0 / 120.0 * DT * functionF_v(Y3_x, Y3_y, Y3_z, i, j, k, time + 80.0 / 120.0 * DT) -
                                                               50.0 / 120.0 * DT * functionF_v(Y2_x, Y2_y, Y2_z, i, j, k, time + 64.0 / 120.0 * DT) -
                                                               40.0 / 120.0 * DT * /*d_Py(i,j,k,time + 80.0/120.0*DT);*/ (Phi_p[i * (NY + 1) * (NZ + 1) + (j + 1) * (NZ + 1) + k] - Phi_p[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k]) / (DY);
#endif
#ifdef DIRICHELET
                grid.v[indexingDirichelety(i, j, k)] = Y3_y[indexingDirichelety(i, j, k)] +
                                                       90.0 / 120.0 * DT * functionF_v(Y3_x, Y3_y, Y3_z, i, j, k, time + 80.0 / 120.0 * DT) -
                                                       50.0 / 120.0 * DT * functionF_v(Y2_x, Y2_y, Y2_z, i, j, k, time + 64.0 / 120.0 * DT) -
                                                       40.0 / 120.0 * DT * /*d_Py(i,j,k,time + 80.0/120.0*DT);*/ (Phi_p[i * (NY + 1) * (NZ + 1) + (j + 1) * (NZ + 1) + k] - Phi_p[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k]) / (DY);
#endif
            }
        }
    }

    for (int i = start; i < NX + 1 + end; i++)
    {
        for (int j = start; j < NY + 1 + end; j++)
        {
            for (int k = start; k < NZ + end; k++)
            {
#ifdef PERIODIC
                grid.w[i * (NY + 1) * NZ + j * NZ + k] = Y3_z[indexingPeriodicz(i, j, k)] +
                                                         90.0 / 120.0 * DT * functionF_w(Y3_x, Y3_y, Y3_z, i, j, k, time + 80.0 / 120.0 * DT) -
                                                         50.0 / 120.0 * DT * functionF_w(Y2_x, Y2_y, Y2_z, i, j, k, time + 64.0 / 120.0 * DT) -
                                                         40.0 / 120.0 * DT * /*d_Pz(i,j,k,time + 80.0/120.0*DT);*/ (Phi_p[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k + 1] - Phi_p[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k]) / (DZ);
#endif
#ifdef DIRICHELET
                grid.w[i * (NY + 1) * (NZ) + j * (NZ) + k] = Y3_z[indexingDiricheletz(i, j, k)] +
                                                             90.0 / 120.0 * DT * functionF_w(Y3_x, Y3_y, Y3_z, i, j, k, time + 80.0 / 120.0 * DT) -
                                                             50.0 / 120.0 * DT * functionF_w(Y2_x, Y2_y, Y2_z, i, j, k, time + 64.0 / 120.0 * DT) -
                                                             40.0 / 120.0 * DT * /*d_Pz(i,j,k,time + 80.0/120.0*DT);*/ (Phi_p[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k + 1] - Phi_p[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k]) / (DZ);
#endif
            }
        }
    }

#ifdef PERIODIC
    for (int i = 0; i < NX; i++)
    {
        for (int j = 0; j < NY; j++)
        {
            for (int k = 0; k < NZ; k++)
            {
                Y2_p[i * NY * NZ + j * NZ + k] = 120.0 / (40.0 * DT) * ((grid.u[indexingPeriodicx(i, j, k)] - grid.u[indexingPeriodicx(i - 1, j, k)]) / (DX) + (grid.v[indexingPeriodicy(i, j, k)] - grid.v[indexingPeriodicy(i, j - 1, k)]) / (DY) + (grid.w[indexingPeriodicz(i, j, k)] - grid.w[indexingPeriodicz(i, j, k - 1)]) / (DZ));
            }
        }
    }
    /////////////////////poisson_solver.solvePoisson(step 3);//////////////////////////////////// -> the solution is stored in Sol_p=Y2_p
    poissonSolver.solveDirichletPoisson(Y2_p, helper);
#endif
#ifdef DIRICHELET
    // boundary.update_boundary(grid.u, grid.v, grid.w, time + DT);
    //   boundary.divergence(grid.u, grid.v, grid.w, Y2_p, time + DT, 64.0);
    for (int i = 0; i < NX + 1; i++)
    {
        for (int j = 0; j < NY + 1; j++)
        {
            for (int k = 0; k < NZ + 1; k++)
            {
                // Y2_p[i * (NY+1) * (NZ+1) + j * (NZ+1) + k] = 120.0 / (64.0 * DT) * ((grid.u[indexingDiricheletx(i, j, k)] - grid.u[indexingDiricheletx(i - 1, j, k)]) / (DX) + (grid.v[indexingDirichelety(i, j, k)] - grid.v[indexingDirichelety(i, j - 1, k)]) / (DY) + (grid.w[indexingDiricheletz(i, j, k)] - grid.w[indexingDiricheletz(i, j, k - 1)]) / (DZ));
                Real diffx = 0.0;
                Real diffy = 0.0;
                Real diffz = 0.0;

                if (i == 0)
                {
                    // continue;
                    diffx = (-8 * boundary.boundary_value_u[0]->value(i - 0.5, j, k, time + DT) + 9 * grid.u[indexingDiricheletx(i, j, k)] - grid.u[indexingDiricheletx(i + 1, j, k)]) / (3 * DX);
                }
                else if (i == NX)
                {
                    // continue;
                    diffx = (8 * boundary.boundary_value_u[0]->value(NX - 0.5, j, k, time + DT) - 9 * grid.u[indexingDiricheletx(i - 1, j, k)] + grid.u[indexingDiricheletx(i - 2, j, k)]) / (3 * DX);
                }
                else
                {
                    diffx = (grid.u[indexingDiricheletx(i, j, k)] - grid.u[indexingDiricheletx(i - 1, j, k)]) / (DX);
                }

                if (j == 0)
                {
                    // continue;
                    diffy = (-8 * boundary.boundary_value_v[0]->value(i, j - 0.5, k, time + DT) + 9 * grid.v[indexingDirichelety(i, j, k)] - grid.v[indexingDirichelety(i, j + 1, k)]) / (3 * DY);
                }
                else if (j == NY)
                {
                    // continue;
                    diffy = (8 * boundary.boundary_value_v[0]->value(i, j - 0.5, k, time + DT) - 9 * grid.v[indexingDirichelety(i, j - 1, k)] + grid.v[indexingDirichelety(i, j - 2, k)]) / (3 * DY);
                }
                else
                {
                    diffy = (grid.v[indexingDirichelety(i, j, k)] - grid.v[indexingDirichelety(i, j - 1, k)]) / (DY);
                }

                if (k == 0)
                {
                    // continue;
                    diffz = (-8 * boundary.boundary_value_w[0]->value(i, j, k - 0.5, time + DT) + 9 * grid.w[indexingDiricheletz(i, j, k)] - grid.w[indexingDiricheletz(i, j, k + 1)]) / (3 * DZ);
                }
                else if (k == NZ)
                {
                    // continue;
                    diffz = (8 * boundary.boundary_value_w[0]->value(i, j, k - 0.5, time + DT) - 9 * grid.w[indexingDiricheletz(i, j, k - 1)] + grid.w[indexingDiricheletz(i, j, k - 2)]) / (3 * DZ);
                }
                else
                {
                    diffz = (grid.w[indexingDiricheletz(i, j, k)] - grid.w[indexingDiricheletz(i, j, k - 1)]) / (DZ);
                }

                Y2_p[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k] = 120.0 / (40.0 * DT) * (diffx + diffy + diffz);
            }
        }
    }
    // boundary.divergence(grid.u, grid.v, grid.w, Y2_p, time + DT, 64.0);

    sum = 0.0;

    for (int i = 0; i < NX + 1; i++)
    {
        for (int j = 0; j < NY + 1; j++)
        {
            for (int k = 0; k < NZ + 1; k++)
            {
                sum += (Y2_p[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k] - 3.0 * std::cos(i * DX) * std::cos(j * DY) * std::cos(k * DZ) * (sin(time + 80.0 * DT / (120.0)) - sin(time + DT))) *
                       (Y2_p[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k] - 3.0 * std::cos(i * DX) * std::cos(j * DY) * std::cos(k * DZ) * (sin(time + 80.0 * DT / (120.0)) - sin(time + DT))) * DX * DY * DZ;
                // Y2_p[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k] = 3.0 * std::cos(i * DX) * std::cos(j * DY) * std::cos(k * DZ) * (sin(time + 80.0 * DT / (120.0)) - sin(time + DT));
            }
        }
    }

    std::cout << "sum = " << sqrt(sum) << std::endl;

    poissonSolver.solveNeumannPoisson(Y2_p);

    sum = 0.0;
    // compute the average value of Y2_p
    for (int i = 0; i < NX + 1; i++)
    {
        for (int j = 0; j < NY + 1; j++)
        {
            for (int k = 0; k < NZ + 1; k++)
            {
                if (Y2_p[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k] > 0.5)
                {
                    // std::cout << "Y2_p[" << i << "][" << j << "][" << k << "] = (2)" << Y2_p[i * (NY+1) * (NZ+1) + j * (NZ+1) + k] << std::endl;
                }
            }
        }
    }

#endif

    for (int i = start; i < NX + end; i++)
    {
        for (int j = start; j < NY + 1 + end; j++)
        {
            for (int k = start; k < NZ + 1 + end; k++)
            {
#ifdef PERIODIC
                grid.u[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k] -= 40.0 * DT / (120.0) * /*(d_Px(i,j,k,time+DT)-d_Px(i,j,k,time+80.0 * DT / (120.0)));*/ (Y2_p[indexingPeriodicp(i + 1, j, k)] - Y2_p[indexingPeriodicp(i, j, k)]) / (DX);
#endif
#ifdef DIRICHELET
                grid.u[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k] -= 40.0 * DT / (120.0) * /*(d_Px(i,j,k,time+DT)-d_Px(i,j,k,time+80.0 * DT / (120.0)));*/ (Y2_p[indexingDiricheletp(i + 1, j, k)] - Y2_p[indexingDiricheletp(i, j, k)]) / (DX);
#endif
            }
        }
    }
    for (int i = start; i < NX + 1 + end; i++)
    {
        for (int j = start; j < NY + end; j++)
        {
            for (int k = start; k < NZ + 1 + end; k++)
            {
#ifdef PERIODIC
                grid.v[i * NY * (NZ + 1) + j * (NZ + 1) + k] -= 40.0 * DT / (120.0) * /*(d_Py(i,j,k,time+DT)-d_Py(i,j,k,time+80.0 * DT / (120.0)));*/ (Y2_p[indexingPeriodicp(i, j + 1, k)] - Y2_p[indexingPeriodicp(i, j, k)]) / (DY);
#endif
#ifdef DIRICHELET
                grid.v[i * NY * (NZ + 1) + j * (NZ + 1) + k] -= 40.0 * DT / (120.0) * /*(d_Py(i,j,k,time+DT)-d_Py(i,j,k,time+80.0 * DT / (120.0)));*/ (Y2_p[indexingDiricheletp(i, j + 1, k)] - Y2_p[indexingDiricheletp(i, j, k)]) / (DY);
#endif
            }
        }
    }
    for (int i = start; i < NX + 1 + end; i++)
    {
        for (int j = start; j < NY + 1 + end; j++)
        {
            for (int k = start; k < NZ + end; k++)
            {
#ifdef PERIODIC
                grid.w[i * (NY + 1) * NZ + j * NZ + k] -= 40.0 * DT / (120.0) * /*(d_Pz(i,j,k,time+DT)-d_Pz(i,j,k,time+80.0 * DT / (120.0)));*/ (Y2_p[indexingPeriodicp(i, j, k + 1)] - Y2_p[indexingPeriodicp(i, j, k)]) / (DZ);
#endif
#ifdef DIRICHELET
                grid.w[i * (NY + 1) * NZ + j * NZ + k] -= 40.0 * DT / (120.0) * /*(d_Pz(i,j,k,time+DT)-d_Pz(i,j,k,time+80.0 * DT / (120.0)));*/ (Y2_p[indexingDiricheletp(i, j, k + 1)] - Y2_p[indexingDiricheletp(i, j, k)]) / (DZ);
#endif
            }
        }
    }
    for (int i = 0; i < NX + 1; i++)
    {
        for (int j = 0; j < NY + 1; j++)
        {
            for (int k = 0; k < NZ + 1; k++)
            {
#ifdef PERIODIC
                grid.p[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k] = Y2_p[indexingPeriodicp(i, j, k)] + Phi_p[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k];
#endif
#ifdef DIRICHELET
                grid.p[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k] = Y2_p[indexingDiricheletp(i, j, k)] + Phi_p[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k];
#endif
            }
        }
    }

#ifdef DIRICHELET
    pressionCorrection(grid.p);
    boundary.update_boundary(grid.u, grid.v, grid.w, time + DT);
#endif
}

Real IcoNS::functionF_u(const std::array<Real, NX *(NY + 1) * (NZ + 1)> &u, const std::array<Real, (NX + 1) * NY *(NZ + 1)> &v, const std::array<Real, (NX + 1) * (NY + 1) * NZ> &w, int i, int j, int k, Real t)
{
#ifdef PERIODIC
    return -(u[indexingPeriodicx(i, j, k)] * (u[indexingPeriodicx(i + 1, j, k)] - u[indexingPeriodicx(i - 1, j, k)]) / (2.0 * DX) +
             (v[indexingPeriodicy(i, j, k)] + v[indexingPeriodicy(i + 1, j, k)] + v[indexingPeriodicy(i, j - 1, k)] + v[indexingPeriodicy(i + 1, j - 1, k)]) / 4.0 * (u[indexingPeriodicx(i, j + 1, k)] - u[indexingPeriodicx(i, j - 1, k)]) / (2.0 * DY) +
             (w[indexingPeriodicz(i, j, k)] + w[indexingPeriodicz(i + 1, j, k)] + w[indexingPeriodicz(i, j, k - 1)] + w[indexingPeriodicz(i + 1, j, k - 1)]) / 4.0 * (u[indexingPeriodicx(i, j, k + 1)] - u[indexingPeriodicx(i, j, k - 1)]) / (2.0 * DZ)) +
           1 / RE * ((u[indexingPeriodicx(i + 1, j, k)] - 2 * u[indexingPeriodicx(i, j, k)] + u[indexingPeriodicx(i - 1, j, k)]) / (DX * DX) + (u[indexingPeriodicx(i, j + 1, k)] - 2 * u[indexingPeriodicx(i, j, k)] + u[indexingPeriodicx(i, j - 1, k)]) / (DY * DY) + (u[indexingPeriodicx(i, j, k + 1)] - 2 * u[indexingPeriodicx(i, j, k)] + u[indexingPeriodicx(i, j, k - 1)]) / (DZ * DZ)) +
           functionG_u(i, j, k, t);
#endif
#ifdef DIRICHELET
    /* return -(u[indexingDiricheletx(i, j, k)] * (u[indexingDiricheletx(i + 1, j, k)] - u[indexingDiricheletx(i - 1, j, k)]) / (2.0 * DX) +
             (v[indexingDirichelety(i, j, k)] + v[indexingDirichelety(i + 1, j, k)] + v[indexingDirichelety(i, j - 1, k)] + v[indexingDirichelety(i + 1, j - 1, k)]) / 4.0 * (u[indexingDiricheletx(i, j + 1, k)] - u[indexingDiricheletx(i, j - 1, k)]) / (2.0 * DY) +
             (w[indexingDiricheletz(i, j, k)] + w[indexingDiricheletz(i + 1, j, k)] + w[indexingDiricheletz(i, j, k - 1)] + w[indexingDiricheletz(i + 1, j, k - 1)]) / 4.0 * (u[indexingDiricheletx(i, j, k + 1)] - u[indexingDiricheletx(i, j, k - 1)]) / (2.0 * DZ)) +
           1 / RE * ((u[indexingDiricheletx(i + 1, j, k)] - 2 * u[indexingDiricheletx(i, j, k)] + u[indexingDiricheletx(i - 1, j, k)]) / (DX * DX) + (u[indexingDiricheletx(i, j + 1, k)] - 2 * u[indexingDiricheletx(i, j, k)] + u[indexingDiricheletx(i, j - 1, k)]) / (DY * DY) + (u[indexingDiricheletx(i, j, k + 1)] - 2 * u[indexingDiricheletx(i, j, k)] + u[indexingDiricheletx(i, j, k - 1)]) / (DZ * DZ)) +
           functionG_u(i, j, k, t); */
    Real diffx, diffy, diffz, u_value, v_value, w_value, diffsecx, diffsecy, diffsecz;

    u_value = u[indexingDiricheletx(i, j, k)];

    if (j == 0 || j == NY)
    {
        v_value = boundary.boundary_value_v[0]->value(i, j - 0.5, k, t);
    }
    else
    {
        v_value = (v[indexingDirichelety(i, j, k)] + v[indexingDirichelety(i + 1, j, k)] + v[indexingDirichelety(i, j - 1, k)] + v[indexingDirichelety(i + 1, j - 1, k)]) / 4.0;
    }

    if (k == 0 || k == NZ)
    {
        w_value = boundary.boundary_value_w[0]->value(i, j, k - 0.5, t);
    }
    else
    {
        w_value = (w[indexingDiricheletz(i, j, k)] + w[indexingDiricheletz(i + 1, j, k)] + w[indexingDiricheletz(i, j, k - 1)] + w[indexingDiricheletz(i + 1, j, k - 1)]) / 4.0;
    }

    if (i == 0)
    {
        diffx = (-4 * boundary.boundary_value_u[0]->value(i - 0.5, j, k, t) + 3 * u[indexingDiricheletx(i, j, k)] + u[indexingDiricheletx(i + 1, j, k)]) / (3.0 * DX);
    }
    else if (i == NX - 1)
    {
        diffx = (4 * boundary.boundary_value_u[0]->value(NX - 0.5, j, k, t) - 3 * u[indexingDiricheletx(i, j, k)] - u[indexingDiricheletx(i - 1, j, k)]) / (3.0 * DX);
    }
    else
    {
        diffx = (u[indexingDiricheletx(i + 1, j, k)] - u[indexingDiricheletx(i - 1, j, k)]) / (2.0 * DX);
    }

    if (j == 0)
    {
        diffy = (-3 * u[indexingDiricheletx(i, j, k)] + 4 * u[indexingDiricheletx(i, j + 1, k)] - u[indexingDiricheletx(i, j + 2, k)]) / (2 * DY);
    }
    else if (j == NY)
    {
        diffy = (3 * u[indexingDiricheletx(i, j, k)] - 4 * u[indexingDiricheletx(i, j - 1, k)] + u[indexingDiricheletx(i, j - 2, k)]) / (2 * DY);
    }
    else
    {
        diffy = (u[indexingDiricheletx(i, j + 1, k)] - u[indexingDiricheletx(i, j - 1, k)]) / (2.0 * DY);
    }

    if (k == 0)
    {
        diffz = (-3 * u[indexingDiricheletx(i, j, k)] + 4 * u[indexingDiricheletx(i, j, k + 1)] - u[indexingDiricheletx(i, j, k + 2)]) / (2 * DZ);
    }
    else if (k == NZ)
    {
        diffz = (3 * u[indexingDiricheletx(i, j, k)] - 4 * u[indexingDiricheletx(i, j, k - 1)] + u[indexingDiricheletx(i, j, k - 2)]) / (2 * DZ);
    }
    else
    {
        diffz = (u[indexingDiricheletx(i, j, k + 1)] - u[indexingDiricheletx(i, j, k - 1)]) / (2.0 * DZ);
    }

    if (i == 0)
    {
        diffsecx = (16 * boundary.boundary_value_u[0]->value(i - 0.5, j, k, t) - 25 * u[indexingDiricheletx(i, j, k)] + 10 * u[indexingDiricheletx(i + 1, j, k)] - u[indexingDiricheletx(i + 2, j, k)]) / (5.0 * DX * DX);
    }
    else if (i == NX - 1)
    {
        diffsecx = (-16 * boundary.boundary_value_u[0]->value(i + 0.5, j, k, t) + 25 * u[indexingDiricheletx(i, j, k)] - 10 * u[indexingDiricheletx(i - 1, j, k)] + u[indexingDiricheletx(i - 2, j, k)]) / (5.0 * DX * DX);
    }
    else
    {
        diffsecx = (u[indexingDiricheletx(i + 1, j, k)] - 2 * u[indexingDiricheletx(i, j, k)] + u[indexingDiricheletx(i - 1, j, k)]) / (DX * DX);
    }

    if (j == 0)
    {
        diffsecy = (2 * u[indexingDiricheletx(i, j, k)] - 5 * u[indexingDiricheletx(i, j + 1, k)] + 4 * u[indexingDiricheletx(i, j + 2, k)] - u[indexingDiricheletx(i, j + 3, k)]) / (DY * DY);
    }
    else if (j == NY)
    {
        diffsecy = (-2 * u[indexingDiricheletx(i, j, k)] + 5 * u[indexingDiricheletx(i, j - 1, k)] - 4 * u[indexingDiricheletx(i, j - 2, k)] + u[indexingDiricheletx(i, j - 3, k)]) / (DY * DY);
    }
    else
    {
        diffsecy = (u[indexingDiricheletx(i, j + 1, k)] - 2 * u[indexingDiricheletx(i, j, k)] + u[indexingDiricheletx(i, j - 1, k)]) / (DY * DY);
    }

    if (k == 0)
    {
        diffsecz = (2 * u[indexingDiricheletx(i, j, k)] - 5 * u[indexingDiricheletx(i, j, k + 1)] + 4 * u[indexingDiricheletx(i, j, k + 2)] - u[indexingDiricheletx(i, j, k + 3)]) / (DZ * DZ);
    }
    else if (k == NZ)
    {
        diffsecz = (-2 * u[indexingDiricheletx(i, j, k)] + 5 * u[indexingDiricheletx(i, j, k - 1)] - 4 * u[indexingDiricheletx(i, j, k - 2)] + u[indexingDiricheletx(i, j, k - 3)]) / (DZ * DZ);
    }
    else
    {
        diffsecz = (u[indexingDiricheletx(i, j, k + 1)] - 2 * u[indexingDiricheletx(i, j, k)] + u[indexingDiricheletx(i, j, k - 1)]) / (DZ * DZ);
    }

    return -1 * (u_value * diffx + v_value * diffy + w_value * diffz) + 1 / RE * (diffsecx + diffsecy + diffsecz) + functionG_u(i, j, k, t);
#endif
}

Real IcoNS::functionF_v(const std::array<Real, NX *(NY + 1) * (NZ + 1)> &u, const std::array<Real, (NX + 1) * NY *(NZ + 1)> &v, const std::array<Real, (NX + 1) * (NY + 1) * NZ> &w, int i, int j, int k, Real t)
{
#ifdef PERIODIC
    return -((u[indexingPeriodicx(i, j, k)] + u[indexingPeriodicx(i, j + 1, k)] + u[indexingPeriodicx(i - 1, j, k)] + u[indexingPeriodicx(i - 1, j + 1, k)]) / 4.0 * (v[indexingPeriodicy(i + 1, j, k)] - v[indexingPeriodicy(i - 1, j, k)]) / (2.0 * DX) +
             v[indexingPeriodicy(i, j, k)] * (v[indexingPeriodicy(i, j + 1, k)] - v[indexingPeriodicy(i, j - 1, k)]) / (2.0 * DY) +
             (w[indexingPeriodicz(i, j, k)] + w[indexingPeriodicz(i, j + 1, k)] + w[indexingPeriodicz(i, j, k - 1)] + w[indexingPeriodicz(i, j + 1, k - 1)]) / 4.0 * (v[indexingPeriodicy(i, j, k + 1)] - v[indexingPeriodicy(i, j, k - 1)]) / (2.0 * DZ)) +
           1 / RE * ((v[indexingPeriodicy(i + 1, j, k)] - 2 * v[indexingPeriodicy(i, j, k)] + v[indexingPeriodicy(i - 1, j, k)]) / (DX * DX) + (v[indexingPeriodicy(i, j + 1, k)] - 2 * v[indexingPeriodicy(i, j, k)] + v[indexingPeriodicy(i, j - 1, k)]) / (DY * DY) + (v[indexingPeriodicy(i, j, k + 1)] - 2 * v[indexingPeriodicy(i, j, k)] + v[indexingPeriodicy(i, j, k - 1)]) / (DZ * DZ)) +
           functionG_v(i, j, k, t);
#endif
#ifdef DIRICHELET
    /* return -((u[indexingDiricheletx(i, j, k)] + u[indexingDiricheletx(i, j + 1, k)] + u[indexingDiricheletx(i - 1, j, k)] + u[indexingDiricheletx(i - 1, j + 1, k)]) / 4.0 * (v[indexingDirichelety(i + 1, j, k)] - v[indexingDirichelety(i - 1, j, k)]) / (2.0 * DX) +
             v[indexingDirichelety(i, j, k)] * (v[indexingDirichelety(i, j + 1, k)] - v[indexingDirichelety(i, j - 1, k)]) / (2.0 * DY) +
             (w[indexingDiricheletz(i, j, k)] + w[indexingDiricheletz(i, j + 1, k)] + w[indexingDiricheletz(i, j, k - 1)] + w[indexingDiricheletz(i, j + 1, k - 1)]) / 4.0 * (v[indexingDirichelety(i, j, k + 1)] - v[indexingDirichelety(i, j, k - 1)]) / (2.0 * DZ)) +
           1 / RE * ((v[indexingDirichelety(i + 1, j, k)] - 2 * v[indexingDirichelety(i, j, k)] + v[indexingDirichelety(i - 1, j, k)]) / (DX * DX) + (v[indexingDirichelety(i, j + 1, k)] - 2 * v[indexingDirichelety(i, j, k)] + v[indexingDirichelety(i, j - 1, k)]) / (DY * DY) + (v[indexingDirichelety(i, j, k + 1)] - 2 * v[indexingDirichelety(i, j, k)] + v[indexingDirichelety(i, j, k - 1)]) / (DZ * DZ)) +
           functionG_v(i, j, k, t); */
    Real diffx, diffy, diffz, u_value, v_value, w_value, diffsecx, diffsecy, diffsecz;
    v_value = v[indexingDirichelety(i, j, k)];

    if (i == 0 || i == NX)
    {
        u_value = boundary.boundary_value_u[0]->value(i - 0.5, j, k, t);
    }
    else
    {
        u_value = (u[indexingDiricheletx(i, j, k)] + u[indexingDiricheletx(i, j + 1, k)] + u[indexingDiricheletx(i - 1, j, k)] + u[indexingDiricheletx(i - 1, j + 1, k)]) / 4.0;
    }

    if (k == 0 || k == NZ)
    {
        w_value = boundary.boundary_value_w[0]->value(i, j, k - 0.5, t);
    }
    else
    {
        w_value = (w[indexingDiricheletz(i, j, k)] + w[indexingDiricheletz(i, j + 1, k)] + w[indexingDiricheletz(i, j, k - 1)] + w[indexingDiricheletz(i, j + 1, k - 1)]) / 4.0;
    }

    if (i == 0)
    {
        diffx = (-3 * v[indexingDirichelety(i, j, k)] + 4 * v[indexingDirichelety(i + 1, j, k)] - v[indexingDirichelety(i + 2, j, k)]) / (2 * DX);
    }
    else if (i == NX)
    {
        diffx = (3 * v[indexingDirichelety(i, j, k)] - 4 * v[indexingDirichelety(i - 1, j, k)] + v[indexingDirichelety(i - 2, j, k)]) / (2 * DX);
    }
    else
    {
        diffx = (v[indexingDirichelety(i + 1, j, k)] - v[indexingDirichelety(i - 1, j, k)]) / (2.0 * DX);
    }

    if (j == 0)
    {
        diffy = (-4 * boundary.boundary_value_v[0]->value(i, j - 0.5, k, t) + 3 * v[indexingDirichelety(i, j, k)] + v[indexingDirichelety(i, j + 1, k)]) / (3.0 * DY);
    }
    else if (j == NY - 1)
    {
        diffy = (4 * boundary.boundary_value_v[0]->value(i, NY - 0.5, k, t) - 3 * v[indexingDirichelety(i, j, k)] - v[indexingDirichelety(i, j - 1, k)]) / (3.0 * DY);
    }
    else
    {
        diffy = (v[indexingDirichelety(i, j + 1, k)] - v[indexingDirichelety(i, j - 1, k)]) / (2.0 * DY);
    }

    if (k == 0)
    {
        diffz = (-3 * v[indexingDirichelety(i, j, k)] + 4 * v[indexingDirichelety(i, j, k + 1)] - v[indexingDirichelety(i, j, k + 2)]) / (2 * DZ);
    }
    else if (k == NZ)
    {
        diffz = (3 * v[indexingDirichelety(i, j, k)] - 4 * v[indexingDirichelety(i, j, k - 1)] + v[indexingDirichelety(i, j, k - 2)]) / (2 * DZ);
    }
    else
    {
        diffz = (v[indexingDirichelety(i, j, k + 1)] - v[indexingDirichelety(i, j, k - 1)]) / (2.0 * DZ);
    }

    if (j == 0)
    {
        diffsecy = (16 * boundary.boundary_value_v[0]->value(i, j - 0.5, k, t) - 25 * v[indexingDirichelety(i, j, k)] + 10 * v[indexingDirichelety(i, j + 1, k)] - v[indexingDirichelety(i, j + 2, k)]) / (5.0 * DY * DY);
    }
    else if (j == NY - 1)
    {
        diffsecy = (-16 * boundary.boundary_value_v[0]->value(i, NY - 0.5, k, t) + 25 * v[indexingDirichelety(i, j, k)] - 10 * v[indexingDirichelety(i, j - 1, k)] + v[indexingDirichelety(i, j - 2, k)]) / (5.0 * DY * DY);
    }
    else
    {
        diffsecy = (v[indexingDirichelety(i, j + 1, k)] - 2 * v[indexingDirichelety(i, j, k)] + v[indexingDirichelety(i, j - 1, k)]) / (DY * DY);
    }

    if (i == 0)
    {
        diffsecx = (2 * v[indexingDirichelety(i, j, k)] - 5 * v[indexingDirichelety(i + 1, j, k)] + 4 * v[indexingDirichelety(i + 2, j, k)] - v[indexingDirichelety(i + 3, j, k)]) / (DX * DX);
    }
    else if (i == NX)
    {
        diffsecx = (-2 * v[indexingDirichelety(i, j, k)] + 5 * v[indexingDirichelety(i - 1, j, k)] - 4 * v[indexingDirichelety(i - 2, j, k)] + v[indexingDirichelety(i - 3, j, k)]) / (DX * DX);
    }
    else
    {
        diffsecx = (v[indexingDirichelety(i + 1, j, k)] - 2 * v[indexingDirichelety(i, j, k)] + v[indexingDirichelety(i - 1, j, k)]) / (DX * DX);
    }

    if (k == 0)
    {
        diffsecz = (2 * v[indexingDirichelety(i, j, k)] - 5 * v[indexingDirichelety(i, j, k + 1)] + 4 * v[indexingDirichelety(i, j, k + 2)] - v[indexingDirichelety(i, j, k + 3)]) / (DZ * DZ);
    }
    else if (k == NZ)
    {
        diffsecz = (-2 * v[indexingDirichelety(i, j, k)] + 5 * v[indexingDirichelety(i, j, k - 1)] - 4 * v[indexingDirichelety(i, j, k - 2)] + v[indexingDirichelety(i, j, k - 3)]) / (DZ * DZ);
    }
    else
    {
        diffsecz = (v[indexingDirichelety(i, j, k + 1)] - 2 * v[indexingDirichelety(i, j, k)] + v[indexingDirichelety(i, j, k - 1)]) / (DZ * DZ);
    }

    return -1 * (u_value * diffx + v_value * diffy + w_value * diffz) + 1 / RE * (diffsecx + diffsecy + diffsecz) + functionG_v(i, j, k, t);
#endif
}

Real IcoNS::functionF_w(const std::array<Real, NX *(NY + 1) * (NZ + 1)> &u, const std::array<Real, (NX + 1) * NY *(NZ + 1)> &v, const std::array<Real, (NX + 1) * (NY + 1) * NZ> &w, int i, int j, int k, Real t)

{
#ifdef PERIODIC
    return -((u[indexingPeriodicx(i, j, k)] + u[indexingPeriodicx(i, j, k + 1)] + u[indexingPeriodicx(i - 1, j, k)] + u[indexingPeriodicx(i - 1, j, k + 1)]) / 4.0 * (w[indexingPeriodicz(i + 1, j, k)] - w[indexingPeriodicz(i - 1, j, k)]) / (2.0 * DX) +
             (v[indexingPeriodicy(i, j, k)] + v[indexingPeriodicy(i, j - 1, k)] + v[indexingPeriodicy(i, j, k + 1)] + v[indexingPeriodicy(i, j - 1, k + 1)]) / 4.0 * (w[indexingPeriodicz(i, j + 1, k)] - w[indexingPeriodicz(i, j - 1, k)]) / (2.0 * DY) +
             w[indexingPeriodicz(i, j, k)] * (w[indexingPeriodicz(i, j, k + 1)] - w[indexingPeriodicz(i, j, k - 1)]) / (2.0 * DZ)) +
           1 / RE * ((w[indexingPeriodicz(i + 1, j, k)] - 2 * w[indexingPeriodicz(i, j, k)] + w[indexingPeriodicz(i - 1, j, k)]) / (DX * DX) + (w[indexingPeriodicz(i, j + 1, k)] - 2 * w[indexingPeriodicz(i, j, k)] + w[indexingPeriodicz(i, j - 1, k)]) / (DY * DY) + (w[indexingPeriodicz(i, j, k + 1)] - 2 * w[indexingPeriodicz(i, j, k)] + w[indexingPeriodicz(i, j, k - 1)]) / (DZ * DZ)) +
           functionG_w(i, j, k, t);
#endif
#ifdef DIRICHELET
    /* return -((u[indexingDiricheletx(i, j, k)] + u[indexingDiricheletx(i, j, k + 1)] + u[indexingDiricheletx(i - 1, j, k)] + u[indexingDiricheletx(i - 1, j, k + 1)]) / 4.0 * (w[indexingDiricheletz(i + 1, j, k)] - w[indexingDiricheletz(i - 1, j, k)]) / (2.0 * DX) +
             (v[indexingDirichelety(i, j, k)] + v[indexingDirichelety(i, j - 1, k)] + v[indexingDirichelety(i, j, k + 1)] + v[indexingDirichelety(i, j - 1, k + 1)]) / 4.0 * (w[indexingDiricheletz(i, j + 1, k)] - w[indexingDiricheletz(i, j - 1, k)]) / (2.0 * DY) +
             w[indexingDiricheletz(i, j, k)] * (w[indexingDiricheletz(i, j, k + 1)] - w[indexingDiricheletz(i, j, k - 1)]) / (2.0 * DZ)) +
           1 / RE * ((w[indexingDiricheletz(i + 1, j, k)] - 2 * w[indexingDiricheletz(i, j, k)] + w[indexingDiricheletz(i - 1, j, k)]) / (DX * DX) + (w[indexingDiricheletz(i, j + 1, k)] - 2 * w[indexingDiricheletz(i, j, k)] + w[indexingDiricheletz(i, j - 1, k)]) / (DY * DY) + (w[indexingDiricheletz(i, j, k + 1)] - 2 * w[indexingDiricheletz(i, j, k)] + w[indexingDiricheletz(i, j, k - 1)]) / (DZ * DZ)) +
           functionG_w(i, j, k, t); */
    Real diffx, diffy, diffz, u_value, v_value, w_value, diffsecx, diffsecy, diffsecz;
    w_value = w[indexingDiricheletz(i, j, k)];
    if (i == 0 || i == NX)
    {
        u_value = boundary.boundary_value_u[0]->value(i - 0.5, j, k, t);
    }
    else
    {
        u_value = (u[indexingDiricheletx(i, j, k)] + u[indexingDiricheletx(i, j, k + 1)] + u[indexingDiricheletx(i - 1, j, k)] + u[indexingDiricheletx(i - 1, j, k + 1)]) / 4.0;
    }
    if (j == 0 || j == NY)
    {
        v_value = boundary.boundary_value_v[0]->value(i, j - 0.5, k, t);
    }
    else
    {
        v_value = (v[indexingDirichelety(i, j, k)] + v[indexingDirichelety(i, j - 1, k)] + v[indexingDirichelety(i, j, k + 1)] + v[indexingDirichelety(i, j - 1, k + 1)]) / 4.0;
    }

    if (i == 0)
    {
        diffx = (-3 * w[indexingDiricheletz(i, j, k)] + 4 * w[indexingDiricheletz(i + 1, j, k)] - w[indexingDiricheletz(i + 2, j, k)]) / (2 * DX);
    }
    else if (i == NX)
    {
        diffx = (3 * w[indexingDiricheletz(i, j, k)] - 4 * w[indexingDiricheletz(i - 1, j, k)] + w[indexingDiricheletz(i - 2, j, k)]) / (2 * DX);
    }
    else
    {
        diffx = (w[indexingDiricheletz(i + 1, j, k)] - w[indexingDiricheletz(i - 1, j, k)]) / (2.0 * DX);
    }

    if (j == 0)
    {
        diffy = (-3 * w[indexingDiricheletz(i, j, k)] + 4 * w[indexingDiricheletz(i, j + 1, k)] - w[indexingDiricheletz(i, j + 2, k)]) / (2 * DY);
    }
    else if (j == NY)
    {
        diffy = (3 * w[indexingDiricheletz(i, j, k)] - 4 * w[indexingDiricheletz(i, j - 1, k)] + w[indexingDiricheletz(i, j - 2, k)]) / (2 * DY);
    }
    else
    {
        diffy = (w[indexingDiricheletz(i, j + 1, k)] - w[indexingDiricheletz(i, j - 1, k)]) / (2.0 * DY);
    }

    if (k == 0)
    {
        diffz = (-4 * boundary.boundary_value_w[0]->value(i, j, k - 0.5, t) + 3 * w[indexingDiricheletz(i, j, k)] + w[indexingDiricheletz(i, j, k + 1)]) / (3.0 * DZ);
    }
    else if (k == NZ - 1)
    {
        diffz = (4 * boundary.boundary_value_w[0]->value(i, j, NZ - 0.5, t) - 3 * w[indexingDiricheletz(i, j, k)] - w[indexingDiricheletz(i, j, k - 1)]) / (3.0 * DZ);
    }
    else
    {
        diffz = (w[indexingDiricheletz(i, j, k + 1)] - w[indexingDiricheletz(i, j, k - 1)]) / (2.0 * DZ);
    }

    if (k == 0)
    {
        diffsecz = (16 * boundary.boundary_value_w[0]->value(i, j, k - 0.5, t) - 25 * w[indexingDiricheletz(i, j, k)] + 10 * w[indexingDiricheletz(i, j, k + 1)] - w[indexingDiricheletz(i, j, k + 2)]) / (5.0 * DZ * DZ);
    }
    else if (k == NZ - 1)
    {
        diffsecz = (-16 * boundary.boundary_value_w[0]->value(i, j, NZ - 0.5, t) + 25 * w[indexingDiricheletz(i, j, k)] - 10 * w[indexingDiricheletz(i, j, k - 1)] + w[indexingDiricheletz(i, j, k - 2)]) / (5.0 * DZ * DZ);
    }
    else
    {
        diffsecz = (w[indexingDiricheletz(i, j, k + 1)] - 2 * w[indexingDiricheletz(i, j, k)] + w[indexingDiricheletz(i, j, k - 1)]) / (DZ * DZ);
    }

    if (i == 0)
    {
        diffsecx = (2 * w[indexingDiricheletz(i, j, k)] - 5 * w[indexingDiricheletz(i + 1, j, k)] + 4 * w[indexingDiricheletz(i + 2, j, k)] - w[indexingDiricheletz(i + 3, j, k)]) / (DX * DX);
    }
    else if (i == NX)
    {
        diffsecx = (-2 * w[indexingDiricheletz(i, j, k)] + 5 * w[indexingDiricheletz(i - 1, j, k)] - 4 * w[indexingDiricheletz(i - 2, j, k)] + w[indexingDiricheletz(i - 3, j, k)]) / (DX * DX);
    }
    else
    {
        diffsecx = (w[indexingDiricheletz(i + 1, j, k)] - 2 * w[indexingDiricheletz(i, j, k)] + w[indexingDiricheletz(i - 1, j, k)]) / (DX * DX);
    }

    if (j == 0)
    {
        diffsecy = (2 * w[indexingDiricheletz(i, j, k)] - 5 * w[indexingDiricheletz(i, j + 1, k)] + 4 * w[indexingDiricheletz(i, j + 2, k)] - w[indexingDiricheletz(i, j + 3, k)]) / (DY * DY);
    }
    else if (j == NY)
    {
        diffsecy = (-2 * w[indexingDiricheletz(i, j, k)] + 5 * w[indexingDiricheletz(i, j - 1, k)] - 4 * w[indexingDiricheletz(i, j - 2, k)] + w[indexingDiricheletz(i, j - 3, k)]) / (DY * DY);
    }
    else
    {
        diffsecy = (w[indexingDiricheletz(i, j + 1, k)] - 2 * w[indexingDiricheletz(i, j, k)] + w[indexingDiricheletz(i, j - 1, k)]) / (DY * DY);
    }

    return -1 * (u_value * diffx + v_value * diffy + w_value * diffz) + 1 / RE * (diffsecx + diffsecy + diffsecz) + functionG_w(i, j, k, t);
#endif
}

// w[((i + 1) % (NX) + NX) * (NY + 1) * NZ + (j) % (NY)*NZ + (k) % (NZ)] - 2.0 * w[lw] + w[((i - 1) % (NX) + NX) * (NY + 1) * NZ + (j) % (NY)*NZ + (k) % (NZ)]) / (DX * DX);
// v[((i + 1) % (NX) + NX) * NY * (NZ + 1) + (j) % (NY) * (NZ + 1) + (k) % (NZ)] - 2.0 * v[lv] + v[((i - 1) % (NX) + NX) * NY * (NZ + 1) + (j) % (NY) * (NZ + 1) + (k) % (NZ)]) / (DX * DX);

/* Real IcoNS::functionF_u(const std::array<Real, NX *(NY + 1) * (NZ + 1)> &u, const std::array<Real, (NX + 1) * NY *(NZ + 1)> &v, const std::array<Real, (NX + 1) * (NY + 1) * NZ> &w, int i, int j, int k, Real t)
{
    int lu = i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k;
    int lv = i * NY * (NZ + 1) + j * (NZ + 1) + k;
    int lw = i * (NY + 1) * NZ + j * NZ + k;

    return -(u[lu] * (u[lu + (NY + 1) * (NZ + 1)] - u[lu - (NY + 1) * (NZ + 1)]) / (2.0 * DX) +
             (v[lv] + v[lv + NY * (NZ + 1)] + v[lv - (NZ + 1)] + v[lv + NY * (NZ + 1) - (NZ + 1)]) / 4.0 * (u[lu + (NZ + 1)] - u[lu - (NZ + 1)]) / (2.0 * DY) +
             (w[lw] + w[lw + (NY + 1) * NZ] + w[lw - 1] + w[lw + (NY + 1) * NZ - 1]) / 4.0 * (u[lu + 1] - u[lu - 1]) / (2.0 * DZ)) +
           1 / RE * ((u[lu + (NY + 1) * (NZ + 1)] - 2 * u[lu] + u[lu - (NY + 1) * (NZ + 1)]) / (DX * DX) + (u[lu + (NZ + 1)] - 2 * u[lu] + u[lu - (NZ + 1)]) / (DY * DY) + (u[lu + 1] - 2 * u[lu] + u[lu - 1]) / (DZ * DZ)) +
             functionG_u(i, j, k, t);
}

Real IcoNS::functionF_v(const std::array<Real, NX *(NY + 1) * (NZ + 1)> &u, const std::array<Real, (NX + 1) * NY *(NZ + 1)> &v, const std::array<Real, (NX + 1) * (NY + 1) * NZ> &w, int i, int j, int k, Real t)
{
    int lu = i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k;
    int lv = i * NY * (NZ + 1) + j * (NZ + 1) + k;
    int lw = i * (NY + 1) * NZ + j * NZ + k;

    return -((u[lu] + u[lu + (NZ + 1)] + u[lu - (NY + 1) * (NZ + 1)] + u[lu - (NY + 1) * (NZ + 1) + (NZ + 1)]) / 4.0 * ((v[lv + NY * (NZ + 1)] - v[lv - NY * (NZ + 1)]) / (2.0 * DX)) +
             v[lv] * (v[lv + (NZ + 1)] - v[lv - (NZ + 1)]) / (2.0 * DY) +
             (w[lw] + w[lw - 1] + w[lw + (NY + 1)] + w[lw + (NY + 1) - 1]) / 4.0 * (v[lv + 1] - v[lv - 1]) / (2.0 * DZ)) +
           (1.0 / RE) * ((v[lv + NY * (NZ + 1)] - 2.0 * v[lv] + v[lv - NY * (NZ + 1)]) / (DX * DX) +
                         (v[lv + (NZ + 1)] - 2.0 * v[lv] + v[lv - (NZ + 1)]) / (DY * DY) +
                         (v[lv + 1] - 2.0 * v[lv] + v[lv - 1]) / (DZ * DZ)) +
           functionG_v(i, j, k, t);
}

Real IcoNS::functionF_w(const std::array<Real, NX *(NY + 1) * (NZ + 1)> &u, const std::array<Real, (NX + 1) * NY *(NZ + 1)> &v, const std::array<Real, (NX + 1) * (NY + 1) * NZ> &w, int i, int j, int k, Real t)

{
    int lu = i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k;
    int lv = i * NY * (NZ + 1) + j * (NZ + 1) + k;
    int lw = i * (NY + 1) * NZ + j * NZ + k;

    return -((u[lu] + u[lu - (NY + 1) * (NZ + 1)] + u[lu + 1] + u[lu - (NZ + 1) * (NY + 1) + 1]) / 4.0 * (w[lw + (NY + 1) * NZ] - w[lw - (NY + 1) * NZ]) / (2.0 * DX) +
             (v[lv + 1] + v[lv - (NZ + 1) + 1] + v[lv] + v[lv - (NZ + 1)]) / 4.0 * (w[lw + NZ] - w[lw - NZ]) / (2.0 * DY) +
             w[lw] * (w[lw + 1] - w[lw - 1]) / (2.0 * DZ)) +
           (1.0 / RE) * ((w[lw + (NY + 1) * NZ] - 2.0 * w[lw] + w[lw - (NY + 1) * NZ]) / (DX * DX) +
                         (w[lw + NZ] - 2.0 * w[lw] + w[lw - NZ]) / (DY * DY) +
                         (w[lw + 1] - 2.0 * w[lw] + w[lw - 1]) / (DZ * DZ)) +
           functionG_w(i, j, k, t);
} */

Real IcoNS::functionG_u(int i, int j, int k, Real t)
{
    Real x = i * DX + DX / 2;
    Real y = j * DY;
    Real z = k * DZ;
    return 0.0;
    return std::sin(x) * std::cos(y) * std::sin(z) * std::cos(t) +
           std::sin(x) * std::cos(x) * std::cos(y) * std::cos(y) * std::sin(z) * std::sin(z) * std::sin(t) * std::sin(t) -
           std::sin(x) * std::cos(x) * std::sin(y) * std::sin(y) * std::sin(z) * std::sin(z) * std::sin(t) * std::sin(t) + 2.0 * std::sin(x) * std::cos(x) * std::cos(y) * std::cos(y) * std::cos(z) * std::cos(z) * std::sin(t) * std::sin(t) +
           3.0 / RE * std::sin(x) * std::cos(y) * std::sin(z) * std::sin(t) - std::sin(x) * std::cos(y) * std::cos(z) * std::sin(t);
    /*     return std::sin(z)*std::cos(t) + std::sin(y)*std::sin(t)*std::cos(z)*std::sin(t)+1.0 / RE *std::sin(z)*std::sin(t);
     */
    // return std::cos(t);
}

Real IcoNS::functionG_v(int i, int j, int k, Real t)
{
    Real x = i * DX;
    Real y = j * DY + DY / 2;
    Real z = k * DZ;
    return 0.0;
    return std::cos(x) * std::sin(y) * std::sin(z) * std::cos(t) -
           std::sin(x) * std::sin(x) * std::sin(y) * std::cos(y) * std::sin(z) * std::sin(z) * std::sin(t) * std::sin(t) +
           std::cos(x) * std::cos(x) * std::sin(y) * std::cos(y) * std::sin(z) * std::sin(z) * std::sin(t) * std::sin(t) +
           2.0 * std::cos(x) * std::cos(x) * std::sin(y) * std::cos(y) * std::cos(z) * std::cos(z) * std::sin(t) * std::sin(t) +
           3.0 / RE * std::cos(x) * std::sin(y) * std::sin(z) * std::sin(t) - std::cos(x) * std::sin(y) * std::cos(z) * std::sin(t);
    /*     return std::sin(x)*std::cos(t) + std::sin(z)*std::sin(t)*std::cos(x)*std::sin(t)+1.0 / RE *std::sin(x)*std::sin(t);
     */
    // return std::cos(t);
}

Real IcoNS::functionG_w(int i, int j, int k, Real t)
{
    Real x = i * DX;
    Real y = j * DY;
    Real z = k * DZ + DZ / 2;
    return 0.0;
    return 2.0 * std::cos(x) * std::cos(y) * std::cos(z) * std::cos(t) -
           2.0 * std::sin(x) * std::sin(x) * std::cos(y) * std::cos(y) * std::sin(z) * std::cos(z) * std::sin(t) * std::sin(t) -
           2.0 * std::cos(x) * std::cos(x) * std::sin(y) * std::sin(y) * std::sin(z) * std::cos(z) * std::sin(t) * std::sin(t) -
           4.0 * std::cos(x) * std::cos(x) * std::cos(y) * std::cos(y) * std::sin(z) * std::cos(z) * std::sin(t) * std::sin(t) +
           6.0 / RE * std::cos(x) * std::cos(y) * std::cos(z) * std::sin(t) - std::cos(x) * std::cos(y) * std::sin(z) * std::sin(t);
    /*     return std::sin(y)*std::cos(t) + std::sin(x)*std::sin(t)*std::cos(y)*std::sin(t)+1.0 / RE *std::sin(y)*std::sin(t);
     */
    // return std::cos(t);
}

Real IcoNS::d_Px(int i, int j, int k, Real t)
{
    return -std::sin((i + 0.5) * DX) * std::cos(j * DY) * std::sin(k * DZ) * std::sin(t);
}

Real IcoNS::d_Py(int i, int j, int k, Real t)
{
    return -std::cos(i * DX) * std::sin((j + 0.5) * DY) * std::sin(k * DZ) * std::sin(t);
}

Real IcoNS::d_Pz(int i, int j, int k, Real t)
{
    return std::cos(i * DX) * std::cos(j * DY) * std::cos((k + 0.5) * DZ) * std::sin(t);
}

#ifdef DIRICHELET
void IcoNS::pressionCorrection(std::array<Real, (NX + 1) * (NY + 1) * (NZ + 1)> &p)
{
    // LEFT FACE
    for (int j = 1; j < NY; j++)
    {
        for (int k = 1; k < NZ; k++)
        {
            p[indexingDiricheletp(0, j, k)] = (4 / 3) * p[indexingDiricheletp(1, j, k)] - (1 / 3) * p[indexingDiricheletp(2, j, k)];
        }
    }
    // RIGHT FACE
    for (int j = 1; j < NY; j++)
    {
        for (int k = 1; k < NZ; k++)
        {
            p[indexingDiricheletp(NX, j, k)] = (4 / 3) * p[indexingDiricheletp(NX - 1, j, k)] - (1 / 3) * p[indexingDiricheletp(NX - 2, j, k)];
        }
    }
    // FRONT FACE
    for (int i = 1; i < NX; i++)
    {
        for (int k = 1; k < NZ; k++)
        {
            p[indexingDiricheletp(i, 0, k)] = (4 / 3) * p[indexingDiricheletp(i, 1, k)] - (1 / 3) * p[indexingDiricheletp(i, 2, k)];
        }
    }
    // BACK FACE
    for (int i = 1; i < NX; i++)
    {
        for (int k = 1; k < NZ; k++)
        {
            p[indexingDiricheletp(i, NY, k)] = (4 / 3) * p[indexingDiricheletp(i, NY - 1, k)] - (1 / 3) * p[indexingDiricheletp(i, NY - 2, k)];
        }
    }
    // LOWER FACE
    for (int i = 1; i < NX; i++)
    {
        for (int j = 1; j < NY; j++)
        {
            p[indexingDiricheletp(i, j, 0)] = (4 / 3) * p[indexingDiricheletp(i, j, 1)] - (1 / 3) * p[indexingDiricheletp(i, j, 2)];
        }
    }
    // UPPER FACE
    for (int i = 1; i < NX; i++)
    {
        for (int j = 1; j < NY; j++)
        {
            p[indexingDiricheletp(i, j, NZ)] = (4 / 3) * p[indexingDiricheletp(i, j, NZ - 1)] - (1 / 3) * p[indexingDiricheletp(i, j, NZ - 2)];
        }
    }
    // 4 X EDGES
    for (int i = 1; i < NX; i++)
    {
        // LOWER FRONT EDGE
        p[indexingDiricheletp(i, 0, 0)] = (4 / 3) * p[indexingDiricheletp(i, 1, 0)] - (1 / 3) * p[indexingDiricheletp(i, 2, 0)];
        // UPPER FRONT EDGE
        p[indexingDiricheletp(i, 0, NZ)] = (4 / 3) * p[indexingDiricheletp(i, 1, NZ)] - (1 / 3) * p[indexingDiricheletp(i, 2, NZ)];
        // LOWER BACK EDGE
        p[indexingDiricheletp(i, NY, 0)] = (4 / 3) * p[indexingDiricheletp(i, NY - 1, 0)] - (1 / 3) * p[indexingDiricheletp(i, NY - 2, 0)];
        // UPPER BACK EDGE
        p[indexingDiricheletp(i, NY, NZ)] = (4 / 3) * p[indexingDiricheletp(i, NY - 1, NZ)] - (1 / 3) * p[indexingDiricheletp(i, NY - 2, NZ)];
    }
    // 4 Y EDGES
    for (int j = 1; j < NY; j++)
    {
        // LOWER LEFT EDGE
        p[indexingDiricheletp(0, j, 0)] = (4 / 3) * p[indexingDiricheletp(1, j, 0)] - (1 / 3) * p[indexingDiricheletp(2, j, 0)];
        // UPPER LEFT EDGE
        p[indexingDiricheletp(0, j, NZ)] = (4 / 3) * p[indexingDiricheletp(1, j, NZ)] - (1 / 3) * p[indexingDiricheletp(2, j, NZ)];
    } // break to exploit locality
    for (int j = 1; j < NY; j++)
    {
        // LOWER RIGHT EDGE
        p[indexingDiricheletp(NX, j, 0)] = (4 / 3) * p[indexingDiricheletp(NX - 1, j, 0)] - (1 / 3) * p[indexingDiricheletp(NX - 2, j, 0)];
        // UPPER RIGHT EDGE
        p[indexingDiricheletp(NX, j, NZ)] = (4 / 3) * p[indexingDiricheletp(NX - 1, j, NZ)] - (1 / 3) * p[indexingDiricheletp(NX - 2, j, NZ)];
    }
    // 4 Z EDGES
    for (int k = 1; k < NZ; k++)
    {
        // FRONT LEFT EDGE
        p[indexingDiricheletp(0, 0, k)] = (4 / 3) * p[indexingDiricheletp(1, 0, k)] - (1 / 3) * p[indexingDiricheletp(2, 0, k)];
    } // break to exploit locality
    for (int k = 1; k < NZ; k++)
    {
        // BACK LEFT EDGE
        p[indexingDiricheletp(0, NY, k)] = (4 / 3) * p[indexingDiricheletp(1, NY, k)] - (1 / 3) * p[indexingDiricheletp(2, NY, k)];
    } // break to exploit locality
    for (int k = 1; k < NZ; k++)
    {
        // FRONT RIGHT EDGE
        p[indexingDiricheletp(NX, 0, k)] = (4 / 3) * p[indexingDiricheletp(NX - 1, 0, k)] - (1 / 3) * p[indexingDiricheletp(NX - 2, 0, k)];
    } // break to exploit locality
    for (int k = 1; k < NZ; k++)
    {
        // BACK RIGHT EDGE
        p[indexingDiricheletp(NX, NY, k)] = (4 / 3) * p[indexingDiricheletp(NX - 1, NY, k)] - (1 / 3) * p[indexingDiricheletp(NX - 2, NY, k)];
    }
    // 8 VERTICES
    // LOWER FRONT LEFT VERTEX (0,0,0)
    p[0] = (4 / 3) * p[1] - (1 / 3) * p[2];
    // LOWER FRONT RIGHT VERTEX (1,0,0)
    p[indexingDiricheletp(NX, 0, 0)] = (4 / 3) * p[indexingDiricheletp(NX - 1, 0, 0)] - (1 / 3) * p[indexingDiricheletp(NX - 2, 0, 0)];
    // LOWER BACK LEFT VERTEX (0,1,0)
    p[indexingDiricheletp(0, NY, 0)] = (4 / 3) * p[indexingDiricheletp(1, NY, 0)] - (1 / 3) * p[indexingDiricheletp(2, NY, 0)];
    // UPPER FRONT LEFT VERTEX (0,0,1)
    p[indexingDiricheletp(0, 0, NZ)] = (4 / 3) * p[indexingDiricheletp(1, 0, NZ)] - (1 / 3) * p[indexingDiricheletp(2, 0, NZ)];
    // UPPER BACK LEFT VERTEX (0,1,1)
    p[indexingDiricheletp(0, NY, NZ)] = (4 / 3) * p[indexingDiricheletp(1, NY, NZ)] - (1 / 3) * p[indexingDiricheletp(2, NY, NZ)];
    // UPPER FRONT RIGHT VERTEX (1,0,1)
    p[indexingDiricheletp(NX, 0, NZ)] = (4 / 3) * p[indexingDiricheletp(NX - 1, 0, NZ)] - (1 / 3) * p[indexingDiricheletp(NX - 2, 0, NZ)];
    // LOWER BACK RIGHT VERTEX (1,1,0)
    p[indexingDiricheletp(NX, NY, 0)] = (4 / 3) * p[indexingDiricheletp(NX - 1, NY, 0)] - (1 / 3) * p[indexingDiricheletp(NX - 2, NY, 0)];
    // UPPER BACK RIGHT VERTEX (1,1,1)
    p[indexingDiricheletp(NX, NY, NZ)] = (4 / 3) * p[indexingDiricheletp(NX - 1, NY, NZ)] - (1 / 3) * p[indexingDiricheletp(NX - 2, NY, NZ)];

    // compute the average pressure and subtract it from the pressure field
    Real p_avg = 0.0;
    for (int i = 0; i < NX + 1; i++)
    {
        for (int j = 0; j < NY + 1; j++)
        {
            for (int k = 0; k < NZ + 1; k++)
            {
                p_avg += p[indexingDiricheletp(i, j, k)];
            }
        }
    }
    p_avg /= (NX + 1) * (NY + 1) * (NZ + 1);
    for (int i = 0; i < NX + 1; i++)
    {
        for (int j = 0; j < NY + 1; j++)
        {
            for (int k = 0; k < NZ + 1; k++)
            {
                p[indexingDiricheletp(i, j, k)] -= p_avg;
            }
        }
    }
}
#endif