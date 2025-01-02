#include "core.hpp"
#include "poissonSolver.hpp"


void IcoNS::solve_time_step(Real time)
{
    //TODO: there was a part with start and end different for periodic and diricazz, i didn't bring it here since lbx variables should act the same

    //TODO: check, adapt the trues with BX BY AND BZ
    PoissonSolver poissonSolver(true,true,true);


    for (int i = 1 + lbx; i < newDimX_x - 1 - rbx; i++)
    {
        for (int j = 1 + lby; j < newDimY_x - 1 - rby; j++)
        {
            for (int k = 1; k < dim_z - 1; k++)

            {
                //TODO: 2deco per la pressione
                Y2_x[indexingDiricheletx(i, j, k)] = grid.u[indexingDiricheletx(i, j, k)] + 64.0 / 120.0 * DT * functionF_u(grid.u, grid.v, grid.w, i, j, k, time) -
                                                                   64.0 / 120.0 * DT * /*d_Px(i,j,k,time);*/(grid.p[(i+1) * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k] - grid.p[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k]) / (DX);
            }
        }
    }


    for (int i = 1 + lbx; i < newDimX_y - 1 - rbx; i++)
    {
        for (int j = 1 + lby; j < newDimY_y - 1 - rby; j++)
        {
            for (int k = 1; k < dim_z - 1; k++)

            {
                Y2_y[indexingDirichelety(i, j, k)] = grid.v[indexingDirichelety(i, j, k)] + 64.0 / 120.0 * DT * functionF_v(grid.u, grid.v, grid.w, i, j, k, time) -
                                                            64.0 / 120.0 * DT * /*d_Py(i,j,k,time);*/(grid.p[i * (NY + 1) * (NZ + 1) + (j+1) * (NZ + 1) + k] - grid.p[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k]) / (DY);
            }
        }
    }


    for (int i = 1 + lbx; i < newDimX_z - 1 - rbx; i++)
    {
        for (int j = 1 + lby; j < newDimY_z - 1 - rby; j++)
        {
            for (int k = 1; k < dim_z_z - 1; k++)
            {
                Y2_z[indexingDiricheletz(i, j, k)] = grid.w[indexingDiricheletz(i, j, k)] + 64.0 / 120.0 * DT * functionF_w(grid.u, grid.v, grid.w, i, j, k, time) -
                                                       64.0 / 120.0 * DT * /*d_Pz(i,j,k,time);*/(grid.p[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k+1] - grid.p[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k]) / (DZ);
            }
        }
    }


//TODO: 2deco parallelization
// dimensions accessible with c2d-> zSize[0] zSize[1] zSize[2] etc.
#ifdef PERIODIC
    for (int i = 0; i < zSize[0]; i++)
    {
        for (int j = 0; j < zSize[1]; j++)
        {
            for (int k = 0; k < zSize[2]; k++)
            {
                Y2_p[i * zSize[1] * zSize[2] + j * zsize[2] + k] = 120.0 / (64.0 * DT) * ((Y2_x[indexingPeriodicx(i, j, k)] - Y2_x[indexingPeriodicx(i - 1, j, k)]) / (DX) + (Y2_y[indexingPeriodicy(i, j, k)] - Y2_y[indexingPeriodicy(i, j - 1, k)]) / (DY) + (Y2_z[indexingPeriodicz(i, j, k)] - Y2_z[indexingPeriodicz(i, j, k - 1)]) / (DZ));
            }
        }
    }
    /////////////////////poisson_solver.solvePoisson(step 1);//////////////////////////////////// -> the solution is stored in Sol_p=Y2_p
    poissonSolver.solveDirichletPoisson(Y2_p,helper);
#endif

#ifdef DIRICHELET
    //TODO: da vedere! exchange data?
    boundary.update_boundary(Y2_x, Y2_y, Y2_z, time + 64.0 / 120.0 * DT);
    
    for (int i = 1; i < NX; i++)
    {
        for (int j = 1; j < NY; j++)
        {
            for (int k = 1; k < NZ; k++)
            {
                Y2_p[i * (NY+1) * (NZ+1) + j * (NZ+1) + k] = 120.0 / (64.0 * DT) * ((Y2_x[indexingDiricheletx(i, j, k)] - Y2_x[indexingDiricheletx(i - 1, j, k)]) / (DX) + (Y2_y[indexingDirichelety(i, j, k)] - Y2_y[indexingDirichelety(i, j - 1, k)]) / (DY) + (Y2_z[indexingDiricheletz(i, j, k)] - Y2_z[indexingDiricheletz(i, j, k - 1)]) / (DZ));
            }
        }
    }
    //TODO: paura, exchange data?, check
    boundary.divergence(Y2_x, Y2_y, Y2_z, Y2_p, time + 64.0 / 120.0 * DT, 64.0);

    poissonSolver.solveNeumannPoisson(Y2_p);
#endif
    //TODO: adapt to the code now
    boundary.update_boundary(Y2_x, Y2_y, Y2_z, time + 64.0 / 120.0 * DT);
    MPI_Barrier(cart_comm);
    exchangeData(Y2_x, newDimX_x, newDimY_x, dim_z, MPI_face_x_x, MPI_face_y_x,0,1);
    exchangeData(Y2_y, newDimX_y, newDimY_y, dim_z, MPI_face_x_y, MPI_face_y_y,1,0);
    exchangeData(Y2_z, newDimX_z, newDimY_z, dim_z_z, MPI_face_x_z, MPI_face_y_z,1,1);


//TODO: correction of Y with pressure, check indexes should not be a problem, ifdef probably removable
 for (int i = 1; i < newDimX_x - 1; i++)
    {
        for (int j = 1; j < newDimY_x - 1; j++)
        {
            for (int k = 0; k < dim_z; k++)
            {
                Y2_x[indexingDiricheletx(i, j, k)] = Y2_x[indexingDiricheletx(i, j, k)] - 64.0 * DT / (120.0) * /*(d_Px(i,j,k,time+64.0 * DT / (120.0))-d_Px(i,j,k,time));*/(Y2_p[indexingDiricheletp(i + 1, j, k)] - Y2_p[indexingDiricheletp(i, j, k)]) / (DX);
            }
        }
    }

    for (int i = 1; i < newDimX_y - 1; i++)
    {
        for (int j = 1; j < newDimY_y - 1; j++)
        {
            for (int k = 0; k < dim_z; k++)
            {
                Y2_y[indexingDirichelety(i, j, k)] = Y2_y[indexingDirichelety(i, j, k)] - 64.0 * DT / (120.0) * /*(d_Py(i,j,k,time+64.0 * DT / (120.0))-d_Py(i,j,k,time));*/(Y2_p[indexingDiricheletp(i, j + 1, k)] - Y2_p[indexingDiricheletp(i, j, k)]) / (DY);
            }
        }
    }

    for (int i = 1; i < newDimX_z - 1; i++)
    {
        for (int j = 1; j < newDimY_z - 1; j++)
        {
            for (int k = 0; k < dim_z_z; k++)
            {
                Y2_z[indexingDiricheletz(i, j, k)] = Y2_z[indexingDiricheletz(i, j, k)] - 64.0 * DT / (120.0) * /*(d_Pz(i,j,k,time+64.0 * DT / (120.0))-d_Pz(i,j,k,time));*/(Y2_p[indexingDiricheletp(i, j, k + 1)] - Y2_p[indexingDiricheletp(i, j, k)]) / (DZ);
            }
        }
    }

    //TODO:?
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


    //TODO: should be same concept as y2s
    for (int i = 1 + lbx; i < newDimX_x - 1 - rbx; i++)
    {
        for (int j = 1 + lby; j < newDimY_x - 1 - rby; j++)
        {
            for (int k = 1; k < dim_z - 1; k++)

            {
                Y3_x[indexingDiricheletx(i, j, k)] = Y2_x[indexingDiricheletx(i, j, k)] +
                                                     50.0 / 120.0 * DT * functionF_u(Y2_x, Y2_y, Y2_z, i, j, k, time + 64.0 / 120.0 * DT) -
                                                     34.0 / 120.0 * DT * functionF_u(grid.u, grid.v, grid.w, i, j, k, time) -
                                                                   16.0 / 120.0 * DT * /*d_Px(i,j,k,time + 64.0/120.0*DT);*/(Phi_p[(i+1) * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k] - Phi_p[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k]) / (DX);
            }
        }
    }
    //TODO:: still same
    for (int i = 1 + lbx; i < newDimX_y - 1 - rbx; i++)
    {
        for (int j = 1 + lby; j < newDimY_y - 1 - rby; j++)
        {
            for (int k = 1; k < dim_z-1; k++)
            {
                Y3_y[indexingDirichelety(i, j, k)] = Y2_y[indexingDirichelety(i, j, k)] +
                                                     50.0 / 120.0 * DT * functionF_v(Y2_x, Y2_y, Y2_z, i, j, k, time + 64.0 / 120.0 * DT) -
                                                     34.0 / 120.0 * DT * functionF_v(grid.u, grid.v, grid.w, i, j, k, time) -
                                                             16.0 / 120.0 * DT * /*d_Py(i,j,k,time + 64.0/120.0*DT);*/(Phi_p[i * (NY + 1) * (NZ + 1) + (j+1) * (NZ + 1) + k] - Phi_p[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k]) / (DY);
            }
        }
    }

//TODO: still same
    for (int i = 1 + lbx; i < newDimX_z - 1 - rbx; i++)
    {
        for (int j = 1 + lby; j < newDimY_z - 1 - rby; j++)
        {
            for (int k = 1; k < dim_z_z - 1; k++)

            {
                Y3_z[indexingDiricheletz(i, j, k)] = Y2_z[indexingDiricheletz(i, j, k)] +
                                                     50.0 / 120.0 * DT * functionF_w(Y2_x, Y2_y, Y2_z, i, j, k, time + 64.0 / 120.0 * DT) -
                                                     34.0 / 120.0 * DT * functionF_w(grid.u, grid.v, grid.w, i, j, k, time) -
                                                       16.0 / 120.0 * DT * /*d_Pz(i,j,k,time + 64.0/120.0*DT);*/(Phi_p[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k+1] - Phi_p[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k]) / (DZ);
            }
        }
    }

//TODO: same pressure part as y2
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
    poissonSolver.solveDirichletPoisson(Y2_p,helper);
#endif
#ifdef DIRICHELET
    boundary.update_boundary(Y3_x, Y3_y, Y3_z, time + 80.0 / 120.0 * DT);

    for (int i = 1; i < NX; i++)
    {
        for (int j = 1; j < NY; j++)
        {
            for (int k = 1; k < NZ; k++)
            {
                Y2_p[i * (NY+1) * (NZ+1) + j * (NZ+1) + k] = 120.0 / (16.0 * DT) * ((Y3_x[indexingDiricheletx(i, j, k)] - Y3_x[indexingDiricheletx(i - 1, j, k)]) / (DX) + (Y3_y[indexingDirichelety(i, j, k)] - Y3_y[indexingDirichelety(i, j - 1, k)]) / (DY) + (Y3_z[indexingDiricheletz(i, j, k)] - Y3_z[indexingDiricheletz(i, j, k - 1)]) / (DZ));
            }
        }
    }
    boundary.divergence(Y3_x, Y3_y, Y3_z, Y2_p, time + 80.0 / 120.0 * DT, 16.0);

    poissonSolver.solveNeumannPoisson(Y2_p);
#endif

//TODO: same correction as y2, same treatment
    for (int i = 1; i < newDimX_x - 1; i++)
    {
        for (int j = 1; j < newDimY_x - 1; j++)
        {
            for (int k = 0; k < dim_z; k++)
            {
                Y3_x[indexingDiricheletx(i, j, k)] = Y3_x[indexingDiricheletx(i, j, k)] - 16.0 * DT / (120.0) * /*(d_Px(i,j,k,time+80.0 * DT / (120.0))-d_Px(i,j,k,time+64.0 * DT / (120.0)));*/(Y2_p[indexingDiricheletp(i + 1, j, k)] - Y2_p[indexingDiricheletp(i, j, k)]) / (DX);
            }
        }
    }

    for (int i = 1; i < newDimX_y - 1; i++)
    {
        for (int j = 1; j < newDimY_y - 1; j++)
        {
            for (int k = 0; k < dim_z; k++)
            {
                Y3_y[indexingDirichelety(i, j, k)] = Y3_y[indexingDirichelety(i, j, k)] - 16.0 * DT / (120.0) * /*(d_Py(i,j,k,time+80.0 * DT / (120.0))-d_Py(i,j,k,time+64.0 * DT / (120.0)));*/(Y2_p[indexingDiricheletp(i, j + 1, k)] - Y2_p[indexingDiricheletp(i, j, k)]) / (DY);
            }
        }
    }

    for (int i = 1; i < newDimX_z - 1; i++)
    {
        for (int j = 1; j < newDimY_z - 1; j++)
        {
            for (int k = 0; k < dim_z_z; k++)
            {
                Y3_z[indexingDiricheletz(i, j, k)] = Y3_z[indexingDiricheletz(i, j, k)] - 16.0 * DT / (120.0) * /*(d_Pz(i,j,k,time+80.0 * DT / (120.0))-d_Pz(i,j,k,time+64.0 * DT / (120.0)));*/(Y2_p[indexingDiricheletp(i, j, k + 1)] - Y2_p[indexingDiricheletp(i, j, k)]) / (DZ);
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

    //TODO: check
    boundary.update_boundary(Y3_x, Y3_y, Y3_z, time + 80.0 / 120.0 * DT);
    MPI_Barrier(cart_comm);
    exchangeData(Y3_x, newDimX_x, newDimY_x, dim_z, MPI_face_x_x, MPI_face_y_x,0,1);
    exchangeData(Y3_y, newDimX_y, newDimY_y, dim_z, MPI_face_x_y, MPI_face_y_y,1,0);
    exchangeData(Y3_z, newDimX_z, newDimY_z, dim_z_z, MPI_face_x_z, MPI_face_y_z,1,1);
    MPI_Barrier(cart_comm);

    //TODO: same as earlier
    for (int i = 1 + lbx; i < newDimX_x - 1 - rbx; i++)
    {
        for (int j = 1 + lby; j < newDimY_x - 1 - rby; j++)
        {
            for (int k = 1; k < dim_z-1; k++)

            {
                grid.u[indexingDiricheletx(i, j, k)] = Y3_x[indexingDiricheletx(i, j, k)] +
                                                       90.0 / 120.0 * DT * functionF_u(Y3_x, Y3_y, Y3_z, i, j, k, time + 80.0 / 120.0 * DT) -
                                                       50.0 / 120.0 * DT * functionF_u(Y2_x, Y2_y, Y2_z, i, j, k, time + 64.0 / 120.0 * DT) -
                                                                     40.0 / 120.0 * DT * /*d_Px(i,j,k,time + 80.0/120.0*DT);*/(Phi_p[(i+1) * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k] - Phi_p[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k]) / (DX);
            }
        }
    }

//TODO: same
    for (int i = 1 + lbx; i < newDimX_y - 1 - rbx; i++)
    {
        for (int j = 1 + lby; j < newDimY_y - 1 - rby; j++)
        {
            for (int k = 1; k < dim_z-1; k++)
            {
                grid.v[indexingDirichelety(i, j, k)] = Y3_y[indexingDirichelety(i, j, k)] +
                                                       90.0 / 120.0 * DT * functionF_v(Y3_x, Y3_y, Y3_z, i, j, k, time + 80.0 / 120.0 * DT) -
                                                       50.0 / 120.0 * DT * functionF_v(Y2_x, Y2_y, Y2_z, i, j, k, time + 64.0 / 120.0 * DT) -
                                                               40.0 / 120.0 * DT * /*d_Py(i,j,k,time + 80.0/120.0*DT);*/(Phi_p[i * (NY + 1) * (NZ + 1) + (j+1) * (NZ + 1) + k] - Phi_p[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k]) / (DY);
            }
        }
    }

    for (int i = 1 + lbx; i < newDimX_z - 1 - rbx; i++)
    {
        for (int j = 1 + lby; j < newDimY_z - 1 - rby; j++)
        {
            for (int k = 1; k < dim_z_z - 1; k++)
            {
                grid.w[indexingDiricheletz(i, j, k)] = Y3_z[indexingDiricheletz(i, j, k)] +
                                                             90.0 / 120.0 * DT * functionF_w(Y3_x, Y3_y, Y3_z, i, j, k, time + 80.0 / 120.0 * DT) -
                                                             50.0 / 120.0 * DT * functionF_w(Y2_x, Y2_y, Y2_z, i, j, k, time + 64.0 / 120.0 * DT) -
                                                         40.0 / 120.0 * DT * /*d_Pz(i,j,k,time + 80.0/120.0*DT);*/(Phi_p[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k+1] - Phi_p[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k]) / (DZ);
           }
        }
    }

    //TODO: same work with pressure
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
    poissonSolver.solveDirichletPoisson(Y2_p,helper);
#endif
#ifdef DIRICHELET
    boundary.update_boundary(grid.u, grid.v, grid.w, time + DT);

    for (int i = 1; i < NX; i++)
    {
        for (int j = 1; j < NY; j++)
        {
            for (int k = 1; k < NZ; k++)
            {
                Y2_p[i * (NY+1) * (NZ+1) + j * (NZ+1) + k] = 120.0 / (40.0 * DT) * ((grid.u[indexingDiricheletx(i, j, k)] - grid.u[indexingDiricheletx(i - 1, j, k)]) / (DX) + (grid.v[indexingDirichelety(i, j, k)] - grid.v[indexingDirichelety(i, j - 1, k)]) / (DY) + (grid.w[indexingDiricheletz(i, j, k)] - grid.w[indexingDiricheletz(i, j, k - 1)]) / (DZ));
            }
        }
    }
    boundary.divergence(grid.u, grid.v, grid.w, Y2_p, time + DT, 40.0);

    poissonSolver.solveNeumannPoisson(Y2_p);
#endif

//TODO: same correction as before to adapt
    for(int i = 1; i < newDimX_x - 1; i++)
    {
        for(int j = 1; j < newDimY_x - 1; j++)
        {
            for(int k = 0; k < dim_z; k++)
            {
                grid.u[indexingDiricheletx(i,j,k)] -= 40.0 * DT / (120.0) * /*(d_Px(i,j,k,time+DT)-d_Px(i,j,k,time+80.0 * DT / (120.0)));*/(Y2_p[indexingDiricheletp(i + 1, j, k)] - Y2_p[indexingDiricheletp(i, j, k)]) / (DX); 
            }
        }
    }
    for(int i = 1; i < newDimX_y - 1; i++)
    {
        for(int j = 1; j < newDimY_y - 1; j++)
        {
            for(int k = 0; k < dim_z; k++)
            {
                grid.v[indexingDirichelety(i,j,k)] -= 40.0 * DT / (120.0) * /*(d_Px(i,j,k,time+DT)-d_Px(i,j,k,time+80.0 * DT / (120.0)));*/(Y2_p[indexingDiricheletp(i, j + 1, k)] - Y2_p[indexingDiricheletp(i, j, k)]) / (DY); 
            }
        }
    }
    for(int i = 1; i < newDimX_z - 1; i++)
    {
        for(int j = 1; j < newDimY_z - 1; j++)
        {
            for(int k = 0; k < dim_z_z; k++)
            {
                grid.w[indexingDiricheletz(i,j,k)] -= 40.0 * DT / (120.0) * /*(d_Px(i,j,k,time+DT)-d_Px(i,j,k,time+80.0 * DT / (120.0)));*/(Y2_p[indexingDiricheletp(i, j, k + 1)] - Y2_p[indexingDiricheletp(i, j, k)]) / (DZ); 
            }
        }
    }
    for(int i = 0; i < NX + 1; i++)
    {
        for(int j = 0; j < NY + 1; j++)
        {
            for(int k = 0; k < NZ + 1; k++)
            {
#ifdef PERIODIC
                grid.p[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k] = Y2_p[indexingPeriodicp(i, j, k)]  + Phi_p[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k]; 
#endif
#ifdef DIRICHELET
                grid.p[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k] = Y2_p[indexingDiricheletp(i, j, k)]  + Phi_p[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k]; 
#endif
            }
        }
    }
}

//TODO: maybe change everything with the indexing function(should not be mandatory though)
Real IcoNS::functionF_u(const std::vector<Real> &u, const std::vector<Real> &v, const std::vector<Real> &w, int i, int j, int k, Real t)
{
    int lu = i * newDimY_x * dim_z + j * dim_z + k;
    int lv = i * newDimY_y * dim_z + j * dim_z + k;
    int lw = i * newDimY_z * dim_z_z + j * dim_z_z + k;

    return -(u[lu] * (u[lu + newDimY_x * dim_z] - u[lu - newDimY_x * dim_z]) / (2.0 * DX) +
             (v[lv] + v[lv + newDimY_y * dim_z] + v[lv - dim_z] + v[lv + newDimY_y * dim_z - dim_z]) / 4.0 * (u[lu + dim_z] - u[lu - dim_z]) / (2.0 * DY) +
             (w[lw] + w[lw + newDimY_z * dim_z_z] + w[lw - 1] + w[lw + newDimY_z * dim_z_z - 1]) / 4.0 * (u[lu + 1] - u[lu - 1]) / (2.0 * DZ)) +
           1 / RE * ((u[lu + newDimY_x * dim_z] - 2 * u[lu] + u[lu - newDimY_x * dim_z]) / (DX * DX) + (u[lu + dim_z] - 2 * u[lu] + u[lu - dim_z]) / (DY * DY) + (u[lu + 1] - 2 * u[lu] + u[lu - 1]) / (DZ * DZ)) +
           functionG_u(i-1 + coords[0] * other_dim_x_x, j-1 + coords[1] * other_dim_y_x, k, t);
}

Real IcoNS::functionF_v(const std::vector<Real> &u, const std::vector<Real> &v, const std::vector<Real> &w, int i, int j, int k, Real t)
{
    int lu = i * newDimY_x * dim_z + j * dim_z + k;
    int lv = i * newDimY_y * dim_z + j * dim_z + k;
    int lw = i * newDimY_z * dim_z_z + j * dim_z_z + k;
    
    return -((u[lu] + u[lu + dim_z] + u[lu - newDimY_x * dim_z] + u[lu - newDimY_x * dim_z + dim_z]) / 4.0 * ((v[lv + newDimY_y * dim_z] - v[lv - newDimY_y * dim_z]) / (2.0 * DX)) +
             v[lv] * (v[lv + dim_z] - v[lv - dim_z]) / (2.0 * DY) +
             (w[lw] + w[lw - 1] + w[lw + dim_z_z] + w[lw + dim_z_z - 1]) / 4.0 * (v[lv + 1] - v[lv - 1]) / (2.0 * DZ)) +
           (1.0 / RE) * ((v[lv + newDimY_y * dim_z] - 2.0 * v[lv] + v[lv - newDimY_y * dim_z]) / (DX * DX) +
                         (v[lv + dim_z] - 2.0 * v[lv] + v[lv - dim_z]) / (DY * DY) +
                         (v[lv + 1] - 2.0 * v[lv] + v[lv - 1]) / (DZ * DZ)) +
           functionG_v(i-1 + coords[0] * other_dim_x_y, j-1 + coords[1] * other_dim_y_y, k, t);
}

Real IcoNS::functionF_w(const std::vector<Real> &u, const std::vector<Real> &v, const std::vector<Real> &w, int i, int j, int k, Real t)

{
    int lu = i * newDimY_x * dim_z + j * dim_z + k;
    int lv = i * newDimY_y * dim_z + j * dim_z + k;
    int lw = i * newDimY_z * dim_z_z + j * dim_z_z + k;

    return -((u[lu] + u[lu - newDimY_x * dim_z] + u[lu + 1] + u[lu - dim_z * newDimY_x + 1]) / 4.0 * (w[lw + newDimY_z * dim_z_z] - w[lw - newDimY_z * dim_z_z]) / (2.0 * DX) +
             (v[lv + 1] + v[lv - dim_z + 1] + v[lv] + v[lv - dim_z]) / 4.0 * (w[lw + dim_z_z] - w[lw - dim_z_z]) / (2.0 * DY) +
             w[lw] * (w[lw + 1] - w[lw - 1]) / (2.0 * DZ)) +
           (1.0 / RE) * ((w[lw + newDimY_z * dim_z_z] - 2.0 * w[lw] + w[lw - newDimY_z * dim_z_z]) / (DX * DX) +
                         (w[lw + dim_z_z] - 2.0 * w[lw] + w[lw - dim_z_z]) / (DY * DY) +
                         (w[lw + 1] - 2.0 * w[lw] + w[lw - 1]) / (DZ * DZ)) +
           functionG_w(i-1 + coords[0] * other_dim_x_z, j-1 + coords[1] * other_dim_y_z, k, t);

}

//TODO: added a term, check just in case
Real IcoNS::functionG_u(int i, int j, int k, Real t)
{
    Real x = i * DX + DX / 2;
    Real y = j * DY;
    Real z = k * DZ;

    return std::sin(x) * std::cos(y) * std::sin(z) * std::cos(t) +
           std::sin(x) * std::cos(x) * std::cos(y) * std::cos(y) * std::sin(z) * std::sin(z) * std::sin(t) * std::sin(t) -
           std::sin(x) * std::cos(x) * std::sin(y) * std::sin(y) * std::sin(z) * std::sin(z) * std::sin(t) * std::sin(t) +
           2 * std::sin(x) * std::cos(x) * std::cos(y) * std::cos(y) * std::cos(z) * std::cos(z) * std::sin(t) * std::sin(t) +
           3.0 / RE * std::sin(x) * std::cos(y) * std::sin(z) * std::sin(t) - std::sin(x) * std::cos(y) * std::sin(z) * std::sin(t);
}

Real IcoNS::functionG_v(int i, int j, int k, Real t)
{
    Real x = i * DX;
    Real y = j * DY + DY / 2;
    Real z = k * DZ;

    return std::cos(x) * std::sin(y) * std::sin(z) * std::cos(t) -
           std::sin(x) * std::sin(x) * std::sin(y) * std::cos(y) * std::sin(z) * std::sin(z) * std::sin(t) * std::sin(t) +
           std::cos(x) * std::cos(x) * std::sin(y) * std::cos(y) * std::sin(z) * std::sin(z) * std::sin(t) * std::sin(t) +
           2.0 * std::cos(x) * std::cos(x) * std::sin(y) * std::cos(y) * std::cos(z) * std::cos(z) * std::sin(t) * std::sin(t) +
           3.0 / RE * std::cos(x) * std::sin(y) * std::sin(z) * std::sin(t) - std::cos(x) * std::sin(y) * std::sin(z) * std::sin(t);
}

Real IcoNS::functionG_w(int i, int j, int k, Real t)
{
    Real x = i * DX;
    Real y = j * DY;
    Real z = k * DZ + DZ / 2;

    return 2 * std::cos(x) * std::cos(y) * std::cos(z) * std::cos(t) -
           2 * std::sin(x) * std::sin(x) * std::cos(y) * std::cos(y) * std::sin(z) * std::cos(z) * std::sin(t) * std::sin(t) -
           2 * std::cos(x) * std::cos(x) * std::sin(y) * std::sin(y) * std::sin(z) * std::cos(z) * std::sin(t) * std::sin(t) -
           4.0 * std::cos(x) * std::cos(x) * std::cos(y) * std::cos(y) * std::sin(z) * std::cos(z) * std::sin(t) * std::sin(t) +
           6.0 / RE * std::cos(x) * std::cos(y) * std::cos(z) * std::sin(t) + std::cos(x) * std::cos(y) * std::cos(z) * std::sin(t);
}
//TODO: exact derivative missing, don't know if needed