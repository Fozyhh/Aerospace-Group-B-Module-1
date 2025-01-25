#include "core.hpp"
#include <fstream>

void IcoNS::solve_time_step(Real time)
{

    // 1) pressure point exchange
    copyPressureToHalo(grid.p, halo_p);
    // MPI_Barrier(cart_comm);
    exchangeData(halo_p, (xSize[2] + 2), (xSize[1] + 2), xSize[0], MPI_face_x_p, MPI_face_y_p, 1, 1);

    // Calculate Y2
    for (int i = 1 + lbx; i < newDimX_x - 1 - rbx; i++)
    {
        for (int j = 1 + lby; j < newDimY_x - 1 - rby; j++)
        {
            for (int k = lbz; k < dim_z - rbz; k++)
            {
                Y2_x[getx(i, j, k)] = grid.u[getx(i, j, k)] +
                                      64.0 / 120.0 * DT * functionF_u(grid.u, grid.v, grid.w, i, j, k, time) -
                                      64.0 / 120.0 * DT * (halo_p[getHaloP(i + 1 - resx, j, k)] - halo_p[getHaloP(i - resx, j, k)]) / (DX);
            }
        }
    }

    for (int i = 1 + lbx; i < newDimX_y - 1 - rbx; i++)
    {
        for (int j = 1 + lby; j < newDimY_y - 1 - rby; j++)
        {
            for (int k = lbz; k < dim_z - rbz; k++)
            {
                Y2_y[gety(i, j, k)] = grid.v[gety(i, j, k)] +
                                      64.0 / 120.0 * DT * functionF_v(grid.u, grid.v, grid.w, i, j, k, time) -
                                      64.0 / 120.0 * DT * (halo_p[getHaloP(i, j + 1 - resy, k)] - halo_p[getHaloP(i, j - resy, k)]) / (DY);
            }
        }
    }

    for (int i = 1 + lbx; i < newDimX_z - 1 - rbx; i++)
    {
        for (int j = 1 + lby; j < newDimY_z - 1 - rby; j++)
        {
            for (int k = lbz; k < dim_z_z - rbz; k++)
            {
                Y2_z[getz(i, j, k)] = grid.w[getz(i, j, k)] +
                                      64.0 / 120.0 * DT * functionF_w(grid.u, grid.v, grid.w, i, j, k, time) -
                                      64.0 / 120.0 * DT * (halo_p[getHaloP(i, j, k + 1)] - halo_p[getHaloP(i, j, k)]) / (DZ);
            }
        }
    }

    // 3) Update boundaries
    boundary.update_boundary(Y2_x, Y2_y, Y2_z, time + 64.0 / 120.0 * DT);
    // MPI_Barrier(cart_comm);
    // Update Halos
    exchangeData(Y2_x, newDimX_x, newDimY_x, dim_z, MPI_face_x_x, MPI_face_y_x, 0, 1);
    exchangeData(Y2_y, newDimX_y, newDimY_y, dim_z, MPI_face_x_y, MPI_face_y_y, 1, 0);
    exchangeData(Y2_z, newDimX_z, newDimY_z, dim_z_z, MPI_face_x_z, MPI_face_y_z, 1, 1);

    for (int i = 1; i < xSize[2] + 1; i++)
    {
        for (int j = 1; j < xSize[1] + 1; j++)
        {
            for (int k = 0; k < xSize[0]; k++)
            {
                if ((lbx && i == 1) || (lby && j == 1) || (lbz && k == 0) || (rbx && i == xSize[2]) || (rby && j == xSize[1]) || (rbz && k == xSize[0] - 1))
                {
                    Y2_p[getp(i - 1, j - 1, k)] = 0.0;
                }
                else
                {
                    Y2_p[getp(i - 1, j - 1, k)] = 120.0 / (64.0 * DT) * ((Y2_x[getx(i + resx, j, k)] - Y2_x[getx(i - 1 + resx, j, k)]) / (DX) + (Y2_y[gety(i, j + resy, k)] - Y2_y[gety(i, j - 1 + resy, k)]) / (DY) + (Y2_z[getz(i, j, k)] - Y2_z[getz(i, j, k - 1)]) / (DZ));
                }
            }
        }
    }

    // Solve for Pressure
    poissonSolver->solveNeumannPoisson(Y2_p);
    // MPI_Barrier(cart_comm);
    copyPressureToHalo(Y2_p, halo_p);
    // MPI_Barrier(cart_comm);
    exchangeData(halo_p, (xSize[2] + 2), (xSize[1] + 2), xSize[0], MPI_face_x_p, MPI_face_y_p, 1, 1);

    // Update Velocities
    for (int i = 1; i < newDimX_x - 1; i++)
    {
        for (int j = 1; j < newDimY_x - 1; j++)
        {
            for (int k = 0; k < dim_z; k++)
            {
                Y2_x[getx(i, j, k)] = Y2_x[getx(i, j, k)] -
                                      64.0 * DT / (120.0) * (halo_p[getHaloP(i + 1 - resx, j, k)] - halo_p[getHaloP(i - resx, j, k)]) / (DX);
            }
        }
    }

    for (int i = 1; i < newDimX_y - 1; i++)
    {
        for (int j = 1; j < newDimY_y - 1; j++)
        {
            for (int k = 0; k < dim_z; k++)
            {
                Y2_y[gety(i, j, k)] = Y2_y[gety(i, j, k)] -
                                      64.0 * DT / (120.0) * (halo_p[getHaloP(i, j + 1 - resy, k)] - halo_p[getHaloP(i, j - resy, k)]) / (DY);
            }
        }
    }

    for (int i = 1; i < newDimX_z - 1; i++)
    {
        for (int j = 1; j < newDimY_z - 1; j++)
        {
            for (int k = 0; k < dim_z_z; k++)
            {
                Y2_z[getz(i, j, k)] = Y2_z[getz(i, j, k)] -
                                      64.0 * DT / (120.0) * (halo_p[getHaloP(i, j, k + 1)] - halo_p[getHaloP(i, j, k)]) / (DZ);
            }
        }
    }

    for (int i = 0; i < xSize[2]; i++)
    {
        for (int j = 0; j < xSize[1]; j++)
        {
            for (int k = 0; k < xSize[0]; k++)
            {
                Phi_p[getp(i, j, k)] = Y2_p[getp(i, j, k)] + grid.p[getp(i, j, k)];
            }
        }
    }

    // MPI_Barrier(cart_comm);
    exchangeData(Y2_x, newDimX_x, newDimY_x, dim_z, MPI_face_x_x, MPI_face_y_x, 0, 1);
    exchangeData(Y2_y, newDimX_y, newDimY_y, dim_z, MPI_face_x_y, MPI_face_y_y, 1, 0);
    exchangeData(Y2_z, newDimX_z, newDimY_z, dim_z_z, MPI_face_x_z, MPI_face_y_z, 1, 1);

    // 3) Phi_p exchange
    copyPressureToHalo(Phi_p, halo_phi);
    // MPI_Barrier(cart_comm);
    exchangeData(halo_phi, (xSize[2] + 2), (xSize[1] + 2), xSize[0], MPI_face_x_p, MPI_face_y_p, 1, 1);

    for (int i = 1 + lbx; i < newDimX_x - 1 - rbx; i++)
    {
        for (int j = 1 + lby; j < newDimY_x - 1 - rby; j++)
        {
            for (int k = lbz; k < dim_z - rbz; k++)

            {
                Y3_x[getx(i, j, k)] = Y2_x[getx(i, j, k)] +
                                      50.0 / 120.0 * DT * functionF_u(Y2_x, Y2_y, Y2_z, i, j, k, time + 64.0 / 120.0 * DT) -
                                      34.0 / 120.0 * DT * functionF_u(grid.u, grid.v, grid.w, i, j, k, time) -
                                      16.0 / 120.0 * DT * (halo_phi[getHaloP(i + 1 - resx, j, k)] - halo_phi[getHaloP(i - resx, j, k)]) / (DX);
            }
        }
    }

    for (int i = 1 + lbx; i < newDimX_y - 1 - rbx; i++)
    {
        for (int j = 1 + lby; j < newDimY_y - 1 - rby; j++)
        {
            for (int k = lbz; k < dim_z - rbz; k++)
            {
                Y3_y[gety(i, j, k)] = Y2_y[gety(i, j, k)] +
                                      50.0 / 120.0 * DT * functionF_v(Y2_x, Y2_y, Y2_z, i, j, k, time + 64.0 / 120.0 * DT) -
                                      34.0 / 120.0 * DT * functionF_v(grid.u, grid.v, grid.w, i, j, k, time) -
                                      16.0 / 120.0 * DT * (halo_phi[getHaloP(i, j + 1 - resy, k)] - halo_phi[getHaloP(i, j - resy, k)]) / (DY);
            }
        }
    }

    for (int i = 1 + lbx; i < newDimX_z - 1 - rbx; i++)
    {
        for (int j = 1 + lby; j < newDimY_z - 1 - rby; j++)
        {
            for (int k = lbz; k < dim_z_z - rbz; k++)

            {
                Y3_z[getz(i, j, k)] = Y2_z[getz(i, j, k)] +
                                      50.0 / 120.0 * DT * functionF_w(Y2_x, Y2_y, Y2_z, i, j, k, time + 64.0 / 120.0 * DT) -
                                      34.0 / 120.0 * DT * functionF_w(grid.u, grid.v, grid.w, i, j, k, time) -
                                      16.0 / 120.0 * DT * (halo_phi[getHaloP(i, j, k + 1)] - halo_phi[getHaloP(i, j, k)]) / (DZ);
            }
        }
    }

    boundary.update_boundary(Y3_x, Y3_y, Y3_z, time + 80.0 / 120.0 * DT);
    // MPI_Barrier(cart_comm);
    exchangeData(Y3_x, newDimX_x, newDimY_x, dim_z, MPI_face_x_x, MPI_face_y_x, 0, 1);
    exchangeData(Y3_y, newDimX_y, newDimY_y, dim_z, MPI_face_x_y, MPI_face_y_y, 1, 0);
    exchangeData(Y3_z, newDimX_z, newDimY_z, dim_z_z, MPI_face_x_z, MPI_face_y_z, 1, 1);

    for (int i = 1; i < xSize[2] + 1; i++)
    {
        for (int j = 1; j < xSize[1] + 1; j++)
        {
            for (int k = 0; k < xSize[0]; k++)
            {
                if ((lbx && i == 1) || (lby && j == 1) || (lbz && k == 0) || (rbx && i == xSize[2]) || (rby && j == xSize[1]) || (rbz && k == xSize[0] - 1))
                {
                    Y2_p[getp(i - 1, j - 1, k)] = 0.0;
                }
                else
                {
                    Y2_p[getp(i - 1, j - 1, k)] = 120.0 / (16.0 * DT) * ((Y3_x[getx(i + resx, j, k)] - Y3_x[getx(i - 1 + resx, j, k)]) / (DX) + (Y3_y[gety(i, j + resy, k)] - Y3_y[gety(i, j + resy - 1, k)]) / (DY) + (Y3_z[getz(i, j, k)] - Y3_z[getz(i, j, k - 1)]) / (DZ));
                }
            }
        }
    }

    poissonSolver->solveNeumannPoisson(Y2_p);
    // MPI_Barrier(cart_comm);
    // 3) y2_p exchange
    copyPressureToHalo(Y2_p, halo_p);
    // MPI_Barrier(cart_comm);
    exchangeData(halo_p, (xSize[2] + 2), (xSize[1] + 2), xSize[0], MPI_face_x_p, MPI_face_y_p, 1, 1);

    for (int i = 1; i < newDimX_x - 1; i++)
    {
        for (int j = 1; j < newDimY_x - 1; j++)
        {
            for (int k = 0; k < dim_z; k++)
            {
                Y3_x[getx(i, j, k)] = Y3_x[getx(i, j, k)] -
                                      16.0 * DT / (120.0) * (halo_p[getHaloP(i + 1 - resx, j, k)] - halo_p[getHaloP(i - resx, j, k)]) / (DX);
            }
        }
    }

    for (int i = 1; i < newDimX_y - 1; i++)
    {
        for (int j = 1; j < newDimY_y - 1; j++)
        {
            for (int k = 0; k < dim_z; k++)
            {
                Y3_y[gety(i, j, k)] = Y3_y[gety(i, j, k)] -
                                      16.0 * DT / (120.0) * (halo_p[getHaloP(i, j + 1 - resy, k)] - halo_p[getHaloP(i, j - resy, k)]) / (DY);
            }
        }
    }

    for (int i = 1; i < newDimX_z - 1; i++)
    {
        for (int j = 1; j < newDimY_z - 1; j++)
        {
            for (int k = 0; k < dim_z_z; k++)
            {
                Y3_z[getz(i, j, k)] = Y3_z[getz(i, j, k)] -
                                      16.0 * DT / (120.0) * (halo_p[getHaloP(i, j, k + 1)] - halo_p[getHaloP(i, j, k)]) / (DZ);
            }
        }
    }

    for (int i = 0; i < xSize[2]; i++)
    {
        for (int j = 0; j < xSize[1]; j++)
        {
            for (int k = 0; k < xSize[0]; k++)
            {
                Phi_p[getp(i, j, k)] = Y2_p[getp(i, j, k)] + Phi_p[getp(i, j, k)]; // Phi_p=phi^3
            }
        }
    }

    // MPI_Barrier(cart_comm);
    exchangeData(Y3_x, newDimX_x, newDimY_x, dim_z, MPI_face_x_x, MPI_face_y_x, 0, 1);
    exchangeData(Y3_y, newDimX_y, newDimY_y, dim_z, MPI_face_x_y, MPI_face_y_y, 1, 0);
    exchangeData(Y3_z, newDimX_z, newDimY_z, dim_z_z, MPI_face_x_z, MPI_face_y_z, 1, 1);

    // 4) Phi_p exchange
    copyPressureToHalo(Phi_p, halo_phi);
    // MPI_Barrier(cart_comm);
    exchangeData(halo_phi, (xSize[2] + 2), (xSize[1] + 2), xSize[0], MPI_face_x_p, MPI_face_y_p, 1, 1);

    for (int i = 1 + lbx; i < newDimX_x - 1 - rbx; i++)
    {
        for (int j = 1 + lby; j < newDimY_x - 1 - rby; j++)
        {
            for (int k = lbz; k < dim_z - rbz; k++)

            {
                grid.u[getx(i, j, k)] = Y3_x[getx(i, j, k)] +
                                        90.0 / 120.0 * DT * functionF_u(Y3_x, Y3_y, Y3_z, i, j, k, time + 80.0 / 120.0 * DT) -
                                        50.0 / 120.0 * DT * functionF_u(Y2_x, Y2_y, Y2_z, i, j, k, time + 64.0 / 120.0 * DT) -
                                        40.0 / 120.0 * DT * (halo_phi[getHaloP(i + 1 - resx, j, k)] - halo_phi[getHaloP(i - resx, j, k)]) / (DX);
            }
        }
    }

    for (int i = 1 + lbx; i < newDimX_y - 1 - rbx; i++)
    {
        for (int j = 1 + lby; j < newDimY_y - 1 - rby; j++)
        {
            for (int k = lbz; k < dim_z - rbz; k++)
            {
                grid.v[gety(i, j, k)] = Y3_y[gety(i, j, k)] +
                                        90.0 / 120.0 * DT * functionF_v(Y3_x, Y3_y, Y3_z, i, j, k, time + 80.0 / 120.0 * DT) -
                                        50.0 / 120.0 * DT * functionF_v(Y2_x, Y2_y, Y2_z, i, j, k, time + 64.0 / 120.0 * DT) -
                                        40.0 / 120.0 * DT * (halo_phi[getHaloP(i, j + 1 - resy, k)] - halo_phi[getHaloP(i, j - resy, k)]) / (DY);
            }
        }
    }

    for (int i = 1 + lbx; i < newDimX_z - 1 - rbx; i++)
    {
        for (int j = 1 + lby; j < newDimY_z - 1 - rby; j++)
        {
            for (int k = lbz; k < dim_z_z - rbz; k++)
            {
                grid.w[getz(i, j, k)] = Y3_z[getz(i, j, k)] +
                                        90.0 / 120.0 * DT * functionF_w(Y3_x, Y3_y, Y3_z, i, j, k, time + 80.0 / 120.0 * DT) -
                                        50.0 / 120.0 * DT * functionF_w(Y2_x, Y2_y, Y2_z, i, j, k, time + 64.0 / 120.0 * DT) -
                                        40.0 / 120.0 * DT * (halo_phi[getHaloP(i, j, k + 1)] - halo_phi[getHaloP(i, j, k)]) / (DZ);
            }
        }
    }

    boundary.update_boundary(grid.u, grid.v, grid.w, time + DT);
    // MPI_Barrier(cart_comm);
    exchangeData(grid.u, newDimX_x, newDimY_x, dim_z, MPI_face_x_x, MPI_face_y_x, 0, 1);
    exchangeData(grid.v, newDimX_y, newDimY_y, dim_z, MPI_face_x_y, MPI_face_y_y, 1, 0);
    exchangeData(grid.w, newDimX_z, newDimY_z, dim_z_z, MPI_face_x_z, MPI_face_y_z, 1, 1);

    for (int i = 1; i < xSize[2] + 1; i++)
    {
        for (int j = 1; j < xSize[1] + 1; j++)
        {
            for (int k = 0; k < xSize[0]; k++)
            {
                if ((lbx && i == 1) || (lby && j == 1) || (lbz && k == 0) || (rbx && i == xSize[2]) || (rby && j == xSize[1]) || (rbz && k == xSize[0] - 1))
                {
                    Y2_p[getp(i - 1, j - 1, k)] = 0.0;
                }
                else
                {
                    Y2_p[getp(i - 1, j - 1, k)] = 120.0 / (40.0 * DT) * ((grid.u[getx(i + resx, j, k)] - grid.u[getx(i + resx - 1, j, k)]) / (DX) + (grid.v[gety(i, j + resy, k)] - grid.v[gety(i, j + resy - 1, k)]) / (DY) + (grid.w[getz(i, j, k)] - grid.w[getz(i, j, k - 1)]) / (DZ));
                }
            }
        }
    }

    poissonSolver->solveNeumannPoisson(Y2_p);
    // MPI_Barrier(cart_comm);
    copyPressureToHalo(Y2_p, halo_p);
    // MPI_Barrier(cart_comm);
    exchangeData(halo_p, (xSize[2] + 2), (xSize[1] + 2), xSize[0], MPI_face_x_p, MPI_face_y_p, 1, 1);

    for (int i = 1; i < newDimX_x - 1; i++)
    {
        for (int j = 1; j < newDimY_x - 1; j++)
        {
            for (int k = 0; k < dim_z; k++)
            {
                grid.u[getx(i, j, k)] -= 40.0 * DT / (120.0) * (halo_p[getHaloP(i + 1 - resx, j, k)] - halo_p[getHaloP(i - resx, j, k)]) / (DX);
            }
        }
    }
    for (int i = 1; i < newDimX_y - 1; i++)
    {
        for (int j = 1; j < newDimY_y - 1; j++)
        {
            for (int k = 0; k < dim_z; k++)
            {
                grid.v[gety(i, j, k)] -= 40.0 * DT / (120.0) * (halo_p[getHaloP(i, j + 1 - resy, k)] - halo_p[getHaloP(i, j - resy, k)]) / (DY);
            }
        }
    }
    for (int i = 1; i < newDimX_z - 1; i++)
    {
        for (int j = 1; j < newDimY_z - 1; j++)
        {
            for (int k = 0; k < dim_z_z; k++)
            {
                grid.w[getz(i, j, k)] -= 40.0 * DT / (120.0) * (halo_p[getHaloP(i, j, k + 1)] - halo_p[getHaloP(i, j, k)]) / (DZ);
            }
        }
    }
    for (int i = 0; i < xSize[2]; i++)
    {
        for (int j = 0; j < xSize[1]; j++)
        {
            for (int k = 0; k < xSize[0]; k++)
            {
                grid.p[getp(i, j, k)] = Y2_p[getp(i, j, k)] + Phi_p[getp(i, j, k)];
            }
        }
    }

    boundary.update_boundary(grid.u, grid.v, grid.w, time);
    // MPI_Barrier(cart_comm);
    exchangeData(grid.u, newDimX_x, newDimY_x, dim_z, MPI_face_x_x, MPI_face_y_x, 0, 1);
    exchangeData(grid.v, newDimX_y, newDimY_y, dim_z, MPI_face_x_y, MPI_face_y_y, 1, 0);
    exchangeData(grid.w, newDimX_z, newDimY_z, dim_z_z, MPI_face_x_z, MPI_face_y_z, 1, 1);
    // MPI_Barrier(cart_comm);
}

Real IcoNS::functionF_u(const std::vector<Real> &u, const std::vector<Real> &v, const std::vector<Real> &w, int i, int j, int k, Real t)
{
    Real value = -(u[getx(i, j, k)] * (u[getx(i + 1, j, k)] - u[getx(i - 1, j, k)]) / (2.0 * DX) +
                   (v[gety(i - resx, j + resy, k)] + v[gety(i + 1 - resx, j + resy, k)] + v[gety(i - resx, j - 1 + resy, k)] + v[gety(i + 1 - resx, j - 1 + resy, k)]) / 4.0 * (u[getx(i, j + 1, k)] - u[getx(i, j - 1, k)]) / (2.0 * DY) +
                   (w[getz(i - resx, j - resy, k)] + w[getz(i + 1 - resx, j - resy, k)] + w[getz(i - resx, j - resy, k - 1)] + w[getz(i + 1 - resx, j - resy, k - 1)]) / 4.0 * (u[getx(i, j, k + 1)] - u[getx(i, j, k - 1)]) / (2.0 * DZ)) +
                 1 / RE * ((u[getx(i + 1, j, k)] - 2 * u[getx(i, j, k)] + u[getx(i - 1, j, k)]) / (DX * DX) + (u[getx(i, j + 1, k)] - 2 * u[getx(i, j, k)] + u[getx(i, j - 1, k)]) / (DY * DY) + (u[getx(i, j, k + 1)] - 2 * u[getx(i, j, k)] + u[getx(i, j, k - 1)]) / (DZ * DZ));

    // if (testCase == 0)
    // {
    //     value += functionG_u(i - 1 + offset_x_x, j - 1 + offset_y_x, k, t);
    // }
    return value;
}

Real IcoNS::functionF_v(const std::vector<Real> &u, const std::vector<Real> &v, const std::vector<Real> &w, int i, int j, int k, Real t)
{
    Real value = -((u[getx(i + resx, j - resy, k)] + u[getx(i + resx, j - resy + 1, k)] + u[getx(i + resx - 1, j - resy, k)] + u[getx(i + resx - 1, j - resy + 1, k)]) / 4.0 * (v[gety(i + 1, j, k)] - v[gety(i - 1, j, k)]) / (2.0 * DX) +
                   v[gety(i, j, k)] * (v[gety(i, j + 1, k)] - v[gety(i, j - 1, k)]) / (2.0 * DY) +
                   (w[getz(i, j - resy, k)] + w[getz(i, j - resy + 1, k)] + w[getz(i, j - resy, k - 1)] + w[getz(i, j - resy + 1, k - 1)]) / 4.0 * (v[gety(i, j, k + 1)] - v[gety(i, j, k - 1)]) / (2.0 * DZ)) +
                 1 / RE * ((v[gety(i + 1, j, k)] - 2 * v[gety(i, j, k)] + v[gety(i - 1, j, k)]) / (DX * DX) + (v[gety(i, j + 1, k)] - 2 * v[gety(i, j, k)] + v[gety(i, j - 1, k)]) / (DY * DY) + (v[gety(i, j, k + 1)] - 2 * v[gety(i, j, k)] + v[gety(i, j, k - 1)]) / (DZ * DZ));

    // if (testCase == 0)
    // {
    //     value += functionG_v(i - 1 + offset_x_y, j - 1 + offset_y_y, k, t);
    // }

    return value;
}

Real IcoNS::functionF_w(const std::vector<Real> &u, const std::vector<Real> &v, const std::vector<Real> &w, int i, int j, int k, Real t)

{
    Real value = -((u[getx(i + resx, j, k)] + u[getx(i + resx, j, k + 1)] + u[getx(i + resx - 1, j, k)] + u[getx(i + resx - 1, j, k + 1)]) / 4.0 * (w[getz(i + 1, j, k)] - w[getz(i - 1, j, k)]) / (2.0 * DX) +
                   (v[gety(i, j + resy, k)] + v[gety(i, j + resy - 1, k)] + v[gety(i, j + resy, k + 1)] + v[gety(i, j + resy - 1, k + 1)]) / 4.0 * (w[getz(i, j + 1, k)] - w[getz(i, j - 1, k)]) / (2.0 * DY) +
                   w[getz(i, j, k)] * (w[getz(i, j, k + 1)] - w[getz(i, j, k - 1)]) / (2.0 * DZ)) +
                 1 / RE * ((w[getz(i + 1, j, k)] - 2 * w[getz(i, j, k)] + w[getz(i - 1, j, k)]) / (DX * DX) + (w[getz(i, j + 1, k)] - 2 * w[getz(i, j, k)] + w[getz(i, j - 1, k)]) / (DY * DY) + (w[getz(i, j, k + 1)] - 2 * w[getz(i, j, k)] + w[getz(i, j, k - 1)]) / (DZ * DZ));
    // if (testCase == 0)
    // {
    //     value += functionG_w(i - 1 + offset_x_z, j - 1 + offset_y_z, k, t);
    // }
    return value;
}

Real IcoNS::functionG_u(int i, int j, int k, Real t)
{
    Real x = i * DX + DX / 2;
    Real y = j * DY;
    Real z = k * DZ;

    return std::sin(x) * std::cos(y) * std::sin(z) * std::cos(t) +
           std::sin(x) * std::cos(x) * std::cos(y) * std::cos(y) * std::sin(z) * std::sin(z) * std::sin(t) * std::sin(t) -
           std::sin(x) * std::cos(x) * std::sin(y) * std::sin(y) * std::sin(z) * std::sin(z) * std::sin(t) * std::sin(t) +
           2 * std::sin(x) * std::cos(x) * std::cos(y) * std::cos(y) * std::cos(z) * std::cos(z) * std::sin(t) * std::sin(t) +
           3.0 / RE * std::sin(x) * std::cos(y) * std::sin(z) * std::sin(t) - std::sin(x) * std::cos(y) * std::cos(z) * std::sin(t);
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
           3.0 / RE * std::cos(x) * std::sin(y) * std::sin(z) * std::sin(t) - std::cos(x) * std::sin(y) * std::cos(z) * std::sin(t);
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
           6.0 / RE * std::cos(x) * std::cos(y) * std::cos(z) * std::sin(t) - std::cos(x) * std::cos(y) * std::sin(z) * std::sin(t);
}
