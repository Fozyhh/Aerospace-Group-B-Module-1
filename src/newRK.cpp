#include "core.hpp"

void IcoNS::solve_time_step(Real time)
{
    int miss = 0;
    for (int i = 0; i < NX; i++)
    {
        for (int j = 0; j < NY + 1; j++)
        {
            for (int k = 0; k < NZ + 1; k++)
            {
                if (j == 0 && k == 0)
                {
                    if (grid.u[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k] != grid.u[i * (NY + 1) * (NZ + 1) + NY * (NZ + 1) + NZ])
                    {
                        miss++;
                        //std::cout << "Difference: " << grid.u[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k] << " - " << grid.u[i * (NY + 1) * (NZ + 1) + NY * (NZ + 1) + NZ] << std::endl;
                    }
                }
            }
        }
    }

    std::cout << "Miss Start: " << miss << std::endl;
    for (int i = 0; i < NX; i++)
    {
        for (int j = 0; j < NY + 1; j++)
        {
            for (int k = 0; k < NZ + 1; k++)
            {
                Y2_x[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k] = grid.u[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k] +
                                                                   64.0 / 120.0 * DT * functionF_u(grid.u, grid.v, grid.w, i % NX, j % NY, k % NZ, time); /* -
                                                                         64.0/120.0*DT*(grid.p[(i+1) * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k]-grid.p[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k])/DX; */
            }
        }
    }

    for (int i = 0; i < NX + 1; i++)
    {
        for (int j = 0; j < NY; j++)
        {
            for (int k = 0; k < NZ + 1; k++)
            {
                Y2_y[i * NY * (NZ + 1) + j * (NZ + 1) + k] = grid.v[i * NY * (NZ + 1) + j * (NZ + 1) + k] +
                                                             64.0 / 120.0 * DT * functionF_v(grid.u, grid.v, grid.w, i % NX, j % NY, k % NZ, time); /* -
                                                                   64.0/120.0*DT*(grid.p[i * (NY + 1) * (NZ + 1) + (j+1) * (NZ + 1) + k]-grid.p[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k])/DY; */
            }
        }
    }

    for (int i = 0; i < NX + 1; i++)
    {
        for (int j = 0; j < NY + 1; j++)
        {
            for (int k = 0; k < NZ; k++)
            {
                Y2_z[i * (NY + 1) * NZ + j * NZ + k] = grid.w[i * (NY + 1) * NZ + j * NZ + k] +
                                                       64.0 / 120.0 * DT * functionF_w(grid.u, grid.v, grid.w, i % NX, j % NY, k % NZ, time); /* -
                                                             64.0/120.0*DT*(grid.p[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k+1]-grid.p[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k])/DZ; */
            }
        }
    }

    /* for (int i = 0; i < NX + 1; i++)
    {
        for (int j = 0; j < NY + 1; j++)
        {
            for (int k = 0; k < NZ + 1; k++)
            {
                if (int(i == 0) || int(i == NX) || int(j == 0) || int(j == NY) || int(k == 0) || int(k == NZ))
                {
                    Y2_x[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k] = ;
                }
            }
        }
    } */

    // boundary.update_boundary(Y2_x, Y2_y, Y2_z, time + 64.0 / 120.0 * DT);
    miss = 0;
    for (int i = 0; i < NX; i++)
    {
        for (int j = 0; j < NY + 1; j++)
        {
            for (int k = 0; k < NZ + 1; k++)
            {
                if (j == 0 && k == 0)
                {
                    if (Y2_x[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k] != Y2_x[i * (NY + 1) * (NZ + 1) + NY * (NZ + 1) + NZ])
                    {
                        miss++;
                        //std::cout << "Difference: " << Y2_x[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k] << " - " << Y2_x[i * (NY + 1) * (NZ + 1) + NY * (NZ + 1) + NZ] << std::endl;
                    }
                }
            }
        }
    }

    std::cout << "Miss first step: " << miss << std::endl;

    for (int i = 0; i < NX; i++)
    {
        for (int j = 0; j < NY + 1; j++)
        {
            for (int k = 0; k < NZ + 1; k++)
            {
                Y3_x[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k] = Y2_x[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k] +
                                                                   50.0 / 120.0 * DT * functionF_u(Y2_x, Y2_y, Y2_z, i % NX, j % NY, k % NZ, time + 64.0 / 120.0 * DT) -
                                                                   34.0 / 120.0 * DT * functionF_u(grid.u, grid.v, grid.w, i % NX, j % NY, k % NZ, time);
            }
        }
    }

    for (int i = 0; i < NX + 1; i++)
    {
        for (int j = 0; j < NY; j++)
        {
            for (int k = 0; k < NZ + 1; k++)
            {
                Y3_y[i * NY * (NZ + 1) + j * (NZ + 1) + k] = Y2_y[i * NY * (NZ + 1) + j * (NZ + 1) + k] +
                                                             50.0 / 120.0 * DT * functionF_v(Y2_x, Y2_y, Y2_z, i % NX, j % NY, k % NZ, time + 64.0 / 120.0 * DT) -
                                                             34.0 / 120.0 * DT * functionF_v(grid.u, grid.v, grid.w, i % NX, j % NY, k % NZ, time);
            }
        }
    }

    for (int i = 0; i < NX + 1; i++)
    {
        for (int j = 0; j < NY + 1; j++)
        {
            for (int k = 0; k < NZ; k++)
            {
                Y3_z[i * (NY + 1) * NZ + j * NZ + k] = Y2_z[i * (NY + 1) * NZ + j * NZ + k] +
                                                       50.0 / 120.0 * DT * functionF_w(Y2_x, Y2_y, Y2_z, i % NX, j % NY, k % NZ, time + 64.0 / 120.0 * DT) -
                                                       34.0 / 120.0 * DT * functionF_w(grid.u, grid.v, grid.w, i % NX, j % NY, k % NZ, time);
            }
        }
    }

    miss = 0;
    for (int i = 0; i < NX; i++)
    {
        for (int j = 0; j < NY + 1; j++)
        {
            for (int k = 0; k < NZ + 1; k++)
            {
                if (j == 0 && k == 0)
                {
                    if (Y3_x[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k] != Y3_x[i * (NY + 1) * (NZ + 1) + NY * (NZ + 1) + NZ])
                    {
                        miss++;
                        //std::cout << "Difference: " << Y3_x[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k] << " - " << Y3_x[i * (NY + 1) * (NZ + 1) + NY * (NZ + 1) + NZ] << std::endl;
                    }
                }
            }
        }
    }

    std::cout << "Miss Second step: " << miss << std::endl;

    // boundary.update_boundary(Y3_x, Y3_y, Y3_z, time + 80.0 / 120.0 * DT);

    for (int i = 0; i < NX; i++)
    {
        for (int j = 0; j < NY + 1; j++)
        {
            for (int k = 0; k < NZ + 1; k++)
            {
                grid.u[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k] = Y3_x[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k] +
                                                                     90.0 / 120.0 * DT * functionF_u(Y3_x, Y3_y, Y3_z, i % NX, j % NY, k % NZ, time + 80.0 / 120.0 * DT) -
                                                                     50.0 / 120.0 * DT * functionF_u(Y2_x, Y2_y, Y2_z, i % NX, j % NY, k % NZ, time + 64.0 / 120.0 * DT);
            }
        }
    }

    for (int i = 0; i < NX + 1; i++)
    {
        for (int j = 0; j < NY; j++)
        {
            for (int k = 0; k < NZ + 1; k++)
            {
                grid.v[i * NY * (NZ + 1) + j * (NZ + 1) + k] = Y3_y[i * NY * (NZ + 1) + j * (NZ + 1) + k] +
                                                               90.0 / 120.0 * DT * functionF_v(Y3_x, Y3_y, Y3_z, i % NX, j % NY, k % NZ, time + 80.0 / 120.0 * DT) -
                                                               50.0 / 120.0 * DT * functionF_v(Y2_x, Y2_y, Y2_z, i % NX, j % NY, k % NZ, time + 64.0 / 120.0 * DT);
            }
        }
    }

    for (int i = 0; i < NX + 1; i++)
    {
        for (int j = 0; j < NY + 1; j++)
        {
            for (int k = 0; k < NZ; k++)
            {
                grid.w[i * (NY + 1) * NZ + j * NZ + k] = Y3_z[i * (NY + 1) * NZ + j * NZ + k] +
                                                         90.0 / 120.0 * DT * functionF_w(Y3_x, Y3_y, Y3_z, i % NX, j % NY, k % NZ, time + 80.0 / 120.0 * DT) -
                                                         50.0 / 120.0 * DT * functionF_w(Y2_x, Y2_y, Y2_z, i % NX, j % NY, k % NZ, time + 64.0 / 120.0 * DT);
            }
        }
    }

    miss = 0;
    for (int i = 0; i < NX; i++)
    {
        for (int j = 0; j < NY + 1; j++)
        {
            for (int k = 0; k < NZ + 1; k++)
            {
                if (j == 0 && k == 0)
                {
                    if (grid.u[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k] != grid.u[i * (NY + 1) * (NZ + 1) + NY * (NZ + 1) + NZ])
                    {
                        miss++;
                    }
                }
            }
        }
    }

    std::cout << "Miss Final: " << miss << std::endl;
}

Real IcoNS::functionF_u(const std::array<Real, NX *(NY + 1) * (NZ + 1)> &u, const std::array<Real, (NX + 1) * NY *(NZ + 1)> &v, const std::array<Real, (NX + 1) * (NY + 1) * NZ> &w, int i, int j, int k, Real t)
{
    int lu = ((i) % (NX)) * (NY + 1) * (NZ + 1) + (j) % (NY) * (NZ + 1) + (k) % (NZ);
    int lv = ((i) % (NX)) * NY * (NZ + 1) + (j) % (NY) * (NZ + 1) + (k) % (NZ);
    int lw = ((i) % (NX)) * (NY + 1) * NZ + (j) % (NY)*NZ + (k) % (NZ);

    return -(u[lu] * (u[((i + 1 + NX) % (NX)) * (NY + 1) * (NZ + 1) + (j) % (NY) * (NZ + 1) + (k) % (NZ)] - u[((i - 1+ NX) % (NX) ) * (NY + 1) * (NZ + 1) + (j) % (NY) * (NZ + 1) + (k) % (NZ)]) / (2.0 * DX) +
             (v[lv] + v[((i + 1 + NX) % (NX)) * NY * (NZ + 1) + (j) % (NY) * (NZ + 1) + (k) % (NZ)] + v[(i) % (NX)*NY * (NZ + 1) + ((j - 1+ NY) % (NY) ) * (NZ + 1) + (k) % (NZ)] + v[((i + 1 + NX) % (NX)) * NY * (NZ + 1) + ((j - 1 + NY) % (NY)) * (NZ + 1) + (k) % (NZ)]) / 4.0 * (u[(i) % (NX) * (NY + 1) * (NZ + 1) + ((j + 1+ NY) % (NY) ) * (NZ + 1) + (k) % (NZ)] - u[(i) % (NX) * (NY + 1) * (NZ + 1) + ((j - 1+ NY) % (NY) ) * (NZ + 1) + (k) % (NZ)]) / (2.0 * DY) +
             (w[lw] + w[((i + 1 + NX) % (NX)) * (NY + 1) * NZ + (j) % (NY)*NZ + (k) % (NZ)] + w[(i) % (NX) * (NY + 1) * NZ + (j) % (NY)*NZ + ((k - 1+ NZ) % (NZ) )] + w[((i + 1 + NX) % (NX)) * (NY + 1) * NZ + (j) % (NY)*NZ + ((k - 1+ NZ) % (NZ) )]) / 4.0 * (u[(i) % (NX) * (NY + 1) * (NZ + 1) + (j) % (NY) * (NZ + 1) + ((k + 1+ NZ) % (NZ) )] - u[(i) % (NX) * (NY + 1) * (NZ + 1) + (j) % (NY) * (NZ + 1) + ((k - 1+ NZ) % (NZ) )]) / (2.0 * DZ)) +
           1 / RE * ((u[((i + 1 + NX) % (NX)) * (NY + 1) * (NZ + 1) + (j) % (NY) * (NZ + 1) + (k) % (NZ)] - 2 * u[lu] + u[((i - 1+ NX) % (NX) ) * (NY + 1) * (NZ + 1) + (j) % (NY) * (NZ + 1) + (k) % (NZ)]) / (DX * DX) + (u[(i) % (NX) * (NY + 1) * (NZ + 1) + ((j + 1+ NY) % (NY) ) * (NZ + 1) + (k) % (NZ)] - 2 * u[lu] + u[(i) % (NX) * (NY + 1) * (NZ + 1) + ((j - 1+ NY) % (NY) ) * (NZ + 1) + (k) % (NZ)]) / (DY * DY) + (u[(i) % (NX) * (NY + 1) * (NZ + 1) + (j) % (NY) * (NZ + 1) + ((k + 1+ NZ) % (NZ) )] - 2 * u[lu] + u[(i) % (NX) * (NY + 1) * (NZ + 1) + (j) % (NY) * (NZ + 1) + ((k - 1+ NZ) % (NZ) )]) / (DZ * DZ)) +
           functionG_u(i, j, k, t);
}

Real IcoNS::functionF_v(const std::array<Real, NX *(NY + 1) * (NZ + 1)> &u, const std::array<Real, (NX + 1) * NY *(NZ + 1)> &v, const std::array<Real, (NX + 1) * (NY + 1) * NZ> &w, int i, int j, int k, Real t)
{
    int lu = (i) % (NX) * (NY + 1) * (NZ + 1) + (j) % (NY) * (NZ + 1) + (k) % (NZ);
    int lv = (i) % (NX)*NY * (NZ + 1) + (j) % (NY) * (NZ + 1) + (k) % (NZ);
    int lw = (i) % (NX) * (NY + 1) * NZ + (j) % (NY)*NZ + (k) % (NZ);


    return -((u[lu] + u[(i) % (NX) * (NY + 1) * (NZ + 1) + ((j + 1+ NY) % (NY) ) * (NZ + 1) + (k) % (NZ)] + u[((i - 1+ NX) % (NX) ) * (NY + 1) * (NZ + 1) + (j) % (NY) * (NZ + 1) + (k) % (NZ)] + u[((i - 1+ NX) % (NX) ) * (NY + 1) * (NZ + 1) + ((j + 1+ NY) % (NY) ) * (NZ + 1) + (k) % (NZ)]) / 4.0 * ((v[((i + 1+ NX) % (NX) ) * NY * (NZ + 1) + (j) % (NY) * (NZ + 1) + (k) % (NZ)] - v[((i - 1+ NX) % (NX) ) * NY * (NZ + 1) + (j) % (NY) * (NZ + 1) + (k) % (NZ)]) / (2.0 * DX)) +
             v[lv] * (v[(i) % (NX)*NY * (NZ + 1) + ((j + 1+ NY) % (NY) ) * (NZ + 1) + (k) % (NZ)] - v[(i) % (NX)*NY * (NZ + 1) + ((j - 1+ NY) % (NY) ) * (NZ + 1) + (k) % (NZ)]) / (2.0 * DY) +
             (w[lw] + w[(i) % (NX) * (NY + 1) * NZ + (j) % (NY)*NZ + ((k - 1+ NZ) % (NZ) )] + w[(i) % (NX) * (NY + 1) * NZ + ((j + 1+ NY) % (NY) ) * NZ + (k) % (NZ)] + w[(i) % (NX) * (NY + 1) * NZ + ((j + 1+ NY) % (NY) ) * NZ + ((k - 1+ NZ) % (NZ) )]) / 4.0 * (v[(i) % (NX)*NY * (NZ + 1) + (j) % (NY) * (NZ + 1) + ((k + 1+ NZ) % (NZ) )] - v[(i) % (NX)*NY * (NZ + 1) + (j) % (NY) * (NZ + 1) + ((k - 1+ NZ) % (NZ) )]) / (2.0 * DZ)) +
           (1.0 / RE) * ((v[((i + 1+ NX) % (NX) ) * NY * (NZ + 1) + (j) % (NY) * (NZ + 1) + (k) % (NZ)] - 2.0 * v[lv] + v[((i - 1+ NX) % (NX) ) * NY * (NZ + 1) + (j) % (NY) * (NZ + 1) + (k) % (NZ)]) / (DX * DX) +
                         (v[(i) % (NX)*NY * (NZ + 1) + ((j + 1 + NY) % (NY)) * (NZ + 1) + (k) % (NZ)] - 2.0 * v[lv] + v[(i) % (NX)*NY * (NZ + 1) + ((j - 1+ NY) % (NY) ) * (NZ + 1) + (k) % (NZ)]) / (DY * DY) +
                         (v[(i) % (NX)*NY * (NZ + 1) + (j) % (NY) * (NZ + 1) + ((k + 1+ NZ) % (NZ) )] - 2.0 * v[lv] + v[(i) % (NX)*NY * (NZ + 1) + (j) % (NY) * (NZ + 1) + ((k - 1+ NZ) % (NZ) )]) / (DZ * DZ)) +
           functionG_v(i, j, k, t);
}

Real IcoNS::functionF_w(const std::array<Real, NX *(NY + 1) * (NZ + 1)> &u, const std::array<Real, (NX + 1) * NY *(NZ + 1)> &v, const std::array<Real, (NX + 1) * (NY + 1) * NZ> &w, int i, int j, int k, Real t)

{
    int lu = (i) % (NX) * (NY + 1) * (NZ + 1) + (j) % (NY) * (NZ + 1) + (k) % (NZ);
    int lv = (i) % (NX)*NY * (NZ + 1) + (j) % (NY) * (NZ + 1) + (k) % (NZ);
    int lw = (i) % (NX) * (NY + 1) * NZ + (j) % (NY)*NZ + (k) % (NZ);
    

    return -((u[lu] + u[((i - 1+ NX) % (NX) ) * (NY + 1) * (NZ + 1) + (j) % (NY) * (NZ + 1) + (k) % (NZ)] + u[(i) % (NX) * (NY + 1) * (NZ + 1) + (j) % (NY) * (NZ + 1) + ((k + 1 + NZ) % (NZ))] + u[((i - 1+ NX) % (NX) ) * (NY + 1) * (NZ + 1) + (j) % (NY) * (NZ + 1) + ((k + 1+ NZ) % (NZ) )]) / 4.0 * (w[((i + 1+ NX) % (NX) ) * (NY + 1) * NZ + (j) % (NY)*NZ + (k) % (NZ)] - w[((i - 1+ NX) % (NX) ) * (NY + 1) * NZ + (j) % (NY)*NZ + (k) % (NZ)]) / (2.0 * DX) +
             (v[(i) % (NX)*NY * (NZ + 1) + (j) % (NY) * (NZ + 1) + ((k + 1 + NZ) % (NZ))] + v[(i) % (NX)*NY * (NZ + 1) + ((j - 1+ NY) % (NY) ) * (NZ + 1) + ((k + 1 + NZ) % (NZ))] + v[lv] + v[(i) % (NX)*NY * (NZ + 1) + ((j - 1+ NY) % (NY) ) * (NZ + 1) + (k) % (NZ)]) / 4.0 * (w[(i) % (NX) * (NY + 1) * NZ + ((j + 1+ NY) % (NY) ) * NZ + (k) % (NZ)] - w[(i) % (NX) * (NY + 1) * NZ + ((j - 1+ NY) % (NY) ) * NZ + (k) % (NZ)]) / (2.0 * DY) +
        w[lw] * (w[(i) % (NX) * (NY + 1) * NZ + (j) % (NY)*NZ + ((k + 1+ NZ) % (NZ) )] - w[(i) % (NX) * (NY + 1) * NZ + (j) % (NY)*NZ + ((k - 1+ NZ) % (NZ) )]) / (2.0 * DZ)) +
        (1.0 / RE) * ((w[((i + 1+NX) % (NX)) * (NY + 1) * NZ + (j) % (NY)*NZ + (k) % (NZ)] - 2.0 * w[lw] + w[((i - 1+ NX) % (NX) ) * (NY + 1) * NZ + (j) % (NY)*NZ + (k) % (NZ)]) / (DX * DX) +
                      (w[(i) % (NX) * (NY + 1) * NZ + ((j + 1+ NY) % (NY) ) * NZ + (k) % (NZ)] - 2.0 * w[lw] + w[(i) % (NX) * (NY + 1) * NZ + ((j - 1+ NY) % (NY) ) * NZ + (k) % (NZ)]) / (DY * DY) +
                      (w[(i) % (NX) * (NY + 1) * NZ + (j) % (NY)*NZ + ((k + 1+ NZ) % (NZ) )] - 2.0 * w[lw] + w[(i) % (NX) * (NY + 1) * NZ + (j) % (NY)*NZ + ((k - 1+ NZ) % (NZ) )]) / (DZ * DZ)) +
        functionG_w(i, j, k, t);
}

//w[((i + 1) % (NX) + NX) * (NY + 1) * NZ + (j) % (NY)*NZ + (k) % (NZ)] - 2.0 * w[lw] + w[((i - 1) % (NX) + NX) * (NY + 1) * NZ + (j) % (NY)*NZ + (k) % (NZ)]) / (DX * DX);
//v[((i + 1) % (NX) + NX) * NY * (NZ + 1) + (j) % (NY) * (NZ + 1) + (k) % (NZ)] - 2.0 * v[lv] + v[((i - 1) % (NX) + NX) * NY * (NZ + 1) + (j) % (NY) * (NZ + 1) + (k) % (NZ)]) / (DX * DX);

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
    return std::sin(x) * std::cos(y) * std::sin(z) * std::cos(t) +
           std::sin(x) * std::cos(x) * std::cos(y) * std::cos(y) * std::sin(z) * std::sin(z) * std::sin(t) * std::sin(t) -
           std::sin(x) * std::cos(x) * std::sin(y) * std::sin(y) * std::sin(z) * std::sin(z) * std::sin(t) * std::sin(t) + 2.0 * std::sin(x) * std::cos(x) * std::cos(y) * std::cos(y) * std::cos(z) * std::cos(z) * std::sin(t) * std::sin(t) +
           3.0 / RE * std::sin(x) * std::cos(y) * std::sin(z) * std::sin(t);
    //return std::sin(z)*std::cos(t) + std::sin(y)*std::sin(t)*std::cos(z)*std::sin(t)+1.0 / RE *std::sin(z)*std::sin(t);
    //return std::cos(t);
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
           3.0 / RE * std::cos(x) * std::sin(y) * std::sin(z) * std::sin(t);
    //return std::sin(x)*std::cos(t) + std::sin(z)*std::sin(t)*std::cos(x)*std::sin(t)+1.0 / RE *std::sin(x)*std::sin(t);
    //return std::cos(t);
}

Real IcoNS::functionG_w(int i, int j, int k, Real t)
{
    Real x = i * DX;
    Real y = j * DY;
    Real z = k * DZ + DZ / 2;
    return 2.0 * std::cos(x) * std::cos(y) * std::cos(z) * std::cos(t) - 
           2.0 * std::sin(x) * std::sin(x) * std::cos(y) * std::cos(y) * std::sin(z) * std::cos(z) * std::sin(t) * std::sin(t) -
           2.0 * std::cos(x) * std::cos(x) * std::sin(y) * std::sin(y) * std::sin(z) * std::cos(z) * std::sin(t) * std::sin(t) -
           4.0 * std::cos(x) * std::cos(x) * std::cos(y) * std::cos(y) * std::sin(z) * std::cos(z) * std::sin(t) * std::sin(t) + 
           6.0 / RE * std::cos(x) * std::cos(y) * std::cos(z) * std::sin(t);
    //return std::sin(y)*std::cos(t) + std::sin(x)*std::sin(t)*std::cos(y)*std::sin(t)+1.0 / RE *std::sin(y)*std::sin(t);
    //return std::cos(t);
}
