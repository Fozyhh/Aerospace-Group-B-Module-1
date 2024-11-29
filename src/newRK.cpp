#include "core.hpp"

void IcoNS::solve_time_step(Real time)
{
    for (int i = 1; i < NX - 1; i++)

    {
        for (int j = 1; j < NY; j++)
        {
            for (int k = 1; k < NZ; k++)
            {
                Y2_x[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k] = grid.u[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k] +
                                                                   64.0 / 120.0 * DT * functionF_u(grid.u, grid.v, grid.w, i, j, k, time);
            }
        }
    }

    for (int i = 1; i < NX; i++)
    {
        for (int j = 1; j < NY - 1; j++)
        {
            for (int k = 1; k < NZ; k++)
            {
                Y2_y[i * NY * (NZ + 1) + j * (NZ + 1) + k] = grid.v[i * NY * (NZ + 1) + j * (NZ + 1) + k] +
                                                             64.0 / 120.0 * DT * functionF_v(grid.u, grid.v, grid.w, i, j, k, time);
            }
        }
    }

    for (int i = 1; i < NX; i++)
    {
        for (int j = 1; j < NY; j++)
        {
            for (int k = 1; k < NZ - 1; k++)
            {
                Y2_z[i * (NY + 1) * NZ + j * NZ + k] = grid.w[i * (NY + 1) * NZ + j * NZ + k] +
                                                       64.0 / 120.0 * DT * functionF_w(grid.u, grid.v, grid.w, i, j, k, time);
            }
        }
    }

    boundary.update_boundary(Y2_x, Y2_y, Y2_z, time + 64.0 / 120.0 * DT);

    for (int i = 1; i < NX - 1; i++)

    {
        for (int j = 1; j < NY; j++)
        {
            for (int k = 1; k < NZ; k++)
            {
                Y3_x[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k] = Y2_x[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k] +
                                                                   50.0 / 120.0 * DT * functionF_u(Y2_x, Y2_y, Y2_z, i, j, k, time + 64.0 / 120.0 * DT) -
                                                                   34.0 / 120.0 * DT * functionF_u(grid.u, grid.v, grid.w, i, j, k, time);
            }
        }
    }

    for (int i = 1; i < NX; i++)
    {
        for (int j = 1; j < NY - 1; j++)
        {
            for (int k = 1; k < NZ; k++)
            {
                Y3_y[i * NY * (NZ + 1) + j * (NZ + 1) + k] = Y2_y[i * NY * (NZ + 1) + j * (NZ + 1) + k] +
                                                             50.0 / 120.0 * DT * functionF_v(Y2_x, Y2_y, Y2_z, i, j, k, time + 64.0 / 120.0 * DT) -
                                                             34.0 / 120.0 * DT * functionF_v(grid.u, grid.v, grid.w, i, j, k, time);
            }
        }
    }

    for (int i = 1; i < NX; i++)
    {
        for (int j = 1; j < NY; j++)
        {
            for (int k = 1; k < NZ - 1; k++)
            {
                Y3_z[i * (NY + 1) * NZ + j * NZ + k] = Y2_z[i * (NY + 1) * NZ + j * NZ + k] +
                                                       50.0 / 120.0 * DT * functionF_w(Y2_x, Y2_y, Y2_z, i, j, k, time + 64.0 / 120.0 * DT) -
                                                       34.0 / 120.0 * DT * functionF_w(grid.u, grid.v, grid.w, i, j, k, time);
            }
        }
    }

    boundary.update_boundary(Y3_x, Y3_y, Y3_z, time + 80.0 / 120.0 * DT);

    for (int i = 1; i < NX - 1; i++)

    {
        for (int j = 1; j < NY; j++)
        {
            for (int k = 1; k < NZ; k++)
            {
                grid.u[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k] = Y3_x[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k] +
                                                                     90.0 / 120.0 * DT * functionF_u(Y3_x, Y3_y, Y3_z, i, j, k, time + 80.0 / 120.0 * DT) -
                                                                     50.0 / 120.0 * DT * functionF_u(Y2_x, Y2_y, Y2_z, i, j, k, time + 64.0 / 120.0 * DT);
            }
        }
    }

    for (int i = 1; i < NX; i++)
    {
        for (int j = 1; j < NY - 1; j++)
        {
            for (int k = 1; k < NZ; k++)
            {
                grid.v[i * NY * (NZ + 1) + j * (NZ + 1) + k] = Y3_y[i * NY * (NZ + 1) + j * (NZ + 1) + k] +
                                                               90.0 / 120.0 * DT * functionF_v(Y3_x, Y3_y, Y3_z, i, j, k, time + 80.0 / 120.0 * DT) -
                                                               50.0 / 120.0 * DT * functionF_v(Y2_x, Y2_y, Y2_z, i, j, k, time + 64.0 / 120.0 * DT);
            }
        }
    }

    for (int i = 1; i < NX; i++)
    {
        for (int j = 1; j < NY; j++)
        {
            for (int k = 1; k < NZ - 1; k++)
            {
                grid.w[i * (NY + 1) * NZ + j * NZ + k] = Y3_z[i * (NY + 1) * NZ + j * NZ + k] +
                                                         90.0 / 120.0 * DT * functionF_w(Y3_x, Y3_y, Y3_z, i, j, k, time + 80.0 / 120.0 * DT) -
                                                         50.0 / 120.0 * DT * functionF_w(Y2_x, Y2_y, Y2_z, i, j, k, time + 64.0 / 120.0 * DT);
            }
        }
    }
}

Real IcoNS::functionF_u(const std::array<Real, NX *(NY + 1) * (NZ + 1)> &u, const std::array<Real, (NX + 1) * NY *(NZ + 1)> &v, const std::array<Real, (NX + 1) * (NY + 1) * NZ> &w, int i, int j, int k, Real t)
{
    int lu = ((i)) * (NY + 1) * (NZ + 1) + (j) * (NZ + 1) + (k);
    int lv = ((i)) * NY * (NZ + 1) + (j) * (NZ + 1) + (k);
    int lw = ((i)) * (NY + 1) * NZ + (j)*NZ + (k);

    return -(u[lu] * (u[(i + 1) * (NY + 1) * (NZ + 1) + (j) * (NZ + 1) + (k)] - u[(i - 1) * (NY + 1) * (NZ + 1) + (j) * (NZ + 1) + (k)]) / (2.0 * DX) +
             (v[lv] + v[(i + 1) * NY * (NZ + 1) + (j) * (NZ + 1) + (k)] + v[(i)*NY * (NZ + 1) + (j - 1) * (NZ + 1) + (k)] + v[(i + 1) * NY * (NZ + 1) + (j - 1) * (NZ + 1) + (k)]) / 4.0 * (u[(i) * (NY + 1) * (NZ + 1) + (j + 1) * (NZ + 1) + (k)] - u[(i) * (NY + 1) * (NZ + 1) + (j - 1) * (NZ + 1) + (k)]) / (2.0 * DY) +
             (w[lw] + w[(i + 1) * (NY + 1) * NZ + (j)*NZ + (k)] + w[(i) * (NY + 1) * NZ + (j)*NZ + (k - 1)] + w[(i + 1) * (NY + 1) * NZ + (j)*NZ + (k - 1)]) / 4.0 * (u[(i) * (NY + 1) * (NZ + 1) + (j) * (NZ + 1) + (k + 1)] - u[(i) * (NY + 1) * (NZ + 1) + (j) * (NZ + 1) + (k - 1)]) / (2.0 * DZ)) +
           1 / RE * ((u[(i + 1) * (NY + 1) * (NZ + 1) + (j) * (NZ + 1) + (k)] - 2 * u[lu] + u[(i - 1) * (NY + 1) * (NZ + 1) + (j) * (NZ + 1) + (k)]) / (DX * DX) + (u[(i) * (NY + 1) * (NZ + 1) + (j + 1) * (NZ + 1) + (k)] - 2 * u[lu] + u[(i) * (NY + 1) * (NZ + 1) + (j - 1) * (NZ + 1) + (k)]) / (DY * DY) + (u[(i) * (NY + 1) * (NZ + 1) + (j) * (NZ + 1) + (k + 1)] - 2 * u[lu] + u[(i) * (NY + 1) * (NZ + 1) + (j) * (NZ + 1) + (k - 1)]) / (DZ * DZ)) +
           functionG_u(i, j, k, t);
}

Real IcoNS::functionF_v(const std::array<Real, NX *(NY + 1) * (NZ + 1)> &u, const std::array<Real, (NX + 1) * NY *(NZ + 1)> &v, const std::array<Real, (NX + 1) * (NY + 1) * NZ> &w, int i, int j, int k, Real t)
{
    int lu = (i) * (NY + 1) * (NZ + 1) + (j) * (NZ + 1) + (k);
    int lv = (i)*NY * (NZ + 1) + (j) * (NZ + 1) + (k);
    int lw = (i) * (NY + 1) * NZ + (j)*NZ + (k);

    return -((u[lu] + u[(i) * (NY + 1) * (NZ + 1) + (j + 1) * (NZ + 1) + (k)] + u[(i - 1) * (NY + 1) * (NZ + 1) + (j) * (NZ + 1) + (k)] + u[(i - 1) * (NY + 1) * (NZ + 1) + (j + 1) * (NZ + 1) + (k)]) / 4.0 * ((v[(i + 1) * NY * (NZ + 1) + (j) * (NZ + 1) + (k)] - v[(i - 1) * NY * (NZ + 1) + (j) * (NZ + 1) + (k)]) / (2.0 * DX)) +
             v[lv] * (v[(i)*NY * (NZ + 1) + (j + 1) * (NZ + 1) + (k)] - v[(i)*NY * (NZ + 1) + (j - 1) * (NZ + 1) + (k)]) / (2.0 * DY) +
             (w[lw] + w[(i) * (NY + 1) * NZ + (j)*NZ + (k - 1)] + w[(i) * (NY + 1) * NZ + (j + 1) * NZ + (k)] + w[(i) * (NY + 1) * NZ + (j + 1) * NZ + (k - 1)]) / 4.0 * (v[(i)*NY * (NZ + 1) + (j) * (NZ + 1) + (k + 1)] - v[(i)*NY * (NZ + 1) + (j) * (NZ + 1) + (k - 1)]) / (2.0 * DZ)) +
           (1.0 / RE) * ((v[(i + 1) * NY * (NZ + 1) + (j) * (NZ + 1) + (k)] - 2.0 * v[lv] + v[(i - 1) * NY * (NZ + 1) + (j) * (NZ + 1) + (k)]) / (DX * DX) +
                         (v[(i)*NY * (NZ + 1) + (j + 1) * (NZ + 1) + (k)] - 2.0 * v[lv] + v[(i)*NY * (NZ + 1) + (j - 1) * (NZ + 1) + (k)]) / (DY * DY) +
                         (v[(i)*NY * (NZ + 1) + (j) * (NZ + 1) + (k + 1)] - 2.0 * v[lv] + v[(i)*NY * (NZ + 1) + (j) * (NZ + 1) + (k - 1)]) / (DZ * DZ)) +
           functionG_v(i, j, k, t);
}

Real IcoNS::functionF_w(const std::array<Real, NX *(NY + 1) * (NZ + 1)> &u, const std::array<Real, (NX + 1) * NY *(NZ + 1)> &v, const std::array<Real, (NX + 1) * (NY + 1) * NZ> &w, int i, int j, int k, Real t)
{
    int lu = (i) * (NY + 1) * (NZ + 1) + (j) * (NZ + 1) + (k);
    int lv = (i)*NY * (NZ + 1) + (j) * (NZ + 1) + (k);
    int lw = (i) * (NY + 1) * NZ + (j)*NZ + (k);

    return -((u[lu] + u[(i - 1) * (NY + 1) * (NZ + 1) + (j) * (NZ + 1) + (k)] + u[(i) * (NY + 1) * (NZ + 1) + (j) * (NZ + 1) + ((k + 1))] + u[((i - 1)) * (NY + 1) * (NZ + 1) + (j) * (NZ + 1) + (k + 1)]) / 4.0 * (w[(i + 1) * (NY + 1) * NZ + (j)*NZ + (k)] - w[(i - 1) * (NY + 1) * NZ + (j)*NZ + (k)]) / (2.0 * DX) +
             (v[(i)*NY * (NZ + 1) + (j) * (NZ + 1) + (k + 1)] + v[(i)*NY * (NZ + 1) + (j - 1) * (NZ + 1) + ((k + 1))] + v[lv] + v[(i)*NY * (NZ + 1) + (j - 1) * (NZ + 1) + (k)]) / 4.0 * (w[(i) * (NY + 1) * NZ + (j + 1) * NZ + (k)] - w[(i) * (NY + 1) * NZ + (j - 1) * NZ + (k)]) / (2.0 * DY) +
             w[lw] * (w[(i) * (NY + 1) * NZ + (j)*NZ + (k + 1)] - w[(i) * (NY + 1) * NZ + (j)*NZ + (k - 1)]) / (2.0 * DZ)) +
           (1.0 / RE) * ((w[(i + 1) * (NY + 1) * NZ + (j)*NZ + (k)] - 2.0 * w[lw] + w[(i - 1) * (NY + 1) * NZ + (j)*NZ + (k)]) / (DX * DX) +
                         (w[(i) * (NY + 1) * NZ + (j + 1) * NZ + (k)] - 2.0 * w[lw] + w[(i) * (NY + 1) * NZ + (j - 1) * NZ + (k)]) / (DY * DY) +
                         (w[(i) * (NY + 1) * NZ + (j)*NZ + (k + 1)] - 2.0 * w[lw] + w[(i) * (NY + 1) * NZ + (j)*NZ + (k - 1)]) / (DZ * DZ)) +
           functionG_w(i, j, k, t);
}

Real IcoNS::functionG_u(int i, int j, int k, Real t)
{
    Real x = i * DX + DX / 2;
    Real y = j * DY;
    Real z = k * DZ;
    return std::sin(x) * std::cos(y) * std::sin(z) * std::cos(t) +
           std::sin(x) * std::cos(x) * std::cos(y) * std::cos(y) * std::sin(z) * std::sin(z) * std::sin(t) * std::sin(t) -
           std::sin(x) * std::cos(x) * std::sin(y) * std::sin(y) * std::sin(z) * std::sin(z) * std::sin(t) * std::sin(t) + 2 * std::sin(x) * std::cos(x) * std::cos(y) * std::cos(y) * std::cos(z) * std::cos(z) * std::sin(t) * std::sin(t) +
           3.0 / RE * std::sin(x) * std::cos(y) * std::sin(z) * std::sin(t);
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
}

Real IcoNS::functionG_w(int i, int j, int k, Real t)
{
    Real x = i * DX;
    Real y = j * DY;
    Real z = k * DZ + DZ / 2;
    return 2 * std::cos(x) * std::cos(y) * std::cos(z) * std::cos(t) - 2 * std::sin(x) * std::sin(x) * std::cos(y) * std::cos(y) * std::sin(z) * std::cos(z) * std::sin(t) * std::sin(t) -
           2 * std::cos(x) * std::cos(x) * std::sin(y) * std::sin(y) * std::sin(z) * std::cos(z) * std::sin(t) * std::sin(t) -
           4.0 * std::cos(x) * std::cos(x) * std::cos(y) * std::cos(y) * std::sin(z) * std::cos(z) * std::sin(t) * std::sin(t) + 6.0 / RE * std::cos(x) * std::cos(y) * std::cos(z) * std::sin(t);
}
