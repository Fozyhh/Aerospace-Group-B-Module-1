#include "core.hpp"

void IcoNS::solve_time_step(double time)
{
    for (size_t i = 1; i < NX - 1; i++)
    {
        for (size_t j = 1; j < NY; j++)
        {
            for (size_t k = 1; k < NZ; k++)
            {
                Y2_x[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k] = grid.u[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k] +
                                                                   64.0 / 120.0 * DT * functionF_u(grid.u, grid.v, grid.w, i, j, k, time);
            }
        }
    }

    for (size_t i = 1; i < NX; i++)
    {
        for (size_t j = 1; j < NY - 1; j++)
        {
            for (size_t k = 1; k < NZ; k++)
            {
                Y2_y[i * NY * (NZ + 1) + j * (NZ + 1) + k] = grid.v[i * NY * (NZ + 1) + j * (NZ + 1) + k] +
                                                             64.0 / 120.0 * DT * functionF_v(grid.u, grid.v, grid.w, i, j, k, time);
            }
        }
    }

    for (size_t i = 1; i < NX; i++)
    {
        for (size_t j = 1; j < NY; j++)
        {
            for (size_t k = 1; k < NZ - 1; k++)
            {
                Y2_z[i * (NY + 1) * NZ + j * NZ + k] = grid.w[i * (NY + 1) * NZ + j * NZ + k] +
                                                       64.0 / 120.0 * DT * functionF_w(grid.u, grid.v, grid.w, i, j, k, time);
            }
        }
    }

    boundary.update_boundary(Y2_x, Y2_y, Y2_z, time + 64.0 / 120.0 * DT);

    for (size_t i = 1; i < NX - 1; i++)
    {
        for (size_t j = 1; j < NY; j++)
        {
            for (size_t k = 1; k < NZ; k++)
            {
                Y3_x[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k] = Y2_x[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k] +
                                                                   50.0 / 120.0 * DT * functionF_u(Y2_x, Y2_y, Y2_z, i, j, k, time + 64.0 / 120.0 * DT) -
                                                                   34.0 / 120.0 * DT * functionF_u(grid.u, grid.v, grid.w, i, j, k, time);
            }
        }
    }

    for (size_t i = 1; i < NX; i++)
    {
        for (size_t j = 1; j < NY - 1; j++)
        {
            for (size_t k = 1; k < NZ; k++)
            {
                Y3_y[i * NY * (NZ + 1) + j * (NZ + 1) + k] = Y2_y[i * NY * (NZ + 1) + j * (NZ + 1) + k] +
                                                             50.0 / 120.0 * DT * functionF_v(Y2_x, Y2_y, Y2_z, i, j, k, time + 64.0 / 120.0 * DT) -
                                                             34.0 / 120.0 * DT * functionF_v(grid.u, grid.v, grid.w, i, j, k, time);
            }
        }
    }

    for (size_t i = 1; i < NX; i++)
    {
        for (size_t j = 1; j < NY; j++)
        {
            for (size_t k = 1; k < NZ - 1; k++)
            {
                Y3_z[i * (NY + 1) * NZ + j * NZ + k] = Y2_z[i * (NY + 1) * NZ + j * NZ + k] +
                                                       50.0 / 120.0 * DT * functionF_w(Y2_x, Y2_y, Y2_z, i, j, k, time + 64.0 / 120.0 * DT) -
                                                       34.0 / 120.0 * DT * functionF_w(grid.u, grid.v, grid.w, i, j, k, time);
            }
        }
    }

    boundary.update_boundary(Y3_x, Y3_y, Y3_z, time + 80.0 / 120.0 * DT);

    for (size_t i = 1; i < NX - 1; i++)
    {
        for (size_t j = 1; j < NY; j++)
        {
            for (size_t k = 1; k < NZ; k++)
            {
                grid.u[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k] = Y3_x[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k] +
                                                                     90.0 / 120.0 * DT * functionF_u(Y3_x, Y3_y, Y3_z, i, j, k, time + 80.0 / 120.0 * DT) -
                                                                     50.0 / 120.0 * DT * functionF_u(Y2_x, Y2_y, Y2_z, i, j, k, time + 64.0 / 120.0 * DT);
            }
        }
    }

    for (size_t i = 1; i < NX; i++)
    {
        for (size_t j = 1; j < NY - 1; j++)
        {
            for (size_t k = 1; k < NZ; k++)
            {
                grid.v[i * NY * (NZ + 1) + j * (NZ + 1) + k] = Y3_y[i * NY * (NZ + 1) + j * (NZ + 1) + k] +
                                                               90.0 / 120.0 * DT * functionF_v(Y3_x, Y3_y, Y3_z, i, j, k, time + 80.0 / 120.0 * DT) -
                                                               50.0 / 120.0 * DT * functionF_v(Y2_x, Y2_y, Y2_z, i, j, k, time + 64.0 / 120.0 * DT);
            }
        }
    }

    for (size_t i = 1; i < NX; i++)
    {
        for (size_t j = 1; j < NY; j++)
        {
            for (size_t k = 1; k < NZ - 1; k++)
            {
                grid.w[i * (NY + 1) * NZ + j * NZ + k] = Y3_z[i * (NY + 1) * NZ + j * NZ + k] +
                                                         90.0 / 120.0 * DT * functionF_w(Y3_x, Y3_y, Y3_z, i, j, k, time + 80.0 / 120.0 * DT) -
                                                         50.0 / 120.0 * DT * functionF_w(Y2_x, Y2_y, Y2_z, i, j, k, time + 64.0 / 120.0 * DT);
            }
        }
    }
}

double IcoNS::functionF_u(const std::array<double, NX *(NY + 1) * (NZ + 1)> &u, const std::array<double, (NX + 1) * NY *(NZ + 1)> &v, const std::array<double, (NX + 1) * (NY + 1) * NZ> &w, size_t i, size_t j, size_t k, double t)
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

double IcoNS::functionF_v(const std::array<double, NX *(NY + 1) * (NZ + 1)> &u, const std::array<double, (NX + 1) * NY *(NZ + 1)> &v, const std::array<double, (NX + 1) * (NY + 1) * NZ> &w, size_t i, size_t j, size_t k, double t)
{
    size_t lu = i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k;
    size_t lv = i * NY * (NZ + 1) + j * (NZ + 1) + k;
    size_t lw = i * (NY + 1) * NZ + j * NZ + k;

    return -((u[lu] + u[lu + (NZ + 1)] + u[lu - (NY + 1) * (NZ + 1)] + u[lu - (NY + 1) * (NZ + 1) + (NZ + 1)]) / 4.0 * ((v[lv + NY * (NZ + 1)] - v[lv - NY * (NZ + 1)]) / (2.0 * DX)) +
             v[lv] * (v[lv + (NZ + 1)] - v[lv - (NZ + 1)]) / (2.0 * DY) +
             (w[lw] + w[lw - 1] + w[lw + (NY + 1)] + w[lw + (NY + 1) - 1]) / 4.0 * /* */ (v[lv + 1] - v[lv - 1]) / (2.0 * DZ)) +
           (1.0 / RE) * ((v[lv + NY * (NZ + 1)] - 2.0 * v[lv] + v[lv - NY * (NZ + 1)]) / (DX * DX) +
                         (v[lv + (NZ + 1)] - 2.0 * v[lv] + v[lv - (NZ + 1)]) / (DY * DY) +
                         (v[lv + 1] - 2.0 * v[lv] + v[lv - 1]) / (DZ * DZ)) +
           functionG_v(i, j, k, t);
}

double IcoNS::functionF_w(const std::array<double, NX *(NY + 1) * (NZ + 1)> &u, const std::array<double, (NX + 1) * NY *(NZ + 1)> &v, const std::array<double, (NX + 1) * (NY + 1) * NZ> &w, size_t i, size_t j, size_t k, double t)
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
}

double IcoNS::functionG_u(size_t i, size_t j, size_t k, double t)
{
    double x = i * DX + DX / 2;
    double y = j * DY;
    double z = k * DZ;
    return std::sin(x) * std::cos(y) * std::sin(z) * std::cos(t) +
           std::sin(x) * std::cos(x) * std::cos(y) * std::cos(y) * std::sin(z) * std::sin(z) * std::sin(t) * std::sin(t) -
           std::sin(x) * std::cos(x) * std::sin(y) * std::sin(y) * std::sin(z) * std::sin(z) * std::sin(t) * std::sin(t) + 2 * std::sin(x) * std::cos(x) * std::cos(y) * std::cos(y) * std::cos(z) * std::cos(z) * std::sin(t) * std::sin(t) +
           3.0 / RE * std::sin(x) * std::cos(y) * std::sin(z) * std::sin(t);
}

double IcoNS::functionG_v(size_t i, size_t j, size_t k, double t)
{
    double x = i * DX;
    double y = j * DY + DY / 2;
    double z = k * DZ;
    return std::cos(x) * std::sin(y) * std::sin(z) * std::cos(t) -
           std::sin(x) * std::sin(x) * std::sin(y) * std::cos(y) * std::sin(z) * std::sin(z) * std::sin(t) * std::sin(t) +
           std::cos(x) * std::cos(x) * std::sin(y) * std::cos(y) * std::sin(z) * std::sin(z) * std::sin(t) * std::sin(t) +
           2.0 * std::cos(x) * std::cos(x) * std::sin(y) * std::cos(y) * std::cos(z) * std::cos(z) * std::sin(t) * std::sin(t) +
           3.0 / RE * std::cos(x) * std::sin(y) * std::sin(z) * std::sin(t);
}

double IcoNS::functionG_w(size_t i, size_t j, size_t k, double t)
{
    double x = i * DX;
    double y = j * DY;
    double z = k * DZ + DZ / 2;
    return 2 * std::cos(x) * std::cos(y) * std::cos(z) * std::cos(t) - 2 * std::sin(x) * std::sin(x) * std::cos(y) * std::cos(y) * std::sin(z) * std::cos(z) * std::sin(t) * std::sin(t) -
           2 * std::cos(x) * std::cos(x) * std::sin(y) * std::sin(y) * std::sin(z) * std::cos(z) * std::sin(t) * std::sin(t) -
           4.0 * std::cos(x) * std::cos(x) * std::cos(y) * std::cos(y) * std::sin(z) * std::cos(z) * std::sin(t) * std::sin(t) + 6.0 / RE * std::cos(x) * std::cos(y) * std::cos(z) * std::sin(t);
}
