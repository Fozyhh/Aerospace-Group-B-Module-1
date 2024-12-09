#include "core.hpp"

//TODO: PARALLELIZATION SHOULD BE COMPLETE

void IcoNS::solve_time_step(Real time)
{
    for (size_t i = 1 + lbx; i < newDimX_x - 1 - rbx; i++)

    {
        for (size_t j = 1 + lby; j < newDimY_x - 1 - rby; j++)
        {
            for (size_t k = 1; k < NZ; k++)
            {
                Y2_x[i * newDimY_x * dim_z + j * dim_z + k] = grid_loc_x[i * newDimY_x * (NZ + 1) + j * (NZ + 1) + k] +
                                                                   64.0 / 120.0 * DT * functionF_u(grid_loc_x, grid_loc_y, grid_loc_z, i, j, k, time);
            }
        }
    }

    for (size_t i = 1 + lbx; i < newDimX_y - 1 - rbx; i++)
    {
        for (size_t j = 1 + lbx; j < newDimY_y - 1 - rby; j++)
        {
            for (size_t k = 1; k < NZ; k++)
            {
                Y2_y[i * newDimY_y * (NZ + 1) + j * (NZ + 1) + k] = grid_loc_y[i * newDimY_y * (NZ + 1) + j * (NZ + 1) + k] +
                                                             64.0 / 120.0 * DT * functionF_v(grid_loc_x, grid_loc_y, grid_loc_z, i, j, k, time);
            }
        }
    }

    for (size_t i = 1 + lbx; i < newDimX_z - 1 - rbx; i++)
    {
        for (size_t j = 1 + lby; j < newDimY_z -1 - rby; j++)
        {
            for (size_t k = 1; k < NZ - 1; k++)
            {
                Y2_z[i * newDimY_z * NZ + j * NZ + k] = grid_loc_z[i * newDimY_z * NZ + j * NZ + k] +
                                                       64.0 / 120.0 * DT * functionF_w(grid_loc_x, grid_loc_y, grid_loc_z, i, j, k, time);
            }
        }
    }

    boundary.update_boundary(Y2_x, Y2_y, Y2_z, time + 64.0 / 120.0 * DT);

    for (size_t i = 1 + lbx; i < newDimX_x - 1 - rbx; i++)
    {
        for (size_t j = 1 + lby; j < newDimY_x - 1 - rby; j++)
        {
            for (size_t k = 1; k < NZ; k++)
            {
                Y3_x[i * newDimY_x * (NZ + 1) + j * (NZ + 1) + k] = Y2_x[i * newDimY_x * (NZ + 1) + j * (NZ + 1) + k] +
                                                                   50.0 / 120.0 * DT * functionF_u(Y2_x, Y2_y, Y2_z, i, j, k, time + 64.0 / 120.0 * DT) -
                                                                   34.0 / 120.0 * DT * functionF_u(grid_loc_x, grid_loc_y, grid_loc_z, i, j, k, time);
            }
        }
    }

    for (size_t i = 1 + lbx; i < newDimX_y - 1 - rbx; i++)
    {
        for (size_t j = 1 + lbx; j < newDimY_y - 1 - rby; j++)
        {
            for (size_t k = 1; k < NZ; k++)
            {
                Y3_y[i * newDimY_y * (NZ + 1) + j * (NZ + 1) + k] = Y2_y[i * newDimY_y * (NZ + 1) + j * (NZ + 1) + k] +
                                                             50.0 / 120.0 * DT * functionF_v(Y2_x, Y2_y, Y2_z, i, j, k, time + 64.0 / 120.0 * DT) -
                                                             34.0 / 120.0 * DT * functionF_v(grid_loc_x, grid_loc_y, grid_loc_z, i, j, k, time);
            }
        }
    }

    for (size_t i = 1 + lbx; i < newDimX_z - 1 - rbx; i++)
    {
        for (size_t j = 1 + lby; j < newDimY_z -1 - rby; j++)
        {
            for (size_t k = 1; k < NZ - 1; k++)
            {
                Y3_z[i * newDimY_z * NZ + j * NZ + k] = Y2_z[i * newDimY_z * NZ + j * NZ + k] +
                                                       50.0 / 120.0 * DT * functionF_w(Y2_x, Y2_y, Y2_z, i, j, k, time + 64.0 / 120.0 * DT) -
                                                       34.0 / 120.0 * DT * functionF_w(grid_loc_x, grid_loc_y, grid_loc_z, i, j, k, time);
            }
        }
    }

    boundary.update_boundary(Y3_x, Y3_y, Y3_z, time + 80.0 / 120.0 * DT);

    for (size_t i = 1 + lbx; i < newDimX_x - 1 - rbx; i++)
    {
        for (size_t j = 1 + lby; j < newDimY_x - 1 - rby; j++)
        {
            for (size_t k = 1; k < NZ; k++)
            {
                grid_loc_x[i * newDimY_x * (NZ + 1) + j * (NZ + 1) + k] = Y3_x[i * newDimY_x * (NZ + 1) + j * (NZ + 1) + k] +
                                                                     90.0 / 120.0 * DT * functionF_u(Y3_x, Y3_y, Y3_z, i, j, k, time + 80.0 / 120.0 * DT) -
                                                                     50.0 / 120.0 * DT * functionF_u(Y2_x, Y2_y, Y2_z, i, j, k, time + 64.0 / 120.0 * DT);
            }
        }
    }

    for (size_t i = 1 + lbx; i < newDimX_y - 1 - rbx; i++)
    {
        for (size_t j = 1 + lbx; j < newDimY_y - 1 - rby; j++)
        {
            for (size_t k = 1; k < NZ; k++)
            {
                grid_loc_y[i * newDimY_y * (NZ + 1) + j * (NZ + 1) + k] = Y3_y[i * newDimY_y * (NZ + 1) + j * (NZ + 1) + k] +
                                                               90.0 / 120.0 * DT * functionF_v(Y3_x, Y3_y, Y3_z, i, j, k, time + 80.0 / 120.0 * DT) -
                                                               50.0 / 120.0 * DT * functionF_v(Y2_x, Y2_y, Y2_z, i, j, k, time + 64.0 / 120.0 * DT);
            }
        }
    }

    for (size_t i = 1 + lbx; i < newDimX_z - 1 - rbx; i++)
    {
        for (size_t j = 1 + lby; j < newDimY_z -1 - rby; j++)
        {
            for (size_t k = 1; k < NZ - 1; k++)
            {
                grid_loc_z[i * newDimY_z * NZ + j * NZ + k] = Y3_z[i * newDimY_z * NZ + j * NZ + k] +
                                                         90.0 / 120.0 * DT * functionF_w(Y3_x, Y3_y, Y3_z, i, j, k, time + 80.0 / 120.0 * DT) -
                                                         50.0 / 120.0 * DT * functionF_w(Y2_x, Y2_y, Y2_z, i, j, k, time + 64.0 / 120.0 * DT);
            }
        }
    }
    
}


Real IcoNS::functionF_u(const std::array<Real, newDimX_x*newDimY_x*dim_z> &u, const std::array<Real, newDimX_x*newDimY_x*dim_z> &v, const std::array<Real, newDimX_x*newDimY_x*dim_z_z> &w, size_t i, size_t j, size_t k, Real t)
{
    size_t lu = i * newDimY_x * (NZ + 1) + j * (NZ + 1) + k;
    size_t lv = i * newDimY_y * (NZ + 1) + j * (NZ + 1) + k;
    size_t lw = i * newDimY_z * NZ + j * NZ + k;

    return -(u[lu] * (u[lu + newDimY_x * (NZ + 1)] - u[lu - newDimY_x * (NZ + 1)]) / (2.0 * DX) +
             (v[lv] + v[lv + newDimY_y * (NZ + 1)] + v[lv - (NZ + 1)] + v[lv + newDimY_y * (NZ + 1) - (NZ + 1)]) / 4.0 * (u[lu + (NZ + 1)] - u[lu - (NZ + 1)]) / (2.0 * DY) +
             (w[lw] + w[lw + newDimY_z * NZ] + w[lw - 1] + w[lw + newDimY_z * NZ - 1]) / 4.0 * (u[lu + 1] - u[lu - 1]) / (2.0 * DZ)) +
           1 / RE * ((u[lu + newDimY_x * (NZ + 1)] - 2 * u[lu] + u[lu - newDimY_x * (NZ + 1)]) / (DX * DX) + (u[lu + (NZ + 1)] - 2 * u[lu] + u[lu - (NZ + 1)]) / (DY * DY) + (u[lu + 1] - 2 * u[lu] + u[lu - 1]) / (DZ * DZ)) +
           functionG_u(i, j, k, t);
}

Real IcoNS::functionF_v(const std::array<Real, newDimX_x*newDimY_x*dim_z> &u, const std::array<Real, newDimX_x*newDimY_x*dim_z> &v, const std::array<Real, newDimX_x*newDimY_x*dim_z_z> &w, size_t i, size_t j, size_t k, Real t)
{
    size_t lu = i * newDimY_x * (NZ + 1) + j * (NZ + 1) + k;
    size_t lv = i * newDimY_y * (NZ + 1) + j * (NZ + 1) + k;
    size_t lw = i * newDimY_z * NZ + j * NZ + k;

    return -((u[lu] + u[lu + (NZ + 1)] + u[lu - newDimY_x * (NZ + 1)] + u[lu - newDimY_x * (NZ + 1) + (NZ + 1)]) / 4.0 * ((v[lv + newDimY_y * (NZ + 1)] - v[lv - newDimY_y * (NZ + 1)]) / (2.0 * DX)) +
             v[lv] * (v[lv + (NZ + 1)] - v[lv - (NZ + 1)]) / (2.0 * DY) +
             (w[lw] + w[lw - 1] + w[lw + newDimY_z] + w[lw + newDimY_z - 1]) / 4.0 * /* */ (v[lv + 1] - v[lv - 1]) / (2.0 * DZ)) +
           (1.0 / RE) * ((v[lv + newDimY_y * (NZ + 1)] - 2.0 * v[lv] + v[lv - newDimY_y * (NZ + 1)]) / (DX * DX) +
                         (v[lv + (NZ + 1)] - 2.0 * v[lv] + v[lv - (NZ + 1)]) / (DY * DY) +
                         (v[lv + 1] - 2.0 * v[lv] + v[lv - 1]) / (DZ * DZ)) +
           functionG_v(i, j, k, t);
}

Real IcoNS::functionF_w(const std::array<Real, newDimX_x*newDimY_x*dim_z> &u, const std::array<Real, newDimX_x*newDimY_x*dim_z> &v, const std::array<Real, newDimX_x*newDimY_x*dim_z_z> &w, size_t i, size_t j, size_t k, Real t)

{
    size_t lu = i * newDimY_x * (NZ + 1) + j * (NZ + 1) + k;
    size_t lv = i * newDimY_y * (NZ + 1) + j * (NZ + 1) + k;
    size_t lw = i * newDimY_z * NZ + j * NZ + k;

    return -((u[lu] + u[lu - newDimY_x * (NZ + 1)] + u[lu + 1] + u[lu - (NZ + 1) * newDimY_x + 1]) / 4.0 * (w[lw + newDimY_x * NZ] - w[lw - newDimY_x * NZ]) / (2.0 * DX) +
             (v[lv + 1] + v[lv - (NZ + 1) + 1] + v[lv] + v[lv - (NZ + 1)]) / 4.0 * (w[lw + NZ] - w[lw - NZ]) / (2.0 * DY) +
             w[lw] * (w[lw + 1] - w[lw - 1]) / (2.0 * DZ)) +
           (1.0 / RE) * ((w[lw + newDimY_z * NZ] - 2.0 * w[lw] + w[lw - newDimY_z * NZ]) / (DX * DX) +
                         (w[lw + NZ] - 2.0 * w[lw] + w[lw - NZ]) / (DY * DY) +
                         (w[lw + 1] - 2.0 * w[lw] + w[lw - 1]) / (DZ * DZ)) +
           functionG_w(i, j, k, t);
}


Real IcoNS::functionG_u(size_t i, size_t j, size_t k, Real t)
{
    Real x = i * DX + DX / 2;
    Real y = j * DY;
    Real z = k * DZ;
    return std::sin(x) * std::cos(y) * std::sin(z) * std::cos(t) +
           std::sin(x) * std::cos(x) * std::cos(y) * std::cos(y) * std::sin(z) * std::sin(z) * std::sin(t) * std::sin(t) -
           std::sin(x) * std::cos(x) * std::sin(y) * std::sin(y) * std::sin(z) * std::sin(z) * std::sin(t) * std::sin(t) + 2 * std::sin(x) * std::cos(x) * std::cos(y) * std::cos(y) * std::cos(z) * std::cos(z) * std::sin(t) * std::sin(t) +
           3.0 / RE * std::sin(x) * std::cos(y) * std::sin(z) * std::sin(t);

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
           3.0 / RE * std::cos(x) * std::sin(y) * std::sin(z) * std::sin(t);
}


Real IcoNS::functionG_w(size_t i, size_t j, size_t k, Real t)
{
    Real x = i * DX;
    Real y = j * DY;
    Real z = k * DZ + DZ / 2;
    return 2 * std::cos(x) * std::cos(y) * std::cos(z) * std::cos(t) - 2 * std::sin(x) * std::sin(x) * std::cos(y) * std::cos(y) * std::sin(z) * std::cos(z) * std::sin(t) * std::sin(t) -
           2 * std::cos(x) * std::cos(x) * std::sin(y) * std::sin(y) * std::sin(z) * std::cos(z) * std::sin(t) * std::sin(t) -
           4.0 * std::cos(x) * std::cos(x) * std::cos(y) * std::cos(y) * std::sin(z) * std::cos(z) * std::sin(t) * std::sin(t) + 6.0 / RE * std::cos(x) * std::cos(y) * std::cos(z) * std::sin(t);
}
