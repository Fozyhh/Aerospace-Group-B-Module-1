#include "core.hpp"


void IcoNS::solve_time_step(Real time)
{
    for (int i = 1 + lbx; i < newDimX_x - 1 - rbx; i++)
    {
        for (int j = 1 + lby; j < newDimY_x - 1 - rby; j++)
        {
            for (int k = lbz; k < dim_z - rbz; k++)
            {
                Y2_x[i * newDimY_x * dim_z + j * dim_z + k] = grid_loc_x[i * newDimY_x * dim_z + j * dim_z + k] +
                                                              64.0 / 120.0 * DT * functionF_u(grid_loc_x, grid_loc_y, grid_loc_z, i, j, k, time);
            }
        }
    }

    for (int i = 1 + lbx; i < newDimX_y - 1 - rbx; i++)
    {
        for (int j = 1 + lby; j < newDimY_y - 1 - rby; j++)
        {
            for (int k = lbz; k < dim_z - rbz; k++)
            {
                Y2_y[i * newDimY_y * dim_z + j * dim_z + k] = grid_loc_y[i * newDimY_y * dim_z + j * dim_z + k] +
                                                                    64.0 / 120.0 * DT * functionF_v(grid_loc_x, grid_loc_y, grid_loc_z, i, j, k, time);
            }
        }
    }

    for (int i = 1 + lbx; i < newDimX_z - 1 - rbx; i++)
    {
        for (int j = 1 + lby; j < newDimY_z - 1 - rby; j++)
        {
            for (int k = 1; k < dim_z_z - 1; k++)
            {
                Y2_z[i * newDimY_z * dim_z_z + j * dim_z_z + k] = grid_loc_z[i * newDimY_z * dim_z_z + j * dim_z_z + k] +
                                                        64.0 / 120.0 * DT * functionF_w(grid_loc_x, grid_loc_y, grid_loc_z, i, j, k, time);
            }
        }
    }

    boundary.update_boundary(Y2_x, Y2_y, Y2_z, time + 64.0 / 120.0 * DT);
    MPI_Barrier(cart_comm);
    exchangeData(Y2_x, newDimX_x, newDimY_x, dim_z, MPI_face_x_x, MPI_face_y_x);
    exchangeData(Y2_y, newDimX_y, newDimY_y, dim_z, MPI_face_x_y, MPI_face_y_y);
    exchangeData(Y2_z, newDimX_z, newDimY_z, dim_z_z, MPI_face_x_z, MPI_face_y_z);
    for (int i = 1 + lbx; i < newDimX_x - 1 - rbx; i++)
    {
        for (int j = 1 + lby; j < newDimY_x - 1 - rby; j++)
        {
            for (int k = lbz; k < dim_z - rbz; k++)
            {
                Y3_x[i * newDimY_x * dim_z + j * dim_z + k] = Y2_x[i * newDimY_x * dim_z + j * dim_z + k] +
                                                                    50.0 / 120.0 * DT * functionF_u(Y2_x, Y2_y, Y2_z, i, j, k, time + 64.0 / 120.0 * DT) -
                                                                    34.0 / 120.0 * DT * functionF_u(grid_loc_x, grid_loc_y, grid_loc_z, i, j, k, time);
            }
        }
    }

    for (int i = 1 + lbx; i < newDimX_y - 1 - rbx; i++)
    {
        for (int j = 1 + lby; j < newDimY_y - 1 - rby; j++)
        {
            for (int k = lbz; k < dim_z - rbz; k++)
            {
                Y3_y[i * newDimY_y * dim_z + j * dim_z + k] = Y2_y[i * newDimY_y * dim_z + j * dim_z + k] +
                                                                    50.0 / 120.0 * DT * functionF_v(Y2_x, Y2_y, Y2_z, i, j, k, time + 64.0 / 120.0 * DT) -
                                                                    34.0 / 120.0 * DT * functionF_v(grid_loc_x, grid_loc_y, grid_loc_z, i, j, k, time);
            }
        }
    }

    for (int i = 1 + lbx; i < newDimX_z - 1 - rbx; i++)
    {
        for (int j = 1 + lby; j < newDimY_z - 1 - rby; j++)
        {
            for (int k = lbz; k < dim_z_z - rbz; k++)
            {
                Y3_z[i * newDimY_z * dim_z_z + j * dim_z_z + k] = Y2_z[i * newDimY_z * dim_z_z + j * dim_z_z + k] +
                                                        50.0 / 120.0 * DT * functionF_w(Y2_x, Y2_y, Y2_z, i, j, k, time + 64.0 / 120.0 * DT) -
                                                        34.0 / 120.0 * DT * functionF_w(grid_loc_x, grid_loc_y, grid_loc_z, i, j, k, time);
            }
        }
    }

    boundary.update_boundary(Y3_x, Y3_y, Y3_z, time + 80.0 / 120.0 * DT);
    MPI_Barrier(cart_comm);
    exchangeData(Y3_x, newDimX_x, newDimY_x, dim_z, MPI_face_x_x, MPI_face_y_x);
    exchangeData(Y3_y, newDimX_y, newDimY_y, dim_z, MPI_face_x_y, MPI_face_y_y);
    exchangeData(Y3_z, newDimX_z, newDimY_z, dim_z_z, MPI_face_x_z, MPI_face_y_z);
    MPI_Barrier(cart_comm);
    for (int i = 1 + lbx; i < newDimX_x - 1 - rbx; i++)
    {
        for (int j = 1 + lby; j < newDimY_x - 1 - rby; j++)
        {
            for (int k = lbz; k < dim_z- rbz; k++)
            {
                grid_loc_x[i * newDimY_x * dim_z + j * dim_z + k] = Y3_x[i * newDimY_x * dim_z + j * dim_z + k] +
                                                                          90.0 / 120.0 * DT * functionF_u(Y3_x, Y3_y, Y3_z, i, j, k, time + 80.0 / 120.0 * DT) -
                                                                          50.0 / 120.0 * DT * functionF_u(Y2_x, Y2_y, Y2_z, i, j, k, time + 64.0 / 120.0 * DT);
            }
        }
    }

    for (int i = 1 + lbx; i < newDimX_y - 1 - rbx; i++)
    {
        for (int j = 1 + lby; j < newDimY_y - 1 - rby; j++)
        {
            for (int k = lbz; k < dim_z-rbz; k++)
            {
                grid_loc_y[i * newDimY_y * dim_z + j * dim_z + k] = Y3_y[i * newDimY_y * dim_z + j * dim_z + k] +
                                                                          90.0 / 120.0 * DT * functionF_v(Y3_x, Y3_y, Y3_z, i, j, k, time + 80.0 / 120.0 * DT) -
                                                                          50.0 / 120.0 * DT * functionF_v(Y2_x, Y2_y, Y2_z, i, j, k, time + 64.0 / 120.0 * DT);
            }
        }
    }

    for (int i = 1 + lbx; i < newDimX_z - 1 - rbx; i++)
    {
        for (int j = 1 + lby; j < newDimY_z - 1 - rby; j++)
        {
            for (int k = lbz; k < dim_z_z - rbz; k++)
            {
                grid_loc_z[i * newDimY_z * dim_z_z + j * dim_z_z + k] = Y3_z[i * newDimY_z * dim_z_z + j * dim_z_z + k] +
                                                              90.0 / 120.0 * DT * functionF_w(Y3_x, Y3_y, Y3_z, i, j, k, time + 80.0 / 120.0 * DT) -
                                                              50.0 / 120.0 * DT * functionF_w(Y2_x, Y2_y, Y2_z, i, j, k, time + 64.0 / 120.0 * DT);
            }
        }
    }
}

Real IcoNS::functionF_u(const std::vector<Real> &u, const std::vector<Real> &v, const std::vector<Real> &w, int i, int j, int k, Real t)
{
    int lu = i * newDimY_x * dim_z + j * dim_z + k;
    int lv = i * newDimY_y * dim_z + j * dim_z + k;
    int lw = i * newDimY_z * dim_z_z + j * dim_z_z + k;

    return -(u[lu] * (u[lu + newDimY_x * dim_z] - u[lu - newDimY_x * dim_z]) / (2.0 * DX) +
             (v[lv] + v[lv + newDimY_y * dim_z] + v[lv - dim_z] + v[lv + newDimY_y * dim_z - dim_z]) / 4.0 * (u[lu + dim_z] - u[lu - dim_z]) / (2.0 * DY) +
             (w[lw] + w[lw + newDimY_z * dim_z_z] + w[i * newDimY_z * dim_z_z + j * dim_z_z + ((k-1)%dim_z_z + dim_z_z) % dim_z_z/*lw - 1*/] + w[(i+1) * newDimY_z * dim_z_z + j * dim_z_z + ((k-1)%dim_z_z + dim_z_z) % dim_z_z/*lw + newDimY_z * dim_z_z - 1*/]) / 4.0 * (u[i * newDimY_x * dim_z + j * dim_z + ((k+1)%dim_z)/*lu + 1*/] - u[i * newDimY_x * dim_z + j * dim_z + ((k-1)%dim_z + dim_z) % dim_z/*lu - 1*/]) / (2.0 * DZ)) +
           1 / RE * ((u[lu + newDimY_x * dim_z] - 2 * u[lu] + u[lu - newDimY_x * dim_z]) / (DX * DX) + (u[lu + dim_z] - 2 * u[lu] + u[lu - dim_z]) / (DY * DY) + (u[i * newDimY_x * dim_z + j * dim_z + ((k+1)%dim_z)/*lu + 1*/] - 2 * u[lu] + u[i * newDimY_x * dim_z + j * dim_z + ((k-1)%dim_z + dim_z) % dim_z/*lu - 1*/]) / (DZ * DZ)) +
           functionG_u(i-1 + coords[0] * other_dim_x_x, j-1 + coords[1] * other_dim_y_x, k, t);
}

Real IcoNS::functionF_v(const std::vector<Real> &u, const std::vector<Real> &v, const std::vector<Real> &w, int i, int j, int k, Real t)
{
    int lu = i * newDimY_x * dim_z + j * dim_z + k;
    int lv = i * newDimY_y * dim_z + j * dim_z + k;
    int lw = i * newDimY_z * dim_z_z + j * dim_z_z + k;
    
    return -((u[lu] + u[lu + dim_z] + u[lu - newDimY_x * dim_z] + u[lu - newDimY_x * dim_z + dim_z]) / 4.0 * ((v[lv + newDimY_y * dim_z] - v[lv - newDimY_y * dim_z]) / (2.0 * DX)) +
             v[lv] * (v[lv + dim_z] - v[lv - dim_z]) / (2.0 * DY) +
             (w[lw] + w[i * newDimY_z * dim_z_z + j * dim_z_z + ((k-1)%dim_z_z + dim_z_z) % dim_z_z/*lw - 1*/] + w[lw + dim_z_z] + w[i * newDimY_z * dim_z_z + (j+1) * dim_z_z + ((k-1)%dim_z_z + dim_z_z) % dim_z_z/*lw + dim_z_z - 1*/]) / 4.0 * (v[i * newDimY_y * dim_z + j * dim_z + ((k+1)%dim_z)/*lv + 1*/] - v[i * newDimY_y * dim_z + j * dim_z + ((k-1)%dim_z + dim_z) % dim_z/*lv - 1*/]) / (2.0 * DZ)) +
           (1.0 / RE) * ((v[lv + newDimY_y * dim_z] - 2.0 * v[lv] + v[lv - newDimY_y * dim_z]) / (DX * DX) +
                         (v[lv + dim_z] - 2.0 * v[lv] + v[lv - dim_z]) / (DY * DY) +
                         (v[i * newDimY_y * dim_z + j * dim_z + ((k+1)%dim_z)/*lv + 1*/] - 2.0 * v[lv] + v[i * newDimY_y * dim_z + j * dim_z + ((k-1)%dim_z + dim_z) % dim_z/*lv - 1*/]) / (DZ * DZ)) +
           functionG_v(i-1 + coords[0] * other_dim_x_y, j-1 + coords[1] * other_dim_y_y, k, t);
}

Real IcoNS::functionF_w(const std::vector<Real> &u, const std::vector<Real> &v, const std::vector<Real> &w, int i, int j, int k, Real t)

{
    int lu = i * newDimY_x * dim_z + j * dim_z + k;
    int lv = i * newDimY_y * dim_z + j * dim_z + k;
    int lw = i * newDimY_z * dim_z_z + j * dim_z_z + k;

    return -((u[lu] + u[lu - newDimY_x * dim_z] + u[i * newDimY_x * dim_z + j * dim_z + ((k+1)%dim_z)/*lu + 1*/] + u[(i-1) * newDimY_x * dim_z + j * dim_z + ((k+1)%dim_z)/*lu - dim_z * newDimY_x + 1*/]) / 4.0 * (w[lw + newDimY_z * dim_z_z] - w[lw - newDimY_z * dim_z_z]) / (2.0 * DX) +
             (v[i * newDimY_y * dim_z + j * dim_z + ((k+1)%dim_z)/*lv + 1*/] + v[i * newDimY_y * dim_z + (j-1) * dim_z + ((k+1)%dim_z)/*lv - dim_z + 1*/] + v[lv] + v[lv - dim_z]) / 4.0 * (w[lw + dim_z_z] - w[lw - dim_z_z]) / (2.0 * DY) +
             w[lw] * (w[i * newDimY_z * dim_z_z + j * dim_z_z + ((k+1)%dim_z_z)/*lw + 1*/] - w[i * newDimY_z * dim_z_z + j * dim_z_z + ((k-1)%dim_z_z + dim_z_z) % dim_z_z/*lw - 1*/]) / (2.0 * DZ)) +
           (1.0 / RE) * ((w[lw + newDimY_z * dim_z_z] - 2.0 * w[lw] + w[lw - newDimY_z * dim_z_z]) / (DX * DX) +
                         (w[lw + dim_z_z] - 2.0 * w[lw] + w[lw - dim_z_z]) / (DY * DY) +
                         (w[i * newDimY_z * dim_z_z + j * dim_z_z + ((k+1)%dim_z_z)/*lw + 1*/] - 2.0 * w[lw] + w[i * newDimY_z * dim_z_z + j * dim_z_z + ((k-1)%dim_z_z + dim_z_z) % dim_z_z/*lw - 1*/]) / (DZ * DZ)) +
           functionG_w(i-1 + coords[0] * other_dim_x_z, j-1 + coords[1] * other_dim_y_z, k, t);
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

    return 2 * std::cos(x) * std::cos(y) * std::cos(z) * std::cos(t) -
           2 * std::sin(x) * std::sin(x) * std::cos(y) * std::cos(y) * std::sin(z) * std::cos(z) * std::sin(t) * std::sin(t) -
           2 * std::cos(x) * std::cos(x) * std::sin(y) * std::sin(y) * std::sin(z) * std::cos(z) * std::sin(t) * std::sin(t) -
           4.0 * std::cos(x) * std::cos(x) * std::cos(y) * std::cos(y) * std::sin(z) * std::cos(z) * std::sin(t) * std::sin(t) +
           6.0 / RE * std::cos(x) * std::cos(y) * std::cos(z) * std::sin(t);
}


