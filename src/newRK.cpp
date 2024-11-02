#include "core.hpp"

void IcoNS::solve_time_step(double time)
{
    std::vector<double> Y2_x(nx * (ny+1) * (nz+1));
    std::vector<double> Y2_y((nx+1) * ny * (nz+1));
    std::vector<double> Y2_z((nx+1) * (ny+1) * nz);
    std::vector<double> Y3_x(nx * (ny+1) * (nz+1));
    std::vector<double> Y3_y((nx+1) * ny * (nz+1));
    std::vector<double> Y3_z((nx+1) * (ny+1) * nz);

    for (size_t i = 1; i < nx - 1; i++)
    {
        for (size_t j = 1; j < ny; j++)
        {
            for (size_t k = 1; k < nz; k++)
            {
                Y2_x[i * (ny+1) * (nz+1) + j * (nz+1) + k] = grid.u[i * (ny+1) * (nz+1) + j * (nz+1) + k] +
                            64.0 / 120.0 * dt * functionF_u(grid.u, grid.v, grid.w, i, j, k, time);
            }
        }
    }

    for (size_t i = 1; i < nx; i++)
    {
        for (size_t j = 1; j < ny - 1; j++)
        {
            for (size_t k = 1; k < nz; k++)
            {
                Y2_y[i * ny * (nz+1) + j * (nz+1) + k] = grid.v[i * ny * (nz+1) + j * (nz+1) + k] +
                            64.0 / 120.0 * dt * functionF_v(grid.u, grid.v, grid.w, i, j, k, time);
            }
        }
    }

    for (size_t i = 1; i < nx; i++)
    {
        for (size_t j = 1; j < ny; j++)
        {
            for (size_t k = 1; k < nz - 1; k++)
            {
                Y2_z[i * (ny+1) * nz + j * nz + k] = grid.w[i * (ny+1) * nz + j * nz + k] +
                            64.0 / 120.0 * dt * functionF_w(grid.u, grid.v, grid.w, i, j, k, time);
            }
        }
    }

    boundary.update_boundary(Y2_x, Y2_y, Y2_z, time + 64.0 / 120.0 * dt);

    for (size_t i = 1; i < nx - 1; i++)
    {
        for (size_t j = 1; j < ny; j++)
        {
            for (size_t k = 1; k < nz; k++)
            {
                Y3_x[i * (ny+1) * (nz+1) + j * (nz+1) + k] = Y2_x[i * (ny+1) * (nz+1) + j * (nz+1) + k] +
                    50.0 / 120.0 * dt * functionF_u(Y2_x, Y2_y, Y2_z, i, j, k, time + 64.0 / 120.0 * dt) -
                    34.0 / 120.0 * dt * functionF_u(grid.u, grid.v, grid.w, i, j, k, time);
            }
        }
    }

    for (size_t i = 1; i < nx; i++)
    {
        for (size_t j = 1; j < ny - 1; j++)
        {
            for (size_t k = 1; k < nz; k++)
            {
                Y3_y[i * ny * (nz+1) + j * (nz+1) + k] = Y2_y[i * ny * (nz+1) + j * (nz+1) + k] +
                    50.0 / 120.0 * dt * functionF_v(Y2_x, Y2_y, Y2_z, i, j, k, time + 64.0 / 120.0 * dt) -
                    34.0 / 120.0 * dt * functionF_v(grid.u, grid.v, grid.w, i, j, k, time);
            }
        }
    }

    for (size_t i = 1; i < nx; i++)
    {
        for (size_t j = 1; j < ny; j++)
        {
            for (size_t k = 1; k < nz - 1; k++)
            {
                Y3_z[i * (ny+1) * nz + j * nz + k] = Y2_z[i * (ny+1) * nz + j * nz + k] +
                    50.0 / 120.0 * dt * functionF_w(Y2_x, Y2_y, Y2_z, i, j, k, time + 64.0 / 120.0 * dt) -
                    34.0 / 120.0 * dt * functionF_w(grid.u, grid.v, grid.w, i, j, k, time);
            }
        }
    }

    boundary.update_boundary(Y3_x, Y3_y, Y3_z, time + 80.0 / 120.0 * dt);

    for (size_t i = 1; i < nx - 1; i++)
    {
        for (size_t j = 1; j < ny; j++)
        {
            for (size_t k = 1; k < nz; k++)
            {
                grid.u[i * (ny+1) * (nz+1) + j * (nz+1) + k] = Y3_x[i * (ny+1) * (nz+1) + j * (nz+1) + k] +
                    90.0 / 120.0 * dt * functionF_u(Y3_x, Y3_y, Y3_z, i, j, k, time + 80.0 / 120.0 * dt) -
                    50.0 / 120.0 * dt * functionF_u(Y2_x, Y2_y, Y2_z, i, j, k, time + 64.0 / 120.0 * dt);
            }
        }
    }

    for (size_t i = 1; i < nx; i++)
    {
        for (size_t j = 1; j < ny - 1; j++)
        {
            for (size_t k = 1; k < nz; k++)
            {
                grid.v[i * ny * (nz+1) + j * (nz+1) + k] = Y3_y[i * ny * (nz+1) + j * (nz+1) + k] +
                    90.0 / 120.0 * dt * functionF_v(Y3_x, Y3_y, Y3_z, i, j, k, time + 80.0 / 120.0 * dt) -
                    50.0 / 120.0 * dt * functionF_v(Y2_x, Y2_y, Y2_z, i, j, k, time + 64.0 / 120.0 * dt);
            }
        }
    }

    for (size_t i = 1; i < nx; i++)
    {
        for (size_t j = 1; j < ny; j++)
        {
            for (size_t k = 1; k < nz - 1; k++)
            {
                grid.w[i * (ny+1) * nz + j * nz + k] = Y3_z[i * (ny+1) * nz + j * nz + k] +
                    90.0 / 120.0 * dt * functionF_w(Y3_x, Y3_y, Y3_z, i, j, k, time + 80.0 / 120.0 * dt) -
                    50.0 / 120.0 * dt * functionF_w(Y2_x, Y2_y, Y2_z, i, j, k, time + 64.0 / 120.0 * dt);
            }
        }
    }
}

double IcoNS::functionF_u(const std::vector<double> &u, const std::vector<double> &v,
                                     const std::vector<double> &w, size_t i, size_t j, size_t k, double t)
{
    size_t lu = i * (ny+1) * (nz+1) + j * (nz+1) + k;
    size_t lv = i * ny * (nz+1) + j * (nz+1) + k;
    size_t lw = i * (ny+1) * nz + j * nz + k;

    return -(u[lu] * (u[lu + (ny+1) * (nz+1)] - u[lu - (ny+1) * (nz+1)]) / (2.0 * dx) +
            (v[lv] + v[lv + ny * (nz+1)] + v[lv - (nz+1)] + v[lv + ny * (nz+1) - (nz+1)]) / 4.0 * (u[lu + (nz+1)] - u[lu - (nz+1)]) / (2.0 * dy) +
            (w[lw] + w[lw + (ny+1) * nz] + w[lw - 1] + w[lw + (ny+1) * nz - 1]) / 4.0 * (u[lu + 1] - u[lu - 1]) / (2.0 * dz)) +
            1 / Re * ((u[lu + (ny+1) * (nz+1)] - 2 * u[lu] + u[lu - (ny+1) * (nz+1)]) / (dx * dx) +
                (u[lu + (nz+1)] - 2 * u[lu] + u[lu - (nz+1)]) / (dy * dy) +
                (u[lu + 1] - 2 * u[lu] + u[lu - 1]) / (dz * dz)) +
            functionG_u(i, j, k, t);
}

double IcoNS::functionF_v(const std::vector<double> &u, const std::vector<double> &v,
                                     const std::vector<double> &w, size_t i, size_t j, size_t k, double t)
{
    size_t lu = i * (ny+1) * (nz+1) + j * (nz+1) + k;
    size_t lv = i * ny * (nz+1) + j * (nz+1) + k;
    size_t lw = i * (ny+1) * nz + j * nz + k;

    return -((u[lu] + u[lu + (nz+1)] + u[lu - (ny+1) * (nz+1)] + u[lu - (ny+1) * (nz+1) + (nz+1)]) / 4.0 * (v[lv + ny * (nz+1)] - v[lv - ny * (nz+1)] / (2.0 * dx)) +
            v[lv] * (v[lv + (nz+1)] - v[lv - (nz+1)]) / (2.0 * dy) +
            (w[lw] + w[lw - 1] + w[lw + (ny+1)] + w[lw + (ny+1) - 1]) / 4.0 * /* */ (v[lv + 1] - v[lv - 1]) / (2.0 * dz)) +
            (1.0 / Re) * ((v[lv + ny * (nz + 1)] - 2.0 * v[lv] + v[lv - ny * (nz + 1)]) / (dx * dx) +
                (v[lv + (nz+1)] - 2.0 * v[lv] + v[lv - (nz+1)]) / (dy * dy) +
                (v[lv + 1] - 2.0 * v[lv] + v[lv - 1]) / (dz * dz)) +
            functionG_v(i, j, k, t);
}

double IcoNS::functionF_w(const std::vector<double> &u, const std::vector<double> &v,
                                     const std::vector<double> &w, size_t i, size_t j, size_t k, double t)
{
    size_t lu = i * (ny+1) * (nz+1) + j * (nz+1) + k;
    size_t lv = i * ny * (nz+1) + j * (nz+1) + k;
    size_t lw = i * (ny+1) * nz + j * nz + k;

    return -((u[lu] + u[lu - (ny+1) * (nz+1)] + u[lu + 1] + u[lu - (nz+1) * (ny+1) + 1]) / 4.0 * (w[lw + (ny+1) * nz] - w[lw - (ny+1) * nz]) / (2.0 * dx) +
            (v[lv + 1] + v[lv - (nz+1) + 1] + v[lv] + v[lv - (nz+1)]) / 4.0 * (w[lw + nz] - w[lw - nz]) / (2.0 * dy) +
            w[lw] * (w[lw + 1] - w[lw - 1]) / (2.0 * dz)) +
            (1.0 / Re) * ((w[lw + (ny+1) * nz] - 2.0 * w[lw] + w[lw - (ny+1) * nz]) / (dx * dx) +
                (w[lw + nz] - 2.0 * w[lw] + w[lw - nz]) / (dy * dy) +
                (w[lw + 1] - 2.0 * w[lw] + w[lw - 1]) / (dz * dz)) +
            functionG_w(i, j, k, t);
}


double IcoNS::functionG_u(size_t i, size_t j, size_t k, double t)
{
    double x = i * dx + dx / 2;
    double y = j * dy;
    double z = k * dz;
    return std::sin(x) * std::cos(y) * std::sin(z) * std::cos(t) + std::sin(x) * std::cos(x) * std::cos(y) * std::cos(y) * std::sin(z) * std::sin(z) * std::sin(t) * std::sin(t) - std::sin(x) * std::cos(x) * std::sin(y) * std::sin(y) * std::sin(z) * std::sin(z) * std::sin(t) * std::sin(t) + 2 * std::sin(x) * std::cos(x) * std::cos(y) * std::cos(y) * std::cos(z) * std::cos(z) * std::sin(t) * std::sin(t) + 3.0 / Re * std::sin(x) * std::cos(y) * std::sin(z) * std::sin(t);
}

double IcoNS::functionG_v(size_t i, size_t j, size_t k, double t)
{
    double x = i * dx;
    double y = j * dy + dy / 2;
    double z = k * dz;
    return std::cos(x) * std::sin(y) * std::sin(z) * std::cos(t) - 
           std::sin(x) * std::sin(x) * std::sin(y) * std::cos(y) * std::sin(z) * std::sin(z) * std::sin(t) * std::sin(t) +
           std::cos(x) * std::cos(x) * std::sin(y) * std::cos(y) * std::sin(z) * std::sin(z) * std::sin(t) * std::sin(t) +
           2.0 * std::cos(x) * std::cos(x) * std::sin(y) * std::cos(y) * std::cos(z) * std::cos(z) * std::sin(t) * std::sin(t) +
           3.0 / Re * std::cos(x) * std::sin(y) * std::sin(z) * std::sin(t);
}

    
double IcoNS::functionG_w(size_t i, size_t j, size_t k, double t)
{
    double x = i * dx;
    double y = j * dy;
    double z = k * dz + dz / 2;
    return 2 * std::cos(x) * std::cos(y) * std::cos(z) * std::cos(t) - 2 * std::sin(x) * std::sin(x) * std::cos(y) * std::cos(y) * std::sin(z) * std::cos(z) * std::sin(t) * std::sin(t) -
           2 * std::cos(x) * std::cos(x) * std::sin(y) * std::sin(y) * std::sin(z) * std::cos(z) * std::sin(t) * std::sin(t) -
           4.0 * std::cos(x) * std::cos(x) * std::cos(y) * std::cos(y) * std::sin(z) * std::cos(z) * std::sin(t) * std::sin(t) + 6.0 / Re * std::cos(x) * std::cos(y) * std::cos(z) * std::sin(t);
}


