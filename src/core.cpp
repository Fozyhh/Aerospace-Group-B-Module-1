
#include <cmath>

#include "core.hpp"
#include <math.h>

#include <fstream>
#include <string>
#include <memory>

void IcoNS::preprocessing(/*std::string &input_file*/)

{
    // read the input file.
    for (size_t i = 1; i < nx - 1; i++)
    {
        for (size_t j = 1; j < ny - 1; j++)
        {
            for (size_t k = 1; k < nz - 1; k++)
            {
                grid.u[i * ny * nz + j * nz + k] = 3.0;
                grid.v[i * ny * nz + j * nz + k] = 3.0;
                grid.w[i * ny * nz + j * nz + k] = 3.0;
            }
        }
    }

    for (size_t i = 0; i < nx; i++)
    {
        for (size_t j = 0; j < ny; j++)
        {
            for (size_t k = 0; k < nz; k++)
            {
                if (i == 0 || i == nx || j == 0 || j == ny || k == 0 || k == nz)
                {
                    grid.u[i * ny * nz + j * nz + k] = 1.0;
                    grid.v[i * ny * nz + j * nz + k] = 1.0;
                    grid.w[i * ny * nz + j * nz + k] = 1.0;
                }
            }
        }
        // pressure initialization.
    }

    // boundary
    auto zero = std::make_shared<FunctionZero>();
    auto frontface_u = std::make_shared<Dirichlet>([](double x, double y, double z, double t)
                                                   { return std::sin(x * t); });

    auto frontface_v = std::make_shared<Dirichlet>([](double x, double y, double z, double t)
                                                   { return std::sin(y * t); });

    auto frontface_w = std::make_shared<Dirichlet>([](double x, double y, double z, double t)
                                                   { return std::sin(z * t); });

    std::cout << "vector building" << std::endl;
    // Order: left, right, front, back, lower, upper
    boundary.addFunction(0, zero);
    boundary.addFunction(0, zero);
    boundary.addFunction(0, frontface_u);
    boundary.addFunction(0, frontface_u);
    boundary.addFunction(0, zero);
    boundary.addFunction(0, zero);
    boundary.addFunction(1, zero);
    boundary.addFunction(1, zero);
    boundary.addFunction(1, frontface_v);
    boundary.addFunction(1, frontface_v);
    boundary.addFunction(1, zero);
    boundary.addFunction(1, zero);
    boundary.addFunction(2, zero);
    boundary.addFunction(2, zero);
    boundary.addFunction(2, frontface_w);
    boundary.addFunction(2, frontface_w);
    boundary.addFunction(2, zero);
    boundary.addFunction(2, zero);
}

void IcoNS::solve()
{
    preprocessing();
    double time = 0.0;
    int i = 0;

    std::ofstream error_log("../ressources/error.log");

    while (time < T)
    {
        // apply_boundary_conditions(time);
        boundary.update_boundary(time);
        solve_time_step(time);
        time += dt;
        output();
        //csv file w/ "," delimiter: time step, iter, L2_error
        error_log << time << "," << i << "," << " Error: " << L2_error(time) << std::endl;
        i++;
    }
}

std::vector<double> IcoNS::functionF(const std::vector<double> &u, const std::vector<double> &v,
                                     const std::vector<double> &w, size_t i, size_t j, size_t k, double t)
{
    std::vector<double> f(3);
    // size_t l = i * ny * nz + j * nz + k;   //think this is wrong, we should have one l for each grid

    size_t lu = i * ny * nz + j * nz + k;
    size_t lv = i * (ny - 1) * nz + j * nz + k;
    size_t lw = i * ny * (nz - 1) + j * (nz - 1) + k;

    std::vector<double> g(3);
    g = functionG(i, j, k, t);

    f[0] = -(u[lu] * (u[lu + nz * ny] - u[lu - nz * ny]) / (2.0 * dx) +
             (v[lv] + v[lv + (ny - 1) * nz] + v[lv - nz] + v[lv + nz * (ny - 1) - nz]) / 4.0 * (u[lu + nz] - u[lu - nz]) / (2.0 * dy) +
             (w[lw] + w[lw + ny * (nz - 1)] + w[lw - 1] + w[lw + ny * (nz - 1) - 1]) / 4.0 * (u[lu + 1] - u[lu - 1]) / (2.0 * dz)) +
           1 / Re * ((u[lu + ny * nz] - 2 * u[lu] + u[lu - ny * nz]) / (dx * dx) + (u[lu + nz] - 2 * u[lu] + u[lu - nz]) / (dy * dy) + (u[lu + 1] - 2 * u[lu] + u[lu - 1]) / (dz * dz)) + g[0];
    f[1] = -((u[lu] + u[lu + nz] + u[lu - ny * nz] + u[lu - ny * nz + nz]) / 4.0 * (v[lv + (ny - 1) * nz] - v[lv - (ny - 1) * nz] / (2.0 * dx)) +
             v[lv] * (v[lv + nz] - v[lv - nz]) / (2.0 * dy) +
             (w[lw] + w[lw - 1] + w[lw + ny] + w[lw + ny - 1]) / 4.0 + (v[lv + 1] - v[lv - 1]) / (2.0 * dz)) +
           1 / Re * ((v[lv + nz * (ny - 1)] - 2.0 * v[lv] + v[lv - nz * (ny - 1)]) / (dx * dx) + (v[lv + nz] - 2.0 * v[lv] + v[lv - nz]) / (dy * dy) + (v[lv + 1] - 2.0 * v[lv] + v[lv - 1]) / (dz * dz)) + g[1];
    f[2] = -((u[lu] + u[lu - ny * nz] + u[lu + 1] + u[lu - nz * ny + 1]) / 4.0 * (w[lw + ny * (nz - 1)] - w[lw - ny * (nz - 1)]) / (2.0 * dx) +
             (v[lv + 1] + v[lv - nz + 1] + v[lv] + v[lv - nz]) / 4.0 * (w[lw + (nz - 1)] - w[lw - (nz - 1)]) / (2.0 * dy) +
             w[lw] * (w[lw + 1] - w[lw - 1]) / (2.0 * dz)) +
           1 / Re * ((w[lw + ny * (nz - 1)] - 2.0 * w[lw] + w[lw - ny * (nz - 1)]) / (dx * dx) + (w[lw + (nz - 1)] - 2.0 * w[lw] + w[lw - (nz - 1)]) / (dy * dy) + (w[lw + 1] - 2.0 * w[lw] + w[lw - 1]) / (dz * dz)) + g[2];

    return f;
}

std::vector<double> IcoNS::functionG(size_t i, size_t j, size_t k, double t)
{
    std::vector<double> g(3);
    double x = i * dx + dx / 2;
    double y = j * dy;
    double z = k * dz;
    g[0] = std::sin(x) * std::cos(y) * std::sin(z) * std::cos(t) + std::sin(x) * std::cos(x) * std::cos(y) * std::cos(y) * std::sin(z) * std::sin(z) * std::sin(t) * std::sin(t) - std::sin(x) * std::cos(x) * std::sin(y) * std::sin(y) * std::sin(z) * std::sin(z) * std::sin(t) * std::sin(t) + 2 * std::sin(x) * std::cos(x) * std::cos(y) * std::cos(y) * std::cos(z) * std::cos(z) * std::sin(t) * std::sin(t) + 3.0 / Re * std::sin(x) * std::cos(y) * std::sin(z) * std::sin(t);
    x = i * dx;
    y = j * dy + dy / 2;
    z = k * dz;
    g[1] = std::cos(x) * std::sin(y) * std::sin(z) * std::cos(t) - std::sin(x) * std::sin(x) * std::sin(y) * std::cos(y) * std::sin(z) * std::sin(z) * std::sin(t) * std::sin(t) +
           std::cos(x) * std::cos(x) * std::sin(y) * std::cos(y) * std::sin(z) * std::sin(z) * std::sin(t) * std::sin(t) +
           2.0 * std::cos(x) * std::cos(x) * std::sin(y) * std::cos(y) * std::cos(z) * std::cos(z) * std::sin(t) * std::sin(t) +
           3.0 / Re * std::cos(x) * std::sin(y) * std::sin(z) * std::sin(t);
    x = i * dx;
    y = j * dy;
    z = k * dz + dz / 2;
    g[2] = 2 * std::cos(x) * std::cos(y) * std::cos(z) * std::cos(t) - 2 * std::sin(x) * std::sin(x) * std::cos(y) * std::cos(y) * std::sin(z) * std::cos(z) * std::sin(t) * std::sin(t) -
           2 * std::cos(x) * std::cos(x) * std::sin(y) * std::sin(y) * std::sin(z) * std::cos(z) * std::sin(t) * std::sin(t) -
           4.0 * std::cos(x) * std::cos(x) * std::cos(y) * std::cos(y) * std::sin(z) * std::cos(z) * std::sin(t) * std::sin(t) + 6.0 / Re * std::cos(x) * std::cos(y) * std::cos(z) * std::sin(t);
    return g;
}

void IcoNS::solve_time_step(double time)
{
    std::vector<double> f(3);
    std::vector<double> f_y2(3);
    std::vector<double> f_y3(3);
    std::vector<double> Y2_x((nx - 1) * ny * nz);
    std::vector<double> Y2_y(nx * (ny - 1) * nz);
    std::vector<double> Y2_z(nx * ny * (nz - 1));
    std::vector<double> Y3_x((nx - 1) * ny * nz);
    std::vector<double> Y3_y(nx * (ny - 1) * nz);
    std::vector<double> Y3_z(nx * ny * (nz - 1));

    for (size_t i = 1; i < nx - 2; i++)
    {
        for (size_t j = 1; j < ny - 2; j++)
        {
            for (size_t k = 1; k < nz - 2; k++)
            {
                f = functionF(grid.u, grid.v, grid.w, i, j, k, time);

                Y2_x[i * ny * nz + j * nz + k] = grid.u[i * ny * nz + j * nz + k] +
                                                 64.0 / 120.0 * dt * f[0];
                Y2_y[i * (ny - 1) * nz + j * nz + k] = grid.v[i * (ny - 1) * nz + j * nz + k] +
                                                       64.0 / 120.0 * dt * f[1];
                Y2_z[i * ny * (nz - 1) + j * (nz - 1) + k] = grid.w[i * ny * (nz - 1) + j * (nz - 1) + k] +
                                                             64.0 / 120.0 * dt * f[2];

                f_y2 = functionF(Y2_x, Y2_y, Y2_z, i, j, k, time + 64.0 / 120.0 * dt);

                Y3_x[i * ny * nz + j * nz + k] = Y2_x[i * ny * nz + j * nz + k] +
                                                 50.0 / 120.0 * dt * f_y2[0] -
                                                 34.0 / 120.0 * dt * f[0];
                Y3_y[i * (ny - 1) * nz + j * nz + k] = Y2_y[i * (ny - 1) * nz + j * nz + k] +
                                                       50.0 / 120.0 * dt * f_y2[1] -
                                                       34.0 / 120.0 * dt * f[1];
                Y3_z[i * ny * (nz - 1) + j * (nz - 1) + k] = Y2_z[i * ny * (nz - 1) + j * (nz - 1) + k] +
                                                             50.0 / 120.0 * dt * f_y2[2] -
                                                             34.0 / 120.0 * dt * f[2];

                f_y3 = functionF(Y3_x, Y3_y, Y3_z, i, j, k, time + 80 / 120 * dt);

                grid.u[i * ny * nz + j * nz + k] = Y3_x[i * ny * nz + j * nz + k] +
                                                   90.0 / 120.0 * dt * f_y3[0] -
                                                   50.0 / 120.0 * dt * f_y2[0];
                grid.v[i * (ny - 1) * nz + j * nz + k] = Y3_y[i * (ny - 1) * nz + j * nz + k] +
                                                         90.0 / 120.0 * dt * f_y3[1] -
                                                         50.0 / 120.0 * dt * f_y2[1];
                grid.w[i * ny * (nz - 1) + j * (nz - 1) + k] = Y3_z[i * ny * (nz - 1) + j * (nz - 1) + k] +
                                                               90.0 / 120.0 * dt * f_y3[2] -
                                                               50.0 / 120.0 * dt * f_y2[2];
            }
        }
    }

    // apply runge kutta on faces outside the common cube
    for (size_t i = 1; i < nx - 2; i++)
    {
        for (size_t k = 1; i < nz - 2; k++)
        {
            size_t j = ny - 2;
            f = functionF(grid.u, grid.v, grid.w, i, j, k, time);

            Y2_x[i * ny * nz + j * nz + k] = grid.u[i * ny * nz + j * nz + k] +
                                             64.0 / 120.0 * dt * f[0];
            Y2_y[i * (ny - 1) * nz + j * nz + k] = grid.v[i * (ny - 1) * nz + j * nz + k] +
                                                   64.0 / 120.0 * dt * f[1];
            Y2_z[i * ny * (nz - 1) + j * (nz - 1) + k] = grid.w[i * ny * (nz - 1) + j * (nz - 1) + k] +
                                                         64.0 / 120.0 * dt * f[2];

            f_y2 = functionF(Y2_x, Y2_y, Y2_z, i, j, k, time + 64.0 / 120.0 * dt);

            Y3_x[i * ny * nz + j * nz + k] = Y2_x[i * ny * nz + j * nz + k] +
                                             50.0 / 120.0 * dt * f_y2[0] -
                                             34.0 / 120.0 * dt * f[0];
            Y3_y[i * (ny - 1) * nz + j * nz + k] = Y2_y[i * (ny - 1) * nz + j * nz + k] +
                                                   50.0 / 120.0 * dt * f_y2[1] -
                                                   34.0 / 120.0 * dt * f[1];
            Y3_z[i * ny * (nz - 1) + j * (nz - 1) + k] = Y2_z[i * ny * (nz - 1) + j * (nz - 1) + k] +
                                                         50.0 / 120.0 * dt * f_y2[2] -
                                                         34.0 / 120.0 * dt * f[2];

            f_y3 = functionF(Y3_x, Y3_y, Y3_z, i, j, k, time + 80 / 120 * dt);

            grid.u[i * ny * nz + j * nz + k] = Y3_x[i * ny * nz + j * nz + k] +
                                               90.0 / 120.0 * dt * f_y3[0] -
                                               50.0 / 120.0 * dt * f_y2[0];
            grid.v[i * (ny - 1) * nz + j * nz + k] = Y3_y[i * (ny - 1) * nz + j * nz + k] +
                                                     90.0 / 120.0 * dt * f_y3[1] -
                                                     50.0 / 120.0 * dt * f_y2[1];
            grid.w[i * ny * (nz - 1) + j * (nz - 1) + k] = Y3_z[i * ny * (nz - 1) + j * (nz - 1) + k] +
                                                           90.0 / 120.0 * dt * f_y3[2] -
                                                           50.0 / 120.0 * dt * f_y2[2];
        }
    }

    for (size_t i = 1; i < nx - 2; i++)
    {
        for (size_t j = 1; j < ny - 1; j++)
        {
            size_t k = nz - 2;
            f = functionF(grid.u, grid.v, grid.w, i, j, k, time);

            Y2_x[i * ny * nz + j * nz + k] = grid.u[i * ny * nz + j * nz + k] +
                                             64.0 / 120.0 * dt * f[0];
            Y2_y[i * (ny - 1) * nz + j * nz + k] = grid.v[i * (ny - 1) * nz + j * nz + k] +
                                                   64.0 / 120.0 * dt * f[1];
            Y2_z[i * ny * (nz - 1) + j * (nz - 1) + k] = grid.w[i * ny * (nz - 1) + j * (nz - 1) + k] +
                                                         64.0 / 120.0 * dt * f[2];

            f_y2 = functionF(Y2_x, Y2_y, Y2_z, i, j, k, time + 64.0 / 120.0 * dt);

            Y3_x[i * ny * nz + j * nz + k] = Y2_x[i * ny * nz + j * nz + k] +
                                             50.0 / 120.0 * dt * f_y2[0] -
                                             34.0 / 120.0 * dt * f[0];
            Y3_y[i * (ny - 1) * nz + j * nz + k] = Y2_y[i * (ny - 1) * nz + j * nz + k] +
                                                   50.0 / 120.0 * dt * f_y2[1] -
                                                   34.0 / 120.0 * dt * f[1];
            Y3_z[i * ny * (nz - 1) + j * (nz - 1) + k] = Y2_z[i * ny * (nz - 1) + j * (nz - 1) + k] +
                                                         50.0 / 120.0 * dt * f_y2[2] -
                                                         34.0 / 120.0 * dt * f[2];

            f_y3 = functionF(Y3_x, Y3_y, Y3_z, i, j, k, time + 80 / 120 * dt);

            grid.u[i * ny * nz + j * nz + k] = Y3_x[i * ny * nz + j * nz + k] +
                                               90.0 / 120.0 * dt * f_y3[0] -
                                               50.0 / 120.0 * dt * f_y2[0];
            grid.v[i * (ny - 1) * nz + j * nz + k] = Y3_y[i * (ny - 1) * nz + j * nz + k] +
                                                     90.0 / 120.0 * dt * f_y3[1] -
                                                     50.0 / 120.0 * dt * f_y2[1];
            grid.w[i * ny * (nz - 1) + j * (nz - 1) + k] = Y3_z[i * ny * (nz - 1) + j * (nz - 1) + k] +
                                                           90.0 / 120.0 * dt * f_y3[2] -
                                                           50.0 / 120.0 * dt * f_y2[2];
        }
    }

    for (size_t j = 1; j < ny - 2; j++)
    {
        for (size_t k = 1; k < nz - 1; k++)
        {
            size_t i = nx - 2;
            f = functionF(grid.u, grid.v, grid.w, i, j, k, time);

            Y2_x[i * ny * nz + j * nz + k] = grid.u[i * ny * nz + j * nz + k] +
                                             64.0 / 120.0 * dt * f[0];
            Y2_y[i * (ny - 1) * nz + j * nz + k] = grid.v[i * (ny - 1) * nz + j * nz + k] +
                                                   64.0 / 120.0 * dt * f[1];
            Y2_z[i * ny * (nz - 1) + j * (nz - 1) + k] = grid.w[i * ny * (nz - 1) + j * (nz - 1) + k] +
                                                         64.0 / 120.0 * dt * f[2];

            f_y2 = functionF(Y2_x, Y2_y, Y2_z, i, j, k, time + 64.0 / 120.0 * dt);

            Y3_x[i * ny * nz + j * nz + k] = Y2_x[i * ny * nz + j * nz + k] +
                                             50.0 / 120.0 * dt * f_y2[0] -
                                             34.0 / 120.0 * dt * f[0];
            Y3_y[i * (ny - 1) * nz + j * nz + k] = Y2_y[i * (ny - 1) * nz + j * nz + k] +
                                                   50.0 / 120.0 * dt * f_y2[1] -
                                                   34.0 / 120.0 * dt * f[1];
            Y3_z[i * ny * (nz - 1) + j * (nz - 1) + k] = Y2_z[i * ny * (nz - 1) + j * (nz - 1) + k] +
                                                         50.0 / 120.0 * dt * f_y2[2] -
                                                         34.0 / 120.0 * dt * f[2];

            f_y3 = functionF(Y3_x, Y3_y, Y3_z, i, j, k, time + 80 / 120 * dt);

            grid.u[i * ny * nz + j * nz + k] = Y3_x[i * ny * nz + j * nz + k] +
                                               90.0 / 120.0 * dt * f_y3[0] -
                                               50.0 / 120.0 * dt * f_y2[0];
            grid.v[i * (ny - 1) * nz + j * nz + k] = Y3_y[i * (ny - 1) * nz + j * nz + k] +
                                                     90.0 / 120.0 * dt * f_y3[1] -
                                                     50.0 / 120.0 * dt * f_y2[1];
            grid.w[i * ny * (nz - 1) + j * (nz - 1) + k] = Y3_z[i * ny * (nz - 1) + j * (nz - 1) + k] +
                                                           90.0 / 120.0 * dt * f_y3[2] -
                                                           50.0 / 120.0 * dt * f_y2[2];
        }
    }

    for (size_t i = 1; i < nx - 2; i++)
    {
        for (size_t j = 1; j < ny - 2; j++)
        {
            size_t k = nz - 2;
            f = functionF(grid.u, grid.v, grid.w, i, j, k, time);

            Y2_x[i * ny * nz + j * nz + k] = grid.u[i * ny * nz + j * nz + k] +
                                             64.0 / 120.0 * dt * f[0];
            Y2_y[i * (ny - 1) * nz + j * nz + k] = grid.v[i * (ny - 1) * nz + j * nz + k] +
                                                   64.0 / 120.0 * dt * f[1];
            Y2_z[i * ny * (nz - 1) + j * (nz - 1) + k] = grid.w[i * ny * (nz - 1) + j * (nz - 1) + k] +
                                                         64.0 / 120.0 * dt * f[2];

            f_y2 = functionF(Y2_x, Y2_y, Y2_z, i, j, k, time + 64.0 / 120.0 * dt);

            Y3_x[i * ny * nz + j * nz + k] = Y2_x[i * ny * nz + j * nz + k] +
                                             50.0 / 120.0 * dt * f_y2[0] -
                                             34.0 / 120.0 * dt * f[0];
            Y3_y[i * (ny - 1) * nz + j * nz + k] = Y2_y[i * (ny - 1) * nz + j * nz + k] +
                                                   50.0 / 120.0 * dt * f_y2[1] -
                                                   34.0 / 120.0 * dt * f[1];
            Y3_z[i * ny * (nz - 1) + j * (nz - 1) + k] = Y2_z[i * ny * (nz - 1) + j * (nz - 1) + k] +
                                                         50.0 / 120.0 * dt * f_y2[2] -
                                                         34.0 / 120.0 * dt * f[2];

            f_y3 = functionF(Y3_x, Y3_y, Y3_z, i, j, k, time + 80 / 120 * dt);

            grid.u[i * ny * nz + j * nz + k] = Y3_x[i * ny * nz + j * nz + k] +
                                               90.0 / 120.0 * dt * f_y3[0] -
                                               50.0 / 120.0 * dt * f_y2[0];
            grid.v[i * (ny - 1) * nz + j * nz + k] = Y3_y[i * (ny - 1) * nz + j * nz + k] +
                                                     90.0 / 120.0 * dt * f_y3[1] -
                                                     50.0 / 120.0 * dt * f_y2[1];
            grid.w[i * ny * (nz - 1) + j * (nz - 1) + k] = Y3_z[i * ny * (nz - 1) + j * (nz - 1) + k] +
                                                           90.0 / 120.0 * dt * f_y3[2] -
                                                           50.0 / 120.0 * dt * f_y2[2];
        }
    }

    for (size_t j = 1; j < ny - 2; j++)
    {
        for (size_t k = 1; k < nz - 2; k++)
        {
            size_t i = nx - 2;
            f = functionF(grid.u, grid.v, grid.w, i, j, k, time);

            Y2_x[i * ny * nz + j * nz + k] = grid.u[i * ny * nz + j * nz + k] +
                                             64.0 / 120.0 * dt * f[0];
            Y2_y[i * (ny - 1) * nz + j * nz + k] = grid.v[i * (ny - 1) * nz + j * nz + k] +
                                                   64.0 / 120.0 * dt * f[1];
            Y2_z[i * ny * (nz - 1) + j * (nz - 1) + k] = grid.w[i * ny * (nz - 1) + j * (nz - 1) + k] +
                                                         64.0 / 120.0 * dt * f[2];

            f_y2 = functionF(Y2_x, Y2_y, Y2_z, i, j, k, time + 64.0 / 120.0 * dt);

            Y3_x[i * ny * nz + j * nz + k] = Y2_x[i * ny * nz + j * nz + k] +
                                             50.0 / 120.0 * dt * f_y2[0] -
                                             34.0 / 120.0 * dt * f[0];
            Y3_y[i * (ny - 1) * nz + j * nz + k] = Y2_y[i * (ny - 1) * nz + j * nz + k] +
                                                   50.0 / 120.0 * dt * f_y2[1] -
                                                   34.0 / 120.0 * dt * f[1];
            Y3_z[i * ny * (nz - 1) + j * (nz - 1) + k] = Y2_z[i * ny * (nz - 1) + j * (nz - 1) + k] +
                                                         50.0 / 120.0 * dt * f_y2[2] -
                                                         34.0 / 120.0 * dt * f[2];

            f_y3 = functionF(Y3_x, Y3_y, Y3_z, i, j, k, time + 80 / 120 * dt);

            grid.u[i * ny * nz + j * nz + k] = Y3_x[i * ny * nz + j * nz + k] +
                                               90.0 / 120.0 * dt * f_y3[0] -
                                               50.0 / 120.0 * dt * f_y2[0];
            grid.v[i * (ny - 1) * nz + j * nz + k] = Y3_y[i * (ny - 1) * nz + j * nz + k] +
                                                     90.0 / 120.0 * dt * f_y3[1] -
                                                     50.0 / 120.0 * dt * f_y2[1];
            grid.w[i * ny * (nz - 1) + j * (nz - 1) + k] = Y3_z[i * ny * (nz - 1) + j * (nz - 1) + k] +
                                                           90.0 / 120.0 * dt * f_y3[2] -
                                                           50.0 / 120.0 * dt * f_y2[2];
        }
    }

    for (size_t i = 1; i < nx - 1; i++)
    {
        for (size_t k = 1; k < nz - 2; k++)
        {
            size_t j = ny - 2;
            f = functionF(grid.u, grid.v, grid.w, i, j, k, time);

            Y2_x[i * ny * nz + j * nz + k] = grid.u[i * ny * nz + j * nz + k] +
                                             64.0 / 120.0 * dt * f[0];
            Y2_y[i * (ny - 1) * nz + j * nz + k] = grid.v[i * (ny - 1) * nz + j * nz + k] +
                                                   64.0 / 120.0 * dt * f[1];
            Y2_z[i * ny * (nz - 1) + j * (nz - 1) + k] = grid.w[i * ny * (nz - 1) + j * (nz - 1) + k] +
                                                         64.0 / 120.0 * dt * f[2];

            f_y2 = functionF(Y2_x, Y2_y, Y2_z, i, j, k, time + 64.0 / 120.0 * dt);

            Y3_x[i * ny * nz + j * nz + k] = Y2_x[i * ny * nz + j * nz + k] +
                                             50.0 / 120.0 * dt * f_y2[0] -
                                             34.0 / 120.0 * dt * f[0];
            Y3_y[i * (ny - 1) * nz + j * nz + k] = Y2_y[i * (ny - 1) * nz + j * nz + k] +
                                                   50.0 / 120.0 * dt * f_y2[1] -
                                                   34.0 / 120.0 * dt * f[1];
            Y3_z[i * ny * (nz - 1) + j * (nz - 1) + k] = Y2_z[i * ny * (nz - 1) + j * (nz - 1) + k] +
                                                         50.0 / 120.0 * dt * f_y2[2] -
                                                         34.0 / 120.0 * dt * f[2];

            f_y3 = functionF(Y3_x, Y3_y, Y3_z, i, j, k, time + 80 / 120 * dt);

            grid.u[i * ny * nz + j * nz + k] = Y3_x[i * ny * nz + j * nz + k] +
                                               90.0 / 120.0 * dt * f_y3[0] -
                                               50.0 / 120.0 * dt * f_y2[0];
            grid.v[i * (ny - 1) * nz + j * nz + k] = Y3_y[i * (ny - 1) * nz + j * nz + k] +
                                                     90.0 / 120.0 * dt * f_y3[1] -
                                                     50.0 / 120.0 * dt * f_y2[1];
            grid.w[i * ny * (nz - 1) + j * (nz - 1) + k] = Y3_z[i * ny * (nz - 1) + j * (nz - 1) + k] +
                                                           90.0 / 120.0 * dt * f_y3[2] -
                                                           50.0 / 120.0 * dt * f_y2[2];
        }
    }

    // print the grid values.
    double diff;
    int count = 0;
    for (size_t i = 0; i < nx; i++)
    {
        for (size_t j = 0; j < ny; j++)
        {
            for (size_t k = 0; k < nz; k++)
            {
                diff = std::abs(grid.u[i * ny * nz + j * nz + k] - std::sin(i * dx) * std::cos(j * dy) * std::sin(k * dz) * std::sin(time));
                // diff = grid.u[i * ny * nz + j * nz + k];

                if (diff > 0.01)
                {
                    // std::cout << "diff[" << i << "," << j << "," << k << "] = " << diff << std::endl;
                    count++;
                }
                // std::cout << "u[" << i << "," << j << "," << k << "] = " << grid.u[i * ny * nz + j * nz + k] << std::endl;
                // std::cout << "v[" << i << "," << j << "," << k << "] = " << grid.v[i * ny * nz + j * nz + k] << std::endl;
                // std::cout << "w[" << i << "," << j << "," << k << "] = " << grid.w[i * ny * nz + j * nz + k] << std::endl;
            }
        }
    }
    std::cout << "at time " << time << " there are " << count << " cells whith error > 0.01" << std::endl;
}

void IcoNS::apply_boundary_conditions(double time)
{
    // compute boundary conditions
    for (size_t i = 0; i < nx; i++)
    {
        for (size_t j = 0; j < ny; j++)
        {
            // bottom face -> k = 0 -> sin(k*dz) = 0, cos(k*dz) = 1
            grid.u[i * ny * nz + j * nz] = 0.0;
            grid.v[i * ny * nz + j * nz] = 0.0;
            grid.w[i * ny * nz + j * nz] = 2 * std::cos(i * dx) * std::cos(j * dy) * std::sin(time);
            // top face -> k*dz = lz
            grid.u[i * ny * nz + j * nz + nz - 1] = std::sin(i * dx) * std::cos(j * dy) * std::sin(lz) * std::sin(time);
            grid.v[i * ny * nz + j * nz + nz - 1] = std::cos(i * dx) * std::sin(j * dy) * std::sin(lz) * std::sin(time);
            grid.w[i * ny * nz + j * nz + nz - 1] = 2 * std::cos(i * dx) * std::cos(j * dy) * std::cos(lz) * std::sin(time);
        }
    }

    for (size_t i = 0; i < nx; i++)
    {
        for (size_t k = 0; k < nz; k++)
        {
            // front face -> j = 0 -> sin(j*dy) = 0, cos(j*dy) = 1
            grid.u[i * ny * nz + k] = std::sin(i * dx) * std::sin(k * dz) * std::sin(time);
            grid.v[i * ny * nz + k] = 0.0;
            grid.w[i * ny * nz + k] = 2 * std::cos(i * dx) * std::cos(k * dz) * std::sin(time);
            // back face -> j*dy = ly
            grid.u[i * ny * nz + (ny - 1) * nz + k] = std::sin(i * dx) * std::cos(ly) * std::sin(k * dz) * std::sin(time);
            grid.v[i * ny * nz + (ny - 1) * nz + k] = std::cos(i * dx) * std::sin(ly) * std::sin(k * dz) * std::sin(time);
            grid.w[i * ny * nz + (ny - 1) * nz + k] = 2 * std::cos(i * dx) * std::cos(ly) * std::cos(k * dz) * std::sin(time);
        }
    }

    for (size_t j = 0; j < ny; j++)
    {
        for (size_t k = 0; k < nz; k++)
        {
            // left face -> i = 0 -> sin(i*dx) = 0, cos(i*dx) = 1
            grid.u[j * nz + k] = 0.0;
            grid.v[j * nz + k] = std::sin(j * dy) * std::sin(k * dz) * std::sin(time);
            grid.w[j * nz + k] = 2 * std::cos(j * dy) * std::cos(k * dz) * std::sin(time);
            // right face -> i*dx = lx
            grid.u[(nx - 1) * ny * nz + j * nz + k] = std::sin(lx) * std::cos(j * dy) * std::sin(k * dz) * std::sin(time);
            grid.v[(nx - 1) * ny * nz + j * nz + k] = std::cos(lx) * std::sin(j * dy) * std::sin(k * dz) * std::sin(time);
            grid.w[(nx - 1) * ny * nz + j * nz + k] = 2 * std::cos(lx) * std::cos(j * dy) * std::cos(k * dz) * std::sin(time);
        }
    }

}
double IcoNS::error_comp_X(const double t)
{
    double error = 0.0;

    // first slice (left face)
    {
        error += ((grid.u[0] - exact_solution.value_x(0, 0, 0, t)) *
                  (grid.u[0] - exact_solution.value_x(0, 0, 0, t)) *
                  dx * dy * dz / 8);

        for (size_t k = 1; k < nz - 1; k++)
        {
            error += ((grid.u[k] - exact_solution.value_x(0, 0, k, t)) *
                      (grid.u[k] - exact_solution.value_x(0, 0, k, t)) *
                      dx * dy * dz / 4);
        }

        error += ((grid.u[nz - 1] - exact_solution.value_x(0, 0, nz - 1, t)) *
                  (grid.u[nz - 1] - exact_solution.value_x(0, 0, nz - 1, t)) *
                  dx * dy * dz / 8);

        for (size_t j = 1; j < ny - 1; j++)
        {
            error += ((grid.u[j * nz] - exact_solution.value_x(0, j, 0, t)) *
                      (grid.u[j * nz] - exact_solution.value_x(0, j, 0, t)) *
                      dx * dy * dz / 4);
            for (size_t k = 1; k < nz - 1; k++)
            {
                error += ((grid.u[j * nz + k] - exact_solution.value_x(0, j, k, t)) *
                          (grid.u[j * nz + k] - exact_solution.value_x(0, j, k, t)) *
                          dx * dy * dz / 2);
            }
            error += ((grid.u[j * nz + nz - 1] - exact_solution.value_x(0, j, nz - 1, t)) *
                      (grid.u[j * nz + nz - 1] - exact_solution.value_x(0, j, nz - 1, t)) *
                      dx * dy * dz / 4);
        }

        error += ((grid.u[(ny - 1) * nz] - exact_solution.value_x(0, (ny - 1), 0, t)) *
                  (grid.u[(ny - 1) * nz] - exact_solution.value_x(0, (ny - 1), 0, t)) *
                  dx * dy * dz / 8);

        for (size_t k = 1; k < nz - 1; k++)
        {
            error += ((grid.u[(ny - 1) * nz + k] - exact_solution.value_x(0, (ny - 1), k, t)) *
                      (grid.u[(ny - 1) * nz + k] - exact_solution.value_x(0, (ny - 1), k, t)) *
                      dx * dy * dz / 4);
        }

        error += ((grid.u[(ny - 1) * nz + nz - 1] - exact_solution.value_x(0, (ny - 1), nz - 1, t)) *
                  (grid.u[(ny - 1) * nz + nz - 1] - exact_solution.value_x(0, (ny - 1), nz - 1, t)) *
                  dx * dy * dz / 8);
    }

    // middle slices
    {
        for (size_t i = 1; i < nx - 1; i++)
        {
            error += ((grid.u[i * ny * nz] - exact_solution.value_x(i, 0, 0, t)) *
                      (grid.u[i * ny * nz] - exact_solution.value_x(i, 0, 0, t)) *
                      dx * dy * dz / 4);

            for (size_t k = 1; k < nz - 1; k++)
            {
                error += ((grid.u[i * ny * nz + k] - exact_solution.value_x(i, 0, k, t)) *
                          (grid.u[i * ny * nz + k] - exact_solution.value_x(i, 0, k, t)) *
                          dx * dy * dz / 2);
            }
            error += ((grid.u[i * ny * nz + nz - 1] - exact_solution.value_x(i, 0, nz - 1, t)) *
                      (grid.u[i * ny * nz + nz - 1] - exact_solution.value_x(i, 0, nz - 1, t)) *
                      dx * dy * dz / 4);

            for (size_t j = 1; j < ny - 1; j++)
            {
                error += ((grid.u[i * ny * nz + j * nz] - exact_solution.value_x(i, j, 0, t)) *
                          (grid.u[i * ny * nz + j * nz] - exact_solution.value_x(i, j, 0, t)) *
                          dx * dy * dz / 2);

                for (size_t k = 1; k < nz - 1; k++)
                {
                    error += ((grid.u[i * ny * nz + j * nz + k] - exact_solution.value_x(i, j, k, t)) *
                              (grid.u[i * ny * nz + j * nz + k] - exact_solution.value_x(i, j, k, t)) *
                              dx * dy * dz);
                }

                error += ((grid.u[i * ny * nz + j * nz + nz - 1] - exact_solution.value_x(i, j, nz - 1, t)) *
                          (grid.u[i * ny * nz + j * nz + nz - 1] - exact_solution.value_x(i, j, nz - 1, t)) *
                          dx * dy * dz / 2);
            }

            error += ((grid.u[i * ny * nz + (ny - 1) * nz] - exact_solution.value_x(i, (ny - 1), nz - 1, t)) *
                      (grid.u[i * ny * nz + (ny - 1) * nz] - exact_solution.value_x(i, (ny - 1), nz - 1, t)) *
                      dx * dy * dz / 4);

            for (size_t k = 1; k < nz - 1; k++)
            {
                error += ((grid.u[i * ny * nz + (ny - 1) * nz + k] - exact_solution.value_x(i, (ny - 1), k, t)) *
                          (grid.u[i * ny * nz + (ny - 1) * nz + k] - exact_solution.value_x(i, (ny - 1), k, t)) *
                          dx * dy * dz / 2);
            }

            error += ((grid.u[i * ny * nz + (ny - 1) * nz + nz - 1] - exact_solution.value_x(i, (ny - 1), nz - 1, t)) *
                      (grid.u[i * ny * nz + (ny - 1) * nz + nz - 1] - exact_solution.value_x(i, (ny - 1), nz - 1, t)) *
                      dx * dy * dz / 4);
        }
    }

    // last slice (right face)
    {
        error += ((grid.u[(nx - 1) * ny * nz] - exact_solution.value_x((nx - 1), 0, 0, t)) *
                  (grid.u[(nx - 1) * ny * nz] - exact_solution.value_x((nx - 1), 0, 0, t)) *
                  dx * dy * dz / 8);

        for (size_t k = 1; k < nz - 1; k++)
        {
            error += ((grid.u[(nx - 1) * ny * nz + k] - exact_solution.value_x((nx - 1), 0, k, t)) *
                      (grid.u[(nx - 1) * ny * nz + k] - exact_solution.value_x((nx - 1), 0, k, t)) *
                      dx * dy * dz / 4);
        }

        error += ((grid.u[(nx - 1) * ny * nz + nz - 1] - exact_solution.value_x((nx - 1), 0, nz - 1, t)) *
                  (grid.u[(nx - 1) * ny * nz + nz - 1] - exact_solution.value_x((nx - 1), 0, nz - 1, t)) *
                  dx * dy * dz / 8);

        for (size_t j = 1; j < ny - 1; j++)
        {
            error += ((grid.u[(nx - 1) * ny * nz + j * nz] - exact_solution.value_x((nx - 1), j, 0, t)) *
                      (grid.u[(nx - 1) * ny * nz + j * nz] - exact_solution.value_x((nx - 1), j, 0, t)) *
                      dx * dy * dz / 4);
            for (size_t k = 1; k < nz - 1; k++)
            {
                error += ((grid.u[(nx - 1) * ny * nz + j * nz + k] - exact_solution.value_x((nx - 1), j, k, t)) *
                          (grid.u[(nx - 1) * ny * nz + j * nz + k] - exact_solution.value_x((nx - 1), j, k, t)) *
                          dx * dy * dz / 2);
            }
            error += ((grid.u[(nx - 1) * ny * nz + j * nz + nz - 1] - exact_solution.value_x((nx - 1), j, nz - 1, t)) *
                      (grid.u[(nx - 1) * ny * nz + j * nz + nz - 1] - exact_solution.value_x((nx - 1), j, nz - 1, t)) *
                      dx * dy * dz / 4);
        }

        error += ((grid.u[(nx - 1) * ny * nz + (ny - 1) * nz] - exact_solution.value_x((nx - 1), (ny - 1), 0, t)) *
                  (grid.u[(nx - 1) * ny * nz + (ny - 1) * nz] - exact_solution.value_x((nx - 1), (ny - 1), 0, t)) *
                  dx * dy * dz / 8);

        for (size_t k = 1; k < nz - 1; k++)
        {
            error += ((grid.u[(nx - 1) * ny * nz + (ny - 1) * nz + k] - exact_solution.value_x((nx - 1), (ny - 1), k, t)) *
                      (grid.u[(nx - 1) * ny * nz + (ny - 1) * nz + k] - exact_solution.value_x((nx - 1), (ny - 1), k, t)) *
                      dx * dy * dz / 4);
        }

        error += ((grid.u[(nx - 1) * ny * nz + (ny - 1) * nz + nz - 1] - exact_solution.value_x((nx - 1), (ny - 1), nz - 1, t)) *
                  (grid.u[(nx - 1) * ny * nz + (ny - 1) * nz + nz - 1] - exact_solution.value_x((nx - 1), (ny - 1), nz - 1, t)) *
                  dx * dy * dz / 8);
    }
    return error;
}

double IcoNS::error_comp_Y(const double t)
{
    double error = 0.0;
    // first slice (left face)
    {
        error += ((grid.w[0] - exact_solution.value_y(0, 0, 0, t)) *
                  (grid.w[0] - exact_solution.value_y(0, 0, 0, t)) *
                  dx * dy * dz / 8);

        for (size_t k = 1; k < nz - 1; k++)
        {
            error += ((grid.w[k] - exact_solution.value_y(0, 0, k, t)) *
                      (grid.w[k] - exact_solution.value_y(0, 0, k, t)) *
                      dx * dy * dz / 4);
        }

        error += ((grid.w[nz - 1] - exact_solution.value_y(0, 0, nz - 1, t)) *
                  (grid.w[nz - 1] - exact_solution.value_y(0, 0, nz - 1, t)) *
                  dx * dy * dz / 8);

        for (size_t j = 1; j < ny - 1; j++)
        {
            error += ((grid.w[j * nz] - exact_solution.value_y(0, j, 0, t)) *
                      (grid.w[j * nz] - exact_solution.value_y(0, j, 0, t)) *
                      dx * dy * dz / 4);
            for (size_t k = 1; k < nz - 1; k++)
            {
                error += ((grid.w[j * nz + k] - exact_solution.value_y(0, j, k, t)) *
                          (grid.w[j * nz + k] - exact_solution.value_y(0, j, k, t)) *
                          dx * dy * dz / 2);
            }
            error += ((grid.w[j * nz + nz - 1] - exact_solution.value_y(0, j, nz - 1, t)) *
                      (grid.w[j * nz + nz - 1] - exact_solution.value_y(0, j, nz - 1, t)) *
                      dx * dy * dz / 4);
        }

        error += ((grid.w[(ny - 1) * nz] - exact_solution.value_y(0, (ny - 1), 0, t)) *
                  (grid.w[(ny - 1) * nz] - exact_solution.value_y(0, (ny - 1), 0, t)) *
                  dx * dy * dz / 8);

        for (size_t k = 1; k < nz - 1; k++)
        {
            error += ((grid.w[(ny - 1) * nz + k] - exact_solution.value_y(0, (ny - 1), k, t)) *
                      (grid.w[(ny - 1) * nz + k] - exact_solution.value_y(0, (ny - 1), k, t)) *
                      dx * dy * dz / 4);
        }

        error += ((grid.w[(ny - 1) * nz + nz - 1] - exact_solution.value_y(0, (ny - 1), nz - 1, t)) *
                  (grid.w[(ny - 1) * nz + nz - 1] - exact_solution.value_y(0, (ny - 1), nz - 1, t)) *
                  dx * dy * dz / 8);
    }

    // middle slices
    {
        for (size_t i = 1; i < nx - 1; i++)
        {
            error += ((grid.w[i * ny * nz] - exact_solution.value_y(i, 0, 0, t)) *
                      (grid.w[i * ny * nz] - exact_solution.value_y(i, 0, 0, t)) *
                      dx * dy * dz / 4);

            for (size_t k = 1; k < nz - 1; k++)
            {
                error += ((grid.w[i * ny * nz + k] - exact_solution.value_y(i, 0, k, t)) *
                          (grid.w[i * ny * nz + k] - exact_solution.value_y(i, 0, k, t)) *
                          dx * dy * dz / 2);
            }
            error += ((grid.w[i * ny * nz + nz - 1] - exact_solution.value_y(i, 0, nz - 1, t)) *
                      (grid.w[i * ny * nz + nz - 1] - exact_solution.value_y(i, 0, nz - 1, t)) *
                      dx * dy * dz / 4);

            for (size_t j = 1; j < ny - 1; j++)
            {
                error += ((grid.w[i * ny * nz + j * nz] - exact_solution.value_y(i, j, 0, t)) *
                          (grid.w[i * ny * nz + j * nz] - exact_solution.value_y(i, j, 0, t)) *
                          dx * dy * dz / 2);

                for (size_t k = 1; k < nz - 1; k++)
                {
                    error += ((grid.w[i * ny * nz + j * nz + k] - exact_solution.value_y(i, j, k, t)) *
                              (grid.w[i * ny * nz + j * nz + k] - exact_solution.value_y(i, j, k, t)) *
                              dx * dy * dz);
                }

                error += ((grid.w[i * ny * nz + j * nz + nz - 1] - exact_solution.value_y(i, j, (nz - 1), t)) *
                          (grid.w[i * ny * nz + j * nz + nz - 1] - exact_solution.value_y(i, j, (nz - 1), t)) *
                          dx * dy * dz / 2);
            }

            error += ((grid.w[i * ny * nz + (ny - 1) * nz] - exact_solution.value_y(i, (ny - 1), (nz - 1), t)) *
                      (grid.w[i * ny * nz + (ny - 1) * nz] - exact_solution.value_y(i, (ny - 1), (nz - 1), t)) *
                      dx * dy * dz / 4);

            for (size_t k = 1; k < nz - 1; k++)
            {
                error += ((grid.w[i * ny * nz + (ny - 1) * nz + k] - exact_solution.value_y(i, (ny - 1), k, t)) *
                          (grid.w[i * ny * nz + (ny - 1) * nz + k] - exact_solution.value_y(i, (ny - 1), k, t)) *
                          dx * dy * dz / 2);
            }

            error += ((grid.w[i * ny * nz + (ny - 1) * nz + nz - 1] - exact_solution.value_y(i, (ny - 1), nz - 1, t)) *
                      (grid.w[i * ny * nz + (ny - 1) * nz + nz - 1] - exact_solution.value_y(i, (ny - 1), nz - 1, t)) *
                      dx * dy * dz / 4);
        }
    }

    // last slice (right face)
    {
        error += ((grid.w[(nx - 1) * ny * nz] - exact_solution.value_y((nx - 1), 0, 0, t)) *
                  (grid.w[(nx - 1) * ny * nz] - exact_solution.value_y((nx - 1), 0, 0, t)) *
                  dx * dy * dz / 8);

        for (size_t k = 1; k < nz - 1; k++)
        {
            error += ((grid.w[(nx - 1) * ny * nz + k] - exact_solution.value_y((nx - 1), 0, k, t)) *
                      (grid.w[(nx - 1) * ny * nz + k] - exact_solution.value_y((nx - 1), 0, k, t)) *
                      dx * dy * dz / 4);
        }

        error += ((grid.w[(nx - 1) * ny * nz + nz - 1] - exact_solution.value_y((nx - 1), 0, nz - 1, t)) *
                  (grid.w[(nx - 1) * ny * nz + nz - 1] - exact_solution.value_y((nx - 1), 0, nz - 1, t)) *
                  dx * dy * dz / 8);

        for (size_t j = 1; j < ny - 1; j++)
        {
            error += ((grid.w[(nx - 1) * ny * nz + j * nz] - exact_solution.value_y((nx - 1), j, 0, t)) *
                      (grid.w[(nx - 1) * ny * nz + j * nz] - exact_solution.value_y((nx - 1), j, 0, t)) *
                      dx * dy * dz / 4);
            for (size_t k = 1; k < nz - 1; k++)
            {
                error += ((grid.w[(nx - 1) * ny * nz + j * nz + k] - exact_solution.value_y((nx - 1), j, k, t)) *
                          (grid.w[(nx - 1) * ny * nz + j * nz + k] - exact_solution.value_y((nx - 1), j, k, t)) *
                          dx * dy * dz / 2);
            }
            error += ((grid.w[(nx - 1) * ny * nz + j * nz + nz - 1] - exact_solution.value_y((nx - 1), j, nz - 1, t)) *
                      (grid.w[(nx - 1) * ny * nz + j * nz + nz - 1] - exact_solution.value_y((nx - 1), j, nz - 1, t)) *
                      dx * dy * dz / 4);
        }

        error += ((grid.w[(nx - 1) * ny * nz + (ny - 1) * nz] - exact_solution.value_y((nx - 1), (ny - 1), 0, t)) *
                  (grid.w[(nx - 1) * ny * nz + (ny - 1) * nz] - exact_solution.value_y((nx - 1), (ny - 1), 0, t)) *
                  dx * dy * dz / 8);

        for (size_t k = 1; k < nz - 1; k++)
        {
            error += ((grid.w[(nx - 1) * ny * nz + (ny - 1) * nz + k] - exact_solution.value_y((nx - 1), (ny - 1), k, t)) *
                      (grid.w[(nx - 1) * ny * nz + (ny - 1) * nz + k] - exact_solution.value_y((nx - 1), (ny - 1), k, t)) *
                      dx * dy * dz / 4);
        }

        error += ((grid.w[(nx - 1) * ny * nz + (ny - 1) * nz + nz - 1] - exact_solution.value_y((nx - 1), (ny - 1), nz - 1, t)) *
                  (grid.w[(nx - 1) * ny * nz + (ny - 1) * nz + nz - 1] - exact_solution.value_y((nx - 1), (ny - 1), nz - 1, t)) *
                  dx * dy * dz / 8);
    }

    return error;
}

double IcoNS::error_comp_Z(const double t)
{
    double error = 0.0;

    // first slice (left face)
    {
        error += ((grid.v[0] - exact_solution.value_z(0, 0, 0, t)) *
                  (grid.v[0] - exact_solution.value_z(0, 0, 0, t)) *
                  dx * dy * dz / 8);

        for (size_t k = 1; k < nz - 1; k++)
        {
            error += ((grid.v[k] - exact_solution.value_z(0, 0, k, t)) *
                      (grid.v[k] - exact_solution.value_z(0, 0, k, t)) *
                      dx * dy * dz / 4);
        }

        error += ((grid.v[nz - 1] - exact_solution.value_z(0, 0, nz - 1, t)) *
                  (grid.v[nz - 1] - exact_solution.value_z(0, 0, nz - 1, t)) *
                  dx * dy * dz / 8);

        for (size_t j = 1; j < ny - 1; j++)
        {
            error += ((grid.v[j * nz] - exact_solution.value_z(0, j, 0, t)) *
                      (grid.v[j * nz] - exact_solution.value_z(0, j, 0, t)) *
                      dx * dy * dz / 4);
            for (size_t k = 1; k < nz - 1; k++)
            {
                error += ((grid.v[j * nz + k] - exact_solution.value_z(0, j, k, t)) *
                          (grid.v[j * nz + k] - exact_solution.value_z(0, j, k, t)) *
                          dx * dy * dz / 2);
            }
            error += ((grid.v[j * nz + nz - 1] - exact_solution.value_z(0, j, nz - 1, t)) *
                      (grid.v[j * nz + nz - 1] - exact_solution.value_z(0, j, nz - 1, t)) *
                      dx * dy * dz / 4);
        }

        error += ((grid.v[(ny - 1) * nz] - exact_solution.value_z(0, (ny - 1), 0, t)) *
                  (grid.v[(ny - 1) * nz] - exact_solution.value_z(0, (ny - 1), 0, t)) *
                  dx * dy * dz / 8);

        for (size_t k = 1; k < nz - 1; k++)
        {
            error += ((grid.v[(ny - 1) * nz + k] - exact_solution.value_z(0, (ny - 1), k, t)) *
                      (grid.v[(ny - 1) * nz + k] - exact_solution.value_z(0, (ny - 1), k, t)) *
                      dx * dy * dz / 4);
        }

        error += ((grid.v[(ny - 1) * nz + nz - 1] - exact_solution.value_z(0, (ny - 1), nz - 1, t)) *
                  (grid.v[(ny - 1) * nz + nz - 1] - exact_solution.value_z(0, (ny - 1), nz - 1, t)) *
                  dx * dy * dz / 8);
    }

    // middle slices
    {
        for (size_t i = 1; i < nx - 1; i++)
        {
            error += ((grid.v[i * ny * nz] - exact_solution.value_z(i, 0, 0, t)) *
                      (grid.v[i * ny * nz] - exact_solution.value_z(i, 0, 0, t)) *
                      dx * dy * dz / 4);

            for (size_t k = 1; k < nz - 1; k++)
            {
                error += ((grid.v[i * ny * nz + k] - exact_solution.value_z(i, 0, k, t)) *
                          (grid.v[i * ny * nz + k] - exact_solution.value_z(i, 0, k, t)) *
                          dx * dy * dz / 2);
            }
            error += ((grid.v[i * ny * nz + nz - 1] - exact_solution.value_z(i, 0, nz - 1, t)) *
                      (grid.v[i * ny * nz + nz - 1] - exact_solution.value_z(i, 0, nz - 1, t)) *
                      dx * dy * dz / 4);

            for (size_t j = 1; j < ny - 1; j++)
            {
                error += ((grid.v[i * ny * nz + j * nz] - exact_solution.value_z(i, j, 0, t)) *
                          (grid.v[i * ny * nz + j * nz] - exact_solution.value_z(i, j, 0, t)) *
                          dx * dy * dz / 2);

                for (size_t k = 1; k < nz - 1; k++)
                {
                    error += ((grid.v[i * ny * nz + j * nz + k] - exact_solution.value_z(i, j, k, t)) *
                              (grid.v[i * ny * nz + j * nz + k] - exact_solution.value_z(i, j, k, t)) *
                              dx * dy * dz);
                }

                error += ((grid.v[i * ny * nz + j * nz + nz - 1] - exact_solution.value_z(i, j, nz - 1, t)) *
                          (grid.v[i * ny * nz + j * nz + nz - 1] - exact_solution.value_z(i, j, nz - 1, t)) *
                          dx * dy * dz / 2);
            }

            error += ((grid.v[i * ny * nz + (ny - 1) * nz] - exact_solution.value_z(i, (ny - 1), nz - 1, t)) *
                      (grid.v[i * ny * nz + (ny - 1) * nz] - exact_solution.value_z(i, (ny - 1), nz - 1, t)) *
                      dx * dy * dz / 4);

            for (size_t k = 1; k < nz - 1; k++)
            {
                error += ((grid.v[i * ny * nz + (ny - 1) * nz + k] - exact_solution.value_z(i, (ny - 1), k, t)) *
                          (grid.v[i * ny * nz + (ny - 1) * nz + k] - exact_solution.value_z(i, (ny - 1), k, t)) *
                          dx * dy * dz / 2);
            }

            error += ((grid.v[i * ny * nz + (ny - 1) * nz + nz - 1] - exact_solution.value_z(i, (ny - 1), nz - 1, t)) *
                      (grid.v[i * ny * nz + (ny - 1) * nz + nz - 1] - exact_solution.value_z(i, (ny - 1), nz - 1, t)) *
                      dx * dy * dz / 4);
        }
    }

    // last slice (right face)
    {
        error += ((grid.v[(nx - 1) * ny * nz] - exact_solution.value_z((nx - 1), 0, 0, t)) *
                  (grid.v[(nx - 1) * ny * nz] - exact_solution.value_z((nx - 1), 0, 0, t)) *
                  dx * dy * dz / 8);

        for (size_t k = 1; k < nz - 1; k++)
        {
            error += ((grid.v[(nx - 1) * ny * nz + k] - exact_solution.value_z((nx - 1), 0, k, t)) *
                      (grid.v[(nx - 1) * ny * nz + k] - exact_solution.value_z((nx - 1), 0, k, t)) *
                      dx * dy * dz / 4);
        }

        error += ((grid.v[(nx - 1) * ny * nz + nz - 1] - exact_solution.value_z((nx - 1), 0, nz - 1, t)) *
                  (grid.v[(nx - 1) * ny * nz + nz - 1] - exact_solution.value_z((nx - 1), 0, nz - 1, t)) *
                  dx * dy * dz / 8);

        for (size_t j = 1; j < ny - 1; j++)
        {
            error += ((grid.v[(nx - 1) * ny * nz + j * nz] - exact_solution.value_z((nx - 1), j, 0, t)) *
                      (grid.v[(nx - 1) * ny * nz + j * nz] - exact_solution.value_z((nx - 1), j, 0, t)) *
                      dx * dy * dz / 4);
            for (size_t k = 1; k < nz - 1; k++)
            {
                error += ((grid.v[(nx - 1) * ny * nz + j * nz + k] - exact_solution.value_z((nx - 1), j, k, t)) *
                          (grid.v[(nx - 1) * ny * nz + j * nz + k] - exact_solution.value_z((nx - 1), j, k, t)) *
                          dx * dy * dz / 2);
            }
            error += ((grid.v[(nx - 1) * ny * nz + j * nz + nz - 1] - exact_solution.value_z((nx - 1), j, nz - 1, t)) *
                      (grid.v[(nx - 1) * ny * nz + j * nz + nz - 1] - exact_solution.value_z((nx - 1), j, nz - 1, t)) *
                      dx * dy * dz / 4);
        }

        error += ((grid.v[(nx - 1) * ny * nz + (ny - 1) * nz] - exact_solution.value_z((nx - 1), (ny - 1), 0, t)) *
                  (grid.v[(nx - 1) * ny * nz + (ny - 1) * nz] - exact_solution.value_z((nx - 1), (ny - 1), 0, t)) *
                  dx * dy * dz / 8);

        for (size_t k = 1; k < nz - 1; k++)
        {
            error += ((grid.v[(nx - 1) * ny * nz + (ny - 1) * nz + k] - exact_solution.value_z((nx - 1), (ny - 1), k, t)) *
                      (grid.v[(nx - 1) * ny * nz + (ny - 1) * nz + k] - exact_solution.value_z((nx - 1), (ny - 1), k, t)) *
                      dx * dy * dz / 4);
        }

        error += ((grid.v[(nx - 1) * ny * nz + (ny - 1) * nz + nz - 1] - exact_solution.value_z((nx - 1), (ny - 1), nz - 1, t)) *
                  (grid.v[(nx - 1) * ny * nz + (ny - 1) * nz + nz - 1] - exact_solution.value_z((nx - 1), (ny - 1), nz - 1, t)) *
                  dx * dy * dz / 8);
    }
    return error;
}

double IcoNS::L2_error(const double t)
{
    double error = 0.0;

    error += error_comp_X(t);
    error += error_comp_Y(t);
    error += error_comp_Z(t);

    return sqrt(error);
}

void IcoNS::output()
{
}
