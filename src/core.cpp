
#include <cmath>

#include "core.hpp"
#include <math.h>

#include <fstream>
#include <string>
#include <memory>

void IcoNS::preprocessing(/*std::string &input_file*/)

{
    // read the input file.
    for (size_t i = 1; i < NX - 1; i++)
    {
        for (size_t j = 1; j < NY - 1; j++)
        {
            for (size_t k = 1; k < NZ - 1; k++)
            {
                grid.u[i][j][k] = 3.0;
                grid.v[i][j][k] = 3.0;
                grid.w[i][j][k] = 3.0;
            }
        }
    }

    for (size_t i = 0; i < NX; i++)
    {
        for (size_t j = 0; j < NY; j++)
        {
            for (size_t k = 0; k < NZ; k++)
            {
                if (i == 0 || i == NX || j == 0 || j == NY || k == 0 || k == NZ)
                {
                    grid.u[i][j][k] = 1.0;
                    grid.v[i][j][k] = 1.0;
                    grid.w[i][j][k] = 1.0;
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
        // csv file w/ "," delimiter: time step, iter, L2_error
        std::cout << time << "," << i << "," << " Error: " << L2_error(time) << std::endl;
        i++;
    }
}

std::array<double, 3> IcoNS::functionF(const std::vector<std::vector<std::vector<double>>> &u,
                                       const std::vector<std::vector<std::vector<double>>> &v,
                                       const std::vector<std::vector<std::vector<double>>> &w,
                                       size_t i, size_t j, size_t k, double t)
{
    std::array<double, 3> f;
    size_t l = i * NY * NZ + j * NZ + k;
    std::vector<double> g(3);
    g = functionG(i * dx, j * dy, k * dz, t);

    f[0] = -(u[i][j][k] * (u[i + 1][j][k] - u[i - 1][j][k]) / (2.0 * dx) +
             (v[i + 1][j][k] + v[i][j][k] + v[i][j - 1][k] + v[i + 1][j - 1][k]) / 4.0 * (u[i][j + 1][k] - u[i][j - 1][k]) / (2.0 * dy) +
             (w[i + 1][j][k] + w[i][j][k] + w[i + 1][j][k - 1] + w[i][j][k - 1]) / 4.0 * (u[i][j][k + 1] - u[i][j][k - 1]) / (2.0 * dz)) +
           1 / Re * ((u[i + 1][j][k] - 2 * u[i][j][k] + u[i - 1][j][k]) / (dx * dx) + (u[i][j + 1][k] - 2 * u[i][j][k] + u[i][j - 1][k]) / (dy * dy) + (u[i][j][k + 1] - 2 * u[i][j][k] + u[i][j][k - 1]) / (dz * dz)) + g[0];

    f[1] = -((u[i][j][k] + u[i][j + 1][k] + u[i - 1][j][k] + u[i - 1][j + 1][k]) / 4.0 * (v[i + 1][j][k] - v[i - 1][j][k] / (2.0 * dx)) +
             v[i][j][k] * (v[i][j + 1][k] - v[i][j - 1][k]) / (2.0 * dy) +
             (w[i][j][k] + w[i][j][k - 1] + w[i][j + 1][k] + w[i][j - 1][k]) / 4.0 + (v[i][j][k + 1] - v[i][j][k - 1]) / (2.0 * dz)) +
           1 / Re * ((v[i + 1][j][k] - 2.0 * v[i][j][k] + v[i - 1][j][k]) / (dx * dx) + (v[i][j + 1][k] - 2.0 * v[i][j][k] + v[i][j - 1][k]) / (dy * dy) + (v[i][j][k + 1] - 2.0 * v[i][j][k] + v[i][j][k - 1]) / (dz * dz)) + g[1];

    f[2] = -((u[i][j][k] + u[i - 1][j][k] + u[i][j][k + 1] + u[i - 1][j][k + 1]) / 4.0 * (w[i + 1][j][k] - w[i - 1][j][k]) / (2.0 * dx) +
             (v[i][j][k + 1] + v[i][j - 1][k + 1] + v[i][j][k] + v[i][j - 1][k]) / 4.0 * (w[i][j + 1][k] - w[i][j - 1][k]) / (2.0 * dy) +
             w[i][j][k] * (w[i][j][k + 1] - w[i][j][k - 1]) / (2.0 * dz)) +
           1 / Re * ((w[i + 1][j][k] - 2.0 * w[i][j][k] + w[i - 1][j][k]) / (dx * dx) + (w[i][j + 1][k] - 2.0 * w[i][j][k] + w[i][j - 1][k]) / (dy * dy) + (w[i][j][k + 1] - 2.0 * w[i][j][k] + w[i][j][k - 1]) / (dz * dz)) + g[2];

    return f;
}

std::array<double, 3> IcoNS::functionF(const std::array<std::array<std::array<double, NZ + 1>, NY + 1>, NX> &u,
                                       const std::array<std::array<std::array<double, NZ + 1>, NY>, NX + 1> &v,
                                       const std::array<std::array<std::array<double, NZ>, NY + 1>, NX + 1> &w,
                                       size_t i, size_t j, size_t k, double t)
{
    std::array<double, 3> f;
    size_t l = i * NY * NZ + j * NZ + k;
    std::vector<double> g(3);
    g = functionG(i * dx, j * dy, k * dz, t);

    f[0] = -(u[i][j][k] * (u[i + 1][j][k] - u[i - 1][j][k]) / (2.0 * dx) +
             (v[i + 1][j][k] + v[i][j][k] + v[i][j - 1][k] + v[i + 1][j - 1][k]) / 4.0 * (u[i][j + 1][k] - u[i][j - 1][k]) / (2.0 * dy) +
             (w[i + 1][j][k] + w[i][j][k] + w[i + 1][j][k - 1] + w[i][j][k - 1]) / 4.0 * (u[i][j][k + 1] - u[i][j][k - 1]) / (2.0 * dz)) +
           1 / Re * ((u[i + 1][j][k] - 2 * u[i][j][k] + u[i - 1][j][k]) / (dx * dx) + (u[i][j + 1][k] - 2 * u[i][j][k] + u[i][j - 1][k]) / (dy * dy) + (u[i][j][k + 1] - 2 * u[i][j][k] + u[i][j][k - 1]) / (dz * dz)) + g[0];

    f[1] = -((u[i][j][k] + u[i][j + 1][k] + u[i - 1][j][k] + u[i - 1][j + 1][k]) / 4.0 * (v[i + 1][j][k] - v[i - 1][j][k] / (2.0 * dx)) +
             v[i][j][k] * (v[i][j + 1][k] - v[i][j - 1][k]) / (2.0 * dy) +
             (w[i][j][k] + w[i][j][k - 1] + w[i][j + 1][k] + w[i][j - 1][k]) / 4.0 + (v[i][j][k + 1] - v[i][j][k - 1]) / (2.0 * dz)) +
           1 / Re * ((v[i + 1][j][k] - 2.0 * v[i][j][k] + v[i - 1][j][k]) / (dx * dx) + (v[i][j + 1][k] - 2.0 * v[i][j][k] + v[i][j - 1][k]) / (dy * dy) + (v[i][j][k + 1] - 2.0 * v[i][j][k] + v[i][j][k - 1]) / (dz * dz)) + g[1];

    f[2] = -((u[i][j][k] + u[i - 1][j][k] + u[i][j][k + 1] + u[i - 1][j][k + 1]) / 4.0 * (w[i + 1][j][k] - w[i - 1][j][k]) / (2.0 * dx) +
             (v[i][j][k + 1] + v[i][j - 1][k + 1] + v[i][j][k] + v[i][j - 1][k]) / 4.0 * (w[i][j + 1][k] - w[i][j - 1][k]) / (2.0 * dy) +
             w[i][j][k] * (w[i][j][k + 1] - w[i][j][k - 1]) / (2.0 * dz)) +
           1 / Re * ((w[i + 1][j][k] - 2.0 * w[i][j][k] + w[i - 1][j][k]) / (dx * dx) + (w[i][j + 1][k] - 2.0 * w[i][j][k] + w[i][j - 1][k]) / (dy * dy) + (w[i][j][k + 1] - 2.0 * w[i][j][k] + w[i][j][k - 1]) / (dz * dz)) + g[2];
    return f;
}

std::vector<double> IcoNS::functionG(size_t i, size_t j, size_t k, double t)
{
    std::vector<double> g(3);
    g[0] = std::sin(i) * std::cos(j) * std::sin(k) * std::cos(t) + std::sin(i) * std::cos(i) * std::cos(j) * std::cos(j) * std::sin(k) * std::sin(k) * std::sin(t) * std::sin(t) - std::sin(i) * std::cos(i) * std::sin(j) * std::sin(j) * std::sin(k) * std::sin(k) * std::sin(t) * std::sin(t) + 2 * std::sin(i) * std::cos(i) * std::cos(j) * std::cos(j) * std::cos(k) * std::cos(k) * std::sin(t) * std::sin(t) + 3.0 / Re * std::sin(i) * std::cos(j) * std::sin(k) * std::sin(t);
    g[1] = std::cos(i) * std::sin(j) * std::sin(k) * std::cos(t) - std::sin(i) * std::sin(i) * std::sin(j) * std::cos(j) * std::sin(k) * std::sin(k) * std::sin(t) * std::sin(t) +
           std::cos(i) * std::cos(i) * std::sin(j) * std::cos(j) * std::sin(k) * std::sin(k) * std::sin(t) * std::sin(t) +
           2.0 * std::cos(i) * std::cos(i) * std::sin(j) * std::cos(j) * std::cos(k) * std::cos(k) * std::sin(t) * std::sin(t) +
           3.0 / Re * std::cos(i) * std::sin(j) * std::sin(k) * std::sin(t);
    g[2] = 2 * std::cos(i) * std::cos(j) * std::cos(k) * std::cos(t) - 2 * std::sin(i) * std::sin(i) * std::cos(j) * std::cos(j) * std::sin(k) * std::cos(k) * std::sin(t) * std::sin(t) -
           2 * std::cos(i) * std::cos(i) * std::sin(j) * std::sin(j) * std::sin(k) * std::cos(k) * std::sin(t) * std::sin(t) -
           4.0 * std::cos(i) * std::cos(i) * std::cos(j) * std::cos(j) * std::sin(k) * std::cos(k) * std::sin(t) * std::sin(t) + 6.0 / Re * std::cos(i) * std::cos(j) * std::cos(k) * std::sin(t);

    // std::cout << "g[0] = " << g[0] << std::endl;
    return g;
}

void IcoNS::solve_time_step(double time)
{
    std::array<double, 3> f;
    std::array<double, 3> f_y2;
    std::array<double, 3> f_y3;
    std::vector<std::vector<std::vector<double>>> Y2_x(NX, std::vector<std::vector<double>>(NY + 1, std::vector<double>(NZ + 1)));
    std::vector<std::vector<std::vector<double>>> Y2_y(NX + 1, std::vector<std::vector<double>>(NY, std::vector<double>(NZ + 1)));
    std::vector<std::vector<std::vector<double>>> Y2_z(NX + 1, std::vector<std::vector<double>>(NY + 1, std::vector<double>(NZ)));
    std::vector<std::vector<std::vector<double>>> Y3_x(NX, std::vector<std::vector<double>>(NY + 1, std::vector<double>(NZ + 1)));
    std::vector<std::vector<std::vector<double>>> Y3_y(NX + 1, std::vector<std::vector<double>>(NY, std::vector<double>(NZ + 1)));
    std::vector<std::vector<std::vector<double>>> Y3_z(NX + 1, std::vector<std::vector<double>>(NY + 1, std::vector<double>(NZ)));

    for (size_t i = 1; i < NX - 1; i++)
    {
        for (size_t j = 1; j < NY - 1; j++)
        {
            for (size_t k = 1; k < NZ - 1; k++)
            {
                f = functionF(grid.u, grid.v, grid.w, i, j, k, time);
                // std::cout << "f[0] = " << f[0] << std::endl;
                //  solve the momentum equations -> TODO later also with pressure
                Y2_x[i][j][k] = grid.u[i][j][k] +
                                64.0 / 120.0 * dt * f[0];
                Y2_y[i][j][k] = /* grid.v[i][j][k] + */
                    64.0 / 120.0 * dt * f[1];
                Y2_z[i][j][k] = /* grid.w[i][j][k] + */
                    64.0 / 120.0 * dt * f[2];
                // std::cout << "dt " << dt << std::endl;

                f_y2 = functionF(Y2_x, Y2_y, Y2_z, i, j, k, time + 64.0 / 120.0 * dt);

                Y3_x[i][j][k] = Y2_x[i][j][k] +
                                50.0 / 120.0 * dt * f_y2[0] -
                                34.0 / 120.0 * dt * f[0];
                Y3_y[i][j][k] = Y2_y[i][j][k] +
                                50.0 / 120.0 * dt * f_y2[1] -
                                34.0 / 120.0 * dt * f[1];
                Y3_z[i][j][k] = Y2_z[i][j][k] +
                                50.0 / 120.0 * dt * f_y2[2] -
                                34.0 / 120.0 * dt * f[2];

                // update the grid values.
                f_y3 = functionF(Y3_x, Y3_y, Y3_z, i, j, k, time + 80 / 120 * dt);

                grid.u[i][j][k] = Y3_x[i][j][k] +
                                  90.0 / 120.0 * dt * f_y3[0] -
                                  50.0 / 120.0 * dt * f_y2[0];
                grid.v[i][j][k] = Y3_y[i][j][k] +
                                  90.0 / 120.0 * dt * f_y3[1] -
                                  50.0 / 120.0 * dt * f_y2[1];
                grid.w[i][j][k] = Y3_z[i][j][k] +
                                  90.0 / 120.0 * dt * f_y3[2] -
                                  50.0 / 120.0 * dt * f_y2[2];
                // std::cout << "u[" << i << "," << j << "," << k << "] = " << grid.u[i * NY * NZ + j * NZ + k] << std::endl;
            }
        }
    }
}

void IcoNS::apply_boundary_conditions(double time)
{
    // compute boundary conditions
    for (size_t i = 0; i < NX; i++)
    {
        for (size_t j = 0; j < NY; j++)
        {
            // bottom face -> k = 0 -> sin(k*dz) = 0, cos(k*dz) = 1
            grid.u[i][j][0] = 0.0;
            grid.v[i][j][0] = 0.0;
            grid.w[i][j][0] = 2 * std::cos(i * dx) * std::cos(j * dy) * std::sin(time);
            // top face -> k*dz = lz
            grid.u[i][j][NZ - 1] = std::sin(i * dx) * std::cos(j * dy) * std::sin(lz) * std::sin(time);
            grid.v[i][j][NZ - 1] = std::cos(i * dx) * std::sin(j * dy) * std::sin(lz) * std::sin(time);
            grid.w[i][j][NZ - 1] = 2 * std::cos(i * dx) * std::cos(j * dy) * std::cos(lz) * std::sin(time);
        }
    }

    for (size_t i = 0; i < NX; i++)
    {
        for (size_t k = 0; k < NZ; k++)
        {
            // front face -> j = 0 -> sin(j*dy) = 0, cos(j*dy) = 1
            grid.u[i][0][k] = std::sin(i * dx) * std::sin(k * dz) * std::sin(time);
            grid.v[i][0][k] = 0.0;
            grid.w[i][0][k] = 2 * std::cos(i * dx) * std::cos(k * dz) * std::sin(time);
            // back face -> j*dy = ly
            grid.u[i][NY - 1][k] = std::sin(i * dx) * std::cos(ly) * std::sin(k * dz) * std::sin(time);
            grid.v[i][NY - 1][k] = std::cos(i * dx) * std::sin(ly) * std::sin(k * dz) * std::sin(time);
            grid.w[i][NY - 1][k] = 2 * std::cos(i * dx) * std::cos(ly) * std::cos(k * dz) * std::sin(time);
        }
    }

    for (size_t j = 0; j < NY; j++)
    {
        for (size_t k = 0; k < NZ; k++)
        {
            // left face -> i = 0 -> sin(i*dx) = 0, cos(i*dx) = 1
            grid.u[0][j][k] = 0.0;
            grid.v[0][j][k] = std::sin(j * dy) * std::sin(k * dz) * std::sin(time);
            grid.w[0][j][k] = 2 * std::cos(j * dy) * std::cos(k * dz) * std::sin(time);
            // right face -> i*dx = lx
            grid.u[NX - 1][j][k] = std::sin(lx) * std::cos(j * dy) * std::sin(k * dz) * std::sin(time);
            grid.v[NX - 1][j][k] = std::cos(lx) * std::sin(j * dy) * std::sin(k * dz) * std::sin(time);
            grid.w[NX - 1][j][k] = 2 * std::cos(lx) * std::cos(j * dy) * std::cos(k * dz) * std::sin(time);
        }
    }
}

double IcoNS::L2_error(const double t)
{
    double error = 0.0;

    error += error_comp_X(t);
    error += error_comp_Y(t);
    error += error_comp_Z(t);

    return sqrt(error);
}

double IcoNS::error_comp_X(const double t)
{
    double error = 0.0;

    // first slice (left face)
    {
        error += ((grid.u[0][0][0] - exact_solution.value_x(0, 0, 0, t)) *
                  (grid.u[0][0][0] - exact_solution.value_x(0, 0, 0, t)) *
                  dx * dy * dz / 8);

        for (size_t k = 1; k < NZ; k++)
        {
            error += ((grid.u[0][0][k] - exact_solution.value_x(0, 0, k, t)) *
                      (grid.u[0][0][k] - exact_solution.value_x(0, 0, k, t)) *
                      dx * dy * dz / 4);
        }

        error += ((grid.u[0][0][NZ] - exact_solution.value_x(0, 0, NZ, t)) *
                  (grid.u[0][0][NZ] - exact_solution.value_x(0, 0, NZ, t)) *
                  dx * dy * dz / 8);

        for (size_t j = 1; j < NY; j++)
        {
            error += ((grid.u[0][j][0] - exact_solution.value_x(0, j, 0, t)) *
                      (grid.u[0][j][0] - exact_solution.value_x(0, j, 0, t)) *
                      dx * dy * dz / 4);
            for (size_t k = 1; k < NZ; k++)
            {
                error += ((grid.u[0][j][k] - exact_solution.value_x(0, j, k, t)) *
                          (grid.u[0][j][k] - exact_solution.value_x(0, j, k, t)) *
                          dx * dy * dz / 2);
            }
            error += ((grid.u[0][j][NZ] - exact_solution.value_x(0, j, NZ, t)) *
                      (grid.u[0][j][NZ] - exact_solution.value_x(0, j, NZ, t)) *
                      dx * dy * dz / 4);
        }

        error += ((grid.u[0][NY][0] - exact_solution.value_x(0, NY, 0, t)) *
                  (grid.u[0][NY][0] - exact_solution.value_x(0, NY, 0, t)) *
                  dx * dy * dz / 8);

        for (size_t k = 1; k < NZ; k++)
        {
            error += ((grid.u[0][NY][k] - exact_solution.value_x(0, NY, k, t)) *
                      (grid.u[0][NY][k] - exact_solution.value_x(0, NY, k, t)) *
                      dx * dy * dz / 4);
        }

        error += ((grid.u[0][NY][NZ] - exact_solution.value_x(0, NY, NZ, t)) *
                  (grid.u[0][NY][NZ] - exact_solution.value_x(0, NY, NZ, t)) *
                  dx * dy * dz / 8);
    }

    // middle slices
    {
        for (size_t i = 1; i < NX - 1; i++)
        {
            error += ((grid.u[i][0][0] - exact_solution.value_x(i, 0, 0, t)) *
                      (grid.u[i][0][0] - exact_solution.value_x(i, 0, 0, t)) *
                      dx * dy * dz / 4);

            for (size_t k = 1; k < NZ; k++)
            {
                error += ((grid.u[i][0][k] - exact_solution.value_x(i, 0, k, t)) *
                          (grid.u[i][0][k] - exact_solution.value_x(i, 0, k, t)) *
                          dx * dy * dz / 2);
            }
            error += ((grid.u[i][0][NZ] - exact_solution.value_x(i, 0, NZ, t)) *
                      (grid.u[i][0][NZ] - exact_solution.value_x(i, 0, NZ, t)) *
                      dx * dy * dz / 4);

            for (size_t j = 1; j < NY; j++)
            {
                error += ((grid.u[i][j][0] - exact_solution.value_x(i, j, 0, t)) *
                          (grid.u[i][j][0] - exact_solution.value_x(i, j, 0, t)) *
                          dx * dy * dz / 2);

                for (size_t k = 1; k < NZ; k++)
                {
                    error += ((grid.u[i][j][k] - exact_solution.value_x(i, j, k, t)) *
                              (grid.u[i][j][k] - exact_solution.value_x(i, j, k, t)) *
                              dx * dy * dz);
                }

                error += ((grid.u[i][j][NZ] - exact_solution.value_x(i, j, NZ, t)) *
                          (grid.u[i][j][NZ] - exact_solution.value_x(i, j, NZ, t)) *
                          dx * dy * dz / 2);
            }

            error += ((grid.u[i][NY][NZ] - exact_solution.value_x(i, NY, NZ, t)) *
                      (grid.u[i][NY][NZ] - exact_solution.value_x(i, NY, NZ, t)) *
                      dx * dy * dz / 4);

            for (size_t k = 1; k < NZ; k++)
            {
                error += ((grid.u[i][NY][k] - exact_solution.value_x(i, NY, k, t)) *
                          (grid.u[i][NY][k] - exact_solution.value_x(i, NY, k, t)) *
                          dx * dy * dz / 2);
            }

            error += ((grid.u[i][NY][NZ] - exact_solution.value_x(i, NY, NZ, t)) *
                      (grid.u[i][NY][NZ] - exact_solution.value_x(i, NY, NZ, t)) *
                      dx * dy * dz / 4);
        }
    }

    // last slice (right face)
    {
        error += ((grid.u[NX - 1][0][0] - exact_solution.value_x((NX - 1), 0, 0, t)) *
                  (grid.u[NX - 1][0][0] - exact_solution.value_x((NX - 1), 0, 0, t)) *
                  dx * dy * dz / 8);

        for (size_t k = 1; k < NZ; k++)
        {
            error += ((grid.u[NX - 1][0][k] - exact_solution.value_x((NX - 1), 0, k, t)) *
                      (grid.u[NX - 1][0][k] - exact_solution.value_x((NX - 1), 0, k, t)) *
                      dx * dy * dz / 4);
        }

        error += ((grid.u[NX - 1][0][NZ] - exact_solution.value_x((NX - 1), 0, NZ, t)) *
                  (grid.u[NX - 1][0][NZ] - exact_solution.value_x((NX - 1), 0, NZ, t)) *
                  dx * dy * dz / 8);

        for (size_t j = 1; j < NY; j++)
        {
            error += ((grid.u[NX - 1][j][0] - exact_solution.value_x((NX - 1), j, 0, t)) *
                      (grid.u[NX - 1][j][0] - exact_solution.value_x((NX - 1), j, 0, t)) *
                      dx * dy * dz / 4);
            for (size_t k = 1; k < NZ; k++)
            {
                error += ((grid.u[NX - 1][j][k] - exact_solution.value_x((NX - 1), j, k, t)) *
                          (grid.u[NX - 1][j][k] - exact_solution.value_x((NX - 1), j, k, t)) *
                          dx * dy * dz / 2);
            }
            error += ((grid.u[NX - 1][j][NZ] - exact_solution.value_x((NX - 1), j, NZ, t)) *
                      (grid.u[NX - 1][j][NZ] - exact_solution.value_x((NX - 1), j, NZ, t)) *
                      dx * dy * dz / 4);
        }

        error += ((grid.u[NX - 1][NY][0] - exact_solution.value_x((NX - 1), NY, 0, t)) *
                  (grid.u[NX - 1][NY][0] - exact_solution.value_x((NX - 1), NY, 0, t)) *
                  dx * dy * dz / 8);

        for (size_t k = 1; k < NZ; k++)
        {
            error += ((grid.u[NX - 1][NY][k] - exact_solution.value_x((NX - 1), NY, k, t)) *
                      (grid.u[NX - 1][NY][k] - exact_solution.value_x((NX - 1), NY, k, t)) *
                      dx * dy * dz / 4);
        }

        error += ((grid.u[NX - 1][NY][NZ] - exact_solution.value_x((NX - 1), NY, NZ, t)) *
                  (grid.u[NX - 1][NY][NZ] - exact_solution.value_x((NX - 1), NY, NZ, t)) *
                  dx * dy * dz / 8);
    }
    return error;
}

double IcoNS::error_comp_Y(const double t)
{
    double error = 0.0;

    // first slice (left face)
    {
        error += ((grid.v[0][0][0] - exact_solution.value_y(0, 0, 0, t)) *
                  (grid.v[0][0][0] - exact_solution.value_y(0, 0, 0, t)) *
                  dx * dy * dz / 8);

        for (size_t k = 1; k < NZ; k++)
        {
            error += ((grid.v[0][0][k] - exact_solution.value_y(0, 0, k, t)) *
                      (grid.v[0][0][k] - exact_solution.value_y(0, 0, k, t)) *
                      dx * dy * dz / 4);
        }

        error += ((grid.v[0][0][NZ] - exact_solution.value_y(0, 0, NZ, t)) *
                  (grid.v[0][0][NZ] - exact_solution.value_y(0, 0, NZ, t)) *
                  dx * dy * dz / 8);

        for (size_t j = 1; j < NY - 1; j++)
        {
            error += ((grid.v[0][j][0] - exact_solution.value_y(0, j, 0, t)) *
                      (grid.v[0][j][0] - exact_solution.value_y(0, j, 0, t)) *
                      dx * dy * dz / 4);
            for (size_t k = 1; k < NZ; k++)
            {
                error += ((grid.v[0][j][k] - exact_solution.value_y(0, j, k, t)) *
                          (grid.v[0][j][k] - exact_solution.value_y(0, j, k, t)) *
                          dx * dy * dz / 2);
            }
            error += ((grid.v[0][j][NZ] - exact_solution.value_y(0, j, NZ, t)) *
                      (grid.v[0][j][NZ] - exact_solution.value_y(0, j, NZ, t)) *
                      dx * dy * dz / 4);
        }

        error += ((grid.v[0][NY - 1][0] - exact_solution.value_y(0, (NY - 1), 0, t)) *
                  (grid.v[0][NY - 1][0] - exact_solution.value_y(0, (NY - 1), 0, t)) *
                  dx * dy * dz / 8);

        for (size_t k = 1; k < NZ; k++)
        {
            error += ((grid.v[0][NY - 1][k] - exact_solution.value_y(0, (NY - 1), k, t)) *
                      (grid.v[0][NY - 1][k] - exact_solution.value_y(0, (NY - 1), k, t)) *
                      dx * dy * dz / 4);
        }

        error += ((grid.v[0][NY - 1][NZ] - exact_solution.value_y(0, (NY - 1), NZ, t)) *
                  (grid.v[0][NY - 1][NZ] - exact_solution.value_y(0, (NY - 1), NZ, t)) *
                  dx * dy * dz / 8);
    }

    // middle slices
    {
        for (size_t i = 1; i < NX; i++)
        {
            error += ((grid.v[i][0][0] - exact_solution.value_y(i, 0, 0, t)) *
                      (grid.v[i][0][0] - exact_solution.value_y(i, 0, 0, t)) *
                      dx * dy * dz / 4);

            for (size_t k = 1; k < NZ; k++)
            {
                error += ((grid.v[i][0][k] - exact_solution.value_y(i, 0, k, t)) *
                          (grid.v[i][0][k] - exact_solution.value_y(i, 0, k, t)) *
                          dx * dy * dz / 2);
            }
            error += ((grid.v[i][0][NZ] - exact_solution.value_y(i, 0, NZ, t)) *
                      (grid.v[i][0][NZ] - exact_solution.value_y(i, 0, NZ, t)) *
                      dx * dy * dz / 4);

            for (size_t j = 1; j < NY - 1; j++)
            {
                error += ((grid.v[i][j][0] - exact_solution.value_y(i, j, 0, t)) *
                          (grid.v[i][j][0] - exact_solution.value_y(i, j, 0, t)) *
                          dx * dy * dz / 2);

                for (size_t k = 1; k < NZ; k++)
                {
                    error += ((grid.v[i][j][k] - exact_solution.value_y(i, j, k, t)) *
                              (grid.v[i][j][k] - exact_solution.value_y(i, j, k, t)) *
                              dx * dy * dz);
                }

                error += ((grid.v[i][j][NZ] - exact_solution.value_y(i, j, NZ, t)) *
                          (grid.v[i][j][NZ] - exact_solution.value_y(i, j, NZ, t)) *
                          dx * dy * dz / 2);
            }

            error += ((grid.v[i][NY - 1][NZ] - exact_solution.value_y(i, (NY - 1), NZ, t)) *
                      (grid.v[i][NY - 1][NZ] - exact_solution.value_y(i, (NY - 1), NZ, t)) *
                      dx * dy * dz / 4);

            for (size_t k = 1; k < NZ; k++)
            {
                error += ((grid.v[i][NY - 1][k] - exact_solution.value_y(i, (NY - 1), k, t)) *
                          (grid.v[i][NY - 1][k] - exact_solution.value_y(i, (NY - 1), k, t)) *
                          dx * dy * dz / 2);
            }

            error += ((grid.v[i][NY - 1][NZ] - exact_solution.value_y(i, (NY - 1), NZ, t)) *
                      (grid.v[i][NY - 1][NZ] - exact_solution.value_y(i, (NY - 1), NZ, t)) *
                      dx * dy * dz / 4);
        }
    }

    // last slice (right face)
    {
        error += ((grid.v[NX][0][0] - exact_solution.value_y(NX, 0, 0, t)) *
                  (grid.v[NX][0][0] - exact_solution.value_y(NX, 0, 0, t)) *
                  dx * dy * dz / 8);

        for (size_t k = 1; k < NZ; k++)
        {
            error += ((grid.v[NX][0][k] - exact_solution.value_y(NX, 0, k, t)) *
                      (grid.v[NX][0][k] - exact_solution.value_y(NX, 0, k, t)) *
                      dx * dy * dz / 4);
        }

        error += ((grid.v[NX][0][NZ] - exact_solution.value_y(NX, 0, NZ, t)) *
                  (grid.v[NX][0][NZ] - exact_solution.value_y(NX, 0, NZ, t)) *
                  dx * dy * dz / 8);

        for (size_t j = 1; j < NY - 1; j++)
        {
            error += ((grid.v[NX][j][0] - exact_solution.value_y(NX, j, 0, t)) *
                      (grid.v[NX][j][0] - exact_solution.value_y(NX, j, 0, t)) *
                      dx * dy * dz / 4);
            for (size_t k = 1; k < NZ; k++)
            {
                error += ((grid.v[NX][j][k] - exact_solution.value_y(NX, j, k, t)) *
                          (grid.v[NX][j][k] - exact_solution.value_y(NX, j, k, t)) *
                          dx * dy * dz / 2);
            }
            error += ((grid.v[NX][j][NZ] - exact_solution.value_y(NX, j, NZ, t)) *
                      (grid.v[NX][j][NZ] - exact_solution.value_y(NX, j, NZ, t)) *
                      dx * dy * dz / 4);
        }

        error += ((grid.v[NX][NY - 1][0] - exact_solution.value_y(NX, (NY - 1), 0, t)) *
                  (grid.v[NX][NY - 1][0] - exact_solution.value_y(NX, (NY - 1), 0, t)) *
                  dx * dy * dz / 8);

        for (size_t k = 1; k < NZ; k++)
        {
            error += ((grid.v[NX][NY - 1][k] - exact_solution.value_y(NX, (NY - 1), k, t)) *
                      (grid.v[NX][NY - 1][k] - exact_solution.value_y(NX, (NY - 1), k, t)) *
                      dx * dy * dz / 4);
        }

        error += ((grid.v[NX][NY - 1][NZ] - exact_solution.value_y(NX, (NY - 1), NZ, t)) *
                  (grid.v[NX][NY - 1][NZ] - exact_solution.value_y(NX, (NY - 1), NZ, t)) *
                  dx * dy * dz / 8);
    }
    return error;
}

double IcoNS::error_comp_Z(const double t)
{
    double error = 0.0;
    // first slice (left face)
    {
        error += ((grid.w[0][0][0] - exact_solution.value_z(0, 0, 0, t)) *
                  (grid.w[0][0][0] - exact_solution.value_z(0, 0, 0, t)) *
                  dx * dy * dz / 8);

        for (size_t k = 1; k < NZ - 1; k++)
        {
            error += ((grid.w[0][0][k] - exact_solution.value_z(0, 0, k, t)) *
                      (grid.w[0][0][k] - exact_solution.value_z(0, 0, k, t)) *
                      dx * dy * dz / 4);
        }

        error += ((grid.w[0][0][NZ - 1] - exact_solution.value_z(0, 0, NZ - 1, t)) *
                  (grid.w[0][0][NZ - 1] - exact_solution.value_z(0, 0, NZ - 1, t)) *
                  dx * dy * dz / 8);

        for (size_t j = 1; j < NY; j++)
        {
            error += ((grid.w[0][j][0] - exact_solution.value_z(0, j, 0, t)) *
                      (grid.w[0][j][0] - exact_solution.value_z(0, j, 0, t)) *
                      dx * dy * dz / 4);
            for (size_t k = 1; k < NZ - 1; k++)
            {
                error += ((grid.w[0][j][k] - exact_solution.value_z(0, j, k, t)) *
                          (grid.w[0][j][k] - exact_solution.value_z(0, j, k, t)) *
                          dx * dy * dz / 2);
            }
            error += ((grid.w[0][j][NZ - 1] - exact_solution.value_z(0, j, NZ - 1, t)) *
                      (grid.w[0][j][NZ - 1] - exact_solution.value_z(0, j, NZ - 1, t)) *
                      dx * dy * dz / 4);
        }

        error += ((grid.w[0][NY][0] - exact_solution.value_z(0, NY, 0, t)) *
                  (grid.w[0][NY][0] - exact_solution.value_z(0, NY, 0, t)) *
                  dx * dy * dz / 8);

        for (size_t k = 1; k < NZ - 1; k++)
        {
            error += ((grid.w[0][NY][k] - exact_solution.value_z(0, NY, k, t)) *
                      (grid.w[0][NY][k] - exact_solution.value_z(0, NY, k, t)) *
                      dx * dy * dz / 4);
        }

        error += ((grid.w[0][NY][NZ - 1] - exact_solution.value_z(0, NY, NZ - 1, t)) *
                  (grid.w[0][NY][NZ - 1] - exact_solution.value_z(0, NY, NZ - 1, t)) *
                  dx * dy * dz / 8);
    }

    // middle slices
    {
        for (size_t i = 1; i < NX; i++)
        {
            error += ((grid.w[i][0][0] - exact_solution.value_z(i, 0, 0, t)) *
                      (grid.w[i][0][0] - exact_solution.value_z(i, 0, 0, t)) *
                      dx * dy * dz / 4);

            for (size_t k = 1; k < NZ - 1; k++)
            {
                error += ((grid.w[i][0][k] - exact_solution.value_z(i, 0, k, t)) *
                          (grid.w[i][0][k] - exact_solution.value_z(i, 0, k, t)) *
                          dx * dy * dz / 2);
            }
            error += ((grid.w[i][0][NZ - 1] - exact_solution.value_z(i, 0, NZ - 1, t)) *
                      (grid.w[i][0][NZ - 1] - exact_solution.value_z(i, 0, NZ - 1, t)) *
                      dx * dy * dz / 4);

            for (size_t j = 1; j < NY; j++)
            {
                error += ((grid.w[i][j][0] - exact_solution.value_z(i, j, 0, t)) *
                          (grid.w[i][j][0] - exact_solution.value_z(i, j, 0, t)) *
                          dx * dy * dz / 2);

                for (size_t k = 1; k < NZ - 1; k++)
                {
                    error += ((grid.w[i][j][k] - exact_solution.value_z(i, j, k, t)) *
                              (grid.w[i][j][k] - exact_solution.value_z(i, j, k, t)) *
                              dx * dy * dz);
                }

                error += ((grid.w[i][j][NZ - 1] - exact_solution.value_z(i, j, (NZ - 1), t)) *
                          (grid.w[i][j][NZ - 1] - exact_solution.value_z(i, j, (NZ - 1), t)) *
                          dx * dy * dz / 2);
            }

            error += ((grid.w[i][NY][NZ - 1] - exact_solution.value_z(i, NY, (NZ - 1), t)) *
                      (grid.w[i][NY][NZ - 1] - exact_solution.value_z(i, NY, (NZ - 1), t)) *
                      dx * dy * dz / 4);

            for (size_t k = 1; k < NZ - 1; k++)
            {
                error += ((grid.w[i][NY][k] - exact_solution.value_z(i, NY, k, t)) *
                          (grid.w[i][NY][k] - exact_solution.value_z(i, NY, k, t)) *
                          dx * dy * dz / 2);
            }

            error += ((grid.w[i][NY][NZ - 1] - exact_solution.value_z(i, NY, NZ - 1, t)) *
                      (grid.w[i][NY][NZ - 1] - exact_solution.value_z(i, NY, NZ - 1, t)) *
                      dx * dy * dz / 4);
        }
    }

    // last slice (right face)
    {
        error += ((grid.w[NX][0][0] - exact_solution.value_z(NX, 0, 0, t)) *
                  (grid.w[NX][0][0] - exact_solution.value_z(NX, 0, 0, t)) *
                  dx * dy * dz / 8);

        for (size_t k = 1; k < NZ - 1; k++)
        {
            error += ((grid.w[NX][0][k] - exact_solution.value_z(NX, 0, k, t)) *
                      (grid.w[NX][0][k] - exact_solution.value_z(NX, 0, k, t)) *
                      dx * dy * dz / 4);
        }

        error += ((grid.w[NX][0][NZ - 1] - exact_solution.value_z(NX, 0, NZ - 1, t)) *
                  (grid.w[NX][0][NZ - 1] - exact_solution.value_z(NX, 0, NZ - 1, t)) *
                  dx * dy * dz / 8);

        for (size_t j = 1; j < NY; j++)
        {
            error += ((grid.w[NX][j][0] - exact_solution.value_z(NX, j, 0, t)) *
                      (grid.w[NX][j][0] - exact_solution.value_z(NX, j, 0, t)) *
                      dx * dy * dz / 4);
            for (size_t k = 1; k < NZ - 1; k++)
            {
                error += ((grid.w[NX][j][k] - exact_solution.value_z(NX, j, k, t)) *
                          (grid.w[NX][j][k] - exact_solution.value_z(NX, j, k, t)) *
                          dx * dy * dz / 2);
            }
            error += ((grid.w[NX][j][NZ - 1] - exact_solution.value_z(NX, j, NZ - 1, t)) *
                      (grid.w[NX][j][NZ - 1] - exact_solution.value_z(NX, j, NZ - 1, t)) *
                      dx * dy * dz / 4);
        }

        error += ((grid.w[NX][NY][0] - exact_solution.value_z(NX, NY, 0, t)) *
                  (grid.w[NX][NY][0] - exact_solution.value_z(NX, NY, 0, t)) *
                  dx * dy * dz / 8);

        for (size_t k = 1; k < NZ - 1; k++)
        {
            error += ((grid.w[NX][NY][k] - exact_solution.value_z(NX, NY, k, t)) *
                      (grid.w[NX][NY][k] - exact_solution.value_z(NX, NY, k, t)) *
                      dx * dy * dz / 4);
        }

        error += ((grid.w[NX][NY][NZ - 1] - exact_solution.value_z(NX, NY, NZ - 1, t)) *
                  (grid.w[NX][NY][NZ - 1] - exact_solution.value_z(NX, NY, NZ - 1, t)) *
                  dx * dy * dz / 8);
    }

    return error;
}

void IcoNS::output()
{
}
