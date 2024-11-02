
#include <cmath>

#include "core.hpp"
#include <math.h>

#include <fstream>
#include <string>
#include <memory>

void IcoNS::preprocessing(/*std::string &input_file*/)

{
    // boundary
    auto u_func = std::make_shared<Dirichlet>([&](double x, double y, double z, double t)
                                              { return std::sin((x + 0.5) * dx) * std::cos(y * dy) * std::sin(z * dz) * std::sin(t); });
    auto v_func = std::make_shared<Dirichlet>([&](double x, double y, double z, double t)
                                              { return std::cos(x * dx) * std::sin((y + 0.5) * dy) * std::sin(z * dz) * std::sin(t); });
    auto w_func = std::make_shared<Dirichlet>([&](double x, double y, double z, double t)
                                              { return 2 * (std::cos(x * dx) * std::cos(y * dy) * std::cos((z + 0.5) * dz) * std::sin(t)); });

    for (size_t i = 0; i < 6 /*nfaces*/; i++)
    {
        boundary.addFunction(U, u_func);
        boundary.addFunction(V, v_func);
        boundary.addFunction(W, w_func);
    }
}

void IcoNS::solve()
{
    double time = 0.0;
    int i = 0;

    std::ofstream error_log("../resources/error.log");

    while (time < T)
    {   
        /*Check::Confront(grid,exact_solution,time,U);
        int p;
        std::cin >> p;*/
        boundary.update_boundary(grid.u, grid.v, grid.w, time);

        // csv file w/ "," delimiter: time step, iter, L2_error
        error_log << time << "," << i << "," << L2_error(time) << std::endl;
        solve_time_step(time);
        // output();
        time += dt;
        i++;
        
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
        error += ((grid.u[0] - exact_solution.value_x(0.5, 0, 0, t)) *
                  (grid.u[0] - exact_solution.value_x(0.5, 0, 0, t)) *
                  dx * dy * dz / 8);

        for (size_t k = 1; k < nz; k++)
        {
            error += ((grid.u[k] - exact_solution.value_x(0.5, 0, k, t)) *
                      (grid.u[k] - exact_solution.value_x(0.5, 0, k, t)) *
                      dx * dy * dz / 4);
        }

        error += ((grid.u[nz] - exact_solution.value_x(0.5, 0, nz, t)) *
                  (grid.u[nz] - exact_solution.value_x(0.5, 0, nz, t)) *
                  dx * dy * dz / 8);

        for (size_t j = 1; j < ny; j++)
        {
            error += ((grid.u[j * (nz + 1)] - exact_solution.value_x(0.5, j, 0, t)) *
                      (grid.u[j * (nz + 1)] - exact_solution.value_x(0.5, j, 0, t)) *
                      dx * dy * dz / 4);
            for (size_t k = 1; k < nz; k++)
            {
                error += ((grid.u[j * (nz + 1) + k] - exact_solution.value_x(0.5, j, k, t)) *
                          (grid.u[j * (nz + 1) + k] - exact_solution.value_x(0.5, j, k, t)) *
                          dx * dy * dz / 2);
            }
            error += ((grid.u[j * (nz + 1) + nz] - exact_solution.value_x(0.5, j, nz, t)) *
                      (grid.u[j * (nz + 1) + nz] - exact_solution.value_x(0.5, j, nz, t)) *
                      dx * dy * dz / 4);
        }

        error += ((grid.u[ny * (nz + 1)] - exact_solution.value_x(0.5, ny, 0, t)) *
                  (grid.u[ny * (nz + 1)] - exact_solution.value_x(0.5, ny, 0, t)) *
                  dx * dy * dz / 8);

        for (size_t k = 1; k < nz; k++)
        {
            error += ((grid.u[ny * (nz + 1) + k] - exact_solution.value_x(0.5, ny, k, t)) *
                      (grid.u[ny * (nz + 1) + k] - exact_solution.value_x(0.5, ny, k, t)) *
                      dx * dy * dz / 4);
        }

        error += ((grid.u[ny * (nz + 1) + nz] - exact_solution.value_x(0.5, ny, nz, t)) *
                  (grid.u[ny * (nz + 1) + nz] - exact_solution.value_x(0.5, ny, nz, t)) *
                  dx * dy * dz / 8);
    }

    // middle slices
    {
        for (size_t i = 1; i < nx - 1; i++)
        {
            error += ((grid.u[i * (ny + 1) * (nz + 1)] - exact_solution.value_x(i + 0.5, 0, 0, t)) *
                      (grid.u[i * (ny + 1) * (nz + 1)] - exact_solution.value_x(i + 0.5, 0, 0, t)) *
                      dx * dy * dz / 4);

            for (size_t k = 1; k < nz; k++)
            {
                error += ((grid.u[i * (ny + 1) * (nz + 1) + k] - exact_solution.value_x(i + 0.5, 0, k, t)) *
                          (grid.u[i * (ny + 1) * (nz + 1) + k] - exact_solution.value_x(i + 0.5, 0, k, t)) *
                          dx * dy * dz / 2);
            }
            error += ((grid.u[i * (ny + 1) * (nz + 1) + nz] - exact_solution.value_x(i + 0.5, 0, nz, t)) *
                      (grid.u[i * (ny + 1) * (nz + 1) + nz] - exact_solution.value_x(i + 0.5, 0, nz, t)) *
                      dx * dy * dz / 4);

            for (size_t j = 1; j < ny; j++)
            {
                error += ((grid.u[i * (ny + 1) * (nz + 1) + j * (nz + 1)] - exact_solution.value_x(i + 0.5, j, 0, t)) *
                          (grid.u[i * (ny + 1) * (nz + 1) + j * (nz + 1)] - exact_solution.value_x(i + 0.5, j, 0, t)) *
                          dx * dy * dz / 2);

                for (size_t k = 1; k < nz; k++)
                {
                    error += ((grid.u[i * (ny + 1) * (nz + 1) + j * (nz + 1) + k] - exact_solution.value_x(i + 0.5, j, k, t)) *
                              (grid.u[i * (ny + 1) * (nz + 1) + j * (nz + 1) + k] - exact_solution.value_x(i + 0.5, j, k, t)) *
                              dx * dy * dz);
                }

                error += ((grid.u[i * (ny + 1) * (nz + 1) + j * (nz + 1) + nz] - exact_solution.value_x(i + 0.5, j, nz, t)) *
                          (grid.u[i * (ny + 1) * (nz + 1) + j * (nz + 1) + nz] - exact_solution.value_x(i + 0.5, j, nz, t)) *
                          dx * dy * dz / 2);
            }

            error += ((grid.u[i * (ny + 1) * (nz + 1) + ny * (nz + 1)] - exact_solution.value_x(i + 0.5, ny, 0, t)) *
                      (grid.u[i * (ny + 1) * (nz + 1) + ny * (nz + 1)] - exact_solution.value_x(i + 0.5, ny, 0, t)) *
                      dx * dy * dz / 4);

            for (size_t k = 1; k < nz; k++)
            {
                error += ((grid.u[i * (ny + 1) * (nz + 1) + ny * (nz + 1) + k] - exact_solution.value_x(i + 0.5, ny, k, t)) *
                          (grid.u[i * (ny + 1) * (nz + 1) + ny * (nz + 1) + k] - exact_solution.value_x(i + 0.5, ny, k, t)) *
                          dx * dy * dz / 2);
            }

            error += ((grid.u[i * (ny + 1) * (nz + 1) + ny * (nz + 1) + nz] - exact_solution.value_x(i + 0.5, ny, nz, t)) *
                      (grid.u[i * (ny + 1) * (nz + 1) + ny * (nz + 1) + nz] - exact_solution.value_x(i + 0.5, ny, nz, t)) *
                      dx * dy * dz / 4);
        }
    }

    // last slice (right face)
    {
        error += ((grid.u[(nx - 1) * (ny + 1) * (nz + 1)] - exact_solution.value_x((nx - 0.5), 0, 0, t)) *
                  (grid.u[(nx - 1) * (ny + 1) * (nz + 1)] - exact_solution.value_x((nx - 0.5), 0, 0, t)) *
                  dx * dy * dz / 8);

        for (size_t k = 1; k < nz; k++)
        {
            error += ((grid.u[(nx - 1) * (ny + 1) * (nz + 1) + k] - exact_solution.value_x((nx - 0.5), 0, k, t)) *
                      (grid.u[(nx - 1) * (ny + 1) * (nz + 1) + k] - exact_solution.value_x((nx - 0.5), 0, k, t)) *
                      dx * dy * dz / 4);
        }

        error += ((grid.u[(nx - 1) * (ny + 1) * (nz + 1) + nz] - exact_solution.value_x((nx - 0.5), 0, nz, t)) *
                  (grid.u[(nx - 1) * (ny + 1) * (nz + 1) + nz] - exact_solution.value_x((nx - 0.5), 0, nz, t)) *
                  dx * dy * dz / 8);

        for (size_t j = 1; j < ny; j++)
        {
            error += ((grid.u[(nx - 1) * (ny + 1) * (nz + 1) + j * (nz + 1)] - exact_solution.value_x((nx - 0.5), j, 0, t)) *
                      (grid.u[(nx - 1) * (ny + 1) * (nz + 1) + j * (nz + 1)] - exact_solution.value_x((nx - 0.5), j, 0, t)) *
                      dx * dy * dz / 4);
            for (size_t k = 1; k < nz; k++)
            {
                error += ((grid.u[(nx - 1) * (ny + 1) * (nz + 1) + j * (nz + 1) + k] - exact_solution.value_x((nx - 0.5), j, k, t)) *
                          (grid.u[(nx - 1) * (ny + 1) * (nz + 1) + j * (nz + 1) + k] - exact_solution.value_x((nx - 0.5), j, k, t)) *
                          dx * dy * dz / 2);
            }
            error += ((grid.u[(nx - 1) * (ny + 1) * (nz + 1) + j * (nz + 1) + nz] - exact_solution.value_x((nx - 0.5), j, nz, t)) *
                      (grid.u[(nx - 1) * (ny + 1) * (nz + 1) + j * (nz + 1) + nz] - exact_solution.value_x((nx - 0.5), j, nz, t)) *
                      dx * dy * dz / 4);
        }

        error += ((grid.u[(nx - 1) * (ny + 1) * (nz + 1) + ny * (nz + 1)] - exact_solution.value_x((nx - 0.5), ny, 0, t)) *
                  (grid.u[(nx - 1) * (ny + 1) * (nz + 1) + ny * (nz + 1)] - exact_solution.value_x((nx - 0.5), ny, 0, t)) *
                  dx * dy * dz / 8);

        for (size_t k = 1; k < nz; k++)
        {
            error += ((grid.u[(nx - 1) * (ny + 1) * (nz + 1) + ny * (nz + 1) + k] - exact_solution.value_x((nx - 0.5), ny, k, t)) *
                      (grid.u[(nx - 1) * (ny + 1) * (nz + 1) + ny * (nz + 1) + k] - exact_solution.value_x((nx - 0.5), ny, k, t)) *
                      dx * dy * dz / 4);
        }

        error += ((grid.u[(nx - 1) * (ny + 1) * (nz + 1) + ny * (nz + 1) + nz] - exact_solution.value_x((nx - 0.5), ny, nz, t)) *
                  (grid.u[(nx - 1) * (ny + 1) * (nz + 1) + ny * (nz + 1) + nz] - exact_solution.value_x((nx - 0.5), ny, nz, t)) *
                  dx * dy * dz / 8);
    }
    return error;
}

double IcoNS::error_comp_Y(const double t)
{
    double error = 0.0;

    // first slice (left face)
    {
        error += ((grid.v[0] - exact_solution.value_y(0, 0.5, 0, t)) *
                  (grid.v[0] - exact_solution.value_y(0, 0.5, 0, t)) *
                  dx * dy * dz / 8);

        for (size_t k = 1; k < nz; k++)
        {
            error += ((grid.v[k] - exact_solution.value_y(0, 0.5, k, t)) *
                      (grid.v[k] - exact_solution.value_y(0, 0.5, k, t)) *
                      dx * dy * dz / 4);
        }

        error += ((grid.v[nz] - exact_solution.value_y(0, 0.5, nz, t)) *
                  (grid.v[nz] - exact_solution.value_y(0, 0.5, nz, t)) *
                  dx * dy * dz / 8);

        for (size_t j = 1; j < ny - 1; j++)
        {
            error += ((grid.v[j * (nz + 1)] - exact_solution.value_y(0, j + 0.5, 0, t)) *
                      (grid.v[j * (nz + 1)] - exact_solution.value_y(0, j + 0.5, 0, t)) *
                      dx * dy * dz / 4);
            for (size_t k = 1; k < nz; k++)
            {
                error += ((grid.v[j * (nz + 1) + k] - exact_solution.value_y(0, j + 0.5, k, t)) *
                          (grid.v[j * (nz + 1) + k] - exact_solution.value_y(0, j + 0.5, k, t)) *
                          dx * dy * dz / 2);
            }
            error += ((grid.v[j * (nz + 1) + nz] - exact_solution.value_y(0, j + 0.5, nz, t)) *
                      (grid.v[j * (nz + 1) + nz] - exact_solution.value_y(0, j + 0.5, nz, t)) *
                      dx * dy * dz / 4);
        }

        error += ((grid.v[(ny - 1) * (nz + 1)] - exact_solution.value_y(0, (ny - 0.5), 0, t)) *
                  (grid.v[(ny - 1) * (nz + 1)] - exact_solution.value_y(0, (ny - 0.5), 0, t)) *
                  dx * dy * dz / 8);

        for (size_t k = 1; k < nz; k++)
        {
            error += ((grid.v[(ny - 1) * (nz + 1) + k] - exact_solution.value_y(0, (ny - 0.5), k, t)) *
                      (grid.v[(ny - 1) * (nz + 1) + k] - exact_solution.value_y(0, (ny - 0.5), k, t)) *
                      dx * dy * dz / 4);
        }

        error += ((grid.v[(ny - 1) * (nz + 1) + nz] - exact_solution.value_y(0, (ny - 0.5), nz, t)) *
                  (grid.v[(ny - 1) * (nz + 1) + nz] - exact_solution.value_y(0, (ny - 0.5), nz, t)) *
                  dx * dy * dz / 8);
    }

    // middle slices
    {
        for (size_t i = 1; i < nx; i++)
        {
            error += ((grid.v[i * ny * (nz + 1)] - exact_solution.value_y(i, 0.5, 0, t)) *
                      (grid.v[i * ny * (nz + 1)] - exact_solution.value_y(i, 0.5, 0, t)) *
                      dx * dy * dz / 4);

            for (size_t k = 1; k < nz; k++)
            {
                error += ((grid.v[i * ny * (nz + 1) + k] - exact_solution.value_y(i, 0.5, k, t)) *
                          (grid.v[i * ny * (nz + 1) + k] - exact_solution.value_y(i, 0.5, k, t)) *
                          dx * dy * dz / 2);
            }
            error += ((grid.v[i * ny * (nz + 1) + nz] - exact_solution.value_y(i, 0.5, nz, t)) *
                      (grid.v[i * ny * (nz + 1) + nz] - exact_solution.value_y(i, 0.5, nz, t)) *
                      dx * dy * dz / 4);

            for (size_t j = 1; j < ny - 1; j++)
            {
                error += ((grid.v[i * ny * (nz + 1) + j * (nz + 1)] - exact_solution.value_y(i, j + 0.5, 0, t)) *
                          (grid.v[i * ny * (nz + 1) + j * (nz + 1)] - exact_solution.value_y(i, j + 0.5, 0, t)) *
                          dx * dy * dz / 2);

                for (size_t k = 1; k < nz; k++)
                {
                    error += ((grid.v[i * ny * (nz + 1) + j * (nz + 1) + k] - exact_solution.value_y(i, j + 0.5, k, t)) *
                              (grid.v[i * ny * (nz + 1) + j * (nz + 1) + k] - exact_solution.value_y(i, j + 0.5, k, t)) *
                              dx * dy * dz);
                }

                error += ((grid.v[i * ny * (nz + 1) + j * (nz + 1) + nz] - exact_solution.value_y(i, j + 0.5, nz, t)) *
                          (grid.v[i * ny * (nz + 1) + j * (nz + 1) + nz] - exact_solution.value_y(i, j + 0.5, nz, t)) *
                          dx * dy * dz / 2);
            }

            error += ((grid.v[i * ny * (nz + 1) + (ny - 1) * (nz + 1)] - exact_solution.value_y(i, (ny - 0.5), 0, t)) *
                      (grid.v[i * ny * (nz + 1) + (ny - 1) * (nz + 1)] - exact_solution.value_y(i, (ny - 0.5), 0, t)) *
                      dx * dy * dz / 4);

            for (size_t k = 1; k < nz; k++)
            {
                error += ((grid.v[i * ny * (nz + 1) + (ny - 1) * (nz + 1) + k] - exact_solution.value_y(i, (ny - 0.5), k, t)) *
                          (grid.v[i * ny * (nz + 1) + (ny - 1) * (nz + 1) + k] - exact_solution.value_y(i, (ny - 0.5), k, t)) *
                          dx * dy * dz / 2);
            }

            error += ((grid.v[i * ny * (nz + 1) + (ny - 1) * (nz + 1) + nz] - exact_solution.value_y(i, (ny - 0.5), nz, t)) *
                      (grid.v[i * ny * (nz + 1) + (ny - 1) * (nz + 1) + nz] - exact_solution.value_y(i, (ny - 0.5), nz, t)) *
                      dx * dy * dz / 4);
        }
    }

    // last slice (right face)
    {
        error += ((grid.v[nx * ny * (nz + 1)] - exact_solution.value_y(nx, 0.5, 0, t)) *
                  (grid.v[nx * ny * (nz + 1)] - exact_solution.value_y(nx, 0.5, 0, t)) *
                  dx * dy * dz / 8);

        for (size_t k = 1; k < nz; k++)
        {
            error += ((grid.v[nx * ny * (nz + 1) + k] - exact_solution.value_y(nx, 0.5, k, t)) *
                      (grid.v[nx * ny * (nz + 1) + k] - exact_solution.value_y(nx, 0.5, k, t)) *
                      dx * dy * dz / 4);
        }

        error += ((grid.v[nx * ny * (nz + 1) + nz] - exact_solution.value_y(nx, 0.5, nz, t)) *
                  (grid.v[nx * ny * (nz + 1) + nz] - exact_solution.value_y(nx, 0.5, nz, t)) *
                  dx * dy * dz / 8);

        for (size_t j = 1; j < ny - 1; j++)
        {
            error += ((grid.v[nx * ny * (nz + 1) + j * (nz + 1)] - exact_solution.value_y(nx, j + 0.5, 0, t)) *
                      (grid.v[nx * ny * (nz + 1) + j * (nz + 1)] - exact_solution.value_y(nx, j + 0.5, 0, t)) *
                      dx * dy * dz / 4);
            for (size_t k = 1; k < nz; k++)
            {
                error += ((grid.v[nx * ny * (nz + 1) + j * (nz + 1) + k] - exact_solution.value_y(nx, j + 0.5, k, t)) *
                          (grid.v[nx * ny * (nz + 1) + j * (nz + 1) + k] - exact_solution.value_y(nx, j + 0.5, k, t)) *
                          dx * dy * dz / 2);
            }
            error += ((grid.v[nx * ny * (nz + 1) + j * (nz + 1) + nz] - exact_solution.value_y(nx, j + 0.5, nz, t)) *
                      (grid.v[nx * ny * (nz + 1) + j * (nz + 1) + nz] - exact_solution.value_y(nx, j + 0.5, nz, t)) *
                      dx * dy * dz / 4);
        }

        error += ((grid.v[nx * ny * (nz + 1) + (ny - 1) * (nz + 1)] - exact_solution.value_y(nx, (ny - 0.5), 0, t)) *
                  (grid.v[nx * ny * (nz + 1) + (ny - 1) * (nz + 1)] - exact_solution.value_y(nx, (ny - 0.5), 0, t)) *
                  dx * dy * dz / 8);

        for (size_t k = 1; k < nz; k++)
        {
            error += ((grid.v[nx * ny * (nz + 1) + (ny - 1) * (nz + 1) + k] - exact_solution.value_y(nx, (ny - 0.5), k, t)) *
                      (grid.v[nx * ny * (nz + 1) + (ny - 1) * (nz + 1) + k] - exact_solution.value_y(nx, (ny - 0.5), k, t)) *
                      dx * dy * dz / 4);
        }

        error += ((grid.v[nx * ny * (nz + 1) + (ny - 1) * (nz + 1) + nz] - exact_solution.value_y(nx, (ny - 0.5), nz, t)) *
                  (grid.v[nx * ny * (nz + 1) + (ny - 1) * (nz + 1) + nz] - exact_solution.value_y(nx, (ny - 0.5), nz, t)) *
                  dx * dy * dz / 8);
    }
    return error;
}

double IcoNS::error_comp_Z(const double t)
{
    double error = 0.0;
    // first slice (left face)
    {
        error += ((grid.w[0] - exact_solution.value_z(0, 0, 0.5, t)) *
                  (grid.w[0] - exact_solution.value_z(0, 0, 0.5, t)) *
                  dx * dy * dz / 8);

        for (size_t k = 1; k < nz - 1; k++)
        {
            error += ((grid.w[k] - exact_solution.value_z(0, 0, k + 0.5, t)) *
                      (grid.w[k] - exact_solution.value_z(0, 0, k + 0.5, t)) *
                      dx * dy * dz / 4);
        }

        error += ((grid.w[nz - 1] - exact_solution.value_z(0, 0, nz - 0.5, t)) *
                  (grid.w[nz - 1] - exact_solution.value_z(0, 0, nz - 0.5, t)) *
                  dx * dy * dz / 8);

        for (size_t j = 1; j < ny; j++)
        {
            error += ((grid.w[j * nz] - exact_solution.value_z(0, j, 0.5, t)) *
                      (grid.w[j * nz] - exact_solution.value_z(0, j, 0.5, t)) *
                      dx * dy * dz / 4);
            for (size_t k = 1; k < nz - 1; k++)
            {
                error += ((grid.w[j * nz + k] - exact_solution.value_z(0, j, k + 0.5, t)) *
                          (grid.w[j * nz + k] - exact_solution.value_z(0, j, k + 0.5, t)) *
                          dx * dy * dz / 2);
            }
            error += ((grid.w[j * nz + nz - 1] - exact_solution.value_z(0, j, nz - 0.5, t)) *
                      (grid.w[j * nz + nz - 1] - exact_solution.value_z(0, j, nz - 0.5, t)) *
                      dx * dy * dz / 4);
        }

        error += ((grid.w[ny * nz] - exact_solution.value_z(0, ny, 0.5, t)) *
                  (grid.w[ny * nz] - exact_solution.value_z(0, ny, 0.5, t)) *
                  dx * dy * dz / 8);

        for (size_t k = 1; k < nz - 1; k++)
        {
            error += ((grid.w[ny * nz + k] - exact_solution.value_z(0, ny, k + 0.5, t)) *
                      (grid.w[ny * nz + k] - exact_solution.value_z(0, ny, k + 0.5, t)) *
                      dx * dy * dz / 4);
        }

        error += ((grid.w[ny * nz + nz - 1] - exact_solution.value_z(0, ny, nz - 0.5, t)) *
                  (grid.w[ny * nz + nz - 1] - exact_solution.value_z(0, ny, nz - 0.5, t)) *
                  dx * dy * dz / 8);
    }

    // middle slices
    {
        for (size_t i = 1; i < nx; i++)
        {
            error += ((grid.w[i * (ny + 1) * nz] - exact_solution.value_z(i, 0, 0.5, t)) *
                      (grid.w[i * (ny + 1) * nz] - exact_solution.value_z(i, 0, 0.5, t)) *
                      dx * dy * dz / 4);

            for (size_t k = 1; k < nz - 1; k++)
            {
                error += ((grid.w[i * (ny + 1) * nz + k] - exact_solution.value_z(i, 0, k + 0.5, t)) *
                          (grid.w[i * (ny + 1) * nz + k] - exact_solution.value_z(i, 0, k + 0.5, t)) *
                          dx * dy * dz / 2);
            }
            error += ((grid.w[i * (ny + 1) * nz + nz - 1] - exact_solution.value_z(i, 0, nz - 0.5, t)) *
                      (grid.w[i * (ny + 1) * nz + nz - 1] - exact_solution.value_z(i, 0, nz - 0.5, t)) *
                      dx * dy * dz / 4);

            for (size_t j = 1; j < ny; j++)
            {
                error += ((grid.w[i * (ny + 1) * nz + j * nz] - exact_solution.value_z(i, j, 0.5, t)) *
                          (grid.w[i * (ny + 1) * nz + j * nz] - exact_solution.value_z(i, j, 0.5, t)) *
                          dx * dy * dz / 2);

                for (size_t k = 1; k < nz - 1; k++)
                {
                    error += ((grid.w[i * (ny + 1) * nz + j * nz + k] - exact_solution.value_z(i, j, k + 0.5, t)) *
                              (grid.w[i * (ny + 1) * nz + j * nz + k] - exact_solution.value_z(i, j, k + 0.5, t)) *
                              dx * dy * dz);
                }

                error += ((grid.w[i * (ny + 1) * nz + j * nz + nz - 1] - exact_solution.value_z(i, j, (nz - 0.5), t)) *
                          (grid.w[i * (ny + 1) * nz + j * nz + nz - 1] - exact_solution.value_z(i, j, (nz - 0.5), t)) *
                          dx * dy * dz / 2);
            }

            error += ((grid.w[i * (ny + 1) * nz + ny * nz] - exact_solution.value_z(i, ny, 0.5, t)) *
                      (grid.w[i * (ny + 1) * nz + ny * nz] - exact_solution.value_z(i, ny, 0.5, t)) *
                      dx * dy * dz / 4);

            for (size_t k = 1; k < nz - 1; k++)
            {
                error += ((grid.w[i * (ny + 1) * nz + ny * nz + k] - exact_solution.value_z(i, ny, k + 0.5, t)) *
                          (grid.w[i * (ny + 1) * nz + ny * nz + k] - exact_solution.value_z(i, ny, k + 0.5, t)) *
                          dx * dy * dz / 2);
            }

            error += ((grid.w[i * (ny + 1) * nz + ny * nz + nz - 1] - exact_solution.value_z(i, ny, nz - 0.5, t)) *
                      (grid.w[i * (ny + 1) * nz + ny * nz + nz - 1] - exact_solution.value_z(i, ny, nz - 0.5, t)) *
                      dx * dy * dz / 4);
        }
    }

    // last slice (right face)
    {
        error += ((grid.w[nx * (ny + 1) * nz] - exact_solution.value_z(nx, 0, 0.5, t)) *
                  (grid.w[nx * (ny + 1) * nz] - exact_solution.value_z(nx, 0, 0.5, t)) *
                  dx * dy * dz / 8);

        for (size_t k = 1; k < nz - 1; k++)
        {
            error += ((grid.w[nx * (ny + 1) * nz + k] - exact_solution.value_z(nx, 0, k + 0.5, t)) *
                      (grid.w[nx * (ny + 1) * nz + k] - exact_solution.value_z(nx, 0, k + 0.5, t)) *
                      dx * dy * dz / 4);
        }

        error += ((grid.w[nx * (ny + 1) * nz + nz - 1] - exact_solution.value_z(nx, 0, nz - 0.5, t)) *
                  (grid.w[nx * (ny + 1) * nz + nz - 1] - exact_solution.value_z(nx, 0, nz - 0.5, t)) *
                  dx * dy * dz / 8);

        for (size_t j = 1; j < ny; j++)
        {
            error += ((grid.w[nx * (ny + 1) * nz + j * nz] - exact_solution.value_z(nx, j, 0.5, t)) *
                      (grid.w[nx * (ny + 1) * nz + j * nz] - exact_solution.value_z(nx, j, 0.5, t)) *
                      dx * dy * dz / 4);
            for (size_t k = 1; k < nz - 1; k++)
            {
                error += ((grid.w[nx * (ny + 1) * nz + j * nz + k] - exact_solution.value_z(nx, j, k + 0.5, t)) *
                          (grid.w[nx * (ny + 1) * nz + j * nz + k] - exact_solution.value_z(nx, j, k + 0.5, t)) *
                          dx * dy * dz / 2);
            }
            error += ((grid.w[nx * (ny + 1) * nz + j * nz + nz - 1] - exact_solution.value_z(nx, j, nz - 0.5, t)) *
                      (grid.w[nx * (ny + 1) * nz + j * nz + nz - 1] - exact_solution.value_z(nx, j, nz - 0.5, t)) *
                      dx * dy * dz / 4);
        }

        error += ((grid.w[nx * (ny + 1) * nz + ny * nz] - exact_solution.value_z(nx, ny, 0.5, t)) *
                  (grid.w[nx * (ny + 1) * nz + ny * nz] - exact_solution.value_z(nx, ny, 0.5, t)) *
                  dx * dy * dz / 8);

        for (size_t k = 1; k < nz - 1; k++)
        {
            error += ((grid.w[nx * (ny + 1) * nz + ny * nz + k] - exact_solution.value_z(nx, ny, k + 0.5, t)) *
                      (grid.w[nx * (ny + 1) * nz + ny * nz + k] - exact_solution.value_z(nx, ny, k + 0.5, t)) *
                      dx * dy * dz / 4);
        }

        error += ((grid.w[nx * (ny + 1) * nz + ny * nz + nz - 1] - exact_solution.value_z(nx, ny, nz - 0.5, t)) *
                  (grid.w[nx * (ny + 1) * nz + ny * nz + nz - 1] - exact_solution.value_z(nx, ny, nz - 0.5, t)) *
                  dx * dy * dz / 8);
    }

    return error;
}

void IcoNS::output()
{
}
