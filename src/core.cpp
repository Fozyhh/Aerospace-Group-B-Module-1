#include "../includes/core.hpp"
#include <math.h>
#include <string>

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
}

void IcoNS::solve()
{
    preprocessing();
    double time = 0.0;
    int i = 0;

    while (time < T)
    {
        solve_time_step();
        time += dt;

        output();
        std::cout << "time step: " << i << " ";
        std::cout << " Error: " << L2_error() << std::endl;
        i++;
    }
}

std::vector<double> IcoNS::functionF(const std::vector<double> &u, const std::vector<double> &v,
                                     const std::vector<double> &w, size_t i, size_t j, size_t k)
{
    std::vector<double> f(3);
    size_t l = i * ny * nz + j * nz + k;

    f[0] = -(u[l] * (u[l + nz * ny] - u[l - nz * ny]) / (2.0 * dx) +
             (v[l + nz * ny] + v[l] + v[l - nz] + v[l + nz * ny - nz]) / 4.0 * (u[l + nz] - u[l - nz]) / (2.0 * dy) +
             (w[l + nz * ny] + w[l] + w[l + nz * ny - 1] + w[l - 1]) / 4.0 * (u[l + 1] - u[l - 1]) / (2.0 * dz)) +
           1 / Re * ((u[l + ny * nz] - 2 * u[l] + u[l - ny * nz]) / (dx * dx) + (u[l + nz] - 2 * u[l] + u[l - nz]) / (dy * dy) + (u[l + 1] - 2 * u[l] + u[l - 1]) / (dz * dz));

    f[1] = -((u[l] + u[l + nz] + u[l - ny * nz] + u[l - ny * nz + nz]) / 4.0 * (v[l + ny * nz] - v[l - ny * nz] / (2.0 * dx)) +
             v[l] * (v[l + nz] - v[l - nz]) / (2.0 * dy) +
             (w[l] + w[l - 1] + w[l + nz] + w[l + nz - 1]) / 4.0 + (v[l + 1] - v[l - 1]) / (2.0 * dz)) +
           1 / Re * ((v[l + nz * ny] - 2.0 * v[l] + v[l - nz * ny]) / (dx * dx) + (v[l + nz] - 2.0 * v[l] + v[l - nz]) / (dy * dy) + (v[l + 1] - 2.0 * v[l] + v[l - 1]) / (dz * dz));

    f[2] = -((u[l] + u[l - ny * nz] + u[l + 1] + u[l - nz * ny + 1]) / 4.0 * (w[l + nz * ny] - w[l - nz * ny]) / (2.0 * dx) +
             (v[l + 1] + v[l - nz + 1] + v[l] + v[l - nz]) / 4.0 * (w[l + nz] - w[l - nz]) / (2.0 * dy) +
             w[l] * (w[l + 1] - w[l - 1]) / (2.0 * dz)) +
           1 / Re * ((w[l + nz * ny] - 2.0 * w[l] + w[l - nz * ny]) / (dx * dx) + (w[l + nz] - 2.0 * w[l] + w[l - nz]) / (dy * dy) + (w[l + 1] - 2.0 * w[l] + w[l - 1]) / (dz * dz));

    return f;
}

void IcoNS::solve_time_step()
{
    std::vector<double> f(3);
    std::vector<double> f_y2(3);
    std::vector<double> f_y3(3);
    std::vector<double> Y2_x(nx * ny * nz);
    std::vector<double> Y2_y(nx * ny * nz);
    std::vector<double> Y2_z(nx * ny * nz);
    std::vector<double> Y3_x(nx * ny * nz);
    std::vector<double> Y3_y(nx * ny * nz);
    std::vector<double> Y3_z(nx * ny * nz);

    for (size_t i = 1; i < nx - 1; i++)
    {
        for (size_t j = 1; j < ny - 1; j++)
        {
            for (size_t k = 1; k < nz - 1; k++)
            {
                f = functionF(grid.u, grid.v, grid.w, i, j, k);
                // solve the momentum equations -> TODO later also with pressure
                Y2_x[i * ny * nz + j * nz + k] = grid.u[i * ny * nz + j * nz + k] +
                                                 64 / 120 * dt * f[0];
                Y2_y[i * ny * nz + j * nz + k] = grid.v[i * ny * nz + j * nz + k] +
                                                 64 / 120 * dt * f[1];
                Y2_z[i * ny * nz + j * nz + k] = grid.w[i * ny * nz + j * nz + k] +
                                                 64 / 120 * dt * f[2];

                f_y2 = functionF(Y2_x, Y2_y, Y2_z, i, j, k);

                Y3_x[i * ny * nz + j * nz + k] = Y2_x[i * ny * nz + j * nz + k] +
                                                 50 / 120 * dt * f_y2[0] -
                                                 34 / 120 * dt * f[0];
                Y3_y[i * ny * nz + j * nz + k] = Y2_y[i * ny * nz + j * nz + k] +
                                                 50 / 120 * dt * f_y2[1] -
                                                 34 / 120 * dt * f[1];
                Y3_z[i * ny * nz + j * nz + k] = Y2_z[i * ny * nz + j * nz + k] +
                                                 50 / 120 * dt * f_y2[2] -
                                                 34 / 120 * dt * f[2];

                // update the grid values.
                f_y3 = functionF(Y3_x, Y3_y, Y3_z, i, j, k);

                grid.u[i * ny * nz + j * nz + k] = Y3_x[i * ny * nz + j * nz + k] +
                                                   90 / 120 * dt * f_y3[0] -
                                                   50 / 120 * dt * f_y2[0];
                grid.v[i * ny * nz + j * nz + k] = Y3_y[i * ny * nz + j * nz + k] +
                                                   90 / 120 * dt * f_y3[1] -
                                                   50 / 120 * dt * f_y2[1];
                grid.w[i * ny * nz + j * nz + k] = Y3_z[i * ny * nz + j * nz + k] +
                                                   90 / 120 * dt * f_y3[2] -
                                                   50 / 120 * dt * f_y2[2];
            }
        }
    }
}
/*
double IcoNS::L2_error()
{

    double sum = 0.0;

    std::vector<double> wx(grid.nx, 1.0);
    std::vector<double> wy(grid.ny, 1.0);
    std::vector<double> wz(grid.nz, 1.0);

    wx[0] = 0.5;
    wx[grid.nx - 1] = 0.5;
    wy[0] = 0.5;
    wy[grid.ny - 1] = 0.5;
    wz[0] = 0.5;
    wz[grid.nz - 1] = 0.5;

    size_t l;

    for (size_t i = 0; i < grid.nx; i++)
    {
        for (size_t j = 0; j < grid.ny; j++)
        {
            for (size_t k = 0; k < grid.nz; k++)
            {
                l = i * grid.ny * grid.nz + j * grid.nz + k;
                sum += (wx[i] * wy[j] * wz[k] * grid.u[l] * grid.u[l]) * dx * dy * dz;
            }
        }
    }

    for (size_t i = 0; i < grid.nx; i++)
    {
        for (size_t j = 0; j < grid.ny; j++)
        {
            for (size_t k = 0; k < grid.nz; k++)
            {
                l = i * grid.ny * grid.nz + j * grid.nz + k;
                sum += (wx[i] * wy[j] * wz[k] * grid.v[l] * grid.v[l]) * dy * dx * dz;
            }
        }
    }

    for (size_t i = 0; i < grid.nx; i++)
    {
        for (size_t j = 0; j < grid.ny; j++)
        {
            for (size_t k = 0; k < grid.nz; k++)
            {
                l = i * grid.ny * grid.nz + j * grid.nz + k;
                sum += (wx[i] * wy[j] * wz[k] * grid.w[l] * grid.w[l]) * dz * dx * dz;
            }
        }
    }
    return sqrt(sum);
}
*/

double IcoNS::error_comp_X()
{
    double error = 0.0;

    // first slice (left face)
    {
        error += ((grid.u[0] - exact_solution.value_x(0)) *
                  (grid.u[0] - exact_solution.value_x(0)) *
                  dx * dy * dz / 8);

        for (size_t k = 1; k < nz - 1; k++)
        {
            error += ((grid.u[k] - exact_solution.value_x(k)) *
                      (grid.u[k] - exact_solution.value_x(k)) *
                      dx * dy * dz / 4);
        }

        error += ((grid.u[nz - 1] - exact_solution.value_x(nz - 1)) *
                  (grid.u[nz - 1] - exact_solution.value_x(nz - 1)) *
                  dx * dy * dz / 8);

        for (size_t j = 1; j < ny - 1; j++)
        {
            error += ((grid.u[j * nz] - exact_solution.value_x(j * nz)) *
                      (grid.u[j * nz] - exact_solution.value_x(j * nz)) *
                      dx * dy * dz / 4);
            for (size_t k = 1; k < nz - 1; k++)
            {
                error += ((grid.u[j * nz + k] - exact_solution.value_x(j * nz + k)) *
                          (grid.u[j * nz + k] - exact_solution.value_x(j * nz + k)) *
                          dx * dy * dz / 2);
            }
            error += ((grid.u[j * nz + nz - 1] - exact_solution.value_x(j * nz + nz - 1)) *
                      (grid.u[j * nz + nz - 1] - exact_solution.value_x(j * nz + nz - 1)) *
                      dx * dy * dz / 4);
        }

        error += ((grid.u[(ny - 1) * nz] - exact_solution.value_x((ny - 1) * nz)) *
                  (grid.u[(ny - 1) * nz] - exact_solution.value_x((ny - 1) * nz)) *
                  dx * dy * dz / 8);

        for (size_t k = 1; k < nz - 1; k++)
        {
            error += ((grid.u[(ny - 1) * nz + k] - exact_solution.value_x((ny - 1) * nz + k)) *
                      (grid.u[(ny - 1) * nz + k] - exact_solution.value_x((ny - 1) * nz + k)) *
                      dx * dy * dz / 4);
        }

        error += ((grid.u[(ny - 1) * nz + nz - 1] - exact_solution.value_x((ny - 1) * nz + nz - 1)) *
                  (grid.u[(ny - 1) * nz + nz - 1] - exact_solution.value_x((ny - 1) * nz + nz - 1)) *
                  dx * dy * dz / 8);
    }

    // middle slices
    {
        for (size_t i = 1; i < nx - 1; i++)
        {
            error += ((grid.u[i * ny * nz] - exact_solution.value_x(i * ny * nz)) *
                      (grid.u[i * ny * nz] - exact_solution.value_x(i * ny * nz)) *
                      dx * dy * dz / 4);

            for (size_t k = 1; k < nz - 1; k++)
            {
                error += ((grid.u[i * ny * nz + k] - exact_solution.value_x(i * ny * nz + k)) *
                          (grid.u[i * ny * nz + k] - exact_solution.value_x(i * ny * nz + k)) *
                          dx * dy * dz / 2);
            }
            error += ((grid.u[i * ny * nz + nz - 1] - exact_solution.value_x(i * ny * nz + nz - 1)) *
                      (grid.u[i * ny * nz + nz - 1] - exact_solution.value_x(i * ny * nz + nz - 1)) *
                      dx * dy * dz / 4);

            for (size_t j = 1; j < ny - 1; j++)
            {
                error += ((grid.u[i * ny * nz + j * nz] - exact_solution.value_x(i * ny * nz + j * nz)) *
                          (grid.u[i * ny * nz + j * nz] - exact_solution.value_x(i * ny * nz + j * nz)) *
                          dx * dy * dz / 2);

                for (size_t k = 1; k < nz - 1; k++)
                {
                    error += ((grid.u[i * ny * nz + j * nz + k] - exact_solution.value_x(i * ny * nz + j * nz + k)) *
                              (grid.u[i * ny * nz + j * nz + k] - exact_solution.value_x(i * ny * nz + j * nz + k)) *
                              dx * dy * dz);
                }

                error += ((grid.u[i * ny * nz + j * nz + nz - 1] - exact_solution.value_x(i * ny * nz + j * nz + nz - 1)) *
                          (grid.u[i * ny * nz + j * nz + nz - 1] - exact_solution.value_x(i * ny * nz + j * nz + nz - 1)) *
                          dx * dy * dz / 2);
            }

            error += ((grid.u[i * ny * nz + (ny - 1) * nz] - exact_solution.value_x(i * ny * nz + (ny - 1) * nz)) *
                      (grid.u[i * ny * nz + (ny - 1) * nz] - exact_solution.value_x(i * ny * nz + (ny - 1) * nz)) *
                      dx * dy * dz / 4);

            for (size_t k = 1; k < nz - 1; k++)
            {
                error += ((grid.u[i * ny * nz + (ny - 1) * nz + k] - exact_solution.value_x(i * ny * nz + (ny - 1) * nz + k)) *
                          (grid.u[i * ny * nz + (ny - 1) * nz + k] - exact_solution.value_x(i * ny * nz + (ny - 1) * nz + k)) *
                          dx * dy * dz / 2);
            }

            error += ((grid.u[i * ny * nz + (ny - 1) * nz + nz - 1] - exact_solution.value_x(i * ny * nz + (ny - 1) * nz + nz - 1)) *
                      (grid.u[i * ny * nz + (ny - 1) * nz + nz - 1] - exact_solution.value_x(i * ny * nz + (ny - 1) * nz + nz - 1)) *
                      dx * dy * dz / 4);
        }
    }

    // last slice (right face)
    {
        error += ((grid.u[(nx - 1) * ny * nz] - exact_solution.value_x((nx - 1) * ny * nz)) *
                  (grid.u[(nx - 1) * ny * nz] - exact_solution.value_x((nx - 1) * ny * nz)) *
                  dx * dy * dz / 8);

        for (size_t k = 1; k < nz - 1; k++)
        {
            error += ((grid.u[(nx - 1) * ny * nz + k] - exact_solution.value_x((nx - 1) * ny * nz + k)) *
                      (grid.u[(nx - 1) * ny * nz + k] - exact_solution.value_x((nx - 1) * ny * nz + k)) *
                      dx * dy * dz / 4);
        }

        error += ((grid.u[(nx - 1) * ny * nz + nz - 1] - exact_solution.value_x((nx - 1) * ny * nz + nz - 1)) *
                  (grid.u[(nx - 1) * ny * nz + nz - 1] - exact_solution.value_x((nx - 1) * ny * nz + nz - 1)) *
                  dx * dy * dz / 8);

        for (size_t j = 1; j < ny - 1; j++)
        {
            error += ((grid.u[(nx - 1) * ny * nz + j * nz] - exact_solution.value_x((nx - 1) * ny * nz + j * nz)) *
                      (grid.u[(nx - 1) * ny * nz + j * nz] - exact_solution.value_x((nx - 1) * ny * nz + j * nz)) *
                      dx * dy * dz / 4);
            for (size_t k = 1; k < nz - 1; k++)
            {
                error += ((grid.u[(nx - 1) * ny * nz + j * nz + k] - exact_solution.value_x((nx - 1) * ny * nz + j * nz + k)) *
                          (grid.u[(nx - 1) * ny * nz + j * nz + k] - exact_solution.value_x((nx - 1) * ny * nz + j * nz + k)) *
                          dx * dy * dz / 2);
            }
            error += ((grid.u[(nx - 1) * ny * nz + j * nz + nz - 1] - exact_solution.value_x((nx - 1) * ny * nz + j * nz + nz - 1)) *
                      (grid.u[(nx - 1) * ny * nz + j * nz + nz - 1] - exact_solution.value_x((nx - 1) * ny * nz + j * nz + nz - 1)) *
                      dx * dy * dz / 4);
        }

        error += ((grid.u[(nx - 1) * ny * nz + (ny - 1) * nz] - exact_solution.value_x((nx - 1) * ny * nz + (ny - 1) * nz)) *
                  (grid.u[(nx - 1) * ny * nz + (ny - 1) * nz] - exact_solution.value_x((nx - 1) * ny * nz + (ny - 1) * nz)) *
                  dx * dy * dz / 8);

        for (size_t k = 1; k < nz - 1; k++)
        {
            error += ((grid.u[(nx - 1) * ny * nz + (ny - 1) * nz + k] - exact_solution.value_x((nx - 1) * ny * nz + (ny - 1) * nz + k)) *
                      (grid.u[(nx - 1) * ny * nz + (ny - 1) * nz + k] - exact_solution.value_x((nx - 1) * ny * nz + (ny - 1) * nz + k)) *
                      dx * dy * dz / 4);
        }

        error += ((grid.u[(nx - 1) * ny * nz + (ny - 1) * nz + nz - 1] - exact_solution.value_x((nx - 1) * ny * nz + (ny - 1) * nz + nz - 1)) *
                  (grid.u[(nx - 1) * ny * nz + (ny - 1) * nz + nz - 1] - exact_solution.value_x((nx - 1) * ny * nz + (ny - 1) * nz + nz - 1)) *
                  dx * dy * dz / 8);
    }
    return error;
}

double IcoNS::error_comp_Y()
{
    double error = 0.0;
    // first slice (left face)
    {
        error += ((grid.w[0] - exact_solution.value_z(0)) *
                  (grid.w[0] - exact_solution.value_z(0)) *
                  dx * dy * dz / 8);

        for (size_t k = 1; k < nz - 1; k++)
        {
            error += ((grid.w[k] - exact_solution.value_z(k)) *
                      (grid.w[k] - exact_solution.value_z(k)) *
                      dx * dy * dz / 4);
        }

        error += ((grid.w[nz - 1] - exact_solution.value_z(nz - 1)) *
                  (grid.w[nz - 1] - exact_solution.value_z(nz - 1)) *
                  dx * dy * dz / 8);

        for (size_t j = 1; j < ny - 1; j++)
        {
            error += ((grid.w[j * nz] - exact_solution.value_z(j * nz)) *
                      (grid.w[j * nz] - exact_solution.value_z(j * nz)) *
                      dx * dy * dz / 4);
            for (size_t k = 1; k < nz - 1; k++)
            {
                error += ((grid.w[j * nz + k] - exact_solution.value_z(j * nz + k)) *
                          (grid.w[j * nz + k] - exact_solution.value_z(j * nz + k)) *
                          dx * dy * dz / 2);
            }
            error += ((grid.w[j * nz + nz - 1] - exact_solution.value_z(j * nz + nz - 1)) *
                      (grid.w[j * nz + nz - 1] - exact_solution.value_z(j * nz + nz - 1)) *
                      dx * dy * dz / 4);
        }

        error += ((grid.w[(ny - 1) * nz] - exact_solution.value_z((ny - 1) * nz)) *
                  (grid.w[(ny - 1) * nz] - exact_solution.value_z((ny - 1) * nz)) *
                  dx * dy * dz / 8);

        for (size_t k = 1; k < nz - 1; k++)
        {
            error += ((grid.w[(ny - 1) * nz + k] - exact_solution.value_z((ny - 1) * nz + k)) *
                      (grid.w[(ny - 1) * nz + k] - exact_solution.value_z((ny - 1) * nz + k)) *
                      dx * dy * dz / 4);
        }

        error += ((grid.w[(ny - 1) * nz + nz - 1] - exact_solution.value_z((ny - 1) * nz + nz - 1)) *
                  (grid.w[(ny - 1) * nz + nz - 1] - exact_solution.value_z((ny - 1) * nz + nz - 1)) *
                  dx * dy * dz / 8);
    }

    // middle slices
    {
        for (size_t i = 1; i < nx - 1; i++)
        {
            error += ((grid.w[i * ny * nz] - exact_solution.value_z(i * ny * nz)) *
                      (grid.w[i * ny * nz] - exact_solution.value_z(i * ny * nz)) *
                      dx * dy * dz / 4);

            for (size_t k = 1; k < nz - 1; k++)
            {
                error += ((grid.w[i * ny * nz + k] - exact_solution.value_z(i * ny * nz + k)) *
                          (grid.w[i * ny * nz + k] - exact_solution.value_z(i * ny * nz + k)) *
                          dx * dy * dz / 2);
            }
            error += ((grid.w[i * ny * nz + nz - 1] - exact_solution.value_z(i * ny * nz + nz - 1)) *
                      (grid.w[i * ny * nz + nz - 1] - exact_solution.value_z(i * ny * nz + nz - 1)) *
                      dx * dy * dz / 4);

            for (size_t j = 1; j < ny - 1; j++)
            {
                error += ((grid.w[i * ny * nz + j * nz] - exact_solution.value_z(i * ny * nz + j * nz)) *
                          (grid.w[i * ny * nz + j * nz] - exact_solution.value_z(i * ny * nz + j * nz)) *
                          dx * dy * dz / 2);

                for (size_t k = 1; k < nz - 1; k++)
                {
                    error += ((grid.w[i * ny * nz + j * nz + k] - exact_solution.value_z(i * ny * nz + j * nz + k)) *
                              (grid.w[i * ny * nz + j * nz + k] - exact_solution.value_z(i * ny * nz + j * nz + k)) *
                              dx * dy * dz);
                }

                error += ((grid.w[i * ny * nz + j * nz + nz - 1] - exact_solution.value_z(i * ny * nz + j * nz + nz - 1)) *
                          (grid.w[i * ny * nz + j * nz + nz - 1] - exact_solution.value_z(i * ny * nz + j * nz + nz - 1)) *
                          dx * dy * dz / 2);
            }

            error += ((grid.w[i * ny * nz + (ny - 1) * nz] - exact_solution.value_z(i * ny * nz + (ny - 1) * nz)) *
                      (grid.w[i * ny * nz + (ny - 1) * nz] - exact_solution.value_z(i * ny * nz + (ny - 1) * nz)) *
                      dx * dy * dz / 4);

            for (size_t k = 1; k < nz - 1; k++)
            {
                error += ((grid.w[i * ny * nz + (ny - 1) * nz + k] - exact_solution.value_z(i * ny * nz + (ny - 1) * nz + k)) *
                          (grid.w[i * ny * nz + (ny - 1) * nz + k] - exact_solution.value_z(i * ny * nz + (ny - 1) * nz + k)) *
                          dx * dy * dz / 2);
            }

            error += ((grid.w[i * ny * nz + (ny - 1) * nz + nz - 1] - exact_solution.value_z(i * ny * nz + (ny - 1) * nz + nz - 1)) *
                      (grid.w[i * ny * nz + (ny - 1) * nz + nz - 1] - exact_solution.value_z(i * ny * nz + (ny - 1) * nz + nz - 1)) *
                      dx * dy * dz / 4);
        }
    }

    // last slice (right face)
    {
        error += ((grid.w[(nx - 1) * ny * nz] - exact_solution.value_z((nx - 1) * ny * nz)) *
                  (grid.w[(nx - 1) * ny * nz] - exact_solution.value_z((nx - 1) * ny * nz)) *
                  dx * dy * dz / 8);

        for (size_t k = 1; k < nz - 1; k++)
        {
            error += ((grid.w[(nx - 1) * ny * nz + k] - exact_solution.value_z((nx - 1) * ny * nz + k)) *
                      (grid.w[(nx - 1) * ny * nz + k] - exact_solution.value_z((nx - 1) * ny * nz + k)) *
                      dx * dy * dz / 4);
        }

        error += ((grid.w[(nx - 1) * ny * nz + nz - 1] - exact_solution.value_z((nx - 1) * ny * nz + nz - 1)) *
                  (grid.w[(nx - 1) * ny * nz + nz - 1] - exact_solution.value_z((nx - 1) * ny * nz + nz - 1)) *
                  dx * dy * dz / 8);

        for (size_t j = 1; j < ny - 1; j++)
        {
            error += ((grid.w[(nx - 1) * ny * nz + j * nz] - exact_solution.value_z((nx - 1) * ny * nz + j * nz)) *
                      (grid.w[(nx - 1) * ny * nz + j * nz] - exact_solution.value_z((nx - 1) * ny * nz + j * nz)) *
                      dx * dy * dz / 4);
            for (size_t k = 1; k < nz - 1; k++)
            {
                error += ((grid.w[(nx - 1) * ny * nz + j * nz + k] - exact_solution.value_z((nx - 1) * ny * nz + j * nz + k)) *
                          (grid.w[(nx - 1) * ny * nz + j * nz + k] - exact_solution.value_z((nx - 1) * ny * nz + j * nz + k)) *
                          dx * dy * dz / 2);
            }
            error += ((grid.w[(nx - 1) * ny * nz + j * nz + nz - 1] - exact_solution.value_z((nx - 1) * ny * nz + j * nz + nz - 1)) *
                      (grid.w[(nx - 1) * ny * nz + j * nz + nz - 1] - exact_solution.value_z((nx - 1) * ny * nz + j * nz + nz - 1)) *
                      dx * dy * dz / 4);
        }

        error += ((grid.w[(nx - 1) * ny * nz + (ny - 1) * nz] - exact_solution.value_z((nx - 1) * ny * nz + (ny - 1) * nz)) *
                  (grid.w[(nx - 1) * ny * nz + (ny - 1) * nz] - exact_solution.value_z((nx - 1) * ny * nz + (ny - 1) * nz)) *
                  dx * dy * dz / 8);

        for (size_t k = 1; k < nz - 1; k++)
        {
            error += ((grid.w[(nx - 1) * ny * nz + (ny - 1) * nz + k] - exact_solution.value_z((nx - 1) * ny * nz + (ny - 1) * nz + k)) *
                      (grid.w[(nx - 1) * ny * nz + (ny - 1) * nz + k] - exact_solution.value_z((nx - 1) * ny * nz + (ny - 1) * nz + k)) *
                      dx * dy * dz / 4);
        }

        error += ((grid.w[(nx - 1) * ny * nz + (ny - 1) * nz + nz - 1] - exact_solution.value_z((nx - 1) * ny * nz + (ny - 1) * nz + nz - 1)) *
                  (grid.w[(nx - 1) * ny * nz + (ny - 1) * nz + nz - 1] - exact_solution.value_z((nx - 1) * ny * nz + (ny - 1) * nz + nz - 1)) *
                  dx * dy * dz / 8);
    }

    return error;
}

double IcoNS::error_comp_Z()
{
    double error = 0.0;

    // first slice (left face)
    {
        error += ((grid.v[0] - exact_solution.value_y(0)) *
                  (grid.v[0] - exact_solution.value_y(0)) *
                  dx * dy * dz / 8);

        for (size_t k = 1; k < nz - 1; k++)
        {
            error += ((grid.v[k] - exact_solution.value_y(k)) *
                      (grid.v[k] - exact_solution.value_y(k)) *
                      dx * dy * dz / 4);
        }

        error += ((grid.v[nz - 1] - exact_solution.value_y(nz - 1)) *
                  (grid.v[nz - 1] - exact_solution.value_y(nz - 1)) *
                  dx * dy * dz / 8);

        for (size_t j = 1; j < ny - 1; j++)
        {
            error += ((grid.v[j * nz] - exact_solution.value_y(j * nz)) *
                      (grid.v[j * nz] - exact_solution.value_y(j * nz)) *
                      dx * dy * dz / 4);
            for (size_t k = 1; k < nz - 1; k++)
            {
                error += ((grid.v[j * nz + k] - exact_solution.value_y(j * nz + k)) *
                          (grid.v[j * nz + k] - exact_solution.value_y(j * nz + k)) *
                          dx * dy * dz / 2);
            }
            error += ((grid.v[j * nz + nz - 1] - exact_solution.value_y(j * nz + nz - 1)) *
                      (grid.v[j * nz + nz - 1] - exact_solution.value_y(j * nz + nz - 1)) *
                      dx * dy * dz / 4);
        }

        error += ((grid.v[(ny - 1) * nz] - exact_solution.value_y((ny - 1) * nz)) *
                  (grid.v[(ny - 1) * nz] - exact_solution.value_y((ny - 1) * nz)) *
                  dx * dy * dz / 8);

        for (size_t k = 1; k < nz - 1; k++)
        {
            error += ((grid.v[(ny - 1) * nz + k] - exact_solution.value_y((ny - 1) * nz + k)) *
                      (grid.v[(ny - 1) * nz + k] - exact_solution.value_y((ny - 1) * nz + k)) *
                      dx * dy * dz / 4);
        }

        error += ((grid.v[(ny - 1) * nz + nz - 1] - exact_solution.value_y((ny - 1) * nz + nz - 1)) *
                  (grid.v[(ny - 1) * nz + nz - 1] - exact_solution.value_y((ny - 1) * nz + nz - 1)) *
                  dx * dy * dz / 8);
    }

    // middle slices
    {
        for (size_t i = 1; i < nx - 1; i++)
        {
            error += ((grid.v[i * ny * nz] - exact_solution.value_y(i * ny * nz)) *
                      (grid.v[i * ny * nz] - exact_solution.value_y(i * ny * nz)) *
                      dx * dy * dz / 4);

            for (size_t k = 1; k < nz - 1; k++)
            {
                error += ((grid.v[i * ny * nz + k] - exact_solution.value_y(i * ny * nz + k)) *
                          (grid.v[i * ny * nz + k] - exact_solution.value_y(i * ny * nz + k)) *
                          dx * dy * dz / 2);
            }
            error += ((grid.v[i * ny * nz + nz - 1] - exact_solution.value_y(i * ny * nz + nz - 1)) *
                      (grid.v[i * ny * nz + nz - 1] - exact_solution.value_y(i * ny * nz + nz - 1)) *
                      dx * dy * dz / 4);

            for (size_t j = 1; j < ny - 1; j++)
            {
                error += ((grid.v[i * ny * nz + j * nz] - exact_solution.value_y(i * ny * nz + j * nz)) *
                          (grid.v[i * ny * nz + j * nz] - exact_solution.value_y(i * ny * nz + j * nz)) *
                          dx * dy * dz / 2);

                for (size_t k = 1; k < nz - 1; k++)
                {
                    error += ((grid.v[i * ny * nz + j * nz + k] - exact_solution.value_y(i * ny * nz + j * nz + k)) *
                              (grid.v[i * ny * nz + j * nz + k] - exact_solution.value_y(i * ny * nz + j * nz + k)) *
                              dx * dy * dz);
                }

                error += ((grid.v[i * ny * nz + j * nz + nz - 1] - exact_solution.value_y(i * ny * nz + j * nz + nz - 1)) *
                          (grid.v[i * ny * nz + j * nz + nz - 1] - exact_solution.value_y(i * ny * nz + j * nz + nz - 1)) *
                          dx * dy * dz / 2);
            }

            error += ((grid.v[i * ny * nz + (ny - 1) * nz] - exact_solution.value_y(i * ny * nz + (ny - 1) * nz)) *
                      (grid.v[i * ny * nz + (ny - 1) * nz] - exact_solution.value_y(i * ny * nz + (ny - 1) * nz)) *
                      dx * dy * dz / 4);

            for (size_t k = 1; k < nz - 1; k++)
            {
                error += ((grid.v[i * ny * nz + (ny - 1) * nz + k] - exact_solution.value_y(i * ny * nz + (ny - 1) * nz + k)) *
                          (grid.v[i * ny * nz + (ny - 1) * nz + k] - exact_solution.value_y(i * ny * nz + (ny - 1) * nz + k)) *
                          dx * dy * dz / 2);
            }

            error += ((grid.v[i * ny * nz + (ny - 1) * nz + nz - 1] - exact_solution.value_y(i * ny * nz + (ny - 1) * nz + nz - 1)) *
                      (grid.v[i * ny * nz + (ny - 1) * nz + nz - 1] - exact_solution.value_y(i * ny * nz + (ny - 1) * nz + nz - 1)) *
                      dx * dy * dz / 4);
        }
    }

    // last slice (right face)
    {
        error += ((grid.v[(nx - 1) * ny * nz] - exact_solution.value_y((nx - 1) * ny * nz)) *
                  (grid.v[(nx - 1) * ny * nz] - exact_solution.value_y((nx - 1) * ny * nz)) *
                  dx * dy * dz / 8);

        for (size_t k = 1; k < nz - 1; k++)
        {
            error += ((grid.v[(nx - 1) * ny * nz + k] - exact_solution.value_y((nx - 1) * ny * nz + k)) *
                      (grid.v[(nx - 1) * ny * nz + k] - exact_solution.value_y((nx - 1) * ny * nz + k)) *
                      dx * dy * dz / 4);
        }

        error += ((grid.v[(nx - 1) * ny * nz + nz - 1] - exact_solution.value_y((nx - 1) * ny * nz + nz - 1)) *
                  (grid.v[(nx - 1) * ny * nz + nz - 1] - exact_solution.value_y((nx - 1) * ny * nz + nz - 1)) *
                  dx * dy * dz / 8);

        for (size_t j = 1; j < ny - 1; j++)
        {
            error += ((grid.v[(nx - 1) * ny * nz + j * nz] - exact_solution.value_y((nx - 1) * ny * nz + j * nz)) *
                      (grid.v[(nx - 1) * ny * nz + j * nz] - exact_solution.value_y((nx - 1) * ny * nz + j * nz)) *
                      dx * dy * dz / 4);
            for (size_t k = 1; k < nz - 1; k++)
            {
                error += ((grid.v[(nx - 1) * ny * nz + j * nz + k] - exact_solution.value_y((nx - 1) * ny * nz + j * nz + k)) *
                          (grid.v[(nx - 1) * ny * nz + j * nz + k] - exact_solution.value_y((nx - 1) * ny * nz + j * nz + k)) *
                          dx * dy * dz / 2);
            }
            error += ((grid.v[(nx - 1) * ny * nz + j * nz + nz - 1] - exact_solution.value_y((nx - 1) * ny * nz + j * nz + nz - 1)) *
                      (grid.v[(nx - 1) * ny * nz + j * nz + nz - 1] - exact_solution.value_y((nx - 1) * ny * nz + j * nz + nz - 1)) *
                      dx * dy * dz / 4);
        }

        error += ((grid.v[(nx - 1) * ny * nz + (ny - 1) * nz] - exact_solution.value_y((nx - 1) * ny * nz + (ny - 1) * nz)) *
                  (grid.v[(nx - 1) * ny * nz + (ny - 1) * nz] - exact_solution.value_y((nx - 1) * ny * nz + (ny - 1) * nz)) *
                  dx * dy * dz / 8);

        for (size_t k = 1; k < nz - 1; k++)
        {
            error += ((grid.v[(nx - 1) * ny * nz + (ny - 1) * nz + k] - exact_solution.value_y((nx - 1) * ny * nz + (ny - 1) * nz + k)) *
                      (grid.v[(nx - 1) * ny * nz + (ny - 1) * nz + k] - exact_solution.value_y((nx - 1) * ny * nz + (ny - 1) * nz + k)) *
                      dx * dy * dz / 4);
        }

        error += ((grid.v[(nx - 1) * ny * nz + (ny - 1) * nz + nz - 1] - exact_solution.value_y((nx - 1) * ny * nz + (ny - 1) * nz + nz - 1)) *
                  (grid.v[(nx - 1) * ny * nz + (ny - 1) * nz + nz - 1] - exact_solution.value_y((nx - 1) * ny * nz + (ny - 1) * nz + nz - 1)) *
                  dx * dy * dz / 8);
    }
    return error;
}

double IcoNS::L2_error()
{
    double error = 0.0;

    error += error_comp_X();
    error += error_comp_Y();
    error += error_comp_Z();

    return sqrt(error);
}

void IcoNS::output()
{
}
