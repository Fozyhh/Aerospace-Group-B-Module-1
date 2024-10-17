#include "../includes/core.hpp"
#include <cmath>

void IcoNS::preprocessing(/* std::string &input_file */)
{
    // read the input file.
    for (size_t i = 1; i < nx - 1; i++)
    {
        for (size_t j = 1; j < ny - 1; j++)
        {
            for (size_t k = 1; k < nz - 1; k++)
            {
                grid.u[i * ny * nz + j * nz + k] = 0.0;
                grid.v[i * ny * nz + j * nz + k] = 0.0;
                grid.w[i * ny * nz + j * nz + k] = 0.0;
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
                    grid.u[i * ny * nz + j * nz + k] = 0.0;
                    grid.v[i * ny * nz + j * nz + k] = 0.0;
                    grid.w[i * ny * nz + j * nz + k] = 0.0;
                }
            }
            // pressure initialization.
        }
    }
}

void IcoNS::solve()
{
    double time = 0.0;

    while (time < T)
    {
        solve_time_step(time);
        time += dt;

        // output();
    }
}

std::vector<double> IcoNS::functionF(const std::vector<double> &u, const std::vector<double> &v,
                                     const std::vector<double> &w, size_t i, size_t j, size_t k, double t)
{
    std::vector<double> f(3);
    size_t l = i * ny * nz + j * nz + k;

    std::vector<double> g(3);
    g = functionG(i, j, k, t);

    f[0] = -(u[l] * (u[l + nz * ny] - u[l - nz * ny]) / (2.0 * dx) +
             (v[l + nz * ny] + v[l] + v[l - nz] + v[l + nz * ny - nz]) / 4.0 * (u[l + nz] - u[l - nz]) / (2.0 * dy) +
             (w[l + nz * ny] + w[l] + w[l + nz * ny - 1] + w[l - 1]) / 4.0 * (u[l + 1] - u[l - 1]) / (2.0 * dz)) +
           1 / Re * ((u[l + ny * nz] - 2 * u[l] + u[l - ny * nz]) / (dx * dx) + (u[l + nz] - 2 * u[l] + u[l - nz]) / (dy * dy) + (u[l + 1] - 2 * u[l] + u[l - 1]) / (dz * dz)) - g[0];

    f[1] = -((u[l] + u[l + nz] + u[l - ny * nz] + u[l - ny * nz + nz]) / 4.0 * (v[l + ny * nz] - v[l - ny * nz] / (2.0 * dx)) +
             v[l] * (v[l + nz] - v[l - nz]) / (2.0 * dy) +
             (w[l] + w[l - 1] + w[l + nz] + w[l + nz - 1]) / 4.0 + (v[l + 1] - v[l - 1]) / (2.0 * dz)) +
           1 / Re * ((v[l + nz * ny] - 2.0 * v[l] + v[l - nz * ny]) / (dx * dx) + (v[l + nz] - 2.0 * v[l] + v[l - nz]) / (dy * dy) + (v[l + 1] - 2.0 * v[l] + v[l - 1]) / (dz * dz)) - g[1];

    f[2] = -((u[l] + u[l - ny * nz] + u[l + 1] + u[l - nz * ny + 1]) / 4.0 * (w[l + nz * ny] - w[l - nz * ny]) / (2.0 * dx) +
             (v[l + 1] + v[l - nz + 1] + v[l] + v[l - nz]) / 4.0 * (w[l + nz] - w[l - nz]) / (2.0 * dy) +
             w[l] * (w[l + 1] - w[l - 1]) / (2.0 * dz)) +
           1 / Re * ((w[l + nz * ny] - 2.0 * w[l] + w[l - nz * ny]) / (dx * dx) + (w[l + nz] - 2.0 * w[l] + w[l - nz]) / (dy * dy) + (w[l + 1] - 2.0 * w[l] + w[l - 1]) / (dz * dz)) - g[2];

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

    //std::cout << "g[0] = " << g[0] << std::endl;
    return g;
}

void IcoNS::solve_time_step(double time)
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
                f = functionF(grid.u, grid.v, grid.w, i, j, k, time);
                std::cout << "f[0] = " << f[0] << std::endl;
                // solve the momentum equations -> TODO later also with pressure
                Y2_x[i * ny * nz + j * nz + k] = grid.u[i * ny * nz + j * nz + k] +
                                                 64.0 / 120.0 * dt * f[0];
                Y2_y[i * ny * nz + j * nz + k] = /* grid.v[i * ny * nz + j * nz + k] + */
                                                 64.0 / 120.0 * dt * f[1];
                Y2_z[i * ny * nz + j * nz + k] = /* grid.w[i * ny * nz + j * nz + k] + */
                                                 64.0 / 120.0 * dt * f[2];
                //std::cout << "dt " << dt << std::endl;

                f_y2 = functionF(Y2_x, Y2_y, Y2_z, i, j, k, time + 64.0 / 120.0 * dt);

                Y3_x[i * ny * nz + j * nz + k] = Y2_x[i * ny * nz + j * nz + k] +
                                                 50.0 / 120.0 * dt * f_y2[0] -
                                                 34.0 / 120.0 * dt * f[0];
                Y3_y[i * ny * nz + j * nz + k] = Y2_y[i * ny * nz + j * nz + k] +
                                                 50.0 / 120.0 * dt * f_y2[1] -
                                                 34.0 / 120.0 * dt * f[1];
                Y3_z[i * ny * nz + j * nz + k] = Y2_z[i * ny * nz + j * nz + k] +
                                                 50.0 / 120.0 * dt * f_y2[2] -
                                                 34.0 / 120.0 * dt * f[2];

                // update the grid values.
                f_y3 = functionF(Y3_x, Y3_y, Y3_z, i, j, k, time + 80 / 120 * dt);

                grid.u[i * ny * nz + j * nz + k] = Y3_x[i * ny * nz + j * nz + k] +
                                                   90.0 / 120.0 * dt * f_y3[0] -
                                                   50.0 / 120.0 * dt * f_y2[0];
                grid.v[i * ny * nz + j * nz + k] = Y3_y[i * ny * nz + j * nz + k] +
                                                   90.0 / 120.0 * dt * f_y3[1] -
                                                   50.0 / 120.0 * dt * f_y2[1];
                grid.w[i * ny * nz + j * nz + k] = Y3_z[i * ny * nz + j * nz + k] +
                                                   90.0 / 120.0 * dt * f_y3[2] -
                                                   50.0 / 120.0 * dt * f_y2[2];
                //std::cout << "u[" << i << "," << j << "," << k << "] = " << grid.u[i * ny * nz + j * nz + k] << std::endl;
            }
        }
    }

    //compute boundary conditions
    for (size_t i = 0; i < nx; i++)
    {
        for (size_t j = 0; j < ny; j++)
        {
            for (size_t k = 0; k < nz; k++)
            {
                if(i == 0 || i == nx-1 || j == 0 || j == ny-1 || k == 0 || k == nz-1)
                {
                    grid.u[i * ny * nz + j * nz + k] = std::sin(i) * std::cos(j) * std::sin(k) * std::cos(time+dt);
                    grid.v[i * ny * nz + j * nz + k] = std::cos(i) * std::sin(j) * std::sin(k) * std::cos(time+dt);
                    grid.w[i * ny * nz + j * nz + k] = 2*std::cos(i) * std::cos(j) * std::cos(k) * std::sin(time+dt);
                }
            }
        }
    }

    // print the grid values.
    double diff;
    for (size_t i = 0; i < nx; i++)
    {
        for (size_t j = 0; j < ny; j++)
        {
            for (size_t k = 0; k < nz; k++)
            {
                diff = std::abs(grid.u[i * ny * nz + j * nz + k] - std::sin(i) * std::cos(j) * std::sin(k) * std::cos(time+dt));
                //diff = grid.u[i * ny * nz + j * nz + k];
                std::cout << "diff[" << i << "," << j << "," << k << "] = " << diff << std::endl;
                // std::cout << "u[" << i << "," << j << "," << k << "] = " << grid.u[i * ny * nz + j * nz + k] << std::endl;
                // std::cout << "v[" << i << "," << j << "," << k << "] = " << grid.v[i * ny * nz + j * nz + k] << std::endl;
                // std::cout << "w[" << i << "," << j << "," << k << "] = " << grid.w[i * ny * nz + j * nz + k] << std::endl;
            }
        }
    }
}

void IcoNS::output()
{
    // write the output file.
}