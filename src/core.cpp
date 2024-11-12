
#include <cmath>

#include "core.hpp"
#include <math.h>

#include <fstream>
#include <string>
#include <memory>

void IcoNS::preprocessing(/*std::string &input_file*/)

{
    // boundary
    auto u_func = std::make_shared<Dirichlet>([&](Real x, Real y, Real z, Real t)
                                              { return std::sin((x + 0.5) * DX) * std::cos(y * DY) * std::sin(z * DZ) * std::sin(t); });
    auto v_func = std::make_shared<Dirichlet>([&](Real x, Real y, Real z, Real t)
                                              { return std::cos(x * DX) * std::sin((y + 0.5) * DY) * std::sin(z * DZ) * std::sin(t); });
    auto w_func = std::make_shared<Dirichlet>([&](Real x, Real y, Real z, Real t)
                                              { return 2 * (std::cos(x * DX) * std::cos(y * DY) * std::cos((z + 0.5) * DZ) * std::sin(t)); });

    for (size_t i = 0; i < 6 /*nfaces*/; i++)
    {
        boundary.addFunction(U, u_func);
        boundary.addFunction(V, v_func);
        boundary.addFunction(W, w_func);
    }
}

void IcoNS::solve()
{
    Real time = 0.0;
    int i = 0;

    std::ofstream error_log("../resources/error.log");

    while (time < T)
    {
        /*Check::Confront(grid,exact_solution,time,U);
        int p;
        std::cin >> p;*/
        boundary.update_boundary(grid.u, grid.v, grid.w, time);

        // csv file w/ "," delimiter: time step, iter, L2_error
        std::cout << time << "," << i << "," << L2_error(time) << std::endl;
        solve_time_step(time);
        // output();
        time += DT;
        i++;
    }
}

Real IcoNS::L2_error(const Real t)
{
    Real error = 0.0;

    error += error_comp_X(t);
    error += error_comp_Y(t);
    error += error_comp_Z(t);

    return sqrt(error);
}

Real IcoNS::error_comp_X(const Real t)
{
    Real error = 0.0;

    // first slice (left face)
    {
        error += ((grid.u[0] - exact_solution.value_x(0.5, 0, 0, t)) *
                  (grid.u[0] - exact_solution.value_x(0.5, 0, 0, t)) *
                  DX * DY * DZ / 8);

        for (size_t k = 1; k < NZ; k++)
        {
            error += ((grid.u[k] - exact_solution.value_x(0.5, 0, k, t)) *
                      (grid.u[k] - exact_solution.value_x(0.5, 0, k, t)) *
                      DX * DY * DZ / 4);
        }

        error += ((grid.u[NZ] - exact_solution.value_x(0.5, 0, NZ, t)) *
                  (grid.u[NZ] - exact_solution.value_x(0.5, 0, NZ, t)) *
                  DX * DY * DZ / 8);

        for (size_t j = 1; j < NY; j++)
        {
            error += ((grid.u[j * (NZ + 1)] - exact_solution.value_x(0.5, j, 0, t)) *
                      (grid.u[j * (NZ + 1)] - exact_solution.value_x(0.5, j, 0, t)) *
                      DX * DY * DZ / 4);
            for (size_t k = 1; k < NZ; k++)
            {
                error += ((grid.u[j * (NZ + 1) + k] - exact_solution.value_x(0.5, j, k, t)) *
                          (grid.u[j * (NZ + 1) + k] - exact_solution.value_x(0.5, j, k, t)) *
                          DX * DY * DZ / 2);
            }
            error += ((grid.u[j * (NZ + 1) + NZ] - exact_solution.value_x(0.5, j, NZ, t)) *
                      (grid.u[j * (NZ + 1) + NZ] - exact_solution.value_x(0.5, j, NZ, t)) *
                      DX * DY * DZ / 4);
        }

        error += ((grid.u[NY * (NZ + 1)] - exact_solution.value_x(0.5, NY, 0, t)) *
                  (grid.u[NY * (NZ + 1)] - exact_solution.value_x(0.5, NY, 0, t)) *
                  DX * DY * DZ / 8);

        for (size_t k = 1; k < NZ; k++)
        {
            error += ((grid.u[NY * (NZ + 1) + k] - exact_solution.value_x(0.5, NY, k, t)) *
                      (grid.u[NY * (NZ + 1) + k] - exact_solution.value_x(0.5, NY, k, t)) *
                      DX * DY * DZ / 4);
        }

        error += ((grid.u[NY * (NZ + 1) + NZ] - exact_solution.value_x(0.5, NY, NZ, t)) *
                  (grid.u[NY * (NZ + 1) + NZ] - exact_solution.value_x(0.5, NY, NZ, t)) *
                  DX * DY * DZ / 8);
    }

    // middle slices
    {
        for (size_t i = 1; i < NX - 1; i++)
        {
            error += ((grid.u[i * (NY + 1) * (NZ + 1)] - exact_solution.value_x(i + 0.5, 0, 0, t)) *
                      (grid.u[i * (NY + 1) * (NZ + 1)] - exact_solution.value_x(i + 0.5, 0, 0, t)) *
                      DX * DY * DZ / 4);

            for (size_t k = 1; k < NZ; k++)
            {
                error += ((grid.u[i * (NY + 1) * (NZ + 1) + k] - exact_solution.value_x(i + 0.5, 0, k, t)) *
                          (grid.u[i * (NY + 1) * (NZ + 1) + k] - exact_solution.value_x(i + 0.5, 0, k, t)) *
                          DX * DY * DZ / 2);
            }
            error += ((grid.u[i * (NY + 1) * (NZ + 1) + NZ] - exact_solution.value_x(i + 0.5, 0, NZ, t)) *
                      (grid.u[i * (NY + 1) * (NZ + 1) + NZ] - exact_solution.value_x(i + 0.5, 0, NZ, t)) *
                      DX * DY * DZ / 4);

            for (size_t j = 1; j < NY; j++)
            {
                error += ((grid.u[i * (NY + 1) * (NZ + 1) + j * (NZ + 1)] - exact_solution.value_x(i + 0.5, j, 0, t)) *
                          (grid.u[i * (NY + 1) * (NZ + 1) + j * (NZ + 1)] - exact_solution.value_x(i + 0.5, j, 0, t)) *
                          DX * DY * DZ / 2);

                for (size_t k = 1; k < NZ; k++)
                {
                    error += ((grid.u[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k] - exact_solution.value_x(i + 0.5, j, k, t)) *
                              (grid.u[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k] - exact_solution.value_x(i + 0.5, j, k, t)) *
                              DX * DY * DZ);
                }

                error += ((grid.u[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + NZ] - exact_solution.value_x(i + 0.5, j, NZ, t)) *
                          (grid.u[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + NZ] - exact_solution.value_x(i + 0.5, j, NZ, t)) *
                          DX * DY * DZ / 2);
            }

            error += ((grid.u[i * (NY + 1) * (NZ + 1) + NY * (NZ + 1)] - exact_solution.value_x(i + 0.5, NY, 0, t)) *
                      (grid.u[i * (NY + 1) * (NZ + 1) + NY * (NZ + 1)] - exact_solution.value_x(i + 0.5, NY, 0, t)) *
                      DX * DY * DZ / 4);

            for (size_t k = 1; k < NZ; k++)
            {
                error += ((grid.u[i * (NY + 1) * (NZ + 1) + NY * (NZ + 1) + k] - exact_solution.value_x(i + 0.5, NY, k, t)) *
                          (grid.u[i * (NY + 1) * (NZ + 1) + NY * (NZ + 1) + k] - exact_solution.value_x(i + 0.5, NY, k, t)) *
                          DX * DY * DZ / 2);
            }

            error += ((grid.u[i * (NY + 1) * (NZ + 1) + NY * (NZ + 1) + NZ] - exact_solution.value_x(i + 0.5, NY, NZ, t)) *
                      (grid.u[i * (NY + 1) * (NZ + 1) + NY * (NZ + 1) + NZ] - exact_solution.value_x(i + 0.5, NY, NZ, t)) *
                      DX * DY * DZ / 4);
        }
    }

    // last slice (right face)
    {
        error += ((grid.u[(NX - 1) * (NY + 1) * (NZ + 1)] - exact_solution.value_x((NX - 0.5), 0, 0, t)) *
                  (grid.u[(NX - 1) * (NY + 1) * (NZ + 1)] - exact_solution.value_x((NX - 0.5), 0, 0, t)) *
                  DX * DY * DZ / 8);

        for (size_t k = 1; k < NZ; k++)
        {
            error += ((grid.u[(NX - 1) * (NY + 1) * (NZ + 1) + k] - exact_solution.value_x((NX - 0.5), 0, k, t)) *
                      (grid.u[(NX - 1) * (NY + 1) * (NZ + 1) + k] - exact_solution.value_x((NX - 0.5), 0, k, t)) *
                      DX * DY * DZ / 4);
        }

        error += ((grid.u[(NX - 1) * (NY + 1) * (NZ + 1) + NZ] - exact_solution.value_x((NX - 0.5), 0, NZ, t)) *
                  (grid.u[(NX - 1) * (NY + 1) * (NZ + 1) + NZ] - exact_solution.value_x((NX - 0.5), 0, NZ, t)) *
                  DX * DY * DZ / 8);

        for (size_t j = 1; j < NY; j++)
        {
            error += ((grid.u[(NX - 1) * (NY + 1) * (NZ + 1) + j * (NZ + 1)] - exact_solution.value_x((NX - 0.5), j, 0, t)) *
                      (grid.u[(NX - 1) * (NY + 1) * (NZ + 1) + j * (NZ + 1)] - exact_solution.value_x((NX - 0.5), j, 0, t)) *
                      DX * DY * DZ / 4);
            for (size_t k = 1; k < NZ; k++)
            {
                error += ((grid.u[(NX - 1) * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k] - exact_solution.value_x((NX - 0.5), j, k, t)) *
                          (grid.u[(NX - 1) * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k] - exact_solution.value_x((NX - 0.5), j, k, t)) *
                          DX * DY * DZ / 2);
            }
            error += ((grid.u[(NX - 1) * (NY + 1) * (NZ + 1) + j * (NZ + 1) + NZ] - exact_solution.value_x((NX - 0.5), j, NZ, t)) *
                      (grid.u[(NX - 1) * (NY + 1) * (NZ + 1) + j * (NZ + 1) + NZ] - exact_solution.value_x((NX - 0.5), j, NZ, t)) *
                      DX * DY * DZ / 4);
        }

        error += ((grid.u[(NX - 1) * (NY + 1) * (NZ + 1) + NY * (NZ + 1)] - exact_solution.value_x((NX - 0.5), NY, 0, t)) *
                  (grid.u[(NX - 1) * (NY + 1) * (NZ + 1) + NY * (NZ + 1)] - exact_solution.value_x((NX - 0.5), NY, 0, t)) *
                  DX * DY * DZ / 8);

        for (size_t k = 1; k < NZ; k++)
        {
            error += ((grid.u[(NX - 1) * (NY + 1) * (NZ + 1) + NY * (NZ + 1) + k] - exact_solution.value_x((NX - 0.5), NY, k, t)) *
                      (grid.u[(NX - 1) * (NY + 1) * (NZ + 1) + NY * (NZ + 1) + k] - exact_solution.value_x((NX - 0.5), NY, k, t)) *
                      DX * DY * DZ / 4);
        }

        error += ((grid.u[(NX - 1) * (NY + 1) * (NZ + 1) + NY * (NZ + 1) + NZ] - exact_solution.value_x((NX - 0.5), NY, NZ, t)) *
                  (grid.u[(NX - 1) * (NY + 1) * (NZ + 1) + NY * (NZ + 1) + NZ] - exact_solution.value_x((NX - 0.5), NY, NZ, t)) *
                  DX * DY * DZ / 8);
    }
    return error;
}

Real IcoNS::error_comp_Y(const Real t)
{
    Real error = 0.0;

    // first slice (left face)
    {
        error += ((grid.v[0] - exact_solution.value_y(0, 0.5, 0, t)) *
                  (grid.v[0] - exact_solution.value_y(0, 0.5, 0, t)) *
                  DX * DY * DZ / 8);

        for (size_t k = 1; k < NZ; k++)
        {
            error += ((grid.v[k] - exact_solution.value_y(0, 0.5, k, t)) *
                      (grid.v[k] - exact_solution.value_y(0, 0.5, k, t)) *
                      DX * DY * DZ / 4);
        }

        error += ((grid.v[NZ] - exact_solution.value_y(0, 0.5, NZ, t)) *
                  (grid.v[NZ] - exact_solution.value_y(0, 0.5, NZ, t)) *
                  DX * DY * DZ / 8);

        for (size_t j = 1; j < NY - 1; j++)
        {
            error += ((grid.v[j * (NZ + 1)] - exact_solution.value_y(0, j + 0.5, 0, t)) *
                      (grid.v[j * (NZ + 1)] - exact_solution.value_y(0, j + 0.5, 0, t)) *
                      DX * DY * DZ / 4);
            for (size_t k = 1; k < NZ; k++)
            {
                error += ((grid.v[j * (NZ + 1) + k] - exact_solution.value_y(0, j + 0.5, k, t)) *
                          (grid.v[j * (NZ + 1) + k] - exact_solution.value_y(0, j + 0.5, k, t)) *
                          DX * DY * DZ / 2);
            }
            error += ((grid.v[j * (NZ + 1) + NZ] - exact_solution.value_y(0, j + 0.5, NZ, t)) *
                      (grid.v[j * (NZ + 1) + NZ] - exact_solution.value_y(0, j + 0.5, NZ, t)) *
                      DX * DY * DZ / 4);
        }

        error += ((grid.v[(NY - 1) * (NZ + 1)] - exact_solution.value_y(0, (NY - 0.5), 0, t)) *
                  (grid.v[(NY - 1) * (NZ + 1)] - exact_solution.value_y(0, (NY - 0.5), 0, t)) *
                  DX * DY * DZ / 8);

        for (size_t k = 1; k < NZ; k++)
        {
            error += ((grid.v[(NY - 1) * (NZ + 1) + k] - exact_solution.value_y(0, (NY - 0.5), k, t)) *
                      (grid.v[(NY - 1) * (NZ + 1) + k] - exact_solution.value_y(0, (NY - 0.5), k, t)) *
                      DX * DY * DZ / 4);
        }

        error += ((grid.v[(NY - 1) * (NZ + 1) + NZ] - exact_solution.value_y(0, (NY - 0.5), NZ, t)) *
                  (grid.v[(NY - 1) * (NZ + 1) + NZ] - exact_solution.value_y(0, (NY - 0.5), NZ, t)) *
                  DX * DY * DZ / 8);
    }

    // middle slices
    {
        for (size_t i = 1; i < NX; i++)
        {
            error += ((grid.v[i * NY * (NZ + 1)] - exact_solution.value_y(i, 0.5, 0, t)) *
                      (grid.v[i * NY * (NZ + 1)] - exact_solution.value_y(i, 0.5, 0, t)) *
                      DX * DY * DZ / 4);

            for (size_t k = 1; k < NZ; k++)
            {
                error += ((grid.v[i * NY * (NZ + 1) + k] - exact_solution.value_y(i, 0.5, k, t)) *
                          (grid.v[i * NY * (NZ + 1) + k] - exact_solution.value_y(i, 0.5, k, t)) *
                          DX * DY * DZ / 2);
            }
            error += ((grid.v[i * NY * (NZ + 1) + NZ] - exact_solution.value_y(i, 0.5, NZ, t)) *
                      (grid.v[i * NY * (NZ + 1) + NZ] - exact_solution.value_y(i, 0.5, NZ, t)) *
                      DX * DY * DZ / 4);

            for (size_t j = 1; j < NY - 1; j++)
            {
                error += ((grid.v[i * NY * (NZ + 1) + j * (NZ + 1)] - exact_solution.value_y(i, j + 0.5, 0, t)) *
                          (grid.v[i * NY * (NZ + 1) + j * (NZ + 1)] - exact_solution.value_y(i, j + 0.5, 0, t)) *
                          DX * DY * DZ / 2);

                for (size_t k = 1; k < NZ; k++)
                {
                    error += ((grid.v[i * NY * (NZ + 1) + j * (NZ + 1) + k] - exact_solution.value_y(i, j + 0.5, k, t)) *
                              (grid.v[i * NY * (NZ + 1) + j * (NZ + 1) + k] - exact_solution.value_y(i, j + 0.5, k, t)) *
                              DX * DY * DZ);
                }

                error += ((grid.v[i * NY * (NZ + 1) + j * (NZ + 1) + NZ] - exact_solution.value_y(i, j + 0.5, NZ, t)) *
                          (grid.v[i * NY * (NZ + 1) + j * (NZ + 1) + NZ] - exact_solution.value_y(i, j + 0.5, NZ, t)) *
                          DX * DY * DZ / 2);
            }

            error += ((grid.v[i * NY * (NZ + 1) + (NY - 1) * (NZ + 1)] - exact_solution.value_y(i, (NY - 0.5), 0, t)) *
                      (grid.v[i * NY * (NZ + 1) + (NY - 1) * (NZ + 1)] - exact_solution.value_y(i, (NY - 0.5), 0, t)) *
                      DX * DY * DZ / 4);

            for (size_t k = 1; k < NZ; k++)
            {
                error += ((grid.v[i * NY * (NZ + 1) + (NY - 1) * (NZ + 1) + k] - exact_solution.value_y(i, (NY - 0.5), k, t)) *
                          (grid.v[i * NY * (NZ + 1) + (NY - 1) * (NZ + 1) + k] - exact_solution.value_y(i, (NY - 0.5), k, t)) *
                          DX * DY * DZ / 2);
            }

            error += ((grid.v[i * NY * (NZ + 1) + (NY - 1) * (NZ + 1) + NZ] - exact_solution.value_y(i, (NY - 0.5), NZ, t)) *
                      (grid.v[i * NY * (NZ + 1) + (NY - 1) * (NZ + 1) + NZ] - exact_solution.value_y(i, (NY - 0.5), NZ, t)) *
                      DX * DY * DZ / 4);
        }
    }

    // last slice (right face)
    {
        error += ((grid.v[NX * NY * (NZ + 1)] - exact_solution.value_y(NX, 0.5, 0, t)) *
                  (grid.v[NX * NY * (NZ + 1)] - exact_solution.value_y(NX, 0.5, 0, t)) *
                  DX * DY * DZ / 8);

        for (size_t k = 1; k < NZ; k++)
        {
            error += ((grid.v[NX * NY * (NZ + 1) + k] - exact_solution.value_y(NX, 0.5, k, t)) *
                      (grid.v[NX * NY * (NZ + 1) + k] - exact_solution.value_y(NX, 0.5, k, t)) *
                      DX * DY * DZ / 4);
        }

        error += ((grid.v[NX * NY * (NZ + 1) + NZ] - exact_solution.value_y(NX, 0.5, NZ, t)) *
                  (grid.v[NX * NY * (NZ + 1) + NZ] - exact_solution.value_y(NX, 0.5, NZ, t)) *
                  DX * DY * DZ / 8);

        for (size_t j = 1; j < NY - 1; j++)
        {
            error += ((grid.v[NX * NY * (NZ + 1) + j * (NZ + 1)] - exact_solution.value_y(NX, j + 0.5, 0, t)) *
                      (grid.v[NX * NY * (NZ + 1) + j * (NZ + 1)] - exact_solution.value_y(NX, j + 0.5, 0, t)) *
                      DX * DY * DZ / 4);
            for (size_t k = 1; k < NZ; k++)
            {
                error += ((grid.v[NX * NY * (NZ + 1) + j * (NZ + 1) + k] - exact_solution.value_y(NX, j + 0.5, k, t)) *
                          (grid.v[NX * NY * (NZ + 1) + j * (NZ + 1) + k] - exact_solution.value_y(NX, j + 0.5, k, t)) *
                          DX * DY * DZ / 2);
            }
            error += ((grid.v[NX * NY * (NZ + 1) + j * (NZ + 1) + NZ] - exact_solution.value_y(NX, j + 0.5, NZ, t)) *
                      (grid.v[NX * NY * (NZ + 1) + j * (NZ + 1) + NZ] - exact_solution.value_y(NX, j + 0.5, NZ, t)) *
                      DX * DY * DZ / 4);
        }

        error += ((grid.v[NX * NY * (NZ + 1) + (NY - 1) * (NZ + 1)] - exact_solution.value_y(NX, (NY - 0.5), 0, t)) *
                  (grid.v[NX * NY * (NZ + 1) + (NY - 1) * (NZ + 1)] - exact_solution.value_y(NX, (NY - 0.5), 0, t)) *
                  DX * DY * DZ / 8);

        for (size_t k = 1; k < NZ; k++)
        {
            error += ((grid.v[NX * NY * (NZ + 1) + (NY - 1) * (NZ + 1) + k] - exact_solution.value_y(NX, (NY - 0.5), k, t)) *
                      (grid.v[NX * NY * (NZ + 1) + (NY - 1) * (NZ + 1) + k] - exact_solution.value_y(NX, (NY - 0.5), k, t)) *
                      DX * DY * DZ / 4);
        }

        error += ((grid.v[NX * NY * (NZ + 1) + (NY - 1) * (NZ + 1) + NZ] - exact_solution.value_y(NX, (NY - 0.5), NZ, t)) *
                  (grid.v[NX * NY * (NZ + 1) + (NY - 1) * (NZ + 1) + NZ] - exact_solution.value_y(NX, (NY - 0.5), NZ, t)) *
                  DX * DY * DZ / 8);
    }
    return error;
}

Real IcoNS::error_comp_Z(const Real t)
{
    Real error = 0.0;
    // first slice (left face)
    {
        error += ((grid.w[0] - exact_solution.value_z(0, 0, 0.5, t)) *
                  (grid.w[0] - exact_solution.value_z(0, 0, 0.5, t)) *
                  DX * DY * DZ / 8);

        for (size_t k = 1; k < NZ - 1; k++)
        {
            error += ((grid.w[k] - exact_solution.value_z(0, 0, k + 0.5, t)) *
                      (grid.w[k] - exact_solution.value_z(0, 0, k + 0.5, t)) *
                      DX * DY * DZ / 4);
        }

        error += ((grid.w[NZ - 1] - exact_solution.value_z(0, 0, NZ - 0.5, t)) *
                  (grid.w[NZ - 1] - exact_solution.value_z(0, 0, NZ - 0.5, t)) *
                  DX * DY * DZ / 8);

        for (size_t j = 1; j < NY; j++)
        {
            error += ((grid.w[j * NZ] - exact_solution.value_z(0, j, 0.5, t)) *
                      (grid.w[j * NZ] - exact_solution.value_z(0, j, 0.5, t)) *
                      DX * DY * DZ / 4);
            for (size_t k = 1; k < NZ - 1; k++)
            {
                error += ((grid.w[j * NZ + k] - exact_solution.value_z(0, j, k + 0.5, t)) *
                          (grid.w[j * NZ + k] - exact_solution.value_z(0, j, k + 0.5, t)) *
                          DX * DY * DZ / 2);
            }
            error += ((grid.w[j * NZ + NZ - 1] - exact_solution.value_z(0, j, NZ - 0.5, t)) *
                      (grid.w[j * NZ + NZ - 1] - exact_solution.value_z(0, j, NZ - 0.5, t)) *
                      DX * DY * DZ / 4);
        }

        error += ((grid.w[NY * NZ] - exact_solution.value_z(0, NY, 0.5, t)) *
                  (grid.w[NY * NZ] - exact_solution.value_z(0, NY, 0.5, t)) *
                  DX * DY * DZ / 8);

        for (size_t k = 1; k < NZ - 1; k++)
        {
            error += ((grid.w[NY * NZ + k] - exact_solution.value_z(0, NY, k + 0.5, t)) *
                      (grid.w[NY * NZ + k] - exact_solution.value_z(0, NY, k + 0.5, t)) *
                      DX * DY * DZ / 4);
        }

        error += ((grid.w[NY * NZ + NZ - 1] - exact_solution.value_z(0, NY, NZ - 0.5, t)) *
                  (grid.w[NY * NZ + NZ - 1] - exact_solution.value_z(0, NY, NZ - 0.5, t)) *
                  DX * DY * DZ / 8);
    }

    // middle slices
    {
        for (size_t i = 1; i < NX; i++)
        {
            error += ((grid.w[i * (NY + 1) * NZ] - exact_solution.value_z(i, 0, 0.5, t)) *
                      (grid.w[i * (NY + 1) * NZ] - exact_solution.value_z(i, 0, 0.5, t)) *
                      DX * DY * DZ / 4);

            for (size_t k = 1; k < NZ - 1; k++)
            {
                error += ((grid.w[i * (NY + 1) * NZ + k] - exact_solution.value_z(i, 0, k + 0.5, t)) *
                          (grid.w[i * (NY + 1) * NZ + k] - exact_solution.value_z(i, 0, k + 0.5, t)) *
                          DX * DY * DZ / 2);
            }
            error += ((grid.w[i * (NY + 1) * NZ + NZ - 1] - exact_solution.value_z(i, 0, NZ - 0.5, t)) *
                      (grid.w[i * (NY + 1) * NZ + NZ - 1] - exact_solution.value_z(i, 0, NZ - 0.5, t)) *
                      DX * DY * DZ / 4);

            for (size_t j = 1; j < NY; j++)
            {
                error += ((grid.w[i * (NY + 1) * NZ + j * NZ] - exact_solution.value_z(i, j, 0.5, t)) *
                          (grid.w[i * (NY + 1) * NZ + j * NZ] - exact_solution.value_z(i, j, 0.5, t)) *
                          DX * DY * DZ / 2);

                for (size_t k = 1; k < NZ - 1; k++)
                {
                    error += ((grid.w[i * (NY + 1) * NZ + j * NZ + k] - exact_solution.value_z(i, j, k + 0.5, t)) *
                              (grid.w[i * (NY + 1) * NZ + j * NZ + k] - exact_solution.value_z(i, j, k + 0.5, t)) *
                              DX * DY * DZ);
                }

                error += ((grid.w[i * (NY + 1) * NZ + j * NZ + NZ - 1] - exact_solution.value_z(i, j, (NZ - 0.5), t)) *
                          (grid.w[i * (NY + 1) * NZ + j * NZ + NZ - 1] - exact_solution.value_z(i, j, (NZ - 0.5), t)) *
                          DX * DY * DZ / 2);
            }

            error += ((grid.w[i * (NY + 1) * NZ + NY * NZ] - exact_solution.value_z(i, NY, 0.5, t)) *
                      (grid.w[i * (NY + 1) * NZ + NY * NZ] - exact_solution.value_z(i, NY, 0.5, t)) *
                      DX * DY * DZ / 4);

            for (size_t k = 1; k < NZ - 1; k++)
            {
                error += ((grid.w[i * (NY + 1) * NZ + NY * NZ + k] - exact_solution.value_z(i, NY, k + 0.5, t)) *
                          (grid.w[i * (NY + 1) * NZ + NY * NZ + k] - exact_solution.value_z(i, NY, k + 0.5, t)) *
                          DX * DY * DZ / 2);
            }

            error += ((grid.w[i * (NY + 1) * NZ + NY * NZ + NZ - 1] - exact_solution.value_z(i, NY, NZ - 0.5, t)) *
                      (grid.w[i * (NY + 1) * NZ + NY * NZ + NZ - 1] - exact_solution.value_z(i, NY, NZ - 0.5, t)) *
                      DX * DY * DZ / 4);
        }
    }

    // last slice (right face)
    {
        error += ((grid.w[NX * (NY + 1) * NZ] - exact_solution.value_z(NX, 0, 0.5, t)) *
                  (grid.w[NX * (NY + 1) * NZ] - exact_solution.value_z(NX, 0, 0.5, t)) *
                  DX * DY * DZ / 8);

        for (size_t k = 1; k < NZ - 1; k++)
        {
            error += ((grid.w[NX * (NY + 1) * NZ + k] - exact_solution.value_z(NX, 0, k + 0.5, t)) *
                      (grid.w[NX * (NY + 1) * NZ + k] - exact_solution.value_z(NX, 0, k + 0.5, t)) *
                      DX * DY * DZ / 4);
        }

        error += ((grid.w[NX * (NY + 1) * NZ + NZ - 1] - exact_solution.value_z(NX, 0, NZ - 0.5, t)) *
                  (grid.w[NX * (NY + 1) * NZ + NZ - 1] - exact_solution.value_z(NX, 0, NZ - 0.5, t)) *
                  DX * DY * DZ / 8);

        for (size_t j = 1; j < NY; j++)
        {
            error += ((grid.w[NX * (NY + 1) * NZ + j * NZ] - exact_solution.value_z(NX, j, 0.5, t)) *
                      (grid.w[NX * (NY + 1) * NZ + j * NZ] - exact_solution.value_z(NX, j, 0.5, t)) *
                      DX * DY * DZ / 4);
            for (size_t k = 1; k < NZ - 1; k++)
            {
                error += ((grid.w[NX * (NY + 1) * NZ + j * NZ + k] - exact_solution.value_z(NX, j, k + 0.5, t)) *
                          (grid.w[NX * (NY + 1) * NZ + j * NZ + k] - exact_solution.value_z(NX, j, k + 0.5, t)) *
                          DX * DY * DZ / 2);
            }
            error += ((grid.w[NX * (NY + 1) * NZ + j * NZ + NZ - 1] - exact_solution.value_z(NX, j, NZ - 0.5, t)) *
                      (grid.w[NX * (NY + 1) * NZ + j * NZ + NZ - 1] - exact_solution.value_z(NX, j, NZ - 0.5, t)) *
                      DX * DY * DZ / 4);
        }

        error += ((grid.w[NX * (NY + 1) * NZ + NY * NZ] - exact_solution.value_z(NX, NY, 0.5, t)) *
                  (grid.w[NX * (NY + 1) * NZ + NY * NZ] - exact_solution.value_z(NX, NY, 0.5, t)) *
                  DX * DY * DZ / 8);

        for (size_t k = 1; k < NZ - 1; k++)
        {
            error += ((grid.w[NX * (NY + 1) * NZ + NY * NZ + k] - exact_solution.value_z(NX, NY, k + 0.5, t)) *
                      (grid.w[NX * (NY + 1) * NZ + NY * NZ + k] - exact_solution.value_z(NX, NY, k + 0.5, t)) *
                      DX * DY * DZ / 4);
        }

        error += ((grid.w[NX * (NY + 1) * NZ + NY * NZ + NZ - 1] - exact_solution.value_z(NX, NY, NZ - 0.5, t)) *
                  (grid.w[NX * (NY + 1) * NZ + NY * NZ + NZ - 1] - exact_solution.value_z(NX, NY, NZ - 0.5, t)) *
                  DX * DY * DZ / 8);
    }

    return error;
}

void IcoNS::output()
{
}
