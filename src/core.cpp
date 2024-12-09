
#include <cmath>

#include "core.hpp"
#include <math.h>

#include <fstream>
#include <string>
#include <memory>
#include <mpi.h>

// #define OUTPUT
// #define OUTPUTERROR
#define VERBOSE

#ifdef VERBOSE
#include <chrono>
#endif
void IcoNS::preprocessing(/*std::string &input_file*/)

{
#ifdef VERBOSE
    std::cout << "*************************************************" << std::endl;
    std::cout << "Incompressible Navier-Stokes equation Solver" << std::endl
              << std::endl
              << std::endl;

    std::cout << "Solving for a Mesh of physical dimension (" << LX << "," << LY << "," << LZ << ") meters." << std::endl
              << "Number of partitions: " << NX << " nx, " << NY << " ny, " << NZ << " nz." << std::endl
              << "Dimension of a single cell:(" << DX << "," << DY << "," << DZ << ")." << std::endl
              << "Reynolds number: " << RE << std::endl
              << "Total lenght of simulation: " << T << " seconds, whit a time step of " << DT << " seconds." << std::endl

              << "------------------------------------------------------------" << std::endl
              << std::endl
              << "Reading Initial condition from file: Not implemented yet, setting all to 0." << std::endl
              << "Reading Boundary conditions from file: Not implemented yet, using default ones" << std::endl;
    std::cout << "*************************************************" << std::endl
              << std::endl;
#endif
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

    setParallelization();
    boundary.initializeBoundary(
        dim_x_x, dim_y_x, dim_x_y, dim_y_y,
        dim_x_z, dim_y_z, dim_z, dim_z_z,
        newDimX_x, newDimY_x, newDimX_y, newDimY_y,
        newDimX_z, newDimY_z);
}

void IcoNS::solve()
{
    Real time = 0.0;
    Real error;
    int i = 0;
#ifdef OUTPUTERROR
    Grid ERROR(grid);
#endif
// std::ofstream error_log("../resources/" + error_file);
#ifdef VERBOSE
    std::cout << "Starting solver" << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
#endif

    while (time < T)
    {
        boundary.update_boundary(grid_loc_x, grid_loc_y, grid_loc_z, time);

        MPI_Barrier(cart_comm);
        exchangeData_x(grid_loc_x);
        exchangeData_y(grid_loc_y);
        exchangeData_z(grid_loc_z);

        auto x = L2_error(time); // every processor calculates his error not counting ghosts(and then some sort of reduce?)
        MPI_Barrier(cart_comm);
        error = 0.0;
        MPI_Reduce(&x, &error, 1, MPI_DOUBLE, MPI_SUM, 0, cart_comm);
        if (rank == 0)
        {
            std::cout << " errorx: " << error << std::endl;
        }
        // reduce
        solve_time_step(time); // adapt cycles to skip ghosts
        MPI_Barrier(cart_comm);
        time += DT;
        i++;
    }
/*while (time < T)
{
    /*Check::Confront(grid,exact_solution,time,U);
    int p;
    std::cin >> p;*/
/*
boundary.update_boundary(grid.u, grid.v, grid.w, time);

// csv file w/ "," delimiter: time step, iter, L2_error
std::cout << time << "," << i << "," << L2_error(time) << std::endl;
solve_time_step(time);
// output();
time += DT;
i++;
}*/
// error = L2_error(time);
// error_log << time << "," << i << "," << error << std::endl;
#ifdef VERBOSE
    if(rank==0){
        std::cout << "At time: " << time << "s of " << T << "s the L2 norm of the error is: " << error << std::endl;
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = end - start;
        std::cout << std::endl
                << "Time: " << duration.count() << std::endl;
    }

#endif
}

Real IcoNS::L2_error(const Real t)
{
    Real error = 0.0;

    error += error_comp_X(t);
    // error += error_comp_Y(t);
    // error += error_comp_Z(t);

    // std::cout << error_comp_X(t) << std::endl;
    // std::cout << error_comp_Y(t) << std::endl;
    // std::cout << error_comp_Z(t) << std::endl << std::endl;

    return sqrt(error);
}

// TODO:
Real IcoNS::error_comp_X(const Real t)
{
    Real error = 0.0;
    int offset_x = coords[0] + dim_x_x;
    int offset_y = (PY - 1 - coords[1]) * dim_y_x;
    // first slice (left face)
    if (lbx)
    {
        {
            if (lby)
            {
                error += ((grid_loc_x[0 + (newDimY_x + 1) * dim_z] - exact_solution.value_x(0.5 + offset_x, 0 + offset_y, 0, t)) *
                          (grid_loc_x[0 + (newDimY_x + 1) * dim_z] - exact_solution.value_x(0.5 + offset_x, 0 + offset_y, 0, t)) *
                          DX * DY * DZ / 8);

                for (size_t k = 1; k < NZ; k++)
                {
                    error += ((grid_loc_x[(newDimY_x + 1) * dim_z + k] - exact_solution.value_x(0.5 + offset_x, 0 + offset_y, k, t)) *
                              (grid_loc_x[(newDimY_x + 1) * dim_z + k] - exact_solution.value_x(0.5 + offset_x, 0 + offset_y, k, t)) *
                              DX * DY * DZ / 4);
                }
            }
            error += ((grid_loc_x[(newDimY_x + 1) * dim_z + NZ] - exact_solution.value_x(0.5 + offset_x, 0 + offset_y, NZ, t)) *
                      (grid_loc_x[(newDimY_x + 1) * dim_z + NZ] - exact_solution.value_x(0.5 + offset_x, 0 + offset_y, NZ, t)) *
                      DX * DY * DZ / 8);

            for (size_t j = 1 + lby; j < newDimY_x - 1 - rby; j++)
            {
                error += ((grid_loc_x[(newDimY_x)*dim_z + j * (NZ + 1)] - exact_solution.value_x(0.5 + offset_x, j + offset_y, 0, t)) *
                          (grid_loc_x[(newDimY_x)*dim_z + j * (NZ + 1)] - exact_solution.value_x(0.5 + offset_x, j + offset_y, 0, t)) *
                          DX * DY * DZ / 4);
                for (size_t k = 1; k < NZ; k++)
                {
                    error += ((grid_loc_x[(newDimY_x)*dim_z + j * (NZ + 1) + k] - exact_solution.value_x(0.5 + offset_x, j + offset_y, k, t)) *
                              (grid_loc_x[(newDimY_x)*dim_z + j * (NZ + 1) + k] - exact_solution.value_x(0.5 + offset_x, j + offset_y, k, t)) *
                              DX * DY * DZ / 2);
                }
                error += ((grid_loc_x[(newDimY_x)*dim_z + j * (NZ + 1) + NZ] - exact_solution.value_x(0.5 + offset_x, j + offset_y, NZ, t)) *
                          (grid_loc_x[(newDimY_x)*dim_z + j * (NZ + 1) + NZ] - exact_solution.value_x(0.5 + offset_x, j + offset_y, NZ, t)) *
                          DX * DY * DZ / 4);
            }
            if (rby)
            {
                error += ((grid_loc_x[(newDimY_x)*dim_z + (newDimY_x - 2) * (NZ + 1)] - exact_solution.value_x(0.5 + offset_x, NY + offset_y, 0, t)) *
                          (grid_loc_x[(newDimY_x)*dim_z + (newDimY_x - 2) * (NZ + 1)] - exact_solution.value_x(0.5 + offset_x, NY + offset_y, 0, t)) *
                          DX * DY * DZ / 8);

                for (size_t k = 1; k < NZ; k++)
                {
                    error += ((grid_loc_x[(newDimY_x)*dim_z + (newDimY_x - 2) * (NZ + 1) + k] - exact_solution.value_x(0.5 + offset_x, NY + offset_y, k, t)) *
                              (grid_loc_x[(newDimY_x)*dim_z + (newDimY_x - 2) * (NZ + 1) + k] - exact_solution.value_x(0.5 + offset_x, NY + offset_y, k, t)) *
                              DX * DY * DZ / 4);
                }

                error += ((grid_loc_x[(newDimY_x)*dim_z + (newDimY_x - 2) * (NZ + 1) + NZ] - exact_solution.value_x(0.5 + offset_x, NY + offset_y, NZ, t)) *
                          (grid_loc_x[(newDimY_x)*dim_z + (newDimY_x - 2) * (NZ + 1) + NZ] - exact_solution.value_x(0.5 + offset_x, NY + offset_y, NZ, t)) *
                          DX * DY * DZ / 8);
            }
        }
    }

    // middle slices
    {
        for (size_t i = 1 + lbx; i < newDimX_x - 1 - rbx; i++)
        {
            if (lby)
            {
                error += ((grid_loc_x[i * newDimY_x * (NZ + 1) + dim_z] - exact_solution.value_x(i + 0.5 + offset_x, 0 + offset_y, 0, t)) *
                          (grid_loc_x[i * newDimY_x * (NZ + 1) + dim_z] - exact_solution.value_x(i + 0.5 + offset_x, 0 + offset_y, 0, t)) *
                          DX * DY * DZ / 4);

                for (size_t k = 1; k < NZ; k++)
                {
                    error += ((grid_loc_x[i * newDimY_x * (NZ + 1) + dim_z + k] - exact_solution.value_x(i + 0.5 + offset_x, 0 + offset_y, k, t)) *
                              (grid_loc_x[i * newDimY_x * (NZ + 1) + dim_z + k] - exact_solution.value_x(i + 0.5 + offset_x, 0 + offset_y, k, t)) *
                              DX * DY * DZ / 2);
                }

                error += ((grid_loc_x[i * newDimY_x * (NZ + 1) + dim_z + NZ] - exact_solution.value_x(i + 0.5 + offset_x, 0 + offset_y, NZ, t)) *
                          (grid_loc_x[i * newDimY_x * (NZ + 1) + dim_z + NZ] - exact_solution.value_x(i + 0.5 + offset_x, 0 + offset_y, NZ, t)) *
                          DX * DY * DZ / 4);
            }
            for (size_t j = 1 + lby; j < newDimY_x - 1 - rby; j++)
            {
                error += ((grid_loc_x[i * newDimY_x * (NZ + 1) + j * (NZ + 1)] - exact_solution.value_x(i + 0.5 + offset_x, j + offset_y, 0, t)) *
                          (grid_loc_x[i * newDimY_x * (NZ + 1) + j * (NZ + 1)] - exact_solution.value_x(i + 0.5 + offset_x, j + offset_y, 0, t)) *
                          DX * DY * DZ / 2);

                for (size_t k = 1; k < NZ; k++)
                {
                    error += ((grid_loc_x[i * newDimY_x * (NZ + 1) + j * (NZ + 1) + k] - exact_solution.value_x(i + 0.5 + offset_x, j + offset_y, k, t)) *
                              (grid_loc_x[i * newDimY_x * (NZ + 1) + j * (NZ + 1) + k] - exact_solution.value_x(i + 0.5 + offset_x, j + offset_y, k, t)) *
                              DX * DY * DZ);
                }

                error += ((grid_loc_x[i * newDimY_x * (NZ + 1) + j * (NZ + 1) + NZ] - exact_solution.value_x(i + 0.5 + offset_x, j + offset_y, NZ, t)) *
                          (grid_loc_x[i * newDimY_x * (NZ + 1) + j * (NZ + 1) + NZ] - exact_solution.value_x(i + 0.5 + offset_x, j + offset_y, NZ, t)) *
                          DX * DY * DZ / 2);
            }
            if (rby)
            {
                error += ((grid_loc_x[i * newDimY_x * (NZ + 1) + (newDimY_x - 2) * (NZ + 1)] - exact_solution.value_x(i + 0.5 + offset_x, NY + offset_y, 0, t)) *
                          (grid_loc_x[i * newDimY_x * (NZ + 1) + (newDimY_x - 2) * (NZ + 1)] - exact_solution.value_x(i + 0.5 + offset_x, NY + offset_y, 0, t)) *
                          DX * DY * DZ / 4);

                for (size_t k = 1; k < NZ; k++)
                {
                    error += ((grid_loc_x[i * newDimY_x * (NZ + 1) + (newDimY_x - 2) * (NZ + 1) + k] - exact_solution.value_x(i + 0.5 + offset_x, NY + offset_y, k, t)) *
                              (grid_loc_x[i * newDimY_x * (NZ + 1) + (newDimY_x - 2) * (NZ + 1) + k] - exact_solution.value_x(i + 0.5 + offset_x, NY + offset_y, k, t)) *
                              DX * DY * DZ / 2);
                }

                error += ((grid_loc_x[i * newDimY_x * (NZ + 1) + (newDimY_x - 2) * (NZ + 1) + NZ] - exact_solution.value_x(i + 0.5 + offset_x, NY + offset_y, NZ, t)) *
                          (grid_loc_x[i * newDimY_x * (NZ + 1) + (newDimY_x - 2) * (NZ + 1) + NZ] - exact_solution.value_x(i + 0.5 + offset_x, NY + offset_y, NZ, t)) *
                          DX * DY * DZ / 4);
            }
        }
    }

    // last slice (right face)
    if (rbx)
    {
        {
            if (lby)
            {
                error += ((grid_loc_x[(newDimX_x - 2) * newDimY_x * (NZ + 1) + dim_z] - exact_solution.value_x((NX - 0.5) + offset_x, 0 + offset_y, 0, t)) *
                          (grid_loc_x[(newDimX_x - 2) * newDimY_x * (NZ + 1) + dim_z] - exact_solution.value_x((NX - 0.5) + offset_x, 0 + offset_y, 0, t)) *
                          DX * DY * DZ / 8);

                for (size_t k = 1; k < NZ; k++)
                {
                    error += ((grid_loc_x[(newDimX_x - 2) * newDimY_x * (NZ + 1) + dim_z + k] - exact_solution.value_x((NX - 0.5) + offset_x, 0 + offset_y, k, t)) *
                              (grid_loc_x[(newDimX_x - 2) * newDimY_x * (NZ + 1) + dim_z + k] - exact_solution.value_x((NX - 0.5) + offset_x, 0 + offset_y, k, t)) *
                              DX * DY * DZ / 4);
                }

                error += ((grid_loc_x[(newDimX_x - 2) * newDimY_x * (NZ + 1) + dim_z + NZ] - exact_solution.value_x((NX - 0.5) + offset_x, 0 + offset_y, NZ, t)) *
                          (grid_loc_x[(newDimX_x - 2) * newDimY_x * (NZ + 1) + dim_z + NZ] - exact_solution.value_x((NX - 0.5) + offset_x, 0 + offset_y, NZ, t)) *
                          DX * DY * DZ / 8);
            }
            for (size_t j = 1 + lby; j < newDimY_x - 1 - rby; j++)
            {
                error += ((grid_loc_x[(newDimX_x - 2) * newDimY_x * (NZ + 1) + j * (NZ + 1)] - exact_solution.value_x((NX - 0.5) + offset_x, j + offset_y, 0, t)) *
                          (grid_loc_x[(newDimX_x - 2) * newDimY_x * (NZ + 1) + j * (NZ + 1)] - exact_solution.value_x((NX - 0.5) + offset_x, j + offset_y, 0, t)) *
                          DX * DY * DZ / 4);
                for (size_t k = 1; k < NZ; k++)
                {
                    error += ((grid_loc_x[(newDimX_x - 2) * newDimY_x * (NZ + 1) + j * (NZ + 1) + k] - exact_solution.value_x((NX - 0.5) + offset_x, j + offset_y, k, t)) *
                              (grid_loc_x[(newDimX_x - 2) * newDimY_x * (NZ + 1) + j * (NZ + 1) + k] - exact_solution.value_x((NX - 0.5) + offset_x, j + offset_y, k, t)) *
                              DX * DY * DZ / 2);
                }
                error += ((grid_loc_x[(newDimX_x - 2) * newDimY_x * (NZ + 1) + j * (NZ + 1) + NZ] - exact_solution.value_x((NX - 0.5) + offset_x, j + offset_y, NZ, t)) *
                          (grid_loc_x[(newDimX_x - 2) * newDimY_x * (NZ + 1) + j * (NZ + 1) + NZ] - exact_solution.value_x((NX - 0.5) + offset_x, j + offset_y, NZ, t)) *
                          DX * DY * DZ / 4);
            }
            if (rby)
            {
                error += ((grid_loc_x[(newDimX_x - 2) * newDimY_x * (NZ + 1) + (newDimY_x - 2) * (NZ + 1)] - exact_solution.value_x((NX - 0.5) + offset_x, NY + offset_y, 0, t)) *
                          (grid_loc_x[(newDimX_x - 2) * newDimY_x * (NZ + 1) + (newDimX_x - 2) * (NZ + 1)] - exact_solution.value_x((NX - 0.5) + offset_x, NY + offset_y, 0, t)) *
                          DX * DY * DZ / 8);

                for (size_t k = 1; k < NZ; k++)
                {
                    error += ((grid_loc_x[(newDimX_x - 2) * newDimY_x * (NZ + 1) + (newDimX_x - 2) * (NZ + 1) + k] - exact_solution.value_x((NX - 0.5) + offset_x, NY + offset_y, k, t)) *
                              (grid_loc_x[(newDimX_x - 2) * newDimY_x * (NZ + 1) + (newDimX_x - 2) * (NZ + 1) + k] - exact_solution.value_x((NX - 0.5) + offset_x, NY + offset_y, k, t)) *
                              DX * DY * DZ / 4);
                }

                error += ((grid_loc_x[(newDimX_x - 2) * newDimY_x * (NZ + 1) + (newDimX_x - 2) * (NZ + 1) + NZ] - exact_solution.value_x((NX - 0.5) + offset_x, NY + offset_y, NZ, t)) *
                          (grid_loc_x[(newDimX_x - 2) * newDimY_x * (NZ + 1) + (newDimX_x - 2) * (NZ + 1) + NZ] - exact_solution.value_x((NX - 0.5) + offset_x, NY + offset_y, NZ, t)) *
                          DX * DY * DZ / 8);
            }
        }
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

void IcoNS::setParallelization()
{

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    std::cout << rank << ", " << size << ", " << dims[0] << ", " << dims[1] << ", " << periods[0] << ", " << periods[1] << std::endl;

    // Create a Cartesian topology (2D)
    // MPI_Dims_create(size, 2, dims);

    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &cart_comm);

    MPI_Cart_coords(cart_comm, rank, 2, coords);

    MPI_Cart_shift(cart_comm, 0, 1, &neighbors[0], &neighbors[2]);

    MPI_Cart_shift(cart_comm, 1, 1, &neighbors[1], &neighbors[3]);

    dim_x_x = NX / PX;
    dim_y_x = (NY + 1) / PY;
    dim_x_y = (NX + 1) / PX;
    dim_y_y = NY / PY;
    dim_x_z = (NX + 1) / PX;
    dim_y_z = (NY + 1) / PY;
    dim_z = NZ + 1;
    dim_z_z = NZ;

    if (NX % PX != 0 && coords[0] == dim_x_x - 1)
        dim_x_x++;
    if ((NY + 1) % PY != 0 && coords[1] == dim_y_x - 1)
    {
        dim_y_x++;
    }

    if ((NX + 1) % PX != 0 && coords[0] == dim_x_y - 1)
        dim_x_y++;
    if ((NY) % PY != 0 && coords[1] == dim_y_y - 1)
    {
        dim_y_y++;
    }

    if ((NX + 1) % PX != 0 && coords[0] == dim_x_z - 1)
        dim_x_z++;
    if ((NY + 1) % PY != 0 && coords[1] == dim_y_z - 1)
    {
        dim_y_z++;
    }
    newDimX_x = dim_x_x + 2;
    newDimY_x = dim_y_x + 2;
    newDimX_y = dim_x_y + 2;
    newDimY_y = dim_y_y + 2;
    newDimX_z = dim_x_z + 2;
    newDimY_z = dim_y_z + 2;

    grid_loc_x.resize(newDimX_x * newDimY_x * (NZ + 1), 0.0);
    grid_loc_y.resize(newDimX_y * newDimY_y * (NZ + 1), 0.0);
    grid_loc_z.resize(newDimX_z * newDimY_z * (NZ), 0.0);
    Y2_x.resize(newDimX_x * newDimY_x * (NZ + 1), 0.0);
    Y2_y.resize(newDimX_y * newDimY_y * (NZ + 1), 0.0);
    Y2_z.resize(newDimX_z * newDimY_z * (NZ), 0.0);
    Y3_x.resize(newDimX_x * newDimY_x * (NZ + 1), 0.0);
    Y3_y.resize(newDimX_y * newDimY_y * (NZ + 1), 0.0);
    Y3_z.resize(newDimX_z * newDimY_z * (NZ), 0.0);

    if (coords[0] == 0)
        lbx++;

    if (coords[0] == PX - 1)
        rbx++;

    if (coords[1] == PY - 1)
        lby++;

    if (coords[1] == 0)
        rby++;

    boundary.setBoundaryOffsets(lbx, rbx, lby, rby);
    boundary.setCoords(coords);

    newDimX_x = (dim_x_x + 2);
    newDimY_x = (dim_y_x + 2);
    newDimX_y = (dim_x_y + 2);
    newDimY_y = (dim_y_y + 2);
    newDimX_z = (dim_x_z + 2);
    newDimY_z = (dim_y_z + 2);

    /*
    glob_address_x_x = (i-1) +coords[0]*dim_x_x;
    glob_address_y_x =(j-1)+ coords[1]*dim_y_x;


    glob_address_x_y = (i-1) +coords[0]*dim_x_y;
    glob_address_y_y =(j-1)+ coords[1]*dim_y_y;
    */

    MPI_Type_vector(dim_x_x, dim_z, (newDimY_x)*dim_z, MPI_INT, &MPI_face_x_x);
    MPI_Type_commit(&MPI_face_x_x);

    MPI_Type_vector(1, dim_z * newDimY_x, 0, MPI_INT, &MPI_face_y_x);
    MPI_Type_commit(&MPI_face_y_x);

    MPI_Type_vector(dim_x_y, dim_z, (newDimY_y)*dim_z, MPI_INT, &MPI_face_x_y);
    MPI_Type_commit(&MPI_face_x_y);

    MPI_Type_vector(1, dim_z * newDimY_y, 0, MPI_INT, &MPI_face_y_y);
    MPI_Type_commit(&MPI_face_y_y);

    MPI_Type_vector(dim_x_z, dim_z_z, (newDimY_z)*dim_z_z, MPI_INT, &MPI_face_x_z);
    MPI_Type_commit(&MPI_face_x_z);

    MPI_Type_vector(1, dim_z_z * newDimY_z, 0, MPI_INT, &MPI_face_y_z);
    MPI_Type_commit(&MPI_face_y_z);
}

void IcoNS::exchangeData_x(std::vector<Real> &grid_loc)
{
    if (neighbors[0] != -2 && neighbors[1] != -2 && neighbors[2] != -2 && neighbors[3] != -2)
    {

        // (x, y-1) <- (x, y)
        MPI_Isend(&grid_loc[newDimY_x * dim_z], 1, MPI_face_y_x, neighbors[0], rank, cart_comm, &req1);
        MPI_Irecv(&grid_loc[(dim_z)*newDimY_x * (newDimX_x - 1)], 1, MPI_face_y_x, neighbors[2], neighbors[2], cart_comm, &req1);

        // (x,y) -> (x, y+1)
        MPI_Isend(&grid_loc[newDimY_x * dim_z * (newDimX_x - 2)], 1, MPI_face_y_x, neighbors[2], rank, cart_comm, &req2);
        MPI_Irecv(&grid_loc[0], 1, MPI_face_y_x, neighbors[0], neighbors[0], cart_comm, &req2);

        // (x-1, y)
        //   ^
        //   |
        // (x, y)
        MPI_Isend(&grid_loc[dim_z * newDimY_x + dim_z], 1, MPI_face_x_x, neighbors[1], rank, cart_comm, &req3);
        MPI_Irecv(&grid_loc[dim_z * newDimY_x + (newDimY_x - 1) * dim_z], 1, MPI_face_x_x, neighbors[3], neighbors[3], cart_comm, &req3);

        // (x, y)
        //   |
        //   V
        // (x+1, y)
        MPI_Isend(&grid_loc[dim_z * newDimY_x + (newDimY_x - 2) * dim_z], 1, MPI_face_x_x, neighbors[3], rank, cart_comm, &req4);
        MPI_Irecv(&grid_loc[dim_z * newDimY_x], 1, MPI_face_x_x, neighbors[1], neighbors[1], cart_comm, &req4);
        MPI_Barrier(cart_comm);
    }
}

void IcoNS::exchangeData_y(std::vector<Real> &grid_loc)
{
    if (neighbors[0] != -2 && neighbors[1] != -2 && neighbors[2] != -2 && neighbors[3] != -2)
    {

        // (x, y-1) <- (x, y)
        MPI_Isend(&grid_loc[newDimY_y * dim_z], 1, MPI_face_y_y, neighbors[0], rank, cart_comm, &req1);
        MPI_Irecv(&grid_loc[(dim_z)*newDimY_y * (newDimX_y - 1)], 1, MPI_face_y_y, neighbors[2], neighbors[2], cart_comm, &req1);

        // (x,y) -> (x, y+1)
        MPI_Isend(&grid_loc[newDimY_y * dim_z * (newDimX_y - 2)], 1, MPI_face_y_y, neighbors[2], rank, cart_comm, &req2);
        MPI_Irecv(&grid_loc[0], 1, MPI_face_y_y, neighbors[0], neighbors[0], cart_comm, &req2);

        // (x-1, y)
        //   ^
        //   |
        // (x, y)
        MPI_Isend(&grid_loc[dim_z * newDimY_y + dim_z], 1, MPI_face_x_y, neighbors[1], rank, cart_comm, &req3);
        MPI_Irecv(&grid_loc[dim_z * newDimY_y + (newDimY_y - 1) * dim_z], 1, MPI_face_x_y, neighbors[3], neighbors[3], cart_comm, &req3);

        // (x, y)
        //   |
        //   V
        // (x+1, y)
        MPI_Isend(&grid_loc[dim_z * newDimY_y + (newDimY_y - 2) * dim_z], 1, MPI_face_x_y, neighbors[3], rank, cart_comm, &req4);
        MPI_Irecv(&grid_loc[dim_z * newDimY_y], 1, MPI_face_x_y, neighbors[1], neighbors[1], cart_comm, &req4);
        MPI_Barrier(cart_comm);
    }
}

void IcoNS::exchangeData_z(std::vector<Real> &grid_loc)
{
    if (neighbors[0] != -2 && neighbors[1] != -2 && neighbors[2] != -2 && neighbors[3] != -2)
    {

        // (x, y-1) <- (x, y)
        MPI_Isend(&grid_loc[newDimY_z * dim_z], 1, MPI_face_y_z, neighbors[0], rank, cart_comm, &req1);
        MPI_Irecv(&grid_loc[(dim_z)*newDimY_z * (newDimX_z - 1)], 1, MPI_face_y_z, neighbors[2], neighbors[2], cart_comm, &req1);

        // (x,y) -> (x, y+1)
        MPI_Isend(&grid_loc[newDimY_z * dim_z * (newDimX_z - 2)], 1, MPI_face_y_z, neighbors[2], rank, cart_comm, &req2);
        MPI_Irecv(&grid_loc[0], 1, MPI_face_y_z, neighbors[0], neighbors[0], cart_comm, &req2);

        // (x-1, y)
        //   ^
        //   |
        // (x, y)
        MPI_Isend(&grid_loc[dim_z * newDimY_z + dim_z], 1, MPI_face_x_z, neighbors[1], rank, cart_comm, &req3);
        MPI_Irecv(&grid_loc[dim_z * newDimY_z + (newDimY_z - 1) * dim_z], 1, MPI_face_x_z, neighbors[3], neighbors[3], cart_comm, &req3);

        // (x, y)
        //   |
        //   V
        // (x+1, y)
        MPI_Isend(&grid_loc[dim_z * newDimY_z + (newDimY_z - 2) * dim_z], 1, MPI_face_x_z, neighbors[3], rank, cart_comm, &req4);
        MPI_Irecv(&grid_loc[dim_z * newDimY_z], 1, MPI_face_x_z, neighbors[1], neighbors[1], cart_comm, &req4);
        MPI_Barrier(cart_comm);
    }
}
// void IcoNS::output()
// {
// }
