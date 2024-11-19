
#include <cmath>

#include "core.hpp"
#include <math.h>

#include <fstream>
#include <string>
#include <memory>

//#define OUTPUT
#define OUTPUTERROR
//#define VERBOSE
#ifdef VERBOSE
    #include <chrono>
#endif
void IcoNS::preprocessing(/*std::string &input_file*/)

{
    #ifdef VERBOSE
        std::cout << "*************************************************" << std::endl;
        std::cout << "Incompressible Navier-Stokes equation Solver" << std::endl << std::endl << std::endl;

        std::cout << "Solving for a Mesh of physical dimension (" << lx << "," << ly << "," << lz <<") meters." << std::endl
        << "Number of partitions: " << nx << " nx, " << ny << " ny, "<< nz << " nz." << std::endl
        << "Dimension of a single cell:(" << dx <<"," << dy << "," << dz <<")." <<std::endl
        << "Reynolds number: " << Re << std::endl
        << "Total lenght of simulation: " << T << " seconds, whit a time step of " << dt << " seconds." << std::endl

        << "------------------------------------------------------------" << std::endl << std::endl
        <<"Reading Initial condition from file: Not implemented yet, setting all to 0." << std::endl
        <<"Reading Boundary conditions from file: Not implemented yet, using default ones" <<std::endl;        
        std::cout << "*************************************************" << std::endl << std::endl;
    #endif
    // boundary
    auto u_func = std::make_shared<Dirichlet>([&](Real x, Real y, Real z, Real t)
                                              { return std::sin((x + 0.5) * dx) * std::cos(y * dy) * std::sin(z * dz) * std::sin(t); });
    auto v_func = std::make_shared<Dirichlet>([&](Real x, Real y, Real z, Real t)
                                              { return std::cos(x * dx) * std::sin((y + 0.5) * dy) * std::sin(z * dz) * std::sin(t); });
    auto w_func = std::make_shared<Dirichlet>([&](Real x, Real y, Real z, Real t)
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
    Real error = 0.0;
    double time = 0.0;
    int i = 0;
    #ifdef OUTPUTERROR
    Grid ERROR(grid);
    #endif
    std::ofstream error_log("../resources/" + error_file);
    #ifdef VERBOSE
    std::cout << "Starting solver" << std::endl;
    auto start =std::chrono::high_resolution_clock::now();
    #endif
    while (time < T)
    {   
        boundary.update_boundary(grid.u, grid.v, grid.w, time);
        #ifdef OUTPUT
            std::string filenameerroru = "../resources/paraview/u/" + std::to_string(time) +".vtk";
            std::string filenameerrorv = "../resources/paraview/v/" + std::to_string(time) +".vtk";
            std::string filenameerrorw = "../resources/paraview/w/" + std::to_string(time) +".vtk";
            output_u(filenameerroru, grid);
            output_v(filenameerrorv, grid);
            output_w(filenameerrorw, grid);
        #endif
        #ifdef OUTPUTERROR
            if(i%5== 0){
                std::string filenameerroru = "../resources/paraviewerror/u/" + std::to_string(time) +".vtk";
                std::string filenameerrorv = "../resources/paraviewerror/v/" + std::to_string(time) +".vtk";
                std::string filenameerrorw = "../resources/paraviewerror/w/" + std::to_string(time) +".vtk";
                for (size_t i = 0; i < nx; i++)
                {
                    for (size_t j = 0; j < ny+1; j++)
                    {
                        for (size_t k = 0; k < nz +1; k++)
                        {
                            ERROR.u[i*(ny+1)*(nz+1) + j*(nz+1) + k] = std::abs(grid.u[i*(ny+1)*(nz+1) + j*(nz+1) + k] - exact_solution.value_x(i+0.5,j,k,time));
                        }
                        
                    }   
                }
                for (size_t i = 0; i < nx+1; i++)
                {
                    for (size_t j = 0; j < ny; j++)
                    {
                        for (size_t k = 0; k < nz +1; k++)
                        {
                            ERROR.v[i*(ny)*(nz+1) + j*(nz+1) + k] = std::abs(grid.v[i*(ny)*(nz+1) + j*(nz+1) + k] - exact_solution.value_y(i,j+0.5,k,time));
                        }
                        
                    }
                    
                }
                for (size_t i = 0; i < nx + 1; i++)
                {
                    for (size_t j = 0; j < ny + 1; j++)
                    {
                        for (size_t k = 0; k < nz; k++)
                        {
                            ERROR.w[i*(ny+1)*(nz) + j*(nz) + k] = std::abs(grid.w[i*(ny+1)*(nz) + j*(nz) + k] - exact_solution.value_z(i,j,k+0.5,time)  );
                        }
                        
                    }
                    
                }
                
                output_u(filenameerroru, ERROR);
                output_v(filenameerrorv, ERROR);
                output_w(filenameerrorw, ERROR);
            }
        #endif
        
        
        //Check::Confront(grid,exact_solution,time,W); int p;std::cin >> p;
        error = L2_error(time);
        #ifdef VERBOSE
            std::cout << "At time: " << time << "s of " << T << "s the L2 norm of the error is: "<< error << std::endl;
            auto tss =std::chrono::high_resolution_clock::now();
        #endif
        error_log << time << "," << i << "," << error << std::endl;
        solve_time_step(time);
        // output();

        
        time += dt;
        i++;
        #ifdef VERBOSE
        auto tse =std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = tse-tss; 
        std::cout << std::endl << "Time: " << duration.count() << std::endl;
        std::cout << "Time per cell: " << duration.count() / (nx*ny*nz) << std::endl;
        #endif
        
    }
    error = L2_error(time);
    error_log << time << "," << i << "," << error << std::endl;
    #ifdef VERBOSE
        std::cout << "At time: " << time << "s of " << T << "s the L2 norm of the error is: "<< error << std::endl;
        auto end =std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = end-start; 
        std::cout << std::endl << "Time: " << duration.count() << std::endl;
        
    #endif
}

void IcoNS::output_u(const std::string& filename, Grid& print) {
        
        std::filesystem::path filepath(filename);
        std::filesystem::create_directories(filepath.parent_path());
        
        std::ofstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Error opening file for writing: " << filename << std::endl;
            return;
        }

        file << "# vtk DataFile Version 3.0\n";
        file << "3D structured grid\n";
        file << "ASCII\n";
        file << "DATASET STRUCTURED_POINTS\n";
        file << "DIMENSIONS " << print.nx << " " << print.ny + 1 << " " << print.nz + 1 << "\n";
        file << "ORIGIN 0 0 0\n";
        file << "SPACING " << dx << " " << dy << " " << dz << "\n";
        file << "POINT_DATA " << (print.nx) * (print.ny + 1) * (print.nz + 1) << "\n";

        // Stampa delle velocità U
        file << "SCALARS velocity_u double 1\n";
        file << "LOOKUP_TABLE default\n";
        
        for(size_t k = 0 ; k < nz+1 ; k++){
            for (size_t j = 0; j < ny+1; j++)
            {
                for (int i = nx-1; i >= 0; i--)
                {   
                    file << print.u[i*(ny+1)*(nz+1) + j * (nz+1) + k] << "\n";
                }
                
            }
        }

        file.close();
        std::cout << "File " << filename << " scritto con successo." << std::endl;
    }
    void IcoNS::output_v(const std::string& filename, Grid& print) {
        
        std::filesystem::path filepath(filename);
        std::filesystem::create_directories(filepath.parent_path());
        
        std::ofstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Error opening file for writing: " << filename << std::endl;
            return;
        }

        file << "# vtk DataFile Version 3.0\n";
        file << "3D structured grid\n";
        file << "ASCII\n";
        file << "DATASET STRUCTURED_POINTS\n";
        file << "DIMENSIONS " << print.nx +1 << " " << print.ny  << " " << print.nz + 1 << "\n";
        file << "ORIGIN 0 0 0\n";
        file << "SPACING " << dx << " " << dy << " " << dz << "\n";
        file << "POINT_DATA " << (print.nx+1) * (print.ny) * (print.nz + 1) << "\n";

        // Stampa delle velocità V
        file << "SCALARS velocity_v double 1\n";
        file << "LOOKUP_TABLE default\n";
        for(int k = nz ; k >= 0 ; k--)
        {
            for (size_t j = 0; j < ny; j++)
            {
                for (size_t i = 0; i <nx +1; i++)
                {   
                    file << print.v[i*(ny)*(nz+1) + j * (nz+1) + k] << "\n";
                }
                
            }
        }
        // for (size_t i = 0; i < print.v.size(); ++i) {
        //     file << print.v[i] << "\n";
        // }

        file.close();
        std::cout << "File " << filename << " scritto con successo." << std::endl;
    }
    void IcoNS::output_w(const std::string& filename, Grid& print) {
        
        std::filesystem::path filepath(filename);
        std::filesystem::create_directories(filepath.parent_path());
        
        std::ofstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Error opening file for writing: " << filename << std::endl;
            return;
        }

        file << "# vtk DataFile Version 3.0\n";
        file << "3D structured grid\n";
        file << "ASCII\n";
        file << "DATASET STRUCTURED_POINTS\n";
        file << "DIMENSIONS " << print.nx +1 << " " << print.ny + 1 << " " << print.nz  << "\n";
        file << "ORIGIN 0 0 0\n";
        file << "SPACING " << dx << " " << dy << " " << dz << "\n";
        file << "POINT_DATA " << (print.nx+1) * (print.ny + 1) * (print.nz ) << "\n";

        // Stampa delle velocità W
        file << "SCALARS velocity_w double 1\n";
        file << "LOOKUP_TABLE default\n";
        for(int k = nz-1 ; k >= 0 ; k--)
        {
            for (size_t j = 0; j < ny+1; j++)
            {
                for (size_t i = 0; i <nx +1; i++)
                {   
                    file << print.w[i*(ny+1)*(nz) + j * (nz) + k] << "\n";
                }
                
            }
        }

        file.close();
        std::cout << "File " << filename << " scritto con successo." << std::endl;
    }

Real IcoNS::L2_error(const Real t)
{
    Real error = 0.0;

    error += error_comp_X(t);
    error += error_comp_Y(t);
    error += error_comp_Z(t);

    /*std::cout << error_comp_X(t) << std::endl;
    std::cout << error_comp_Y(t) << std::endl;
    std::cout << error_comp_Z(t) << std::endl << std::endl;*/

    return sqrt(error);
}

Real IcoNS::error_comp_X(const Real t)
{
    Real error = 0.0;
    
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

Real IcoNS::error_comp_Y(const Real t)
{
    Real error = 0.0;

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

Real IcoNS::error_comp_Z(const Real t)
{
    Real error = 0.0;
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

// void IcoNS::output()
// {
// }
