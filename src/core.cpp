
#include <cmath>
#include "core.hpp"
#include <math.h>
#include <fstream>
#include <string>
#include <memory>
#include <mpi.h>

// #define OUTPUT
// #define OUTPUTERROR
//#define VERBOSE

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

    for (int i = 0; i < 6 /*nfaces*/; i++)
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

    if (NX % PX != 0 && coords[0] == PX- 1)
    {
        dim_x_x++;
    }
    if ((NY + 1) % PY != 0 && coords[1] == 0)
    {
        dim_y_x++;
    }

    if ((NX + 1) % PX != 0 && coords[0] == PX - 1)
        dim_x_y++;
    if ((NY) % PY != 0 && coords[1] == 0)
    {
        dim_y_y++;
    }

    if ((NX + 1) % PX != 0 && coords[0] == PX - 1)
        dim_x_z++;
    if ((NY + 1) % PY != 0 && coords[1] == 0)
    {
        dim_y_z++;
    }

    //TODO: controllare che la riga in più non porti errori agli address globali
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
    /*
    //Resize only for boundary faces
        bool is_edge1 = coords[0] == PX - 1;
        bool is_edge2 = coords[1] == 0;

        if (is_edge1 == true || is_edge2 == true) {
            grid_loc_x.resize(newDimX_x * newDimY_x * (NZ + 1), 0.0);
            grid_loc_y.resize(newDimX_y * newDimY_y * (NZ + 1), 0.0);
            grid_loc_z.resize(newDimX_z * newDimY_z * (NZ), 0.0);
            Y2_x.resize(newDimX_x * newDimY_x * (NZ + 1), 0.0);
            Y2_y.resize(newDimX_y * newDimY_y * (NZ + 1), 0.0);
            Y2_z.resize(newDimX_z * newDimY_z * (NZ), 0.0);
            Y3_x.resize(newDimX_x * newDimY_x * (NZ + 1), 0.0);
            Y3_y.resize(newDimX_y * newDimY_y * (NZ + 1), 0.0);
            Y3_z.resize(newDimX_z * newDimY_z * (NZ), 0.0);
        }
 */
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
    boundary.setOtherDim(other_dim_x_x,other_dim_y_x,other_dim_x_y,other_dim_y_y,other_dim_x_z,other_dim_y_z);

    MPI_Type_vector(dim_x_x, dim_z, (newDimY_x)*dim_z, MPI_DOUBLE, &MPI_face_x_x);
    MPI_Type_commit(&MPI_face_x_x);

    MPI_Type_vector(1, dim_z * newDimY_x, 0, MPI_DOUBLE, &MPI_face_y_x);
    MPI_Type_commit(&MPI_face_y_x);

    MPI_Type_vector(dim_x_y, dim_z, (newDimY_y)*dim_z, MPI_DOUBLE, &MPI_face_x_y);
    MPI_Type_commit(&MPI_face_x_y);

    MPI_Type_vector(1, dim_z * newDimY_y, 0, MPI_DOUBLE, &MPI_face_y_y);
    MPI_Type_commit(&MPI_face_y_y);

    MPI_Type_vector(dim_x_z, dim_z_z, (newDimY_z)*dim_z_z, MPI_DOUBLE, &MPI_face_x_z);
    MPI_Type_commit(&MPI_face_x_z);

    MPI_Type_vector(1, dim_z_z * newDimY_z, 0, MPI_DOUBLE, &MPI_face_y_z);
    MPI_Type_commit(&MPI_face_y_z);
}

void IcoNS::exchangeData(std::vector<Real> &grid_loc, int newDimX, int newDimY, int dim_z, MPI_Datatype MPI_face_x, MPI_Datatype MPI_face_y)
{
    if (!(/*dirichlet on y && */ coords[1] == 0))
    {
        MPI_Isend(&grid_loc[dim_z * newDimY + (newDimY - 2) * dim_z], 1, MPI_face_x, neighbors[1], rank, cart_comm, &reqs[3]);
        MPI_Irecv(&grid_loc[dim_z * newDimY + (newDimY - 1) * dim_z], 1, MPI_face_x, neighbors[1], neighbors[1], cart_comm, &reqs[3]);
        MPI_Wait(&reqs[3], MPI_SUCCESS);
    }
    if (!(/*dirichlet on y && */ coords[1] == PY - 1))
    {
        MPI_Isend(&grid_loc[dim_z * newDimY + dim_z], 1, MPI_face_x, neighbors[3], rank, cart_comm, &reqs[2]);
        MPI_Irecv(&grid_loc[dim_z * newDimY], 1, MPI_face_x, neighbors[3], neighbors[3], cart_comm, &reqs[2]);
        MPI_Wait(&reqs[2], MPI_SUCCESS);
    }
    // inviare al neighbours 0 dal processore con coords 0 vuol dire inviare al processo dall'altra parte >> inutile se dirichlet
    if (!(/*dirichlet on x axis &&*/ coords[0] == 0))
    {
        MPI_Isend(&grid_loc[(newDimY)*dim_z], 1, MPI_face_y, neighbors[0], rank, cart_comm, &reqs[0]);
        MPI_Irecv(&grid_loc[0], 1, MPI_face_y, neighbors[0], neighbors[0], cart_comm, &reqs[0]);
        MPI_Wait(&reqs[0], MPI_SUCCESS);
    }
    if (!(/*dirichlet on x &&*/ coords[0] == PX - 1))
    {
        MPI_Irecv(&grid_loc[(dim_z)*newDimY * (newDimX - 1)], 1, MPI_face_y, neighbors[2], neighbors[2], cart_comm, &reqs[1]);
        MPI_Isend(&grid_loc[newDimY * dim_z * (newDimX - 2)], 1, MPI_face_y, neighbors[2], rank, cart_comm, &reqs[1]);
        MPI_Wait(&reqs[1], MPI_SUCCESS);
    }
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
    double x=0;
    while (time < T)
    {
        boundary.update_boundary(grid_loc_x, grid_loc_y, grid_loc_z, time);

        MPI_Barrier(cart_comm);
        exchangeData(grid_loc_x, newDimX_x, newDimY_x,dim_z,MPI_face_x_x,MPI_face_y_x);
        exchangeData(grid_loc_y, newDimX_y, newDimY_y,dim_z,MPI_face_x_y,MPI_face_y_y);
        exchangeData(grid_loc_z, newDimX_z, newDimY_z,dim_z_z,MPI_face_x_z,MPI_face_y_z);
        x = L2_error(time); // every processor calculates his error not counting ghosts
        MPI_Barrier(cart_comm);
        error = 0.0;
        MPI_Reduce(&x, &error, 1, MPI_DOUBLE, MPI_SUM, 0, cart_comm);

        if (rank == 0)
        {
            error=sqrt(error);
            std::cout << " error: " << error << std::endl;
        }
        // {
        //     for(int r =0; r<4 ;r++){
        //     if(rank==r){
        //         std::cout << "rsnk: " << rank << std::endl;
        //         for (int in = 0 ; in < newDimX_x ; in++)
        //         {
        //             for(int j=0; j < newDimY_x ;j++){
        //                 for(int k=0; k <dim_z;k++){
        //                     std::cout << grid_loc[100%] Built target main
        //                 }
        //                 std::cout<< std::endl;
        //             }
        //             std::cout<< std::endl;
        //         }

        //     }
        //     MPI_Barrier(cart_comm);
        //     }
        // }
        // reduce
        solve_time_step(time); // adapt cycles to skip ghosts
        MPI_Barrier(cart_comm);
        time += DT;
        i++;
    }

    // write to output file
    output_x();
    //output_y();
    //output_z();

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
    error += error_comp_Y(t);
    error += error_comp_Z(t);

    //std::cout << error_comp_X(t) << std::endl;
    // std::cout << error_comp_Y(t) << std::endl;
    // std::cout << error_comp_Z(t) << std::endl << std::endl;

    return error;//sqrt(error);
}

// TODO:
Real IcoNS::error_comp_X(const Real t)
{
    Real error = 0.0;
    int offset_x = coords[0] * other_dim_x_x -1;
    int offset_y = (PY - 1 - coords[1]) * other_dim_y_x -1;

    //first slice (left face)
    if (lbx)
    {
        {
            if (lby)
            {
                error += ((grid_loc_x[0 + (newDimY_x + 1) * dim_z] - exact_solution.value_x(0.5, 0, 0, t)) *
                          (grid_loc_x[0 + (newDimY_x + 1) * dim_z] - exact_solution.value_x(0.5, 0, 0, t)) *
                          DX * DY * DZ / 8);

                for (int k = 1; k < NZ; k++)
                {
                    error += ((grid_loc_x[(newDimY_x + 1) * dim_z + k] - exact_solution.value_x(0.5, 0, k, t)) *
                              (grid_loc_x[(newDimY_x + 1) * dim_z + k] - exact_solution.value_x(0.5, 0, k, t)) *
                              DX * DY * DZ / 4);
                }

                error += ((grid_loc_x[(newDimY_x + 1) * dim_z + NZ] - exact_solution.value_x(0.5, 0, NZ, t)) *
                            (grid_loc_x[(newDimY_x + 1) * dim_z + NZ] - exact_solution.value_x(0.5, 0, NZ, t)) *
                            DX * DY * DZ / 8);
            }
            for (int j = 1 + lby; j < newDimY_x - 1 - rby; j++)
            {
                error += ((grid_loc_x[(newDimY_x)*dim_z + j * (NZ + 1)] - exact_solution.value_x(0.5, j + offset_y, 0, t)) *
                          (grid_loc_x[(newDimY_x)*dim_z + j * (NZ + 1)] - exact_solution.value_x(0.5, j + offset_y, 0, t)) *
                          DX * DY * DZ / 4);
                for (int k = 1; k < NZ; k++)
                {
                    error += ((grid_loc_x[(newDimY_x)*dim_z + j * (NZ + 1) + k] - exact_solution.value_x(0.5, j + offset_y, k, t)) *
                              (grid_loc_x[(newDimY_x)*dim_z + j * (NZ + 1) + k] - exact_solution.value_x(0.5, j + offset_y, k, t)) *
                              DX * DY * DZ / 2);
                }
                error += ((grid_loc_x[(newDimY_x)*dim_z + j * (NZ + 1) + NZ] - exact_solution.value_x(0.5, j + offset_y, NZ, t)) *
                          (grid_loc_x[(newDimY_x)*dim_z + j * (NZ + 1) + NZ] - exact_solution.value_x(0.5, j + offset_y, NZ, t)) *
                          DX * DY * DZ / 4);
            }
            if (rby)
            {
                error += ((grid_loc_x[(newDimY_x)*dim_z + (newDimY_x - 2) * (NZ + 1)] - exact_solution.value_x(0.5, NY, 0, t)) *
                          (grid_loc_x[(newDimY_x)*dim_z + (newDimY_x - 2) * (NZ + 1)] - exact_solution.value_x(0.5, NY, 0, t)) *
                          DX * DY * DZ / 8);

                for (int k = 1; k < NZ; k++)
                {
                    error += ((grid_loc_x[(newDimY_x)*dim_z + (newDimY_x - 2) * (NZ + 1) + k] - exact_solution.value_x(0.5, NY, k, t)) *
                              (grid_loc_x[(newDimY_x)*dim_z + (newDimY_x - 2) * (NZ + 1) + k] - exact_solution.value_x(0.5, NY, k, t)) *
                              DX * DY * DZ / 4);
                }

                error += ((grid_loc_x[(newDimY_x)*dim_z + (newDimY_x - 2) * (NZ + 1) + NZ] - exact_solution.value_x(0.5, NY, NZ, t)) *
                          (grid_loc_x[(newDimY_x)*dim_z + (newDimY_x - 2) * (NZ + 1) + NZ] - exact_solution.value_x(0.5, NY, NZ, t)) *
                          DX * DY * DZ / 8);
            }
        }
    }

    // middle slices
    {
        for (int i = 1 + lbx; i < newDimX_x - 1 - rbx; i++)
        {
            if (lby)
            {
                error += ((grid_loc_x[i * newDimY_x * (NZ + 1) + dim_z] - exact_solution.value_x(i + 0.5 + offset_x, 0, 0, t)) *
                          (grid_loc_x[i * newDimY_x * (NZ + 1) + dim_z] - exact_solution.value_x(i + 0.5 + offset_x, 0, 0, t)) *
                          DX * DY * DZ / 4);

                for (int k = 1; k < NZ; k++)
                {
                    error += ((grid_loc_x[i * newDimY_x * (NZ + 1) + dim_z + k] - exact_solution.value_x(i + 0.5 + offset_x, 0, k, t)) *
                              (grid_loc_x[i * newDimY_x * (NZ + 1) + dim_z + k] - exact_solution.value_x(i + 0.5 + offset_x, 0, k, t)) *
                              DX * DY * DZ / 2);
                }

                error += ((grid_loc_x[i * newDimY_x * (NZ + 1) + dim_z + NZ] - exact_solution.value_x(i + 0.5 + offset_x, 0, NZ, t)) *
                          (grid_loc_x[i * newDimY_x * (NZ + 1) + dim_z + NZ] - exact_solution.value_x(i + 0.5 + offset_x, 0, NZ, t)) *
                          DX * DY * DZ / 4);
            }
            for (int j = 1 + lby; j < newDimY_x - 1 - rby; j++)
            {
                error += ((grid_loc_x[i * newDimY_x * (NZ + 1) + j * (NZ + 1)] - exact_solution.value_x(i + 0.5 + offset_x, j + offset_y, 0, t)) *
                          (grid_loc_x[i * newDimY_x * (NZ + 1) + j * (NZ + 1)] - exact_solution.value_x(i + 0.5 + offset_x, j + offset_y, 0, t)) *
                          DX * DY * DZ / 2);

                for (int k = 1; k < NZ; k++)
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
                error += ((grid_loc_x[i * newDimY_x * (NZ + 1) + (newDimY_x - 2) * (NZ + 1)] - exact_solution.value_x(i + 0.5 + offset_x, NY, 0, t)) *
                          (grid_loc_x[i * newDimY_x * (NZ + 1) + (newDimY_x - 2) * (NZ + 1)] - exact_solution.value_x(i + 0.5 + offset_x, NY, 0, t)) *
                          DX * DY * DZ / 4);

                for (int k = 1; k < NZ; k++)
                {
                    error += ((grid_loc_x[i * newDimY_x * (NZ + 1) + (newDimY_x - 2) * (NZ + 1) + k] - exact_solution.value_x(i + 0.5 + offset_x, NY, k, t)) *
                              (grid_loc_x[i * newDimY_x * (NZ + 1) + (newDimY_x - 2) * (NZ + 1) + k] - exact_solution.value_x(i + 0.5 + offset_x, NY, k, t)) *
                              DX * DY * DZ / 2);
                }

                error += ((grid_loc_x[i * newDimY_x * (NZ + 1) + (newDimY_x - 2) * (NZ + 1) + NZ] - exact_solution.value_x(i + 0.5 + offset_x, NY, NZ, t)) *
                          (grid_loc_x[i * newDimY_x * (NZ + 1) + (newDimY_x - 2) * (NZ + 1) + NZ] - exact_solution.value_x(i + 0.5 + offset_x, NY, NZ, t)) *
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
                error += ((grid_loc_x[(newDimX_x - 2) * newDimY_x * (NZ + 1) + dim_z] - exact_solution.value_x((NX - 0.5), 0, 0, t)) *
                          (grid_loc_x[(newDimX_x - 2) * newDimY_x * (NZ + 1) + dim_z] - exact_solution.value_x((NX - 0.5), 0, 0, t)) *
                          DX * DY * DZ / 8);

                for (int k = 1; k < NZ; k++)
                {
                    error += ((grid_loc_x[(newDimX_x - 2) * newDimY_x * (NZ + 1) + dim_z + k] - exact_solution.value_x((NX - 0.5), 0, k, t)) *
                              (grid_loc_x[(newDimX_x - 2) * newDimY_x * (NZ + 1) + dim_z + k] - exact_solution.value_x((NX - 0.5), 0, k, t)) *
                              DX * DY * DZ / 4);
                }

                error += ((grid_loc_x[(newDimX_x - 2) * newDimY_x * (NZ + 1) + dim_z + NZ] - exact_solution.value_x((NX - 0.5), 0, NZ, t)) *
                          (grid_loc_x[(newDimX_x - 2) * newDimY_x * (NZ + 1) + dim_z + NZ] - exact_solution.value_x((NX - 0.5), 0, NZ, t)) *
                          DX * DY * DZ / 8);
            }
            for (int j = 1 + lby; j < newDimY_x - 1 - rby; j++)
            {
                error += ((grid_loc_x[(newDimX_x - 2) * newDimY_x * (NZ + 1) + j * (NZ + 1)] - exact_solution.value_x((NX - 0.5), j + offset_y, 0, t)) *
                          (grid_loc_x[(newDimX_x - 2) * newDimY_x * (NZ + 1) + j * (NZ + 1)] - exact_solution.value_x((NX - 0.5), j + offset_y, 0, t)) *
                          DX * DY * DZ / 4);
                for (int k = 1; k < NZ; k++)
                {
                    error += ((grid_loc_x[(newDimX_x - 2) * newDimY_x * (NZ + 1) + j * (NZ + 1) + k] - exact_solution.value_x((NX - 0.5), j + offset_y, k, t)) *
                              (grid_loc_x[(newDimX_x - 2) * newDimY_x * (NZ + 1) + j * (NZ + 1) + k] - exact_solution.value_x((NX - 0.5), j + offset_y, k, t)) *
                              DX * DY * DZ / 2);
                }
                error += ((grid_loc_x[(newDimX_x - 2) * newDimY_x * (NZ + 1) + j * (NZ + 1) + NZ] - exact_solution.value_x((NX - 0.5), j + offset_y, NZ, t)) *
                          (grid_loc_x[(newDimX_x - 2) * newDimY_x * (NZ + 1) + j * (NZ + 1) + NZ] - exact_solution.value_x((NX - 0.5), j + offset_y, NZ, t)) *
                          DX * DY * DZ / 4);
            }
            if (rby)
            {
                error += ((grid_loc_x[(newDimX_x - 2) * newDimY_x * (NZ + 1) + (newDimY_x - 2) * (NZ + 1)] - exact_solution.value_x((NX - 0.5), NY, 0, t)) *
                          (grid_loc_x[(newDimX_x - 2) * newDimY_x * (NZ + 1) + (newDimY_x - 2) * (NZ + 1)] - exact_solution.value_x((NX - 0.5), NY, 0, t)) *
                          DX * DY * DZ / 8);

                for (int k = 1; k < NZ; k++)
                {
                    error += ((grid_loc_x[(newDimX_x - 2) * newDimY_x * (NZ + 1) + (newDimY_x - 2) * (NZ + 1) + k] - exact_solution.value_x((NX - 0.5), NY, k, t)) *
                              (grid_loc_x[(newDimX_x - 2) * newDimY_x * (NZ + 1) + (newDimY_x - 2) * (NZ + 1) + k] - exact_solution.value_x((NX - 0.5), NY, k, t)) *
                              DX * DY * DZ / 4);
                }

                error += ((grid_loc_x[(newDimX_x - 2) * newDimY_x * (NZ + 1) + (newDimY_x - 2) * (NZ + 1) + NZ] - exact_solution.value_x((NX - 0.5), NY, NZ, t)) *
                          (grid_loc_x[(newDimX_x - 2) * newDimY_x * (NZ + 1) + (newDimY_x - 2) * (NZ + 1) + NZ] - exact_solution.value_x((NX - 0.5), NY, NZ, t)) *
                          DX * DY * DZ / 8);
            }
        }
    }

    return error;
}

Real IcoNS::error_comp_Y(const Real t)
{
    Real error = 0.0;
    int offset_x = coords[0] * other_dim_x_y -1;
    int offset_y = (PY - 1 - coords[1]) * other_dim_y_y -1;
    // first slice (left face)
    if(lbx)
    {
        if(lby){
            error += ((grid_loc_y[(newDimY_y+1)*dim_z] - exact_solution.value_y(0, 0.5, 0, t)) *
                    (grid_loc_y[(newDimY_y+1)*dim_z] - exact_solution.value_y(0, 0.5, 0, t)) *
                    DX * DY * DZ / 8);

            for (int k = 1; k < NZ; k++)
            {
                error += ((grid_loc_y[(newDimY_y+1)*dim_z + k] - exact_solution.value_y(0, 0.5, k, t)) *
                        (grid_loc_y[(newDimY_y+1)*dim_z + k] - exact_solution.value_y(0, 0.5, k, t)) *
                        DX * DY * DZ / 4);
            }

            error += ((grid_loc_y[(newDimY_y+1)*dim_z + NZ] - exact_solution.value_y(0, 0.5, NZ, t)) *
                    (grid_loc_y[(newDimY_y+1)*dim_z+ NZ] - exact_solution.value_y(0, 0.5, NZ, t)) *
                    DX * DY * DZ / 8);
        }
        for (int j = 1 + lby; j < newDimY_y - 1 -rby; j++)
        {
            error += ((grid_loc_y[(newDimY_y)*dim_z + j * (NZ + 1)] - exact_solution.value_y(0, j + 0.5 + offset_y, 0, t)) *
                      (grid_loc_y[(newDimY_y)*dim_z + j * (NZ + 1)] - exact_solution.value_y(0, j + 0.5 + offset_y, 0, t)) *
                      DX * DY * DZ / 4);
            for (int k = 1; k < NZ; k++)
            {
                error += ((grid_loc_y[(newDimY_y)*dim_z + j * (NZ + 1) + k] - exact_solution.value_y(0, j + 0.5 + offset_y, k, t)) *
                          (grid_loc_y[(newDimY_y)*dim_z + j * (NZ + 1) + k] - exact_solution.value_y(0, j + 0.5 + offset_y, k, t)) *
                          DX * DY * DZ / 2);
            }
            error += ((grid_loc_y[(newDimY_y)*dim_z + j * (NZ + 1) + NZ] - exact_solution.value_y(0, j + 0.5 + offset_y, NZ, t)) *
                      (grid_loc_y[(newDimY_y)*dim_z + j * (NZ + 1) + NZ] - exact_solution.value_y(0, j + 0.5 + offset_y, NZ, t)) *
                      DX * DY * DZ / 4);
        }
        if(rby){
            error += ((grid_loc_y[(newDimY_y)*dim_z + (newDimY_y-2) * (NZ + 1)] - exact_solution.value_y(0, (NY - 0.5), 0, t)) *
                    (grid_loc_y[(newDimY_y)*dim_z + (newDimY_y-2) * (NZ + 1)] - exact_solution.value_y(0, (NY - 0.5), 0, t)) *
                    DX * DY * DZ / 8);

            for (int k = 1; k < NZ; k++)
            {
                error += ((grid_loc_y[(newDimY_y)*dim_z + (newDimY_y-2) * (NZ + 1) + k] - exact_solution.value_y(0, (NY - 0.5), k, t)) *
                        (grid_loc_y[(newDimY_y)*dim_z + (newDimY_y-2) * (NZ + 1) + k] - exact_solution.value_y(0, (NY - 0.5), k, t)) *
                        DX * DY * DZ / 4);
            }

            error += ((grid_loc_y[(newDimY_y)*dim_z + (newDimY_y-2) * (NZ + 1) + NZ] - exact_solution.value_y(0, (NY - 0.5), NZ, t)) *
                    (grid_loc_y[(newDimY_y)*dim_z + (newDimY_y-2) * (NZ + 1) + NZ] - exact_solution.value_y(0, (NY - 0.5), NZ, t)) *
                    DX * DY * DZ / 8);
        }
    }

    // middle slices
    {
        for (int i = 1 + lbx; i < newDimX_y -1 - rbx; i++)
        {
            if(lby){
                error += ((grid_loc_y[i * newDimY_y * (NZ + 1) + dim_z] - exact_solution.value_y(i + offset_x, 0.5, 0, t)) *
                        (grid_loc_y[i * newDimY_y * (NZ + 1) + dim_z] - exact_solution.value_y(i + offset_x, 0.5, 0, t)) *
                        DX * DY * DZ / 4);

                for (int k = 1; k < NZ; k++)
                {
                    error += ((grid_loc_y[i * newDimY_y * (NZ + 1) + dim_z + k] - exact_solution.value_y(i + offset_x, 0.5, k, t)) *
                            (grid_loc_y[i * newDimY_y * (NZ + 1) + dim_z + k] - exact_solution.value_y(i + offset_x, 0.5, k, t)) *
                            DX * DY * DZ / 2);
                }
                error += ((grid_loc_y[i * newDimY_y * (NZ + 1) + dim_z + NZ] - exact_solution.value_y(i + offset_x, 0.5, NZ, t)) *
                        (grid_loc_y[i * newDimY_y * (NZ + 1) + dim_z + NZ] - exact_solution.value_y(i + offset_x, 0.5, NZ, t)) *
                        DX * DY * DZ / 4);
            }
            for (int j = 1 + lby; j < newDimY_y - 1 - rby; j++)
            {
                error += ((grid_loc_y[i * newDimY_y * (NZ + 1) + j * (NZ + 1)] - exact_solution.value_y(i + offset_x, j + 0.5 + offset_y, 0, t)) *
                          (grid_loc_y[i * newDimY_y * (NZ + 1) + j * (NZ + 1)] - exact_solution.value_y(i + offset_x, j + 0.5 + offset_y, 0, t)) *
                          DX * DY * DZ / 2);

                for (int k = 1; k < NZ; k++)
                {
                    error += ((grid_loc_y[i * newDimY_y * (NZ + 1) + j * (NZ + 1) + k] - exact_solution.value_y(i + offset_x, j + 0.5 + offset_y, k, t)) *
                              (grid_loc_y[i * newDimY_y * (NZ + 1) + j * (NZ + 1) + k] - exact_solution.value_y(i + offset_x, j + 0.5 + offset_y, k, t)) *
                              DX * DY * DZ);
                }

                error += ((grid_loc_y[i * newDimY_y * (NZ + 1) + j * (NZ + 1) + NZ] - exact_solution.value_y(i + offset_x, j + 0.5 + offset_y, NZ, t)) *
                          (grid_loc_y[i * newDimY_y * (NZ + 1) + j * (NZ + 1) + NZ] - exact_solution.value_y(i + offset_x, j + 0.5 + offset_y, NZ, t)) *
                          DX * DY * DZ / 2);
            }
            if(rby){
                error += ((grid_loc_y[i * newDimY_y * (NZ + 1) + (newDimY_y-2) * (NZ + 1)] - exact_solution.value_y(i + offset_x, (NY - 0.5), 0, t)) *
                        (grid_loc_y[i * newDimY_y * (NZ + 1) + (newDimY_y-2) * (NZ + 1)] - exact_solution.value_y(i + offset_x, (NY - 0.5), 0, t)) *
                        DX * DY * DZ / 4);

                for (int k = 1; k < NZ; k++)
                {
                    error += ((grid_loc_y[i * newDimY_y * (NZ + 1) + (newDimY_y-2) * (NZ + 1) + k] - exact_solution.value_y(i + offset_x, (NY - 0.5), k, t)) *
                            (grid_loc_y[i * newDimY_y * (NZ + 1) + (newDimY_y-2) * (NZ + 1) + k] - exact_solution.value_y(i + offset_x, (NY - 0.5), k, t)) *
                            DX * DY * DZ / 2);
                }

                error += ((grid_loc_y[i * newDimY_y * (NZ + 1) + (newDimY_y-2) * (NZ + 1) + NZ] - exact_solution.value_y(i + offset_x, (NY - 0.5), NZ, t)) *
                        (grid_loc_y[i * newDimY_y * (NZ + 1) + (newDimY_y-2) * (NZ + 1) + NZ] - exact_solution.value_y(i + offset_x, (NY - 0.5), NZ, t)) *
                        DX * DY * DZ / 4);
            }
        }
    }

    // last slice (right face)
    if(rbx)
    {
        if(lby){
            error += ((grid_loc_y[(newDimX_y-2) * newDimY_y * (NZ + 1) + dim_z] - exact_solution.value_y(NX, 0.5, 0, t)) *
                    (grid_loc_y[(newDimX_y-2) * newDimY_y * (NZ + 1) + dim_z] - exact_solution.value_y(NX, 0.5, 0, t)) *
                    DX * DY * DZ / 8);

            for (int k = 1; k < NZ; k++)
            {
                error += ((grid_loc_y[(newDimX_y-2) * newDimY_y * (NZ + 1) + dim_z + k] - exact_solution.value_y(NX, 0.5, k, t)) *
                        (grid_loc_y[(newDimX_y-2) * newDimY_y * (NZ + 1) + dim_z + k] - exact_solution.value_y(NX, 0.5, k, t)) *
                        DX * DY * DZ / 4);
            }

            error += ((grid_loc_y[(newDimX_y-2) * newDimY_y * (NZ + 1) + dim_z + NZ] - exact_solution.value_y(NX, 0.5, NZ, t)) *
                    (grid_loc_y[(newDimX_y-2) * newDimY_y * (NZ + 1) + dim_z + NZ] - exact_solution.value_y(NX, 0.5, NZ, t)) *
                    DX * DY * DZ / 8);
        }
        for (int j = 1 + lby; j < newDimY_y - 1 - rby; j++)
        {
            error += ((grid_loc_y[(newDimX_y-2) * newDimY_y * (NZ + 1) + j * (NZ + 1)] - exact_solution.value_y(NX, j + 0.5 + offset_y, 0, t)) *
                      (grid_loc_y[(newDimX_y-2) * newDimY_y * (NZ + 1) + j * (NZ + 1)] - exact_solution.value_y(NX, j + 0.5 + offset_y, 0, t)) *
                      DX * DY * DZ / 4);
            for (int k = 1; k < NZ; k++)
            {
                error += ((grid_loc_y[(newDimX_y-2) * newDimY_y * (NZ + 1) + j * (NZ + 1) + k] - exact_solution.value_y(NX, j + 0.5 + offset_y, k, t)) *
                          (grid_loc_y[(newDimX_y-2) * newDimY_y * (NZ + 1) + j * (NZ + 1) + k] - exact_solution.value_y(NX, j + 0.5 + offset_y, k, t)) *
                          DX * DY * DZ / 2);
            }
            error += ((grid_loc_y[(newDimX_y-2) * newDimY_y * (NZ + 1) + j * (NZ + 1) + NZ] - exact_solution.value_y(NX, j + 0.5 + offset_y, NZ, t)) *
                      (grid_loc_y[(newDimX_y-2) * newDimY_y * (NZ + 1) + j * (NZ + 1) + NZ] - exact_solution.value_y(NX, j + 0.5 + offset_y, NZ, t)) *
                      DX * DY * DZ / 4);
        }
        if(rby){
            error += ((grid_loc_y[(newDimX_y-2) * newDimY_y * (NZ + 1) + (newDimY_y-2) * (NZ + 1)] - exact_solution.value_y(NX, (NY - 0.5), 0, t)) *
                    (grid_loc_y[(newDimX_y-2) * newDimY_y * (NZ + 1) + (newDimY_y-2) * (NZ + 1)] - exact_solution.value_y(NX, (NY - 0.5), 0, t)) *
                    DX * DY * DZ / 8);

            for (int k = 1; k < NZ; k++)
            {
                error += ((grid_loc_y[(newDimX_y-2) * newDimY_y * (NZ + 1) + (newDimY_y-2) * (NZ + 1) + k] - exact_solution.value_y(NX, (NY - 0.5), k, t)) *
                        (grid_loc_y[(newDimX_y-2) * newDimY_y * (NZ + 1) + (newDimY_y-2) * (NZ + 1) + k] - exact_solution.value_y(NX, (NY - 0.5), k, t)) *
                        DX * DY * DZ / 4);
            }

            error += ((grid_loc_y[(newDimX_y-2) * newDimY_y * (NZ + 1) + (newDimY_y-2) * (NZ + 1) + NZ] - exact_solution.value_y(NX, (NY - 0.5), NZ, t)) *
                    (grid_loc_y[(newDimX_y-2) * newDimY_y * (NZ + 1) + (newDimY_y-2) * (NZ + 1) + NZ] - exact_solution.value_y(NX, (NY - 0.5), NZ, t)) *
                    DX * DY * DZ / 8);
        }
    }

    return error;
}

Real IcoNS::error_comp_Z(const Real t)
{
    Real error = 0.0;
    int offset_x = coords[0] * other_dim_x_z -1;
    int offset_y = (PY - 1 - coords[1]) * other_dim_y_z -1;
    // first slice (left face)
    if(lbx)
    {
        if(lby){
            error += ((grid_loc_z[(newDimY_z +1) *dim_z_z] - exact_solution.value_z(0, 0, 0.5, t)) *
                    (grid_loc_z[(newDimY_z +1) *dim_z_z] - exact_solution.value_z(0, 0, 0.5, t)) *
                    DX * DY * DZ / 8);

            for (int k = 1; k < NZ - 1; k++)
            {
                error += ((grid_loc_z[(newDimY_z +1) *dim_z_z + k] - exact_solution.value_z(0, 0, k + 0.5, t)) *
                        (grid_loc_z[(newDimY_z +1) *dim_z_z + k] - exact_solution.value_z(0, 0, k + 0.5, t)) *
                        DX * DY * DZ / 4);
            }

            error += ((grid_loc_z[(newDimY_z +1) *dim_z_z + NZ - 1] - exact_solution.value_z(0, 0, NZ - 0.5, t)) *
                    (grid_loc_z[(newDimY_z +1) *dim_z_z + NZ - 1] - exact_solution.value_z(0, 0, NZ - 0.5, t)) *
                    DX * DY * DZ / 8);
        }
        for (int j = 1 + lby; j < newDimY_z-1-rby; j++)
        {
            error += ((grid_loc_z[(newDimY_z) *dim_z_z + j * NZ] - exact_solution.value_z(0, j + offset_y, 0.5, t)) *
                      (grid_loc_z[(newDimY_z) *dim_z_z + j * NZ] - exact_solution.value_z(0, j + offset_y, 0.5, t)) *
                      DX * DY * DZ / 4);
            for (int k = 1; k < NZ - 1; k++)
            {
                error += ((grid_loc_z[(newDimY_z) *dim_z_z + j * NZ + k] - exact_solution.value_z(0, j + offset_y, k + 0.5, t)) *
                          (grid_loc_z[(newDimY_z) *dim_z_z + j * NZ + k] - exact_solution.value_z(0, j + offset_y, k + 0.5, t)) *
                          DX * DY * DZ / 2);
            }
            error += ((grid_loc_z[(newDimY_z) *dim_z_z + j * NZ + NZ - 1] - exact_solution.value_z(0, j + offset_y, NZ - 0.5, t)) *
                      (grid_loc_z[(newDimY_z) *dim_z_z + j * NZ + NZ - 1] - exact_solution.value_z(0, j + offset_y, NZ - 0.5, t)) *
                      DX * DY * DZ / 4);
        }
        if(rby){
            error += ((grid_loc_z[(newDimY_z) *dim_z_z + (newDimY_z-2)*dim_z_z] - exact_solution.value_z(0, NY, 0.5, t)) *
                    (grid_loc_z[(newDimY_z) *dim_z_z + (newDimY_z-2)*dim_z_z] - exact_solution.value_z(0, NY, 0.5, t)) *
                    DX * DY * DZ / 8);

            for (int k = 1; k < NZ - 1; k++)
            {
                error += ((grid_loc_z[(newDimY_z) *dim_z_z + (newDimY_z-2)*dim_z_z + k] - exact_solution.value_z(0, NY, k + 0.5, t)) *
                        (grid_loc_z[(newDimY_z) *dim_z_z + (newDimY_z-2)*dim_z_z + k] - exact_solution.value_z(0, NY, k + 0.5, t)) *
                        DX * DY * DZ / 4);
            }

            error += ((grid_loc_z[(newDimY_z) *dim_z_z + (newDimY_z-2)*dim_z_z + NZ - 1] - exact_solution.value_z(0, NY, NZ - 0.5, t)) *
                    (grid_loc_z[(newDimY_z) *dim_z_z + (newDimY_z-2)*dim_z_z + NZ - 1] - exact_solution.value_z(0, NY, NZ - 0.5, t)) *
                    DX * DY * DZ / 8);
        }
    }

    // middle slices
    {
        for (int i = 1 + lbx; i < newDimX_z - 1 -rbx; i++)
        {
            if(lby){
                error += ((grid_loc_z[i * newDimY_z * NZ  + dim_z_z] - exact_solution.value_z(i + offset_x, 0, 0.5, t)) *
                        (grid_loc_z[i * newDimY_z * NZ  + dim_z_z] - exact_solution.value_z(i + offset_x, 0, 0.5, t)) *
                        DX * DY * DZ / 4);

                for (int k = 1; k < NZ - 1; k++)
                {
                    error += ((grid_loc_z[i * newDimY_z * NZ  + dim_z_z + k] - exact_solution.value_z(i + offset_x, 0, k + 0.5, t)) *
                            (grid_loc_z[i * newDimY_z * NZ  + dim_z_z + k] - exact_solution.value_z(i + offset_x, 0, k + 0.5, t)) *
                            DX * DY * DZ / 2);
                }
                error += ((grid_loc_z[i * newDimY_z * NZ  + dim_z_z + NZ - 1] - exact_solution.value_z(i + offset_x, 0, NZ - 0.5, t)) *
                        (grid_loc_z[i * newDimY_z * NZ  + dim_z_z + NZ - 1] - exact_solution.value_z(i + offset_x, 0, NZ - 0.5, t)) *
                        DX * DY * DZ / 4);
            }
            for (int j = 1 + lby; j < newDimY_z - 1 -rby; j++)
            {
                error += ((grid_loc_z[i * newDimY_z * NZ + j * NZ] - exact_solution.value_z(i + offset_x, j + offset_y, 0.5, t)) *
                          (grid_loc_z[i * newDimY_z * NZ + j * NZ] - exact_solution.value_z(i + offset_x, j + offset_y, 0.5, t)) *
                          DX * DY * DZ / 2);

                for (int k = 1; k < NZ - 1; k++)
                {
                    error += ((grid_loc_z[i * newDimY_z * NZ + j * NZ + k] - exact_solution.value_z(i + offset_x, j + offset_y, k + 0.5, t)) *
                              (grid_loc_z[i * newDimY_z * NZ + j * NZ + k] - exact_solution.value_z(i + offset_x, j + offset_y, k + 0.5, t)) *
                              DX * DY * DZ);
                }

                error += ((grid_loc_z[i * newDimY_z * NZ + j * NZ + NZ - 1] - exact_solution.value_z(i + offset_x, j + offset_y, (NZ - 0.5), t)) *
                          (grid_loc_z[i * newDimY_z * NZ + j * NZ + NZ - 1] - exact_solution.value_z(i + offset_x, j + offset_y, (NZ - 0.5), t)) *
                          DX * DY * DZ / 2);
            }
            if(rby){
                error += ((grid_loc_z[i * newDimY_z * NZ + (newDimY_z-2) * NZ] - exact_solution.value_z(i + offset_x, NY, 0.5, t)) *
                        (grid_loc_z[i * newDimY_z * NZ + (newDimY_z-2) * NZ] - exact_solution.value_z(i + offset_x, NY, 0.5, t)) *
                        DX * DY * DZ / 4);

                for (int k = 1; k < NZ - 1; k++)
                {
                    error += ((grid_loc_z[i * newDimY_z * NZ + (newDimY_z-2) * NZ + k] - exact_solution.value_z(i + offset_x, NY, k + 0.5, t)) *
                            (grid_loc_z[i * newDimY_z * NZ + (newDimY_z-2) * NZ + k] - exact_solution.value_z(i + offset_x, NY, k + 0.5, t)) *
                            DX * DY * DZ / 2);
                }

                error += ((grid_loc_z[i * newDimY_z * NZ + (newDimY_z-2) * NZ + NZ - 1] - exact_solution.value_z(i + offset_x, NY, NZ - 0.5, t)) *
                        (grid_loc_z[i * newDimY_z * NZ + (newDimY_z-2) * NZ + NZ - 1] - exact_solution.value_z(i + offset_x, NY, NZ - 0.5, t)) *
                        DX * DY * DZ / 4);
            }
        }
    }

    // last slice (right face)
    if(rbx)
    {
        if(lby){
            error += ((grid_loc_z[(newDimX_z-2) * (newDimY_z) * NZ + dim_z_z] - exact_solution.value_z(NX, 0, 0.5, t)) *
                    (grid_loc_z[(newDimX_z-2) * (newDimY_z) * NZ + dim_z_z] - exact_solution.value_z(NX, 0, 0.5, t)) *
                    DX * DY * DZ / 8);

            for (int k = 1; k < NZ - 1; k++)
            {
                error += ((grid_loc_z[(newDimX_z-2) * (newDimY_z) * NZ + dim_z_z + k] - exact_solution.value_z(NX, 0, k + 0.5, t)) *
                        (grid_loc_z[(newDimX_z-2) * (newDimY_z) * NZ + dim_z_z + k] - exact_solution.value_z(NX, 0, k + 0.5, t)) *
                        DX * DY * DZ / 4);
            }

            error += ((grid_loc_z[(newDimX_z-2) * (newDimY_z) * NZ + dim_z_z + NZ - 1] - exact_solution.value_z(NX, 0, NZ - 0.5, t)) *
                    (grid_loc_z[(newDimX_z-2) * (newDimY_z) * NZ + dim_z_z + NZ - 1] - exact_solution.value_z(NX, 0, NZ - 0.5, t)) *
                    DX * DY * DZ / 8);
        }
        for (int j = 1 + lby; j < newDimY_z - 1 - rby; j++)
        {
            error += ((grid_loc_z[(newDimX_z-2) * (newDimY_z) * NZ + j * NZ] - exact_solution.value_z(NX, j + offset_y, 0.5, t)) *
                      (grid_loc_z[(newDimX_z-2) * (newDimY_z) * NZ + j * NZ] - exact_solution.value_z(NX, j + offset_y, 0.5, t)) *
                      DX * DY * DZ / 4);
            for (int k = 1; k < NZ - 1; k++)
            {
                error += ((grid_loc_z[(newDimX_z-2) * (newDimY_z) * NZ + j * NZ + k] - exact_solution.value_z(NX, j + offset_y, k + 0.5, t)) *
                          (grid_loc_z[(newDimX_z-2) * (newDimY_z) * NZ + j * NZ + k] - exact_solution.value_z(NX, j + offset_y, k + 0.5, t)) *
                          DX * DY * DZ / 2);
            }
            error += ((grid_loc_z[(newDimX_z-2) * (newDimY_z) * NZ + j * NZ + NZ - 1] - exact_solution.value_z(NX, j + offset_y, NZ - 0.5, t)) *
                      (grid_loc_z[(newDimX_z-2) * (newDimY_z) * NZ + j * NZ + NZ - 1] - exact_solution.value_z(NX, j + offset_y, NZ - 0.5, t)) *
                      DX * DY * DZ / 4);
        }
        if(rby){
            error += ((grid_loc_z[(newDimX_z-2) * (newDimY_z) * NZ  + (newDimY_z-2) * NZ] - exact_solution.value_z(NX, NY, 0.5, t)) *
                    (grid_loc_z[(newDimX_z-2) * (newDimY_z) * NZ  + (newDimY_z-2) * NZ] - exact_solution.value_z(NX, NY, 0.5, t)) *
                    DX * DY * DZ / 8);

            for (int k = 1; k < NZ - 1; k++)
            {
                error += ((grid_loc_z[(newDimX_z-2) * (newDimY_z) * NZ  + (newDimY_z-2) * NZ + k] - exact_solution.value_z(NX, NY, k + 0.5, t)) *
                        (grid_loc_z[(newDimX_z-2) * (newDimY_z) * NZ  + (newDimY_z-2) * NZ + k] - exact_solution.value_z(NX, NY, k + 0.5, t)) *
                        DX * DY * DZ / 4);
            }

            error += ((grid_loc_z[(newDimX_z-2) * (newDimY_z) * NZ  + (newDimY_z-2) * NZ + NZ - 1] - exact_solution.value_z(NX, NY, NZ - 0.5, t)) *
                    (grid_loc_z[(newDimX_z-2) * (newDimY_z) * NZ  + (newDimY_z-2) * NZ + NZ - 1] - exact_solution.value_z(NX, NY, NZ - 0.5, t)) *
                    DX * DY * DZ / 8);
        }
    }

    return error;
}



void IcoNS::output_x() {
    MPI_File fh;
    MPI_Offset offset = 0;

    // Open the file
    MPI_File_open(cart_comm, "solution_x.vtk",
                  MPI_MODE_CREATE | MPI_MODE_WRONLY,
                  MPI_INFO_NULL, &fh);

    // Write the header (rank 0 only)
    if (rank == 0) {
        std::ostringstream full_header;
        full_header << std::fixed << std::setprecision(6);

        // Header section
        full_header << "# vtk DataFile Version 3.0\n"
                   << "Solution x\n"
                   << "ASCII\n"
                   << "DATASET STRUCTURED_GRID\n"
                   << "DIMENSIONS 1 " << NY + 1 << " " << NZ + 1 << "\n"
                   << "POINTS " << (NY + 1) * (NZ + 1) << " float\n";

        // Points section
        for (int k = 0; k < NZ + 1; k++) {
            for (int j = 0; j < NY + 1; j++) {
                full_header << "0.000000 "
                           << j * DY << " "
                           << k * DZ << "\n";
            }
        }

        // Data header section
        full_header << "\nPOINT_DATA " << (NY + 1) * (NZ + 1) << "\n"
                   << "SCALARS v float\n"
                   << "LOOKUP_TABLE default\n";

        std::string header_str = full_header.str();
        MPI_File_write_at(fh, 0, header_str.c_str(),
                         header_str.size(), MPI_CHAR, MPI_STATUS_IGNORE);
        offset = header_str.size();
    }

    // Broadcast offset to all processes
    MPI_Bcast(&offset, 1, MPI_OFFSET, 0, cart_comm);

    // Write data values
    for (int k = 0; k < newDimY_x; k++) {
        for (int j = 0; j < newDimY_y; j++) {
            int global_k = coords[1] * (newDimX_x - 1) + k;
            int global_j = coords[0] * (newDimX_y - 1) + j;

            if (global_k < NZ + 1 && global_j < NY + 1) {
                std::ostringstream value_str;
                value_str << std::fixed << std::setprecision(6)
                         << grid_loc_y[k * newDimY_y + j] << "\n";

                std::string str = value_str.str();
                MPI_Offset pos = offset +
                               (global_k * (NY + 1) + global_j) *
                               (str.size() * sizeof(char));

                MPI_File_write_at(fh, pos, str.c_str(), str.size(),
                                 MPI_CHAR, MPI_STATUS_IGNORE);
            }
        }
    }

    MPI_File_close(&fh);
}


void IcoNS::output_y(){
    MPI_File fh;
    MPI_Offset headerOffset, rowOffsetX, rowOffsetY, colOffsetX, colOffsetY;

    // Open the file
    MPI_File_open(cart_comm, "solution_y.vtk", MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);

    // Write the header
    std::ostringstream header;
    std::ostringstream scalarsu;
    std::ostringstream scalarsw;
    std::ostringstream scalarsp;
    header << "# vtk DataFile Version 3.0\n"
            << "Solution y\n"
            << "BINARY\n"
            << "DATASET STRUCTURED_GRID\n"
            << "DIMENSIONS " << NX + 1 << " " << 1 << " " << NZ + 1 << "\n"
            << "ORIGIN " << 0 << " " << 0 << " " << 0 << "\n"
            << "POINTS " << (NX + 1) * (NZ + 1) * 1 << " float\n";
    scalarsu << "SCALARS u float 1\nLOOKUP_TABLE default\n";
    scalarsw << "SCALARS w float 1\nLOOKUP_TABLE default\n";
    scalarsp << "SCALARS p float 1\nLOOKUP_TABLE default\n";

    headerOffset = header.str().size();

    // This offset determines the cut of the plane in the z direction
    int zOffset = 0;

    if(rank == 0){
        // Write header to file
        MPI_File_write_at(fh, 0, header.str().c_str(), header.str().size(), MPI_CHAR, MPI_STATUS_IGNORE);
        // Write scalars u header after metadata
        MPI_File_write_at(fh, headerOffset, scalarsu.str().c_str(), scalarsu.str().size(), MPI_CHAR, MPI_STATUS_IGNORE);
        // Write scalars v header after scalars u
        MPI_File_write_at(fh, headerOffset + scalarsu.str().size() + NY, scalarsw.str().c_str(), scalarsw.str().size(), MPI_CHAR, MPI_STATUS_IGNORE);
        // Write scalars p header after scalars v
        MPI_File_write_at(fh, headerOffset + scalarsu.str().size() + 2 * NY + scalarsw.str().size(), scalarsp.str().c_str(), scalarsp.str().size(), MPI_CHAR, MPI_STATUS_IGNORE);
    }

    // Write the data
    // Each processor must write all the rows of its local grid but in a different possition in the line
    // depending on the processor rank
    // go through all the rows of the processor
    colOffsetX = coords[0];
    for (size_t i = 0; i < newDimY_x; i++){
        // Write the whole line of the processor
        rowOffsetX = headerOffset + (rank * other_dim_y_x) * i;
        MPI_File_write_at_all(fh, rowOffsetX + colOffsetX, &grid_loc_x[zOffset], newDimX_x, MPI_DOUBLE, MPI_STATUS_IGNORE);
    }

    colOffsetY = coords[0];
    for (size_t i = 0; i < newDimY_y; i++){
        // Write the whole line of the processor
        rowOffsetY = headerOffset + NX + 1 + (rank * other_dim_y_y) * i;
        MPI_File_write_at_all(fh, rowOffsetY + colOffsetY, &grid_loc_y[zOffset], newDimX_y, MPI_DOUBLE, MPI_STATUS_IGNORE);
    }

    // TODO: Write preassure

}

void IcoNS::output_z()
{
    MPI_File fh;
    MPI_Offset headerOffset, rowOffsetX, rowOffsetY, colOffsetX, colOffsetY;

    // Open the file
    MPI_File_open(cart_comm, "solution_z.vtk", MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);

    // Write the header
    std::ostringstream header;
    std::ostringstream scalarsu;
    std::ostringstream scalarsv;
    std::ostringstream scalarsp;
    header << "# vtk DataFile Version 3.0\n"
            << "Solution\n"
            << "BINARY\n"
            << "DATASET STRUCTURED_GRID\n"
            << "DIMENSIONS " << NX + 1 << " " << NY + 1 << " " << 1 << "\n"
            << "ORIGIN " << 0 << " " << 0 << " " << 0 << "\n"
            << "POINTS " << (NX + 1) * (NY + 1) * 1 << " float\n";
    scalarsu << "SCALARS u float 1\nLOOKUP_TABLE default\n";
    scalarsv << "SCALARS v float 1\nLOOKUP_TABLE default\n";
    scalarsp << "SCALARS p float 1\nLOOKUP_TABLE default\n";

    headerOffset = header.str().size();

    // This offset determines the cut of the plane in the z direction
    int zOffset = 0;

    if(rank == 0){
        // Write header to file
        MPI_File_write_at(fh, 0, header.str().c_str(), header.str().size(), MPI_CHAR, MPI_STATUS_IGNORE);
        // Write scalars u header after metadata
        MPI_File_write_at(fh, headerOffset, scalarsu.str().c_str(), scalarsu.str().size(), MPI_CHAR, MPI_STATUS_IGNORE);
        // Write scalars v header after scalars u
        MPI_File_write_at(fh, headerOffset + scalarsu.str().size() + NY, scalarsv.str().c_str(), scalarsv.str().size(), MPI_CHAR, MPI_STATUS_IGNORE);
        // Write scalars p header after scalars v
        MPI_File_write_at(fh, headerOffset + scalarsu.str().size() + 2 * NY + scalarsv.str().size(), scalarsp.str().c_str(), scalarsp.str().size(), MPI_CHAR, MPI_STATUS_IGNORE);
    }

    // Write the data
    // Each processor must write all the rows of its local grid but in a different possition in the line
    // depending on the processor rank
    // go through all the rows of the processor
    colOffsetX = coords[0];
    for (size_t i = 0; i < newDimY_x; i++){
        // Write the whole line of the processor
        rowOffsetX = headerOffset + (rank * other_dim_y_x) * i;

        MPI_File_write_at_all(fh, rowOffsetX + colOffsetX, &grid_loc_x[zOffset], newDimX_x, MPI_DOUBLE, MPI_STATUS_IGNORE);
    }

    colOffsetY = coords[0];
    for (size_t i = 0; i < newDimY_y; i++){
        // Write the whole line of the processor
        rowOffsetY = headerOffset + NX + 1 + (rank * other_dim_y_y) * i;
        MPI_File_write_at_all(fh, rowOffsetY + colOffsetY, &grid_loc_y[zOffset], newDimX_y, MPI_DOUBLE, MPI_STATUS_IGNORE);
    }

    // TODO: Write preassure

}
