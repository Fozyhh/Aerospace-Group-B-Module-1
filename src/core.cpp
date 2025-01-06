#include "constants.hpp"
#include "core.hpp"
#include "utils.hpp"
#include <cmath>
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

    MPI_Cart_shift(cart_comm, 1, 1, &neighbors[3], &neighbors[1]);
    std::cout << rank << "- neigh: " << neighbors[1] << neighbors[3] << std::endl;

    if (NX % PX != 0 && coords[0] == PX- 1)
    {
        dim_x_x+= NX%PX;
    }
    if ((NY + 1) % PY != 0 && coords[1] == PY-1)
    {
        dim_y_x+= (NY+1)%PY;
    }

    if ((NX + 1) % PX != 0 && coords[0] == PX - 1)
    {
        dim_x_y+= (NX+1)%PX;
    }
    if ((NY) % PY != 0 && coords[1] == PY-1)
    {
        dim_y_y+= NY%PY;
    }

    if ((NX + 1) % PX != 0 && coords[0] == PX - 1)
    {
        dim_x_z+= (NX+1)%PX;
    }
    if ((NY + 1) % PY != 0 && coords[1] == PY-1)
    {
        dim_y_z+= (NY+1)%PY;
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
    if(BX){
        if (coords[0] == 0)
            lbx++;

        if (coords[0] == PX - 1)
            rbx++;
    }
    if(BY){
        if (coords[1] == PY - 1)
            rby++;

        if (coords[1] == 0)
            lby++;
    }

    if (coords[0] == 0)
        firstX=1;

    if (coords[0] == PX - 1)
        lastX=1;

    if (coords[1] == PY - 1)
        lastY=1;

    if (coords[1] == 0)
        firstY=1;

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

void IcoNS::exchangeData(std::vector<Real> &grid_loc, int newDimX, int newDimY, int dim_z, MPI_Datatype MPI_face_x, MPI_Datatype MPI_face_y, int sameX, int sameY)
{
    // if(!dirichletX){
    //     p=0;
    // }
    if (!(BY && coords[1] == 0)){
        MPI_Isend(&grid_loc[dim_z * newDimY + dim_z + dim_z * sameY * firstY], 1, MPI_face_x, neighbors[3], 11, cart_comm, &reqs[2]);
    }

    if (!(BY && coords[1] == PY-1))
    {

        MPI_Irecv(&grid_loc[dim_z * newDimY + (newDimY - 1) * dim_z], 1, MPI_face_x, neighbors[1], 11, cart_comm, &reqs[3]);
        MPI_Wait(&reqs[3], MPI_SUCCESS);
        MPI_Isend(&grid_loc[dim_z * newDimY + (newDimY - 2 - sameY * lastY) * dim_z], 1, MPI_face_x, neighbors[1], 12, cart_comm, &reqs[3]);
    }
    if (!(BY &&  coords[1] == 0))
    {
        MPI_Irecv(&grid_loc[dim_z * newDimY], 1, MPI_face_x, neighbors[3], 12, cart_comm, &reqs[2]);
        MPI_Wait(&reqs[2], MPI_SUCCESS);
    }
    // inviare al neighbours 0 dal processore con coords 0 vuol dire inviare al processo dall'altra parte >> inutile se dirichlet
    if (!(BX && coords[0] == 0))
    {
        MPI_Isend(&grid_loc[(newDimY)*dim_z + (newDimY)*dim_z * sameX * firstX], 1, MPI_face_y, neighbors[0], 10, cart_comm, &reqs[0]);

        //MPI_Wait(&reqs[0], MPI_SUCCESS);
        //MPI_Irecv(&grid_loc[(dim_z)*newDimY * (newDimX - 1)], 1, MPI_face_y, neighbors[0], neighbors[0], cart_comm, &reqs[0]);

    }
    if (!(BX && coords[0] == PX - 1)){
        MPI_Irecv(&grid_loc[(dim_z)*newDimY * (newDimX - 1)], 1, MPI_face_y, neighbors[2], 10, cart_comm, &reqs[1]);
        MPI_Wait(&reqs[1], MPI_SUCCESS);
        MPI_Isend(&grid_loc[newDimY * dim_z * (newDimX - 2 - sameX * lastX)], 1, MPI_face_y, neighbors[2], 9, cart_comm, &reqs[1]);
        //MPI_Irecv(&grid_loc[0], 1, MPI_face_y, neighbors[2], neighbors[2], cart_comm, &reqs[1]);
        //MPI_Wait(&reqs[1], MPI_SUCCESS);
    }
    if (!(BX && coords[0] == 0)){
        MPI_Irecv(&grid_loc[0], 1, MPI_face_y, neighbors[0], 9  , cart_comm, &reqs[0]);
        MPI_Wait(&reqs[0], MPI_SUCCESS);
    }
    //MPI_Waitall(test,&reqs[0],MPI_SUCCESS);

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
        exchangeData(grid_loc_x, newDimX_x, newDimY_x,dim_z,MPI_face_x_x,MPI_face_y_x,0,1);
        exchangeData(grid_loc_y, newDimX_y, newDimY_y,dim_z,MPI_face_x_y,MPI_face_y_y,1,0);
        exchangeData(grid_loc_z, newDimX_z, newDimY_z,dim_z_z,MPI_face_x_z,MPI_face_y_z,1,1);
        x = L2_error(time); // every processor calculates his error not counting ghosts
        MPI_Barrier(cart_comm);
        error = 0.0;
        MPI_Reduce(&x, &error, 1, MPI_DOUBLE, MPI_SUM, 0, cart_comm);

        if (rank == 0)
        {
            error=sqrt(error);
            std::cout << " error: " << error << std::endl;
        }
        /*{
            for(int r =0; r<4 ;r++){
            if(rank==r){
                std::cout << "rsnk: " << rank << std::endl;
                for (int in = 0 ; in < newDimX_x ; in++)
                {
                    for(int j=0; j < newDimY_x ;j++){
                        for(int k=0; k <dim_z;k++){
                            std::cout << grid_loc_x[in * newDimY_x * dim_z + j * dim_z + k] << " ";
                        }
                        std::cout<< std::endl;
                    }
                    std::cout<< std::endl;
                }

            }
            MPI_Barrier(cart_comm);
            }
        }*/
        // reduce
        solve_time_step(time); // adapt cycles to skip ghosts
        MPI_Barrier(cart_comm);
        time += DT;
        i++;
    }
    output();
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

    // std::cout << error_comp_X(t) << std::endl;
    // std::cout << error_comp_Y(t) << std::endl;
    // std::cout << error_comp_Z(t) << std::endl << std::endl;

    return error;//sqrt(error);
}

Real IcoNS::error_comp_X(const Real t)
{
    Real error = 0.0;
    int offset_x = coords[0] * other_dim_x_x -1;
    int offset_y = coords[1] * other_dim_y_x -1;
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
    int offset_y = coords[1] * other_dim_y_y -1;
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
    int offset_y = coords[1] * other_dim_y_z -1;
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

/*
vtk file : 3 slices for x=0, y=0, z=0
    all of it for all vars u,v,w,p
two profile*.dat , containing 3 1D arrays of the solution at time=0 and time=final
*/
void IcoNS::parse_input(const std::string& input_file) {
    std::ifstream file(input_file);
    if (!file.is_open()) {
        std::cerr << "Process " << rank << ": Error - Cannot open file " << input_file << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    std::string line;
    // Skip comments and empty lines until we find RE
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream iss(line);
        if (!(iss >> RE)) continue;
        break;
    }

    // Skip comments and empty lines until we find DT
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream iss(line);
        if (!(iss >> DT)) continue;
        break;
    }

    // Skip comments and empty lines until we find T
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream iss(line);
        if (!(iss >> T)) continue;
        break;
    }

    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;

        std::istringstream iss(line);
        std::string sx, sy, sz;

        if (!(iss >> sx >> sy >> sz)) continue;

        try {
            LX = evaluateExpression(sx);
            LY = evaluateExpression(sy);
            LZ = evaluateExpression(sz);
            break;
        } catch (const std::exception& e) {
            std::cerr << "Error parsing dimensions: " << e.what() << std::endl;
            continue;
        }
    }

    // Skip comments and empty lines until we find grid points
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream iss(line);
        if (!(iss >> NX >> NY >> NZ)) continue;
        break;
    }

    // Skip comments and empty lines until we find process grid
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream iss(line);
        if (!(iss >> PX >> PY >> PZ)) continue;
        break;
    }

    // Skip comments and empty lines until we find boundaries

    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream iss(line);
        int x_boundary, y_boundary, z_boundary;
        if (!(iss >> x_boundary >> y_boundary >> z_boundary)) continue;

        // Convert integers to booleans and assign to extern variables
        BX = static_cast<bool>(x_boundary);
        BY = static_cast<bool>(y_boundary);
        BZ = static_cast<bool>(z_boundary);

        // Update periods array for MPI cart create
        periods[0] = !BX;
        periods[1] = !BY;

        break;
    }

    if (rank == 0) {
        std::cout << "\nConfiguration:\n"
                  << "Reynolds number: " << RE << "\n"
                  << "Time step: " << DT << "\n"
                  << "Total time: " << T << "\n"
                  << "Domain size: " << LX << " x " << LY << " x " << LZ << "\n"
                  << "Grid points: " << NX << " x " << NY << " x " << NZ << "\n"
                  << "Process grid: " << PX << " x " << PY << " x " << PZ << "\n"
                  << "Boundary conditions: " << BX << " " << BY << " " << BZ << "\n";
    }

    // Verify that the process grid matches the number of processes
    if (PX * PY * PZ != size) {
        if (rank == 0) {
            std::cerr << "Error: Number of processes (" << size
                      << ") does not match process grid ("
                      << PX << " x " << PY << " x " << PZ
                      << " = " << (PX * PY * PZ) << ")\n";
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Calculate grid spacing
    DX = LX / NX;
    DY = LY / NY;
    DZ = LZ / NZ;

    // Make sure all processes have read the file before proceeding
    MPI_Barrier(MPI_COMM_WORLD);
}

void IcoNS::output(){
    exchangeData(grid_loc_x, newDimX_x, newDimY_x,dim_z,MPI_face_x_x,MPI_face_y_x,0,1);
    exchangeData(grid_loc_y, newDimX_y, newDimY_y,dim_z,MPI_face_x_y,MPI_face_y_y,1,0);
    exchangeData(grid_loc_z, newDimX_z, newDimY_z,dim_z_z,MPI_face_x_z,MPI_face_y_z,1,1);
    MPI_Barrier(cart_comm);
    output_x();
    output_y();
    output_z(); 
}


void IcoNS::output_x(){
    MPI_File fh;
    MPI_Offset offset = 0;
    const float x_middle = /*SX + */LX / 2;
    if(rank==0)
        std::remove("solution_x.vtk");
    MPI_File_open(cart_comm, "solution_x.vtk", MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);

    //===========================================
    // Header Writing (Rank 0 only)
    //===========================================
    std::ostringstream header1, header2, header3,header4;
    header1 << std::fixed << std::setprecision(6);
    header2 << std::fixed << std::setprecision(6);
    header4 << std::fixed << std::setprecision(6);

    // VTK metadata
    header1 << "# vtk DataFile Version 3.0\n"
                << "Solution x\n"
                << "ASCII\n"
                << "DATASET STRUCTURED_GRID\n"
                << "DIMENSIONS 1 " << NY + 1 << " " << NZ + 1 << "\n"
                << "POINTS " << (NY + 1) * (NZ + 1) << " float\n";

    // Define data format
    header2 << "\nPOINT_DATA " << (NY + 1) * (NZ + 1) << "\n"
                << "SCALARS u float\n"
                << "LOOKUP_TABLE default\n";

    header3     << "SCALARS y float\n"
                << "LOOKUP_TABLE default\n";
    header4     << "SCALARS z float\n"
                << "LOOKUP_TABLE default\n";
    // Write header to file
    // std::string header_str = full_header.str();
    // MPI_File_write_at(fh, 0, header_str.c_str(), header_str.size(), MPI_CHAR, MPI_STATUS_IGNORE);
    // offset = header_str.size();

    // Sync offset across all processes
    // MPI_Bcast(&offset, 1, MPI_OFFSET, 0, cart_comm);

    //===========================================
    // Data Collection
    //===========================================
    std::vector<float> minepoints(PX*PY, 0.0f);
    std::vector<float> minevaluesx(PX*PY, 0.0f);
    std::vector<float> minevaluesy(PX*PY, 0.0f);
    std::vector<float> minevaluesz(PX*PY, 0.0f);
    std::vector<float> globalpoints(PX*PY);
    std::vector<float> globalvaluesx(PX*PY);
    std::vector<float> globalvaluesy(PX*PY);
    std::vector<float> globalvaluesz(PX*PY);

    std::ostringstream points, valuesx,valuesy,valuesz;

    // Find process containing middle slice
    int x_index = static_cast<int>(x_middle);
    int offset_x_x = coords[0] * other_dim_x_x -1;
    int offset_y_x = coords[1] * other_dim_y_x -1;
    int offset_x_y = coords[0] * other_dim_x_y -1;
    int offset_y_y = coords[1] * other_dim_y_y -1;
    int offset_x_z = coords[0] * other_dim_x_z -1;
    int offset_y_z = coords[1] * other_dim_y_z -1;

    if (x_index >= offset_x_x && x_index < dim_x_x + offset_x_x) {

        int local_x_x = x_index - offset_x_x + 1;
        int local_x_y = x_index - offset_x_y + 1;
        int local_x_z = x_index - offset_x_z + 1;

        for(int j = 1; j < newDimY_x - 1; j++){
            for(int k=0; k < dim_z; k++){
                
                // Write grid points coordinate
                points << x_middle << " "
                           << static_cast<float>(j + offset_y_x) * DY << " "
                           << static_cast<float>(k) * DZ << "\n";
                
                //valuesx[rank*(dim_y_x * dim_z + dim_z) + ((j-1) * dim_z + k)] = grid_loc_y[local_x* newDimY_x * dim_z + j * dim_z + k];
                valuesx << grid_loc_x[local_x_x* newDimY_x * dim_z + j * dim_z + k] << "\n";

                if(lby &&j==1){
                    valuesy << (boundary.boundary_value_v[2]->value(x_index,j + offset_y_y-0.5,k,T) + boundary.boundary_value_v[2]->value(x_index + 1,j + offset_y_y-0.5,k,T))/2 << "\n"; 
                }else if(rby && j==newDimY_x-2){
                    valuesy << (boundary.boundary_value_v[3]->value(x_index,j + offset_y_y-0.5,k,T) + boundary.boundary_value_v[3]->value(x_index + 1,j + offset_y_y-0.5,k,T))/2 << "\n"; 
                }else{
                    valuesy << (grid_loc_y[local_x_y*newDimY_y * dim_z + j * dim_z + k] + grid_loc_y[local_x_y*newDimY_y * dim_z + (j-1) * dim_z + k] +
                                grid_loc_y[(local_x_y+1)*newDimY_y * dim_z + j * dim_z + k] + grid_loc_y[(local_x_y+1)*newDimY_y * dim_z + (j-1) * dim_z + k])/4 << "\n";
                }  

                if(k==0){
                    valuesz << boundary.boundary_value_w[4]->value(x_index + 0.5,j + offset_y_z,k - 0.5,T) << "\n"; 
                }else if(k==dim_z -1){
                    valuesz << boundary.boundary_value_w[5]->value(x_index + 0.5,j + offset_y_z,k - 0.5,T) << "\n"; 
                }else{
                    valuesz << (grid_loc_z[local_x_z*newDimY_z * dim_z_z + j * dim_z_z + k] + grid_loc_z[local_x_z*newDimY_z * dim_z_z + j * dim_z_z + k-1] +
                                grid_loc_z[(local_x_z+1)*newDimY_z * dim_z_z + j * dim_z_z + k] + grid_loc_z[(local_x_z+1)*newDimY_z * dim_z_z + j * dim_z_z + k-1])/4 << "\n";
                }              
            }
        }
    }
    minepoints[rank] = points.str().size();
    minevaluesx[rank] = valuesx.str().size();
    minevaluesy[rank] = valuesy.str().size();
    minevaluesz[rank] = valuesz.str().size();
    // Gather data from all processes
    MPI_Allreduce(minepoints.data(), globalpoints.data(), minepoints.size(), MPI_FLOAT, MPI_SUM, cart_comm);
    MPI_Allreduce(minevaluesx.data(), globalvaluesx.data(), minevaluesx.size(), MPI_FLOAT, MPI_SUM, cart_comm);
    MPI_Allreduce(minevaluesy.data(), globalvaluesy.data(), minevaluesy.size(), MPI_FLOAT, MPI_SUM, cart_comm);
    MPI_Allreduce(minevaluesz.data(), globalvaluesz.data(), minevaluesz.size(), MPI_FLOAT, MPI_SUM, cart_comm);

    //===========================================
    // Data Writing (Rank 0 only)
    //===========================================
    int my_points=0, allpoints=0, my_valuesx=0, allvaluesx=0,my_valuesy=0,allvaluesy=0,my_valuesz=0,allvaluesz=0;
    for(int i=0; i < PX*PY; i++){
        allpoints += globalpoints[i];
        allvaluesx += globalvaluesx[i];
        allvaluesy += globalvaluesy[i];
        allvaluesz += globalvaluesz[i];
        if(rank > i){
            my_points += globalpoints[i];
            my_valuesx += globalvaluesx[i];
            my_valuesy += globalvaluesy[i];
            my_valuesz += globalvaluesz[i];
        }
    }
    if (rank == 0) {
        MPI_File_write_at(fh, offset, header1.str().c_str(), header1.str().size(), MPI_CHAR, MPI_STATUS_IGNORE);
    }
    offset += header1.str().size();
    MPI_File_write_at(fh, offset + my_points, points.str().c_str(), points.str().size(), MPI_CHAR, MPI_STATUS_IGNORE);
    offset += allpoints;
    if(rank==0){
        MPI_File_write_at(fh, offset, header2.str().c_str(), header2.str().size(), MPI_CHAR, MPI_STATUS_IGNORE);  
    }
    offset += header2.str().size() ;
    MPI_File_write_at(fh, offset + my_valuesx, valuesx.str().c_str(), valuesx.str().size(), MPI_CHAR, MPI_STATUS_IGNORE);
    offset += allvaluesx;
    if (rank == 0) {
        MPI_File_write_at(fh, offset, header3.str().c_str(), header3.str().size(), MPI_CHAR, MPI_STATUS_IGNORE);
    }
    offset += header3.str().size();
    MPI_File_write_at(fh, offset + my_valuesy, valuesy.str().c_str(), valuesy.str().size(), MPI_CHAR, MPI_STATUS_IGNORE);
    offset += allvaluesy;
    if (rank == 0) {
        MPI_File_write_at(fh, offset, header4.str().c_str(), header4.str().size(), MPI_CHAR, MPI_STATUS_IGNORE);
    }
    offset += header4.str().size();
    MPI_File_write_at(fh, offset + my_valuesz, valuesz.str().c_str(), valuesz.str().size(), MPI_CHAR, MPI_STATUS_IGNORE);
    offset += allvaluesz;

    MPI_File_close(&fh);
}

void IcoNS::output_y(){
    MPI_File fh;
    MPI_Offset offset = 0;
    const float x_middle = /*SY + */LY / 2;
    if(rank==0)
        std::remove("solution_y.vtk");
    MPI_File_open(cart_comm, "solution_y.vtk", MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);

    //===========================================
    // Header Writing (Rank 0 only)
    //===========================================
    std::ostringstream header1, header2, header3,header4;
    header1 << std::fixed << std::setprecision(6);
    header2 << std::fixed << std::setprecision(6);
    header4 << std::fixed << std::setprecision(6);

    // VTK metadata
    header1 << "# vtk DataFile Version 3.0\n"
                << "Solution x\n"
                << "ASCII\n"
                << "DATASET STRUCTURED_GRID\n"
                << "DIMENSIONS 1 " << NY + 1 << " " << NZ + 1 << "\n"
                << "POINTS " << (NY + 1) * (NZ + 1) << " float\n";

    // Define data format
    header2 << "\nPOINT_DATA " << (NY + 1) * (NZ + 1) << "\n"
                << "SCALARS u float\n"
                << "LOOKUP_TABLE default\n";

    header3     << "SCALARS y float\n"
                << "LOOKUP_TABLE default\n";
    header4     << "SCALARS z float\n"
                << "LOOKUP_TABLE default\n";
    // Write header to file
    // std::string header_str = full_header.str();
    // MPI_File_write_at(fh, 0, header_str.c_str(), header_str.size(), MPI_CHAR, MPI_STATUS_IGNORE);
    // offset = header_str.size();

    // Sync offset across all processes
    // MPI_Bcast(&offset, 1, MPI_OFFSET, 0, cart_comm);

    //===========================================
    // Data Collection
    //===========================================
    std::vector<float> minepoints(PX*PY, 0.0f);
    std::vector<float> minevaluesx(PX*PY, 0.0f);
    std::vector<float> minevaluesy(PX*PY, 0.0f);
    std::vector<float> minevaluesz(PX*PY, 0.0f);
    std::vector<float> globalpoints(PX*PY);
    std::vector<float> globalvaluesx(PX*PY);
    std::vector<float> globalvaluesy(PX*PY);
    std::vector<float> globalvaluesz(PX*PY);

    std::ostringstream points, valuesx,valuesy,valuesz;

    // Find process containing middle slice
    int x_index = static_cast<int>(x_middle);
    int offset_x_x = coords[0] * other_dim_x_x -1;
    int offset_y_x = coords[1] * other_dim_y_x -1;
    int offset_x_y = coords[0] * other_dim_x_y -1;
    int offset_y_y = coords[1] * other_dim_y_y -1;
    int offset_x_z = coords[0] * other_dim_x_z -1;
    int offset_y_z = coords[1] * other_dim_y_z -1;

    if (x_index >= offset_x_x && x_index < dim_x_x + offset_x_x) {

        int local_x_x = x_index - offset_x_x + 1;
        int local_x_y = x_index - offset_x_y + 1;
        int local_x_z = x_index - offset_x_z + 1;

        for(int j = 1; j < newDimY_x - 1; j++){
            for(int k=0; k < dim_z; k++){
                
                // Write grid points coordinate
                points << x_middle << " "
                           << static_cast<float>(j + offset_y_x) * DY << " "
                           << static_cast<float>(k) * DZ << "\n";
                
                //valuesx[rank*(dim_y_x * dim_z + dim_z) + ((j-1) * dim_z + k)] = grid_loc_y[local_x* newDimY_x * dim_z + j * dim_z + k];
                valuesx << grid_loc_x[local_x_x* newDimY_x * dim_z + j * dim_z + k] << "\n";

                if(lby &&j==1){
                    valuesy << (boundary.boundary_value_v[2]->value(x_index,j + offset_y_y-0.5,k,T) + boundary.boundary_value_v[2]->value(x_index + 1,j + offset_y_y-0.5,k,T))/2 << "\n"; 
                }else if(rby && j==newDimY_x-2){
                    valuesy << (boundary.boundary_value_v[3]->value(x_index,j + offset_y_y-0.5,k,T) + boundary.boundary_value_v[3]->value(x_index + 1,j + offset_y_y-0.5,k,T))/2 << "\n"; 
                }else{
                    valuesy << (grid_loc_y[local_x_y*newDimY_y * dim_z + j * dim_z + k] + grid_loc_y[local_x_y*newDimY_y * dim_z + (j-1) * dim_z + k] +
                                grid_loc_y[(local_x_y+1)*newDimY_y * dim_z + j * dim_z + k] + grid_loc_y[(local_x_y+1)*newDimY_y * dim_z + (j-1) * dim_z + k])/4 << "\n";
                }  

                if(k==0){
                    valuesz << boundary.boundary_value_w[4]->value(x_index + 0.5,j + offset_y_z,k - 0.5,T) << "\n"; 
                }else if(k==dim_z -1){
                    valuesz << boundary.boundary_value_w[5]->value(x_index + 0.5,j + offset_y_z,k - 0.5,T) << "\n"; 
                }else{
                    valuesz << (grid_loc_z[local_x_z*newDimY_z * dim_z_z + j * dim_z_z + k] + grid_loc_z[local_x_z*newDimY_z * dim_z_z + j * dim_z_z + k-1] +
                                grid_loc_z[(local_x_z+1)*newDimY_z * dim_z_z + j * dim_z_z + k] + grid_loc_z[(local_x_z+1)*newDimY_z * dim_z_z + j * dim_z_z + k-1])/4 << "\n";
                }              
            }
        }
    }
    minepoints[rank] = points.str().size();
    minevaluesx[rank] = valuesx.str().size();
    minevaluesy[rank] = valuesy.str().size();
    minevaluesz[rank] = valuesz.str().size();
    // Gather data from all processes
    MPI_Allreduce(minepoints.data(), globalpoints.data(), minepoints.size(), MPI_FLOAT, MPI_SUM, cart_comm);
    MPI_Allreduce(minevaluesx.data(), globalvaluesx.data(), minevaluesx.size(), MPI_FLOAT, MPI_SUM, cart_comm);
    MPI_Allreduce(minevaluesy.data(), globalvaluesy.data(), minevaluesy.size(), MPI_FLOAT, MPI_SUM, cart_comm);
    MPI_Allreduce(minevaluesz.data(), globalvaluesz.data(), minevaluesz.size(), MPI_FLOAT, MPI_SUM, cart_comm);

    //===========================================
    // Data Writing (Rank 0 only)
    //===========================================
    int my_points=0, allpoints=0, my_valuesx=0, allvaluesx=0,my_valuesy=0,allvaluesy=0,my_valuesz=0,allvaluesz=0;
    for(int i=0; i < PX*PY; i++){
        allpoints += globalpoints[i];
        allvaluesx += globalvaluesx[i];
        allvaluesy += globalvaluesy[i];
        allvaluesz += globalvaluesz[i];
        if(rank > i){
            my_points += globalpoints[i];
            my_valuesx += globalvaluesx[i];
            my_valuesy += globalvaluesy[i];
            my_valuesz += globalvaluesz[i];
        }
    }
    if (rank == 0) {
        MPI_File_write_at(fh, offset, header1.str().c_str(), header1.str().size(), MPI_CHAR, MPI_STATUS_IGNORE);
    }
    offset += header1.str().size();
    MPI_File_write_at(fh, offset + my_points, points.str().c_str(), points.str().size(), MPI_CHAR, MPI_STATUS_IGNORE);
    offset += allpoints;
    if(rank==0){
        MPI_File_write_at(fh, offset, header2.str().c_str(), header2.str().size(), MPI_CHAR, MPI_STATUS_IGNORE);  
    }
    offset += header2.str().size() ;
    MPI_File_write_at(fh, offset + my_valuesx, valuesx.str().c_str(), valuesx.str().size(), MPI_CHAR, MPI_STATUS_IGNORE);
    offset += allvaluesx;
    if (rank == 0) {
        MPI_File_write_at(fh, offset, header3.str().c_str(), header3.str().size(), MPI_CHAR, MPI_STATUS_IGNORE);
    }
    offset += header3.str().size();
    MPI_File_write_at(fh, offset + my_valuesy, valuesy.str().c_str(), valuesy.str().size(), MPI_CHAR, MPI_STATUS_IGNORE);
    offset += allvaluesy;
    if (rank == 0) {
        MPI_File_write_at(fh, offset, header4.str().c_str(), header4.str().size(), MPI_CHAR, MPI_STATUS_IGNORE);
    }
    offset += header4.str().size();
    MPI_File_write_at(fh, offset + my_valuesz, valuesz.str().c_str(), valuesz.str().size(), MPI_CHAR, MPI_STATUS_IGNORE);
    offset += allvaluesz;

    MPI_File_close(&fh);
}
void IcoNS::output_z(){
    MPI_File fh;
    MPI_Offset offset = 0;
    const float x_middle = /*SZ + */LZ / 2;
    if(rank==0)
        std::remove("solution_z.vtk");
    MPI_File_open(cart_comm, "solution_z.vtk", MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);

    //===========================================
    // Header Writing (Rank 0 only)
    //===========================================
    std::ostringstream header1, header2, header3,header4;
    header1 << std::fixed << std::setprecision(6);
    header2 << std::fixed << std::setprecision(6);
    header4 << std::fixed << std::setprecision(6);

    // VTK metadata
    header1 << "# vtk DataFile Version 3.0\n"
                << "Solution x\n"
                << "ASCII\n"
                << "DATASET STRUCTURED_GRID\n"
                << "DIMENSIONS 1 " << NY + 1 << " " << NZ + 1 << "\n"
                << "POINTS " << (NY + 1) * (NZ + 1) << " float\n";

    // Define data format
    header2 << "\nPOINT_DATA " << (NY + 1) * (NZ + 1) << "\n"
                << "SCALARS u float\n"
                << "LOOKUP_TABLE default\n";

    header3     << "SCALARS y float\n"
                << "LOOKUP_TABLE default\n";
    header4     << "SCALARS z float\n"
                << "LOOKUP_TABLE default\n";
    // Write header to file
    // std::string header_str = full_header.str();
    // MPI_File_write_at(fh, 0, header_str.c_str(), header_str.size(), MPI_CHAR, MPI_STATUS_IGNORE);
    // offset = header_str.size();

    // Sync offset across all processes
    // MPI_Bcast(&offset, 1, MPI_OFFSET, 0, cart_comm);

    //===========================================
    // Data Collection
    //===========================================
    std::vector<float> minepoints(PX*PY, 0.0f);
    std::vector<float> minevaluesx(PX*PY, 0.0f);
    std::vector<float> minevaluesy(PX*PY, 0.0f);
    std::vector<float> minevaluesz(PX*PY, 0.0f);
    std::vector<float> globalpoints(PX*PY);
    std::vector<float> globalvaluesx(PX*PY);
    std::vector<float> globalvaluesy(PX*PY);
    std::vector<float> globalvaluesz(PX*PY);

    std::ostringstream points, valuesx,valuesy,valuesz;

    // Find process containing middle slice
    int x_index = static_cast<int>(x_middle);
    int offset_x_x = coords[0] * other_dim_x_x -1;
    int offset_y_x = coords[1] * other_dim_y_x -1;
    int offset_x_y = coords[0] * other_dim_x_y -1;
    int offset_y_y = coords[1] * other_dim_y_y -1;
    int offset_x_z = coords[0] * other_dim_x_z -1;
    int offset_y_z = coords[1] * other_dim_y_z -1;

    if (x_index >= offset_x_x && x_index < dim_x_x + offset_x_x) {

        int local_x_x = x_index - offset_x_x + 1;
        int local_x_y = x_index - offset_x_y + 1;
        int local_x_z = x_index - offset_x_z + 1;

        for(int j = 1; j < newDimY_x - 1; j++){
            for(int k=0; k < dim_z; k++){
                
                // Write grid points coordinate
                points << x_middle << " "
                           << static_cast<float>(j + offset_y_x) * DY << " "
                           << static_cast<float>(k) * DZ << "\n";
                
                //valuesx[rank*(dim_y_x * dim_z + dim_z) + ((j-1) * dim_z + k)] = grid_loc_y[local_x* newDimY_x * dim_z + j * dim_z + k];
                valuesx << grid_loc_x[local_x_x* newDimY_x * dim_z + j * dim_z + k] << "\n";

                if(lby &&j==1){
                    valuesy << (boundary.boundary_value_v[2]->value(x_index,j + offset_y_y-0.5,k,T) + boundary.boundary_value_v[2]->value(x_index + 1,j + offset_y_y-0.5,k,T))/2 << "\n"; 
                }else if(rby && j==newDimY_x-2){
                    valuesy << (boundary.boundary_value_v[3]->value(x_index,j + offset_y_y-0.5,k,T) + boundary.boundary_value_v[3]->value(x_index + 1,j + offset_y_y-0.5,k,T))/2 << "\n"; 
                }else{
                    valuesy << (grid_loc_y[local_x_y*newDimY_y * dim_z + j * dim_z + k] + grid_loc_y[local_x_y*newDimY_y * dim_z + (j-1) * dim_z + k] +
                                grid_loc_y[(local_x_y+1)*newDimY_y * dim_z + j * dim_z + k] + grid_loc_y[(local_x_y+1)*newDimY_y * dim_z + (j-1) * dim_z + k])/4 << "\n";
                }  

                if(k==0){
                    valuesz << boundary.boundary_value_w[4]->value(x_index + 0.5,j + offset_y_z,k - 0.5,T) << "\n"; 
                }else if(k==dim_z -1){
                    valuesz << boundary.boundary_value_w[5]->value(x_index + 0.5,j + offset_y_z,k - 0.5,T) << "\n"; 
                }else{
                    valuesz << (grid_loc_z[local_x_z*newDimY_z * dim_z_z + j * dim_z_z + k] + grid_loc_z[local_x_z*newDimY_z * dim_z_z + j * dim_z_z + k-1] +
                                grid_loc_z[(local_x_z+1)*newDimY_z * dim_z_z + j * dim_z_z + k] + grid_loc_z[(local_x_z+1)*newDimY_z * dim_z_z + j * dim_z_z + k-1])/4 << "\n";
                }              
            }
        }
    }
    minepoints[rank] = points.str().size();
    minevaluesx[rank] = valuesx.str().size();
    minevaluesy[rank] = valuesy.str().size();
    minevaluesz[rank] = valuesz.str().size();
    // Gather data from all processes
    MPI_Allreduce(minepoints.data(), globalpoints.data(), minepoints.size(), MPI_FLOAT, MPI_SUM, cart_comm);
    MPI_Allreduce(minevaluesx.data(), globalvaluesx.data(), minevaluesx.size(), MPI_FLOAT, MPI_SUM, cart_comm);
    MPI_Allreduce(minevaluesy.data(), globalvaluesy.data(), minevaluesy.size(), MPI_FLOAT, MPI_SUM, cart_comm);
    MPI_Allreduce(minevaluesz.data(), globalvaluesz.data(), minevaluesz.size(), MPI_FLOAT, MPI_SUM, cart_comm);

    //===========================================
    // Data Writing (Rank 0 only)
    //===========================================
    int my_points=0, allpoints=0, my_valuesx=0, allvaluesx=0,my_valuesy=0,allvaluesy=0,my_valuesz=0,allvaluesz=0;
    for(int i=0; i < PX*PY; i++){
        allpoints += globalpoints[i];
        allvaluesx += globalvaluesx[i];
        allvaluesy += globalvaluesy[i];
        allvaluesz += globalvaluesz[i];
        if(rank > i){
            my_points += globalpoints[i];
            my_valuesx += globalvaluesx[i];
            my_valuesy += globalvaluesy[i];
            my_valuesz += globalvaluesz[i];
        }
    }
    if (rank == 0) {
        MPI_File_write_at(fh, offset, header1.str().c_str(), header1.str().size(), MPI_CHAR, MPI_STATUS_IGNORE);
    }
    offset += header1.str().size();
    MPI_File_write_at(fh, offset + my_points, points.str().c_str(), points.str().size(), MPI_CHAR, MPI_STATUS_IGNORE);
    offset += allpoints;
    if(rank==0){
        MPI_File_write_at(fh, offset, header2.str().c_str(), header2.str().size(), MPI_CHAR, MPI_STATUS_IGNORE);  
    }
    offset += header2.str().size() ;
    MPI_File_write_at(fh, offset + my_valuesx, valuesx.str().c_str(), valuesx.str().size(), MPI_CHAR, MPI_STATUS_IGNORE);
    offset += allvaluesx;
    if (rank == 0) {
        MPI_File_write_at(fh, offset, header3.str().c_str(), header3.str().size(), MPI_CHAR, MPI_STATUS_IGNORE);
    }
    offset += header3.str().size();
    MPI_File_write_at(fh, offset + my_valuesy, valuesy.str().c_str(), valuesy.str().size(), MPI_CHAR, MPI_STATUS_IGNORE);
    offset += allvaluesy;
    if (rank == 0) {
        MPI_File_write_at(fh, offset, header4.str().c_str(), header4.str().size(), MPI_CHAR, MPI_STATUS_IGNORE);
    }
    offset += header4.str().size();
    MPI_File_write_at(fh, offset + my_valuesz, valuesz.str().c_str(), valuesz.str().size(), MPI_CHAR, MPI_STATUS_IGNORE);
    offset += allvaluesz;

    MPI_File_close(&fh);
}