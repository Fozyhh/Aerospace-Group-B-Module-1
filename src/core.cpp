#include "constants.hpp"
#include "core.hpp"
#include "utils.hpp"
#include <cmath>
#include <math.h>
#include <fstream>
#include <string>
#include <memory>
#include <mpi.h>

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
    setBoundaryConditions();

    setParallelization();
   
    for (int i = 0; i < NX * NY * (NZ/2 + 1); i++)
    {
        helper[i][0] = 0.0;
        helper[i][1] = 0.0;
    }

    boundary.initializeBoundary(
        dim_x_x, dim_y_x, dim_x_y, dim_y_y,
        dim_x_z, dim_y_z, dim_z, dim_z_z,
        newDimX_x, newDimY_x, newDimX_y, newDimY_y,
        newDimX_z, newDimY_z,
        c2d->zSize);
}


void IcoNS::setBoundaryConditions(){
    std::shared_ptr<BoundaryFunction> u_func;
    std::shared_ptr<BoundaryFunction> v_func;
    std::shared_ptr<BoundaryFunction> w_func;

    if(testCase==1){
        u_func = std::make_shared<Dirichlet>([&](Real x, Real y, Real z, Real t)
                                                { return 0; });
        v_func = std::make_shared<Dirichlet>([&](Real x, Real y, Real z, Real t)
                                                { 
                                                    //TODO: da capire u(1,0,0) che cella è
                                                    if(x==1.0 && y==0 && z==0){
                                                        return 1;
                                                    }else{
                                                        return 0;
                                                    } 
                                                });
        w_func = std::make_shared<Dirichlet>([&](Real x, Real y, Real z, Real t)
                                                { return 0; });
    }else if(testCase==2){
        u_func = std::make_shared<Dirichlet>([&](Real x, Real y, Real z, Real t)
                                                { return 0; });
        v_func = std::make_shared<Dirichlet>([&](Real x, Real y, Real z, Real t)
                                                { 
                                                    if(x==-0.5 && y==0 && z==0){
                                                        return 1;
                                                    }else{
                                                        return 0;
                                                    } 
                                                });
        w_func = std::make_shared<Dirichlet>([&](Real x, Real y, Real z, Real t)
                                                { return 0; });
    }else{
        u_func = std::make_shared<Dirichlet>([&](Real x, Real y, Real z, Real t)
                                                { return std::sin((x + 0.5) * DX) * std::cos(y * DY) * std::sin(z * DZ) * std::sin(t); });
        v_func = std::make_shared<Dirichlet>([&](Real x, Real y, Real z, Real t)
                                                { return std::cos(x * DX) * std::sin((y + 0.5) * DY) * std::sin(z * DZ) * std::sin(t); });
        w_func = std::make_shared<Dirichlet>([&](Real x, Real y, Real z, Real t)
                                                { return 2 * (std::cos(x * DX) * std::cos(y * DY) * std::cos((z + 0.5) * DZ) * std::sin(t)); }); 
    }

    for (int i = 0; i < 6 /*nfaces*/; i++)
    {
        boundary.addFunction(U, u_func);
        boundary.addFunction(V, v_func);
        boundary.addFunction(W, w_func);
    }
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

    newDimX_x = dim_x_x + 2;
    newDimY_x = dim_y_x + 2;
    newDimX_y = dim_x_y + 2;
    newDimY_y = dim_y_y + 2;
    newDimX_z = dim_x_z + 2;
    newDimY_z = dim_y_z + 2;

    grid.u.resize(newDimX_x * newDimY_x * (NZ + 1), 0.0);
    grid.v.resize(newDimX_y * newDimY_y * (NZ + 1), 0.0);
    grid.w.resize(newDimX_z * newDimY_z * (NZ), 0.0);
    Y2_x.resize(newDimX_x * newDimY_x * (NZ + 1), 0.0);
    Y2_y.resize(newDimX_y * newDimY_y * (NZ + 1), 0.0);
    Y2_z.resize(newDimX_z * newDimY_z * (NZ), 0.0);

    Y3_x.resize(newDimX_x * newDimY_x * (NZ + 1), 0.0);
    Y3_y.resize(newDimX_y * newDimY_y * (NZ + 1), 0.0);
    Y3_z.resize(newDimX_z * newDimY_z * (NZ), 0.0);

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
    }
    if (!(BX && coords[0] == PX - 1)){
        MPI_Irecv(&grid_loc[(dim_z)*newDimY * (newDimX - 1)], 1, MPI_face_y, neighbors[2], 10, cart_comm, &reqs[1]);
        MPI_Wait(&reqs[1], MPI_SUCCESS);
        MPI_Isend(&grid_loc[newDimY * dim_z * (newDimX - 2 - sameX * lastX)], 1, MPI_face_y, neighbors[2], 9, cart_comm, &reqs[1]);
    }
    if (!(BX && coords[0] == 0)){
        MPI_Irecv(&grid_loc[0], 1, MPI_face_y, neighbors[0], 9  , cart_comm, &reqs[0]);
        MPI_Wait(&reqs[0], MPI_SUCCESS);
    }
}

void IcoNS::solve()
{
    Real time = 0.0;
    Real error;
    int i = 0;

#ifdef VERBOSE
    std::cout << "Starting solver" << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
#endif
    double x=0;
    while (i < Nt)
    {
        boundary.update_boundary(grid.u, grid.v, grid.w, time);

        MPI_Barrier(cart_comm);
        exchangeData(grid.u, newDimX_x, newDimY_x,dim_z,MPI_face_x_x,MPI_face_y_x,0,1);
        exchangeData(grid.v, newDimX_y, newDimY_y,dim_z,MPI_face_x_y,MPI_face_y_y,1,0);
        exchangeData(grid.w, newDimX_z, newDimY_z,dim_z_z,MPI_face_x_z,MPI_face_y_z,1,1);
        
        if(testCase==0){
            L2_error(time);
        }
        
        solve_time_step(time);
        MPI_Barrier(cart_comm);
        time += DT;
        i++;
    }
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

void IcoNS::L2_error(const Real t)
{
    Real error = 0.0, totalError=0.0;

    error += error_comp_X(t);
    error += error_comp_Y(t);
    error += error_comp_Z(t);
    error += error_comp_P(t);

    MPI_Barrier(cart_comm);

    MPI_Reduce(&error, &totalError, 1, MPI_DOUBLE, MPI_SUM, 0, cart_comm);
    if (rank == 0)
    {
        totalError=sqrt(totalError);
        std::cout << " totalError: " << totalError << std::endl;
    }
    // std::cout << error_comp_X(t) << std::endl;
    // std::cout << error_comp_Y(t) << std::endl;
    // std::cout << error_comp_Z(t) << std::endl;
    // std::cout << error_comp_P(t) << std::endl << std::endl;

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
                error += ((grid.u[0 + (newDimY_x + 1) * dim_z] - exact_solution.value_x(0.5, 0, 0, t)) *
                          (grid.u[0 + (newDimY_x + 1) * dim_z] - exact_solution.value_x(0.5, 0, 0, t)) *
                          DX * DY * DZ / 8);

                for (int k = 1; k < NZ; k++)
                {
                    error += ((grid.u[(newDimY_x + 1) * dim_z + k] - exact_solution.value_x(0.5, 0, k, t)) *
                              (grid.u[(newDimY_x + 1) * dim_z + k] - exact_solution.value_x(0.5, 0, k, t)) *
                              DX * DY * DZ / 4);
                }

                error += ((grid.u[(newDimY_x + 1) * dim_z + NZ] - exact_solution.value_x(0.5, 0, NZ, t)) *
                            (grid.u[(newDimY_x + 1) * dim_z + NZ] - exact_solution.value_x(0.5, 0, NZ, t)) *
                            DX * DY * DZ / 8);
            }
            for (int j = 1 + lby; j < newDimY_x - 1 - rby; j++)
            {
                error += ((grid.u[(newDimY_x)*dim_z + j * (NZ + 1)] - exact_solution.value_x(0.5, j + offset_y, 0, t)) *
                          (grid.u[(newDimY_x)*dim_z + j * (NZ + 1)] - exact_solution.value_x(0.5, j + offset_y, 0, t)) *
                          DX * DY * DZ / 4);
                for (int k = 1; k < NZ; k++)
                {
                    error += ((grid.u[(newDimY_x)*dim_z + j * (NZ + 1) + k] - exact_solution.value_x(0.5, j + offset_y, k, t)) *
                              (grid.u[(newDimY_x)*dim_z + j * (NZ + 1) + k] - exact_solution.value_x(0.5, j + offset_y, k, t)) *
                              DX * DY * DZ / 2);
                }
                error += ((grid.u[(newDimY_x)*dim_z + j * (NZ + 1) + NZ] - exact_solution.value_x(0.5, j + offset_y, NZ, t)) *
                          (grid.u[(newDimY_x)*dim_z + j * (NZ + 1) + NZ] - exact_solution.value_x(0.5, j + offset_y, NZ, t)) *
                          DX * DY * DZ / 4);
            }
            if (rby)
            {
                error += ((grid.u[(newDimY_x)*dim_z + (newDimY_x - 2) * (NZ + 1)] - exact_solution.value_x(0.5, NY, 0, t)) *
                          (grid.u[(newDimY_x)*dim_z + (newDimY_x - 2) * (NZ + 1)] - exact_solution.value_x(0.5, NY, 0, t)) *
                          DX * DY * DZ / 8);

                for (int k = 1; k < NZ; k++)
                {
                    error += ((grid.u[(newDimY_x)*dim_z + (newDimY_x - 2) * (NZ + 1) + k] - exact_solution.value_x(0.5, NY, k, t)) *
                              (grid.u[(newDimY_x)*dim_z + (newDimY_x - 2) * (NZ + 1) + k] - exact_solution.value_x(0.5, NY, k, t)) *
                              DX * DY * DZ / 4);
                }

                error += ((grid.u[(newDimY_x)*dim_z + (newDimY_x - 2) * (NZ + 1) + NZ] - exact_solution.value_x(0.5, NY, NZ, t)) *
                          (grid.u[(newDimY_x)*dim_z + (newDimY_x - 2) * (NZ + 1) + NZ] - exact_solution.value_x(0.5, NY, NZ, t)) *
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
                error += ((grid.u[i * newDimY_x * (NZ + 1) + dim_z] - exact_solution.value_x(i + 0.5 + offset_x, 0, 0, t)) *
                          (grid.u[i * newDimY_x * (NZ + 1) + dim_z] - exact_solution.value_x(i + 0.5 + offset_x, 0, 0, t)) *
                          DX * DY * DZ / 4);

                for (int k = 1; k < NZ; k++)
                {
                    error += ((grid.u[i * newDimY_x * (NZ + 1) + dim_z + k] - exact_solution.value_x(i + 0.5 + offset_x, 0, k, t)) *
                              (grid.u[i * newDimY_x * (NZ + 1) + dim_z + k] - exact_solution.value_x(i + 0.5 + offset_x, 0, k, t)) *
                              DX * DY * DZ / 2);
                }

                error += ((grid.u[i * newDimY_x * (NZ + 1) + dim_z + NZ] - exact_solution.value_x(i + 0.5 + offset_x, 0, NZ, t)) *
                          (grid.u[i * newDimY_x * (NZ + 1) + dim_z + NZ] - exact_solution.value_x(i + 0.5 + offset_x, 0, NZ, t)) *
                          DX * DY * DZ / 4);
            }
            for (int j = 1 + lby; j < newDimY_x - 1 - rby; j++)
            {
                error += ((grid.u[i * newDimY_x * (NZ + 1) + j * (NZ + 1)] - exact_solution.value_x(i + 0.5 + offset_x, j + offset_y, 0, t)) *
                          (grid.u[i * newDimY_x * (NZ + 1) + j * (NZ + 1)] - exact_solution.value_x(i + 0.5 + offset_x, j + offset_y, 0, t)) *
                          DX * DY * DZ / 2);

                for (int k = 1; k < NZ; k++)
                {
                    error += ((grid.u[i * newDimY_x * (NZ + 1) + j * (NZ + 1) + k] - exact_solution.value_x(i + 0.5 + offset_x, j + offset_y, k, t)) *
                              (grid.u[i * newDimY_x * (NZ + 1) + j * (NZ + 1) + k] - exact_solution.value_x(i + 0.5 + offset_x, j + offset_y, k, t)) *
                              DX * DY * DZ);
                }

                error += ((grid.u[i * newDimY_x * (NZ + 1) + j * (NZ + 1) + NZ] - exact_solution.value_x(i + 0.5 + offset_x, j + offset_y, NZ, t)) *
                          (grid.u[i * newDimY_x * (NZ + 1) + j * (NZ + 1) + NZ] - exact_solution.value_x(i + 0.5 + offset_x, j + offset_y, NZ, t)) *
                          DX * DY * DZ / 2);
            }
            if (rby)
            {
                error += ((grid.u[i * newDimY_x * (NZ + 1) + (newDimY_x - 2) * (NZ + 1)] - exact_solution.value_x(i + 0.5 + offset_x, NY, 0, t)) *
                          (grid.u[i * newDimY_x * (NZ + 1) + (newDimY_x - 2) * (NZ + 1)] - exact_solution.value_x(i + 0.5 + offset_x, NY, 0, t)) *
                          DX * DY * DZ / 4);

                for (int k = 1; k < NZ; k++)
                {
                    error += ((grid.u[i * newDimY_x * (NZ + 1) + (newDimY_x - 2) * (NZ + 1) + k] - exact_solution.value_x(i + 0.5 + offset_x, NY, k, t)) *
                              (grid.u[i * newDimY_x * (NZ + 1) + (newDimY_x - 2) * (NZ + 1) + k] - exact_solution.value_x(i + 0.5 + offset_x, NY, k, t)) *
                              DX * DY * DZ / 2);
                }

                error += ((grid.u[i * newDimY_x * (NZ + 1) + (newDimY_x - 2) * (NZ + 1) + NZ] - exact_solution.value_x(i + 0.5 + offset_x, NY, NZ, t)) *
                          (grid.u[i * newDimY_x * (NZ + 1) + (newDimY_x - 2) * (NZ + 1) + NZ] - exact_solution.value_x(i + 0.5 + offset_x, NY, NZ, t)) *
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
                error += ((grid.u[(newDimX_x - 2) * newDimY_x * (NZ + 1) + dim_z] - exact_solution.value_x((NX - 0.5), 0, 0, t)) *
                          (grid.u[(newDimX_x - 2) * newDimY_x * (NZ + 1) + dim_z] - exact_solution.value_x((NX - 0.5), 0, 0, t)) *
                          DX * DY * DZ / 8);

                for (int k = 1; k < NZ; k++)
                {
                    error += ((grid.u[(newDimX_x - 2) * newDimY_x * (NZ + 1) + dim_z + k] - exact_solution.value_x((NX - 0.5), 0, k, t)) *
                              (grid.u[(newDimX_x - 2) * newDimY_x * (NZ + 1) + dim_z + k] - exact_solution.value_x((NX - 0.5), 0, k, t)) *
                              DX * DY * DZ / 4);
                }

                error += ((grid.u[(newDimX_x - 2) * newDimY_x * (NZ + 1) + dim_z + NZ] - exact_solution.value_x((NX - 0.5), 0, NZ, t)) *
                          (grid.u[(newDimX_x - 2) * newDimY_x * (NZ + 1) + dim_z + NZ] - exact_solution.value_x((NX - 0.5), 0, NZ, t)) *
                          DX * DY * DZ / 8);
            }
            for (int j = 1 + lby; j < newDimY_x - 1 - rby; j++)
            {
                error += ((grid.u[(newDimX_x - 2) * newDimY_x * (NZ + 1) + j * (NZ + 1)] - exact_solution.value_x((NX - 0.5), j + offset_y, 0, t)) *
                          (grid.u[(newDimX_x - 2) * newDimY_x * (NZ + 1) + j * (NZ + 1)] - exact_solution.value_x((NX - 0.5), j + offset_y, 0, t)) *
                          DX * DY * DZ / 4);
                for (int k = 1; k < NZ; k++)
                {
                    error += ((grid.u[(newDimX_x - 2) * newDimY_x * (NZ + 1) + j * (NZ + 1) + k] - exact_solution.value_x((NX - 0.5), j + offset_y, k, t)) *
                              (grid.u[(newDimX_x - 2) * newDimY_x * (NZ + 1) + j * (NZ + 1) + k] - exact_solution.value_x((NX - 0.5), j + offset_y, k, t)) *
                              DX * DY * DZ / 2);
                }
                error += ((grid.u[(newDimX_x - 2) * newDimY_x * (NZ + 1) + j * (NZ + 1) + NZ] - exact_solution.value_x((NX - 0.5), j + offset_y, NZ, t)) *
                          (grid.u[(newDimX_x - 2) * newDimY_x * (NZ + 1) + j * (NZ + 1) + NZ] - exact_solution.value_x((NX - 0.5), j + offset_y, NZ, t)) *
                          DX * DY * DZ / 4);
            }
            if (rby)
            {
                error += ((grid.u[(newDimX_x - 2) * newDimY_x * (NZ + 1) + (newDimY_x - 2) * (NZ + 1)] - exact_solution.value_x((NX - 0.5), NY, 0, t)) *
                          (grid.u[(newDimX_x - 2) * newDimY_x * (NZ + 1) + (newDimY_x - 2) * (NZ + 1)] - exact_solution.value_x((NX - 0.5), NY, 0, t)) *
                          DX * DY * DZ / 8);

                for (int k = 1; k < NZ; k++)
                {
                    error += ((grid.u[(newDimX_x - 2) * newDimY_x * (NZ + 1) + (newDimY_x - 2) * (NZ + 1) + k] - exact_solution.value_x((NX - 0.5), NY, k, t)) *
                              (grid.u[(newDimX_x - 2) * newDimY_x * (NZ + 1) + (newDimY_x - 2) * (NZ + 1) + k] - exact_solution.value_x((NX - 0.5), NY, k, t)) *
                              DX * DY * DZ / 4);
                }

                error += ((grid.u[(newDimX_x - 2) * newDimY_x * (NZ + 1) + (newDimY_x - 2) * (NZ + 1) + NZ] - exact_solution.value_x((NX - 0.5), NY, NZ, t)) *
                          (grid.u[(newDimX_x - 2) * newDimY_x * (NZ + 1) + (newDimY_x - 2) * (NZ + 1) + NZ] - exact_solution.value_x((NX - 0.5), NY, NZ, t)) *
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
            error += ((grid.v[(newDimY_y+1)*dim_z] - exact_solution.value_y(0, 0.5, 0, t)) *
                    (grid.v[(newDimY_y+1)*dim_z] - exact_solution.value_y(0, 0.5, 0, t)) *
                    DX * DY * DZ / 8);

            for (int k = 1; k < NZ; k++)
            {
                error += ((grid.v[(newDimY_y+1)*dim_z + k] - exact_solution.value_y(0, 0.5, k, t)) *
                        (grid.v[(newDimY_y+1)*dim_z + k] - exact_solution.value_y(0, 0.5, k, t)) *
                        DX * DY * DZ / 4);
            }

            error += ((grid.v[(newDimY_y+1)*dim_z + NZ] - exact_solution.value_y(0, 0.5, NZ, t)) *
                    (grid.v[(newDimY_y+1)*dim_z+ NZ] - exact_solution.value_y(0, 0.5, NZ, t)) *
                    DX * DY * DZ / 8);
        }
        for (int j = 1 + lby; j < newDimY_y - 1 -rby; j++)
        {
            error += ((grid.v[(newDimY_y)*dim_z + j * (NZ + 1)] - exact_solution.value_y(0, j + 0.5 + offset_y, 0, t)) *
                      (grid.v[(newDimY_y)*dim_z + j * (NZ + 1)] - exact_solution.value_y(0, j + 0.5 + offset_y, 0, t)) *
                      DX * DY * DZ / 4);
            for (int k = 1; k < NZ; k++)
            {
                error += ((grid.v[(newDimY_y)*dim_z + j * (NZ + 1) + k] - exact_solution.value_y(0, j + 0.5 + offset_y, k, t)) *
                          (grid.v[(newDimY_y)*dim_z + j * (NZ + 1) + k] - exact_solution.value_y(0, j + 0.5 + offset_y, k, t)) *
                          DX * DY * DZ / 2);
            }
            error += ((grid.v[(newDimY_y)*dim_z + j * (NZ + 1) + NZ] - exact_solution.value_y(0, j + 0.5 + offset_y, NZ, t)) *
                      (grid.v[(newDimY_y)*dim_z + j * (NZ + 1) + NZ] - exact_solution.value_y(0, j + 0.5 + offset_y, NZ, t)) *
                      DX * DY * DZ / 4);
        }
        if(rby){
            error += ((grid.v[(newDimY_y)*dim_z + (newDimY_y-2) * (NZ + 1)] - exact_solution.value_y(0, (NY - 0.5), 0, t)) *
                    (grid.v[(newDimY_y)*dim_z + (newDimY_y-2) * (NZ + 1)] - exact_solution.value_y(0, (NY - 0.5), 0, t)) *
                    DX * DY * DZ / 8);

            for (int k = 1; k < NZ; k++)
            {
                error += ((grid.v[(newDimY_y)*dim_z + (newDimY_y-2) * (NZ + 1) + k] - exact_solution.value_y(0, (NY - 0.5), k, t)) *
                        (grid.v[(newDimY_y)*dim_z + (newDimY_y-2) * (NZ + 1) + k] - exact_solution.value_y(0, (NY - 0.5), k, t)) *
                        DX * DY * DZ / 4);
            }

            error += ((grid.v[(newDimY_y)*dim_z + (newDimY_y-2) * (NZ + 1) + NZ] - exact_solution.value_y(0, (NY - 0.5), NZ, t)) *
                    (grid.v[(newDimY_y)*dim_z + (newDimY_y-2) * (NZ + 1) + NZ] - exact_solution.value_y(0, (NY - 0.5), NZ, t)) *
                    DX * DY * DZ / 8);
        }
    }

    // middle slices
    {
        for (int i = 1 + lbx; i < newDimX_y -1 - rbx; i++)
        {
            if(lby){
                error += ((grid.v[i * newDimY_y * (NZ + 1) + dim_z] - exact_solution.value_y(i + offset_x, 0.5, 0, t)) *
                        (grid.v[i * newDimY_y * (NZ + 1) + dim_z] - exact_solution.value_y(i + offset_x, 0.5, 0, t)) *
                        DX * DY * DZ / 4);

                for (int k = 1; k < NZ; k++)
                {
                    error += ((grid.v[i * newDimY_y * (NZ + 1) + dim_z + k] - exact_solution.value_y(i + offset_x, 0.5, k, t)) *
                            (grid.v[i * newDimY_y * (NZ + 1) + dim_z + k] - exact_solution.value_y(i + offset_x, 0.5, k, t)) *
                            DX * DY * DZ / 2);
                }
                error += ((grid.v[i * newDimY_y * (NZ + 1) + dim_z + NZ] - exact_solution.value_y(i + offset_x, 0.5, NZ, t)) *
                        (grid.v[i * newDimY_y * (NZ + 1) + dim_z + NZ] - exact_solution.value_y(i + offset_x, 0.5, NZ, t)) *
                        DX * DY * DZ / 4);
            }
            for (int j = 1 + lby; j < newDimY_y - 1 - rby; j++)
            {
                error += ((grid.v[i * newDimY_y * (NZ + 1) + j * (NZ + 1)] - exact_solution.value_y(i + offset_x, j + 0.5 + offset_y, 0, t)) *
                          (grid.v[i * newDimY_y * (NZ + 1) + j * (NZ + 1)] - exact_solution.value_y(i + offset_x, j + 0.5 + offset_y, 0, t)) *
                          DX * DY * DZ / 2);

                for (int k = 1; k < NZ; k++)
                {
                    error += ((grid.v[i * newDimY_y * (NZ + 1) + j * (NZ + 1) + k] - exact_solution.value_y(i + offset_x, j + 0.5 + offset_y, k, t)) *
                              (grid.v[i * newDimY_y * (NZ + 1) + j * (NZ + 1) + k] - exact_solution.value_y(i + offset_x, j + 0.5 + offset_y, k, t)) *
                              DX * DY * DZ);
                }

                error += ((grid.v[i * newDimY_y * (NZ + 1) + j * (NZ + 1) + NZ] - exact_solution.value_y(i + offset_x, j + 0.5 + offset_y, NZ, t)) *
                          (grid.v[i * newDimY_y * (NZ + 1) + j * (NZ + 1) + NZ] - exact_solution.value_y(i + offset_x, j + 0.5 + offset_y, NZ, t)) *
                          DX * DY * DZ / 2);
            }
            if(rby){
                error += ((grid.v[i * newDimY_y * (NZ + 1) + (newDimY_y-2) * (NZ + 1)] - exact_solution.value_y(i + offset_x, (NY - 0.5), 0, t)) *
                        (grid.v[i * newDimY_y * (NZ + 1) + (newDimY_y-2) * (NZ + 1)] - exact_solution.value_y(i + offset_x, (NY - 0.5), 0, t)) *
                        DX * DY * DZ / 4);

                for (int k = 1; k < NZ; k++)
                {
                    error += ((grid.v[i * newDimY_y * (NZ + 1) + (newDimY_y-2) * (NZ + 1) + k] - exact_solution.value_y(i + offset_x, (NY - 0.5), k, t)) *
                            (grid.v[i * newDimY_y * (NZ + 1) + (newDimY_y-2) * (NZ + 1) + k] - exact_solution.value_y(i + offset_x, (NY - 0.5), k, t)) *
                            DX * DY * DZ / 2);
                }

                error += ((grid.v[i * newDimY_y * (NZ + 1) + (newDimY_y-2) * (NZ + 1) + NZ] - exact_solution.value_y(i + offset_x, (NY - 0.5), NZ, t)) *
                        (grid.v[i * newDimY_y * (NZ + 1) + (newDimY_y-2) * (NZ + 1) + NZ] - exact_solution.value_y(i + offset_x, (NY - 0.5), NZ, t)) *
                        DX * DY * DZ / 4);
            }
        }
    }

    // last slice (right face)
    if(rbx)
    {
        if(lby){
            error += ((grid.v[(newDimX_y-2) * newDimY_y * (NZ + 1) + dim_z] - exact_solution.value_y(NX, 0.5, 0, t)) *
                    (grid.v[(newDimX_y-2) * newDimY_y * (NZ + 1) + dim_z] - exact_solution.value_y(NX, 0.5, 0, t)) *
                    DX * DY * DZ / 8);

            for (int k = 1; k < NZ; k++)
            {
                error += ((grid.v[(newDimX_y-2) * newDimY_y * (NZ + 1) + dim_z + k] - exact_solution.value_y(NX, 0.5, k, t)) *
                        (grid.v[(newDimX_y-2) * newDimY_y * (NZ + 1) + dim_z + k] - exact_solution.value_y(NX, 0.5, k, t)) *
                        DX * DY * DZ / 4);
            }

            error += ((grid.v[(newDimX_y-2) * newDimY_y * (NZ + 1) + dim_z + NZ] - exact_solution.value_y(NX, 0.5, NZ, t)) *
                    (grid.v[(newDimX_y-2) * newDimY_y * (NZ + 1) + dim_z + NZ] - exact_solution.value_y(NX, 0.5, NZ, t)) *
                    DX * DY * DZ / 8);
        }
        for (int j = 1 + lby; j < newDimY_y - 1 - rby; j++)
        {
            error += ((grid.v[(newDimX_y-2) * newDimY_y * (NZ + 1) + j * (NZ + 1)] - exact_solution.value_y(NX, j + 0.5 + offset_y, 0, t)) *
                      (grid.v[(newDimX_y-2) * newDimY_y * (NZ + 1) + j * (NZ + 1)] - exact_solution.value_y(NX, j + 0.5 + offset_y, 0, t)) *
                      DX * DY * DZ / 4);
            for (int k = 1; k < NZ; k++)
            {
                error += ((grid.v[(newDimX_y-2) * newDimY_y * (NZ + 1) + j * (NZ + 1) + k] - exact_solution.value_y(NX, j + 0.5 + offset_y, k, t)) *
                          (grid.v[(newDimX_y-2) * newDimY_y * (NZ + 1) + j * (NZ + 1) + k] - exact_solution.value_y(NX, j + 0.5 + offset_y, k, t)) *
                          DX * DY * DZ / 2);
            }
            error += ((grid.v[(newDimX_y-2) * newDimY_y * (NZ + 1) + j * (NZ + 1) + NZ] - exact_solution.value_y(NX, j + 0.5 + offset_y, NZ, t)) *
                      (grid.v[(newDimX_y-2) * newDimY_y * (NZ + 1) + j * (NZ + 1) + NZ] - exact_solution.value_y(NX, j + 0.5 + offset_y, NZ, t)) *
                      DX * DY * DZ / 4);
        }
        if(rby){
            error += ((grid.v[(newDimX_y-2) * newDimY_y * (NZ + 1) + (newDimY_y-2) * (NZ + 1)] - exact_solution.value_y(NX, (NY - 0.5), 0, t)) *
                    (grid.v[(newDimX_y-2) * newDimY_y * (NZ + 1) + (newDimY_y-2) * (NZ + 1)] - exact_solution.value_y(NX, (NY - 0.5), 0, t)) *
                    DX * DY * DZ / 8);

            for (int k = 1; k < NZ; k++)
            {
                error += ((grid.v[(newDimX_y-2) * newDimY_y * (NZ + 1) + (newDimY_y-2) * (NZ + 1) + k] - exact_solution.value_y(NX, (NY - 0.5), k, t)) *
                        (grid.v[(newDimX_y-2) * newDimY_y * (NZ + 1) + (newDimY_y-2) * (NZ + 1) + k] - exact_solution.value_y(NX, (NY - 0.5), k, t)) *
                        DX * DY * DZ / 4);
            }

            error += ((grid.v[(newDimX_y-2) * newDimY_y * (NZ + 1) + (newDimY_y-2) * (NZ + 1) + NZ] - exact_solution.value_y(NX, (NY - 0.5), NZ, t)) *
                    (grid.v[(newDimX_y-2) * newDimY_y * (NZ + 1) + (newDimY_y-2) * (NZ + 1) + NZ] - exact_solution.value_y(NX, (NY - 0.5), NZ, t)) *
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
            error += ((grid.w[(newDimY_z +1) *dim_z_z] - exact_solution.value_z(0, 0, 0.5, t)) *
                    (grid.w[(newDimY_z +1) *dim_z_z] - exact_solution.value_z(0, 0, 0.5, t)) *
                    DX * DY * DZ / 8);

            for (int k = 1; k < NZ - 1; k++)
            {
                error += ((grid.w[(newDimY_z +1) *dim_z_z + k] - exact_solution.value_z(0, 0, k + 0.5, t)) *
                        (grid.w[(newDimY_z +1) *dim_z_z + k] - exact_solution.value_z(0, 0, k + 0.5, t)) *
                        DX * DY * DZ / 4);
            }

            error += ((grid.w[(newDimY_z +1) *dim_z_z + NZ - 1] - exact_solution.value_z(0, 0, NZ - 0.5, t)) *
                    (grid.w[(newDimY_z +1) *dim_z_z + NZ - 1] - exact_solution.value_z(0, 0, NZ - 0.5, t)) *
                    DX * DY * DZ / 8);
        }
        for (int j = 1 + lby; j < newDimY_z-1-rby; j++)
        {
            error += ((grid.w[(newDimY_z) *dim_z_z + j * NZ] - exact_solution.value_z(0, j + offset_y, 0.5, t)) *
                      (grid.w[(newDimY_z) *dim_z_z + j * NZ] - exact_solution.value_z(0, j + offset_y, 0.5, t)) *
                      DX * DY * DZ / 4);
            for (int k = 1; k < NZ - 1; k++)
            {
                error += ((grid.w[(newDimY_z) *dim_z_z + j * NZ + k] - exact_solution.value_z(0, j + offset_y, k + 0.5, t)) *
                          (grid.w[(newDimY_z) *dim_z_z + j * NZ + k] - exact_solution.value_z(0, j + offset_y, k + 0.5, t)) *
                          DX * DY * DZ / 2);
            }
            error += ((grid.w[(newDimY_z) *dim_z_z + j * NZ + NZ - 1] - exact_solution.value_z(0, j + offset_y, NZ - 0.5, t)) *
                      (grid.w[(newDimY_z) *dim_z_z + j * NZ + NZ - 1] - exact_solution.value_z(0, j + offset_y, NZ - 0.5, t)) *
                      DX * DY * DZ / 4);
        }
        if(rby){
            error += ((grid.w[(newDimY_z) *dim_z_z + (newDimY_z-2)*dim_z_z] - exact_solution.value_z(0, NY, 0.5, t)) *
                    (grid.w[(newDimY_z) *dim_z_z + (newDimY_z-2)*dim_z_z] - exact_solution.value_z(0, NY, 0.5, t)) *
                    DX * DY * DZ / 8);

            for (int k = 1; k < NZ - 1; k++)
            {
                error += ((grid.w[(newDimY_z) *dim_z_z + (newDimY_z-2)*dim_z_z + k] - exact_solution.value_z(0, NY, k + 0.5, t)) *
                        (grid.w[(newDimY_z) *dim_z_z + (newDimY_z-2)*dim_z_z + k] - exact_solution.value_z(0, NY, k + 0.5, t)) *
                        DX * DY * DZ / 4);
            }

            error += ((grid.w[(newDimY_z) *dim_z_z + (newDimY_z-2)*dim_z_z + NZ - 1] - exact_solution.value_z(0, NY, NZ - 0.5, t)) *
                    (grid.w[(newDimY_z) *dim_z_z + (newDimY_z-2)*dim_z_z + NZ - 1] - exact_solution.value_z(0, NY, NZ - 0.5, t)) *
                    DX * DY * DZ / 8);
        }
    }

    // middle slices
    {
        for (int i = 1 + lbx; i < newDimX_z - 1 -rbx; i++)
        {
            if(lby){
                error += ((grid.w[i * newDimY_z * NZ  + dim_z_z] - exact_solution.value_z(i + offset_x, 0, 0.5, t)) *
                        (grid.w[i * newDimY_z * NZ  + dim_z_z] - exact_solution.value_z(i + offset_x, 0, 0.5, t)) *
                        DX * DY * DZ / 4);

                for (int k = 1; k < NZ - 1; k++)
                {
                    error += ((grid.w[i * newDimY_z * NZ  + dim_z_z + k] - exact_solution.value_z(i + offset_x, 0, k + 0.5, t)) *
                            (grid.w[i * newDimY_z * NZ  + dim_z_z + k] - exact_solution.value_z(i + offset_x, 0, k + 0.5, t)) *
                            DX * DY * DZ / 2);
                }
                error += ((grid.w[i * newDimY_z * NZ  + dim_z_z + NZ - 1] - exact_solution.value_z(i + offset_x, 0, NZ - 0.5, t)) *
                        (grid.w[i * newDimY_z * NZ  + dim_z_z + NZ - 1] - exact_solution.value_z(i + offset_x, 0, NZ - 0.5, t)) *
                        DX * DY * DZ / 4);
            }
            for (int j = 1 + lby; j < newDimY_z - 1 -rby; j++)
            {
                error += ((grid.w[i * newDimY_z * NZ + j * NZ] - exact_solution.value_z(i + offset_x, j + offset_y, 0.5, t)) *
                          (grid.w[i * newDimY_z * NZ + j * NZ] - exact_solution.value_z(i + offset_x, j + offset_y, 0.5, t)) *
                          DX * DY * DZ / 2);

                for (int k = 1; k < NZ - 1; k++)
                {
                    error += ((grid.w[i * newDimY_z * NZ + j * NZ + k] - exact_solution.value_z(i + offset_x, j + offset_y, k + 0.5, t)) *
                              (grid.w[i * newDimY_z * NZ + j * NZ + k] - exact_solution.value_z(i + offset_x, j + offset_y, k + 0.5, t)) *
                              DX * DY * DZ);
                }

                error += ((grid.w[i * newDimY_z * NZ + j * NZ + NZ - 1] - exact_solution.value_z(i + offset_x, j + offset_y, (NZ - 0.5), t)) *
                          (grid.w[i * newDimY_z * NZ + j * NZ + NZ - 1] - exact_solution.value_z(i + offset_x, j + offset_y, (NZ - 0.5), t)) *
                          DX * DY * DZ / 2);
            }
            if(rby){
                error += ((grid.w[i * newDimY_z * NZ + (newDimY_z-2) * NZ] - exact_solution.value_z(i + offset_x, NY, 0.5, t)) *
                        (grid.w[i * newDimY_z * NZ + (newDimY_z-2) * NZ] - exact_solution.value_z(i + offset_x, NY, 0.5, t)) *
                        DX * DY * DZ / 4);

                for (int k = 1; k < NZ - 1; k++)
                {
                    error += ((grid.w[i * newDimY_z * NZ + (newDimY_z-2) * NZ + k] - exact_solution.value_z(i + offset_x, NY, k + 0.5, t)) *
                            (grid.w[i * newDimY_z * NZ + (newDimY_z-2) * NZ + k] - exact_solution.value_z(i + offset_x, NY, k + 0.5, t)) *
                            DX * DY * DZ / 2);
                }

                error += ((grid.w[i * newDimY_z * NZ + (newDimY_z-2) * NZ + NZ - 1] - exact_solution.value_z(i + offset_x, NY, NZ - 0.5, t)) *
                        (grid.w[i * newDimY_z * NZ + (newDimY_z-2) * NZ + NZ - 1] - exact_solution.value_z(i + offset_x, NY, NZ - 0.5, t)) *
                        DX * DY * DZ / 4);
            }
        }
    }

    // last slice (right face)
    if(rbx)
    {
        if(lby){
            error += ((grid.w[(newDimX_z-2) * (newDimY_z) * NZ + dim_z_z] - exact_solution.value_z(NX, 0, 0.5, t)) *
                    (grid.w[(newDimX_z-2) * (newDimY_z) * NZ + dim_z_z] - exact_solution.value_z(NX, 0, 0.5, t)) *
                    DX * DY * DZ / 8);

            for (int k = 1; k < NZ - 1; k++)
            {
                error += ((grid.w[(newDimX_z-2) * (newDimY_z) * NZ + dim_z_z + k] - exact_solution.value_z(NX, 0, k + 0.5, t)) *
                        (grid.w[(newDimX_z-2) * (newDimY_z) * NZ + dim_z_z + k] - exact_solution.value_z(NX, 0, k + 0.5, t)) *
                        DX * DY * DZ / 4);
            }

            error += ((grid.w[(newDimX_z-2) * (newDimY_z) * NZ + dim_z_z + NZ - 1] - exact_solution.value_z(NX, 0, NZ - 0.5, t)) *
                    (grid.w[(newDimX_z-2) * (newDimY_z) * NZ + dim_z_z + NZ - 1] - exact_solution.value_z(NX, 0, NZ - 0.5, t)) *
                    DX * DY * DZ / 8);
        }
        for (int j = 1 + lby; j < newDimY_z - 1 - rby; j++)
        {
            error += ((grid.w[(newDimX_z-2) * (newDimY_z) * NZ + j * NZ] - exact_solution.value_z(NX, j + offset_y, 0.5, t)) *
                      (grid.w[(newDimX_z-2) * (newDimY_z) * NZ + j * NZ] - exact_solution.value_z(NX, j + offset_y, 0.5, t)) *
                      DX * DY * DZ / 4);
            for (int k = 1; k < NZ - 1; k++)
            {
                error += ((grid.w[(newDimX_z-2) * (newDimY_z) * NZ + j * NZ + k] - exact_solution.value_z(NX, j + offset_y, k + 0.5, t)) *
                          (grid.w[(newDimX_z-2) * (newDimY_z) * NZ + j * NZ + k] - exact_solution.value_z(NX, j + offset_y, k + 0.5, t)) *
                          DX * DY * DZ / 2);
            }
            error += ((grid.w[(newDimX_z-2) * (newDimY_z) * NZ + j * NZ + NZ - 1] - exact_solution.value_z(NX, j + offset_y, NZ - 0.5, t)) *
                      (grid.w[(newDimX_z-2) * (newDimY_z) * NZ + j * NZ + NZ - 1] - exact_solution.value_z(NX, j + offset_y, NZ - 0.5, t)) *
                      DX * DY * DZ / 4);
        }
        if(rby){
            error += ((grid.w[(newDimX_z-2) * (newDimY_z) * NZ  + (newDimY_z-2) * NZ] - exact_solution.value_z(NX, NY, 0.5, t)) *
                    (grid.w[(newDimX_z-2) * (newDimY_z) * NZ  + (newDimY_z-2) * NZ] - exact_solution.value_z(NX, NY, 0.5, t)) *
                    DX * DY * DZ / 8);

            for (int k = 1; k < NZ - 1; k++)
            {
                error += ((grid.w[(newDimX_z-2) * (newDimY_z) * NZ  + (newDimY_z-2) * NZ + k] - exact_solution.value_z(NX, NY, k + 0.5, t)) *
                        (grid.w[(newDimX_z-2) * (newDimY_z) * NZ  + (newDimY_z-2) * NZ + k] - exact_solution.value_z(NX, NY, k + 0.5, t)) *
                        DX * DY * DZ / 4);
            }

            error += ((grid.w[(newDimX_z-2) * (newDimY_z) * NZ  + (newDimY_z-2) * NZ + NZ - 1] - exact_solution.value_z(NX, NY, NZ - 0.5, t)) *
                    (grid.w[(newDimX_z-2) * (newDimY_z) * NZ  + (newDimY_z-2) * NZ + NZ - 1] - exact_solution.value_z(NX, NY, NZ - 0.5, t)) *
                    DX * DY * DZ / 8);
        }
    }

    return error;
}

//TODO: parallelization
Real IcoNS::error_comp_P(const Real t)
{
    Real error = 0.0;
    // first slice (left face)
    {
        error += ((grid.p[0] - exact_solution.value_p(0, 0, 0, t)) *
                  (grid.p[0] - exact_solution.value_p(0, 0, 0, t)) *
                  DX * DY * DZ / 8);

        for (int k = 1; k < NZ; k++)
        {
            error += ((grid.p[k] - exact_solution.value_p(0, 0, k, t)) *
                      (grid.p[k] - exact_solution.value_p(0, 0, k, t)) *
                      DX * DY * DZ / 4);
        }

        error += ((grid.p[NZ] - exact_solution.value_p(0, 0, NZ, t)) *
                  (grid.p[NZ] - exact_solution.value_p(0, 0, NZ, t)) *
                  DX * DY * DZ / 8);

        for (int j = 1; j < NY; j++)
        {
            error += ((grid.p[j * (NZ + 1)] - exact_solution.value_p(0, j, 0, t)) *
                      (grid.p[j * (NZ + 1)] - exact_solution.value_p(0, j, 0, t)) *
                      DX * DY * DZ / 4);
            for (int k = 1; k < NZ; k++)
            {
                error += ((grid.p[j * (NZ + 1) + k] - exact_solution.value_p(0, j, k, t)) *
                          (grid.p[j * (NZ + 1) + k] - exact_solution.value_p(0, j, k, t)) *
                          DX * DY * DZ / 2);
            }
            error += ((grid.p[j * (NZ + 1) + NZ] - exact_solution.value_p(0, j, NZ, t)) *
                      (grid.p[j * (NZ + 1) + NZ] - exact_solution.value_p(0, j, NZ, t)) *
                      DX * DY * DZ / 4);
        }

        error += ((grid.p[NY * (NZ + 1)] - exact_solution.value_p(0, NY, 0, t)) *
                  (grid.p[NY * (NZ + 1)] - exact_solution.value_p(0, NY, 0, t)) *
                  DX * DY * DZ / 8);

        for (int k = 1; k < NZ; k++)
        {
            error += ((grid.p[NY * (NZ + 1) + k] - exact_solution.value_p(0, NY, k, t)) *
                      (grid.p[NY * (NZ + 1) + k] - exact_solution.value_p(0, NY, k, t)) *
                      DX * DY * DZ / 4);
        }

        error += ((grid.p[NY * (NZ + 1) + NZ] - exact_solution.value_p(0, NY, NZ, t)) *
                  (grid.p[NY * (NZ + 1) + NZ] - exact_solution.value_p(0, NY, NZ, t)) *
                  DX * DY * DZ / 8);
    }

    // middle slices
    {
        for (int i = 1; i < NX; i++)
        {
            error += ((grid.p[i * (NY + 1) * (NZ + 1)] - exact_solution.value_p(i, 0, 0, t)) *
                      (grid.p[i * (NY + 1) * (NZ + 1)] - exact_solution.value_p(i, 0, 0, t)) *
                      DX * DY * DZ / 4);

            for (int k = 1; k < NZ; k++)
            {
                error += ((grid.p[i * (NY + 1) * (NZ + 1) + k] - exact_solution.value_p(i, 0, k, t)) *
                          (grid.p[i * (NY + 1) * (NZ + 1) + k] - exact_solution.value_p(i, 0, k, t)) *
                          DX * DY * DZ / 2);
            }
            error += ((grid.p[i * (NY + 1) * (NZ + 1) + NZ] - exact_solution.value_p(i, 0, NZ, t)) *
                      (grid.p[i * (NY + 1) * (NZ + 1) + NZ] - exact_solution.value_p(i, 0, NZ, t)) *
                      DX * DY * DZ / 4);

            for (int j = 1; j < NY; j++)
            {
                error += ((grid.p[i * (NY + 1) * (NZ + 1) + j * (NZ + 1)] - exact_solution.value_p(i, j, 0, t)) *
                          (grid.p[i * (NY + 1) * (NZ + 1) + j * (NZ + 1)] - exact_solution.value_p(i, j, 0, t)) *
                          DX * DY * DZ / 2);

                for (int k = 1; k < NZ; k++)
                {
                    error += ((grid.p[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k] - exact_solution.value_p(i, j, k, t)) *
                              (grid.p[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k] - exact_solution.value_p(i, j, k, t)) *
                              DX * DY * DZ);
                }

                error += ((grid.p[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + NZ] - exact_solution.value_p(i, j, NZ, t)) *
                          (grid.p[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + NZ] - exact_solution.value_p(i, j, NZ, t)) *
                          DX * DY * DZ / 2);
            }

            error += ((grid.p[i * (NY + 1) * (NZ + 1) + NY * (NZ + 1)] - exact_solution.value_p(i, NY, 0, t)) *
                      (grid.p[i * (NY + 1) * (NZ + 1) + NY * (NZ + 1)] - exact_solution.value_p(i, NY, 0, t)) *
                      DX * DY * DZ / 4);

            for (int k = 1; k < NZ; k++)
            {
                error += ((grid.p[i * (NY + 1) * (NZ + 1) + NY * (NZ + 1) + k] - exact_solution.value_p(i, NY, k, t)) *
                          (grid.p[i * (NY + 1) * (NZ + 1) + NY * (NZ + 1) + k] - exact_solution.value_p(i, NY, k, t)) *
                          DX * DY * DZ / 2);
            }

            error += ((grid.p[i * (NY + 1) * (NZ + 1) + NY * (NZ + 1) + NZ] - exact_solution.value_p(i, NY, NZ, t)) *
                      (grid.p[i * (NY + 1) * (NZ + 1) + NY * (NZ + 1) + NZ] - exact_solution.value_p(i, NY, NZ, t)) *
                      DX * DY * DZ / 4);
        }
    }

    // last slice (right face)
    {
        error += ((grid.p[NX * (NY + 1) * (NZ + 1)] - exact_solution.value_p(NX, 0, 0, t)) *
                  (grid.p[NX * (NY + 1) * (NZ + 1)] - exact_solution.value_p(NX, 0, 0, t)) *
                  DX * DY * DZ / 8);

        for (int k = 1; k < NZ; k++)
        {
            error += ((grid.p[NX * (NY + 1) * (NZ + 1) + k] - exact_solution.value_p(NX, 0, k, t)) *
                      (grid.p[NX * (NY + 1) * (NZ + 1) + k] - exact_solution.value_p(NX, 0, k, t)) *
                      DX * DY * DZ / 4);
        }

        error += ((grid.p[NX * (NY + 1) * (NZ + 1) + NZ] - exact_solution.value_p(NX, 0, NZ, t)) *
                  (grid.p[NX * (NY + 1) * (NZ + 1) + NZ] - exact_solution.value_p(NX, 0, NZ, t)) *
                  DX * DY * DZ / 8);

        for (int j = 1; j < NY; j++)
        {
            error += ((grid.p[NX * (NY + 1) * (NZ + 1) + j * (NZ + 1)] - exact_solution.value_p(NX, j, 0, t)) *
                      (grid.p[NX * (NY + 1) * (NZ + 1) + j * (NZ + 1)] - exact_solution.value_p(NX, j, 0, t)) *
                      DX * DY * DZ / 4);
            for (int k = 1; k < NZ; k++)
            {
                error += ((grid.p[NX * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k] - exact_solution.value_p(NX, j, k, t)) *
                          (grid.p[NX * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k] - exact_solution.value_p(NX, j, k, t)) *
                          DX * DY * DZ / 2);
            }
            error += ((grid.p[NX * (NY + 1) * (NZ + 1) + j * (NZ + 1) + NZ] - exact_solution.value_p(NX, j, NZ, t)) *
                      (grid.p[NX * (NY + 1) * (NZ + 1) + j * (NZ + 1) + NZ] - exact_solution.value_p(NX, j, NZ, t)) *
                      DX * DY * DZ / 4);
        }

        error += ((grid.p[NX * (NY + 1) * (NZ + 1) + NY * (NZ + 1)] - exact_solution.value_p(NX, NY, 0, t)) *
                  (grid.p[NX * (NY + 1) * (NZ + 1) + NY * (NZ + 1)] - exact_solution.value_p(NX, NY, 0, t)) *
                  DX * DY * DZ / 8);

        for (int k = 1; k < NZ; k++)
        {
            error += ((grid.p[NX * (NY + 1) * (NZ + 1) + NY * (NZ + 1) + k] - exact_solution.value_p(NX, NY, k, t)) *
                      (grid.p[NX * (NY + 1) * (NZ + 1) + NY * (NZ + 1) + k] - exact_solution.value_p(NX, NY, k, t)) *
                      DX * DY * DZ / 4);
        }

        error += ((grid.p[NX * (NY + 1) * (NZ + 1) + NY * (NZ + 1) + NZ] - exact_solution.value_p(NX, NY, NZ, t)) *
                  (grid.p[NX * (NY + 1) * (NZ + 1) + NY * (NZ + 1) + NZ] - exact_solution.value_p(NX, NY, NZ, t)) *
                  DX * DY * DZ / 8);
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

    // Skip comments and empty lines until we find test case number
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream iss(line);
        if (!(iss >> testCase)) continue;
        break;
    }

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

    // Skip comments and empty lines until we find Nt
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream iss(line);
        if (!(iss >> Nt)) continue;
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

    NX -= 1;
    NY -= 1;
    NZ -= 1;

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
                  << "Number of timesteps: " << Nt << "\n"
                  << "Domain size: " << LX << " x " << LY << " x " << LZ << "\n"
                  << "Grid points: " << NX+1 << " x " << NY+1 << " x " << NZ+1 << "\n"
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

void IcoNS::Output(){

    const std::string filename = "out" + std::to_string(testCase) + ".dat";
    MPI_File fh;
    MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
    MPI_Offset offset = coords[1] * zSize[1] * sizeof(double) * 7;
    
    // LINE 1
    if(coords[0] == (PX-1)/2){
        const double xCoord = /*SX+*/LX/2;
        double yCoord = /*SY+*/coords[1] * zSize[1] * DY;
        const double zCoord = /*SZ+*/LZ/2;
        double u, v, w, p;
        for (int i = 0; i < zSize[1]; i++)
        {
            MPI_File_write_at(fh, offset, &xCoord, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);
            offset += sizeof(double);

            MPI_File_write_at(fh, offset, &yCoord, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);
            yCoord += DY;
            offset += sizeof(double);

            MPI_File_write_at(fh, offset, &zCoord, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);
            offset += sizeof(double);

            if(PX%2 == 0){
                u = grid.u[getx(newDimX_x-1, i, (dim_z-1)/2)];
                v = (grid.v[gety(newDimX_y-1, i, (dim_z-1)/2)] +
                    grid.v[gety(newDimX_y-1, i-1, (dim_z-1)/2)] +
                    grid.v[gety(newDimX_y, i, (dim_z-1)/2)] +
                    grid.v[gety(newDimX_y, i-1, (dim_z-1)/2)])/4;
                w = (grid.w[getz(newDimX_z-1, i, (dim_z_z-1)/2)] +
                    grid.w[getz(newDimX_z-1, i, (dim_z_z-1)/2-1)] +
                    grid.w[getz(newDimX_z, i, (dim_z_z-1)/2)] +
                    grid.w[getz(newDimX_z, i, (dim_z_z-1)/2-1)])/4;
                p = grid.p[getp(zSize[0]-1, i, (zSize[2]-1)/2)];
                    // + grid.p[getp(zSize[0], i, (zSize[2]-1)/2)])/2;
            }else{
                u = grid.u[getx((newDimX_x-1)/2, i, (dim_z-1)/2)];
                v = (grid.v[gety((newDimX_y-1)/2, i, (dim_z-1)/2)] +
                    grid.v[gety((newDimX_y-1)/2, i-1, (dim_z-1)/2)] +
                    grid.v[gety((newDimX_y-1)/2+1, i, (dim_z-1)/2)] +
                    grid.v[gety((newDimX_y-1)/2+1, i-1, (dim_z-1)/2)])/4;
                w = (grid.w[getz((newDimX_z-1)/2, i, (dim_z_z-1)/2)] +
                    grid.w[getz((newDimX_z-1)/2+1, i, (dim_z_z-1)/2)] +
                    grid.w[getz((newDimX_z-1)/2, i, (dim_z_z-1)/2-1)] +
                    grid.w[getz((newDimX_z-1)/2+1, i, (dim_z_z-1)/2-1)])/4;
                p = (grid.p[getp((zSize[0]-1)/2, i, (zSize[2]-1)/2)] +
                    grid.p[getp((zSize[0]-1)/2+1, i, (zSize[2]-1)/2)])/2;
            }
            MPI_File_write_at(fh, offset, &u, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);
            offset += sizeof(double);
            MPI_File_write_at(fh, offset, &v, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);
            offset += sizeof(double);
            MPI_File_write_at(fh, offset, &w, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);
            offset += sizeof(double);
            MPI_File_write_at(fh, offset, &p, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);
            offset += sizeof(double);
        }
    }

    // LINE 2
    if(coords[1] == (PY-1)/2){
        offset = PY * zSize[1] * sizeof(double) * 7 + coords[0] * zSize[0] * sizeof(double) * 7;
        double xCoord = /*SX+*/coords[0] * zSize[0] * DX;
        const double yCoord = /*SY+*/LY/2;
        const double zCoord = /*SZ+*/LZ/2;
        double u, v, w, p;
        for (int i = 0; i < zSize[0]; i++)
        {
            MPI_File_write_at(fh, offset, &xCoord, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);
            xCoord += DX;
            offset += sizeof(double);

            MPI_File_write_at(fh, offset, &yCoord, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);
            offset += sizeof(double);

            MPI_File_write_at(fh, offset, &zCoord, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);
            offset += sizeof(double);

            if(PY%2 == 0){
                u = (grid.u[getx(i, newDimY_x-1, (dim_z-1)/2)] +
                    grid.u[getx(i-1, newDimY_x-1, (dim_z-1)/2)] +
                    grid.u[getx(i, newDimY_x, (dim_z-1)/2)] +
                    grid.u[getx(i-1, newDimY_x, (dim_z-1)/2)])/4;
                v = grid.v[gety(i, newDimY_y-1, (dim_z-1)/2)];
                w = (grid.w[getz(i, newDimY_z-1, (dim_z-1)/2)] +
                    grid.w[getz(i, newDimY_z-1, (dim_z-1)/2-1)] +
                    grid.w[getz(i, newDimY_z, (dim_z-1)/2)] +
                    grid.w[getz(i, newDimY_z, (dim_z-1)/2-1)])/4;
                p = grid.p[getp(i, zSize[1]-1, (zSize[2]-1)/2)];
                    // + grid.p[getp(i, zSize[1], (zSize[2]-1)/2)])/2;
            }else{
                u = (grid.u[getx(i, (newDimY_x-1)/2, (dim_z-1)/2)] +
                    grid.u[getx(i-1, (newDimY_x-1)/2, (dim_z-1)/2)] +
                    grid.u[getx(i, (newDimY_x-1)/2+1, (dim_z-1)/2)] +
                    grid.u[getx(i-1, (newDimY_x-1)/2+1, (dim_z-1)/2)])/4;
                v = grid.v[gety(i, (newDimY_y-1)/2, (dim_z-1)/2)];
                w = (grid.w[getz(i, (newDimY_z-1)/2, (dim_z-1)/2)] +
                    grid.w[getz(i, (newDimY_z-1)/2+1, (dim_z-1)/2)] +
                    grid.w[getz(i, (newDimY_z-1)/2, (dim_z-1)/2)-1] +
                    grid.w[getz(i, (newDimY_z-1)/2+1, (dim_z-1)/2-1)-1])/4;
                p = (grid.p[getp(i, (zSize[1]-1)/2, (zSize[2]-1)/2)] +
                    grid.p[getp(i, (zSize[1]-1)/2+1, (zSize[2]-1)/2)])/2;
            }
            MPI_File_write_at(fh, offset, &u, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);
            offset += sizeof(double);
            MPI_File_write_at(fh, offset, &v, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);
            offset += sizeof(double);
            MPI_File_write_at(fh, offset, &w, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);
            offset += sizeof(double);
            MPI_File_write_at(fh, offset, &p, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);
            offset += sizeof(double);
        }
    }

    // LINE 3
    if(testCase == 2){
        if(coords[0] == (PX-1)/2 && coords[1] == (PY-1)/2)
        {
            offset = (PX * zSize[0] + PY * zSize[1]) * sizeof(double) * 7;
            const double xCoord = /*SX+*/LX/2;
            const double yCoord = /*SY+*/LY/2;
            double zCoord = /*SZ+*/0.0;
            
            double u, v, w, p;
            for(int k = 0; k < zSize[2]; k++)
            {
                MPI_File_write_at(fh, offset, &xCoord, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);
                offset += sizeof(double);

                MPI_File_write_at(fh, offset, &yCoord, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);
                offset += sizeof(double);

                MPI_File_write_at(fh, offset, &zCoord, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);
                zCoord += DZ;
                offset += sizeof(double);

                int iX, iY, iZ, iP, jX, jY, jZ, jP;
                if(PX % 2 == 0)
                {
                    iX = newDimX_x-1;
                    iY = newDimX_y-1;
                    iZ = newDimX_z-1;
                    iP = zSize[0]-1;
                }else
                {
                    iX = (newDimX_x-1)/2;
                    iY = (newDimX_y-1)/2;
                    iZ = (newDimX_z-1)/2;
                    iP = (zSize[0]-1)/2;
                }
                if(PY % 2 == 0)
                {
                    jX = newDimY_x-1;
                    jY = newDimY_y-1;
                    jZ = newDimY_z-1;
                    jP = zSize[1]-1;
                }else
                {
                    jX = (newDimY_x-1)/2;
                    jY = (newDimY_y-1)/2;
                    jZ = (newDimY_z-1)/2;
                    jP = (zSize[1]-1)/2;
                }
                u = (grid.u[getx(iX, jX, k)] +
                    grid.u[getx(iX, jX+1, k)] +
                    grid.u[getx(iX, jX, k+1)] +
                    grid.u[getx(iX, jX+1, k+1)])/4;
                v = (grid.v[gety(iY, jY, k)] +
                    grid.v[gety(iY+1, jY, k)] +
                    grid.v[gety(iY, jY, k+1)] +
                    grid.v[gety(iY+1, jY, k+1)])/4;
                w = grid.w[getz(iZ, jZ, k)] +
                    grid.w[getz(iZ+1, jZ, k)] +
                    grid.w[getz(iZ, jZ+1, k)] +
                    grid.w[getz(iZ+1, jZ+1, k)]/4;
                p = grid.p[getp(iP, jP, k)];
                    // + grid.p[getp(iP+1, jP, k)]
                    // + grid.p[getp(iP, jP+1, k)]
                    // + grid.p[getp(iP, jP, k+1)]
                    // + grid.p[getp(iP+1, jP+1, k)]
                    // + grid.p[getp(iP+1, jP, k+1)]
                    // + grid.p[getp(iP, jP+1, k+1)]
                    // + grid.p[getp(iP+1, jP+1, k+1)])/8;

                MPI_File_write_at(fh, offset, &u, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);
                offset += sizeof(double);
                MPI_File_write_at(fh, offset, &v, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);
                offset += sizeof(double);
                MPI_File_write_at(fh, offset, &w, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);
                offset += sizeof(double);
                MPI_File_write_at(fh, offset, &p, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);
                offset += sizeof(double);
            }
        }
    }

    MPI_File_close(&fh);

    if(rank == 0){
        int count = 0;
        std::ifstream input(filename, std::ios::binary);
        std::ofstream output("profile" + std::to_string(testCase) + ".dat");
        output << "Line 1" << std::endl;
        output << "x y z u v w p" << std::endl;
        double value;
        while (input.read(reinterpret_cast<char*>(&value), sizeof(double))) {
            count++;
            output << value << " ";
            if (count % 7 == 0) {
                output << std::endl;
            }
            if (count == ((NY+1)*7)) {
                output << std::endl;
                output << "Line 2" << std::endl;
                output << "x y z u v w p" << std::endl;
            }
            if (count == ((NX+1)*7 + (NY+1)*7)) {
                output << std::endl;
                output << "Line 3" << std::endl;
                output << "x y z u v w p" << std::endl;
            }
        }

    input.close();
    output.close();
    std::remove(filename.c_str());
    }
}