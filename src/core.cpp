#include "core.hpp"

void IcoNS::setParallelization()
{
    dims[0] = PX;
    dims[1] = PY;

    dim_x_x = NX / PX;
    dim_y_x = (NY + 1) / PY;
    dim_x_y = (NX + 1) / PX;
    dim_y_y = NY / PY;
    dim_x_z = (NX + 1) / PX;
    dim_y_z = (NY + 1) / PY;
    dim_z = NZ + 1;
    dim_z_z = NZ;

    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &cart_comm);

    MPI_Cart_coords(cart_comm, rank, 2, coords);

    MPI_Cart_shift(cart_comm, 0, 1, &neighbors[0], &neighbors[2]);

    MPI_Cart_shift(cart_comm, 1, 1, &neighbors[3], &neighbors[1]);

    offset_x_x = coords[0] * dim_x_x + std::max(0, coords[0] - (PX - NX % PX));
    offset_y_x = coords[1] * dim_y_x + std::max(0, coords[1] - (PY - (NY + 1) % PY));

    if ((PX - 1 - coords[0]) < NX % PX)
    {
        dim_x_x++;
        resx = 1;
    }

    if ((PY - 1 - coords[1]) < (NY + 1) % PY)
    {
        dim_y_x++;
    }

    offset_x_y = coords[0] * dim_x_y + std::max(0, coords[0] - (PX - (NX + 1) % PX));
    offset_y_y = coords[1] * dim_y_y + std::max(0, coords[1] - (PY - NY % PY));
    if ((PX - 1 - coords[0]) < (NX + 1) % PX)
    {
        dim_x_y++;
    }
    if ((PY - 1 - coords[1]) < NY % PY)
    {
        dim_y_y++;
        resy = 1;
    }

    offset_x_z = coords[0] * dim_x_z + std::max(0, coords[0] - (PX - (NX + 1) % PX));
    offset_y_z = coords[1] * dim_y_z + std::max(0, coords[1] - (PY - (NY + 1) % PY));
    if ((PX - 1 - coords[0]) < (NX + 1) % PX)
    {
        dim_x_z++;
    }
    if ((PY - 1 - coords[1]) < (NY + 1) % PY)
    {
        dim_y_z++;
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
    y2Grid.u.resize(newDimX_x * newDimY_x * (NZ + 1), 0.0);
    y2Grid.v.resize(newDimX_y * newDimY_y * (NZ + 1), 0.0);
    y2Grid.w.resize(newDimX_z * newDimY_z * (NZ), 0.0);

    y3Grid.u.resize(newDimX_x * newDimY_x * (NZ + 1), 0.0);
    y3Grid.v.resize(newDimX_y * newDimY_y * (NZ + 1), 0.0);
    y3Grid.w.resize(newDimX_z * newDimY_z * (NZ), 0.0);

    grid.p.resize((xSize[2] + 2) * (xSize[1] + 2) * xSize[0], 0.0);
    halo_phi.resize((xSize[2] + 2) * (xSize[1] + 2) * xSize[0], 0.0);

    if (BX)
    {
        if (coords[0] == 0)
            lbx++;

        if (coords[0] == PX - 1)
            rbx++;
    }
    if (BY)
    {
        if (coords[1] == PY - 1)
            rby++;

        if (coords[1] == 0)
            lby++;
    }

    if (BZ)
    {
        lbz = 1;
        rbz = 1;
    }

    if (coords[0] == 0)
        firstX = 1;

    if (coords[0] == PX - 1)
        lastX = 1;

    if (coords[1] == PY - 1)
        lastY = 1;

    if (coords[1] == 0)
        firstY = 1;

    boundary.setBoundaryOffsets(lbx, rbx, lby, rby, lbz, rbz);
    boundary.setCoords(coords);
    boundary.setOffsets(offset_x_x, offset_y_x, offset_x_y, offset_y_y, offset_x_z, offset_y_z);

    MPI_Type_vector(dim_x_x, dim_z, (newDimY_x)*dim_z, MPI_REALL, &MPI_face_x_x);
    MPI_Type_commit(&MPI_face_x_x);

    MPI_Type_vector(1, dim_z * newDimY_x, 0, MPI_REALL, &MPI_face_y_x);
    MPI_Type_commit(&MPI_face_y_x);

    MPI_Type_vector(dim_x_y, dim_z, (newDimY_y)*dim_z, MPI_REALL, &MPI_face_x_y);
    MPI_Type_commit(&MPI_face_x_y);

    MPI_Type_vector(1, dim_z * newDimY_y, 0, MPI_REALL, &MPI_face_y_y);
    MPI_Type_commit(&MPI_face_y_y);

    MPI_Type_vector(dim_x_z, dim_z_z, (newDimY_z)*dim_z_z, MPI_REALL, &MPI_face_x_z);
    MPI_Type_commit(&MPI_face_x_z);

    MPI_Type_vector(1, dim_z_z * newDimY_z, 0, MPI_REALL, &MPI_face_y_z);
    MPI_Type_commit(&MPI_face_y_z);

    MPI_Type_vector(xSize[2], xSize[0], (xSize[1] + 2) * xSize[0], MPI_REALL, &MPI_face_x_p);
    MPI_Type_commit(&MPI_face_x_p);

    MPI_Type_vector(1, xSize[0] * (xSize[1] + 2), 0, MPI_REALL, &MPI_face_y_p);
    MPI_Type_commit(&MPI_face_y_p);
}


void IcoNS::set2Decomp()
{
    c2d = new C2Decomp(NZ+1, NY+1, NX+1, PY, PX, periodss);

    // x-pencil size
    xSize[0] = c2d->xSize[0]; 
    xSize[1] = c2d->xSize[1]; 
    xSize[2] = c2d->xSize[2]; 

    // y-pencil size
    ySize[0] = c2d->ySize[0];
    ySize[1] = c2d->ySize[1];
    ySize[2] = c2d->ySize[2];

    // z-pencil size
    zSize[0] = c2d->zSize[0];
    zSize[1] = c2d->zSize[1];
    zSize[2] = c2d->zSize[2];

    // pencils allocation
    c2d->allocX(Y2_p);
}


void IcoNS::setPoissonSolver()
{
    if (BX && BY && BZ)
    {
        poissonSolver = new NeumannPoissonSolver(c2d);
    }
    else if (BX && BY && !BZ)
    {
        poissonSolver = new PeriodicPoissonSolver(c2d);
    }else
    {
        throw unsupportedBoundaryException();
    }
}


void IcoNS::exchangeData(std::vector<Real> &grid_loc, int newDimX, int newDimY, int dim_z, MPI_Datatype MPI_face_x, MPI_Datatype MPI_face_y, int sameX, int sameY)
{
    if (!(BY && coords[1] == 0))
    {
        MPI_Isend(&grid_loc[dim_z * newDimY + dim_z + dim_z * sameY * lastY], 1, MPI_face_x, neighbors[3], 11, cart_comm, &reqs[2]);
    }

    if (!(BY && coords[1] == PY - 1))
    {

        MPI_Irecv(&grid_loc[dim_z * newDimY + (newDimY - 1) * dim_z], 1, MPI_face_x, neighbors[1], 11, cart_comm, &reqs[3]);
        MPI_Wait(&reqs[3], &status);
        MPI_Isend(&grid_loc[dim_z * newDimY + (newDimY - 2 - sameY * firstY) * dim_z], 1, MPI_face_x, neighbors[1], 12, cart_comm, &reqs[3]);
    }
    if (!(BY && coords[1] == 0))
    {
        MPI_Irecv(&grid_loc[dim_z * newDimY], 1, MPI_face_x, neighbors[3], 12, cart_comm, &reqs[2]);
        MPI_Wait(&reqs[2], &status);
    }

    if (!(BX && coords[0] == 0))
    {
        MPI_Isend(&grid_loc[(newDimY)*dim_z + (newDimY)*dim_z * sameX * firstX], 1, MPI_face_y, neighbors[0], 10, cart_comm, &reqs[0]);
    }
    if (!(BX && coords[0] == PX - 1))
    {
        MPI_Irecv(&grid_loc[(dim_z)*newDimY * (newDimX - 1)], 1, MPI_face_y, neighbors[2], 10, cart_comm, &reqs[1]);
        MPI_Wait(&reqs[1], &status);
        MPI_Isend(&grid_loc[newDimY * dim_z * (newDimX - 2 - sameX * lastX)], 1, MPI_face_y, neighbors[2], 9, cart_comm, &reqs[1]);
    }
    if (!(BX && coords[0] == 0))
    {
        MPI_Irecv(&grid_loc[0], 1, MPI_face_y, neighbors[0], 9, cart_comm, &reqs[0]);
        MPI_Wait(&reqs[0], &status);
    }
}


void IcoNS::copyPressureToHalo(Real *p, std::vector<Real> &halo)
{
    for (int i = 0; i < xSize[2]; i++)
    {
        for (int j = 0; j < xSize[1]; j++)
        {
            for (int k = 0; k < xSize[0]; k++)
            {
                halo[(i + 1) * (xSize[1] + 2) * xSize[0] + (j + 1) * xSize[0] + k] = p[i * xSize[1] * xSize[0] + j * xSize[0] + k];
            }
        }
    }
}


void IcoNS::solve()
{
    if(poissonSolver == nullptr) return;
    Real time = 0.0;
    int i = 0;
    double solve_start_time = MPI_Wtime();

    boundary.update_boundary(grid.u, grid.v, grid.w, time);
    MPI_Barrier(cart_comm);
    exchangeData(grid.u, newDimX_x, newDimY_x, dim_z, MPI_face_x_x, MPI_face_y_x, 0, 1);
    exchangeData(grid.v, newDimX_y, newDimY_y, dim_z, MPI_face_x_y, MPI_face_y_y, 1, 0);
    exchangeData(grid.w, newDimX_z, newDimY_z, dim_z_z, MPI_face_x_z, MPI_face_y_z, 1, 1);

    while (i < Nt)
    {
        #ifdef ERROR
        if(testCase==0){
            L2_error(time);
        }
        #endif
        solve_time_step(time);
        time += DT;
        i++;
    }

    double solve_end_time = MPI_Wtime();
    if (rank == 0)
    {
        std::cout << std::endl;
        std::cout << "Total time without init and output: " << solve_end_time - solve_start_time << " seconds" << std::endl;
    }

    output();
    MPI_Barrier(cart_comm);
    cleaner();
}


void IcoNS::parse_input(const std::string &input_file)
{
    std::ifstream file(input_file);
    if (!file.is_open())
    {
        std::cerr << "Process " << rank << ": Error - Cannot open file " << input_file << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    std::string line;

    // Skip comments and empty lines until we find test case number
    while (std::getline(file, line))
    {
        if (line.empty() || line[0] == '#')
            continue;
        std::istringstream iss(line);
        if (!(iss >> testCase))
            continue;
        break;
    }

    // Skip comments and empty lines until we find RE
    while (std::getline(file, line))
    {
        if (line.empty() || line[0] == '#')
            continue;
        std::istringstream iss(line);
        if (!(iss >> RE))
            continue;
        break;
    }

    // Skip comments and empty lines until we find DT
    while (std::getline(file, line))
    {
        if (line.empty() || line[0] == '#')
            continue;
        std::istringstream iss(line);
        if (!(iss >> DT))
            continue;
        break;
    }

    // Skip comments and empty lines until we find Nt
    while (std::getline(file, line))
    {
        if (line.empty() || line[0] == '#')
            continue;
        std::istringstream iss(line);
        if (!(iss >> Nt))
            continue;
        break;
    }

    // Skip comments and empty lines until we find grid points
    while (std::getline(file, line))
    {
        if (line.empty() || line[0] == '#')
            continue;
        std::istringstream iss(line);
        if (!(iss >> SX >> SY >> SZ))
            continue;
        break;
    }
    while (std::getline(file, line))
    {
        if (line.empty() || line[0] == '#')
            continue;

        std::istringstream iss(line);
        std::string sx, sy, sz;

        if (!(iss >> sx >> sy >> sz))
            continue;

        try
        {
            LX = evaluateExpression(sx);
            LY = evaluateExpression(sy);
            LZ = evaluateExpression(sz);
            break;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Error parsing dimensions: " << e.what() << std::endl;
            continue;
        }
    }

    // Skip comments and empty lines until we find grid points
    while (std::getline(file, line))
    {
        if (line.empty() || line[0] == '#')
            continue;
        std::istringstream iss(line);
        if (!(iss >> NX >> NY >> NZ))
            continue;
        break;
    }

    NX -= 1;
    NY -= 1;
    NZ -= 1;

    // Skip comments and empty lines until we find process grid
    while (std::getline(file, line))
    {
        if (line.empty() || line[0] == '#')
            continue;
        std::istringstream iss(line);
        if (!(iss >> PX >> PY >> PZ))
            continue;
        break;
    }

    // Skip comments and empty lines until we find boundaries

    while (std::getline(file, line))
    {
        if (line.empty() || line[0] == '#')
            continue;
        std::istringstream iss(line);
        int x_boundary, y_boundary, z_boundary;
        if (!(iss >> x_boundary >> y_boundary >> z_boundary))
            continue;

        // Convert integers to booleans and assign to extern variables
        BX = static_cast<bool>(x_boundary);
        BY = static_cast<bool>(y_boundary);
        BZ = static_cast<bool>(z_boundary);

        // Update periods array for MPI cart create
        periods[0] = !BX;
        periods[1] = !BY;

        break;
    }

    if (rank == 0)
    {
        std::cout << "\nConfiguration:\n"
                  << "Reynolds number: " << RE << "\n"
                  << "Time step: " << DT << "\n"
                  << "Number of timesteps: " << Nt << "\n"
                  << "Starting point of the domain: " << SX << "x" << SY << "x" << SZ << "\n"
                  << "Domain size: " << LX << " x " << LY << " x " << LZ << "\n"
                  << "Grid points: " << NX + 1 << " x " << NY + 1 << " x " << NZ + 1 << "\n"
                  << "Process grid: " << PX << " x " << PY << " x " << PZ << "\n"
                  << "Boundary conditions: " << BX << " " << BY << " " << BZ << "\n";
    }
    // Verify that the process grid matches the number of processes
    if (PX * PY * PZ != size)
    {
        if (rank == 0)
        {
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


void IcoNS::output()
{
    exchangeData(grid.p, (xSize[2] + 2), (xSize[1] + 2), xSize[0], MPI_face_x_p, MPI_face_y_p, 1, 1);
    MPI_Barrier(cart_comm);
    output_x();
    output_y();
    output_z();
    output_profile();
}


void IcoNS::L2_error(const Real t)
{
    Real error = 0.0, velocityError = 0.0, pressureError = 0.0, totalPressureError;

    error += error_comp_X(t);
    error += error_comp_Y(t);
    error += error_comp_Z(t);
    pressureError += error_comp_P(t);

    MPI_Barrier(cart_comm);

    MPI_Reduce(&error, &velocityError, 1, MPI_REALL, MPI_SUM, 0, cart_comm);
    MPI_Reduce(&pressureError, &totalPressureError, 1, MPI_REALL, MPI_SUM, 0, cart_comm);

    if (rank == 0)
    {
        velocityError = sqrt(velocityError);
        totalPressureError = sqrt(totalPressureError);
        std::cout << " Time : " << t << std::endl
                  << " Velocity error: " << velocityError << " Pressure error: " << totalPressureError << std::endl
                  << std::endl;
    }
}


void IcoNS::cleaner()
{
    c2d->deallocXYZ(poissonSolver->py);
    c2d->deallocXYZ(poissonSolver->pz);
    c2d->deallocXYZ(Y2_p);
    FFTW_PREFIX(cleanup)();
    delete poissonSolver;
    c2d->decompInfoFinalize();
    delete c2d;
    MPI_Type_free(&MPI_face_x_x);
    MPI_Type_free(&MPI_face_y_x);
    MPI_Type_free(&MPI_face_x_y);
    MPI_Type_free(&MPI_face_y_y);
    MPI_Type_free(&MPI_face_x_z);
    MPI_Type_free(&MPI_face_y_z);
    MPI_Type_free(&MPI_face_x_p);
    MPI_Type_free(&MPI_face_y_p);
}