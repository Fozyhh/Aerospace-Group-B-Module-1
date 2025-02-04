#include "core.hpp"

template <typename T>
T to_big_endian(T value)
{
    T result = 0;
    uint8_t *p = reinterpret_cast<uint8_t *>(&result);
    uint8_t *q = reinterpret_cast<uint8_t *>(&value);

    for (size_t i = 0; i < sizeof(T); ++i)
    {
        p[i] = q[sizeof(T) - i - 1];
    }
    return result;
}

void IcoNS::output_x()
{
    MPI_File fh;
    const Real x_middle_x = NX / 2;
    const Real x_middle = (NX + 1) / 2;

    if (rank == 0)
        std::remove("solution_x.vtk");
    MPI_Barrier(cart_comm);
    MPI_File_open(cart_comm, "solution_x.vtk", MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);

    //===========================================
    // Header Writing (Rank 0 only)
    //===========================================
    std::ostringstream header1, header2, header3, header4, header5, header6;
    header1 << std::fixed << std::setprecision(6);
    header2 << std::fixed << std::setprecision(6);
    header3 << std::fixed << std::setprecision(6);
    header4 << std::fixed << std::setprecision(6);
    header5 << std::fixed << std::setprecision(6);
    header6 << std::fixed << std::setprecision(6);

    // VTK metadata
    header1 << "# vtk DataFile Version 3.0\n"
            << "Solution x\n"
            << "BINARY\n"
            << "DATASET STRUCTURED_GRID\n"
            << "DIMENSIONS 1 " << NY + 1 << " " << NZ + 1 << "\n"
            << "POINTS " << (NY + 1) * (NZ + 1) << " " << STRINGA_REAL << "\n";

    // Define data format
    header2 << "POINT_DATA " << (NY + 1) * (NZ + 1) << "\n"
            << "SCALARS u " << STRINGA_REAL << "\n"
            << "LOOKUP_TABLE default\n";
    header3 << "SCALARS v " << STRINGA_REAL << "\n"
            << "LOOKUP_TABLE default\n";
    header4 << "SCALARS w " << STRINGA_REAL << "\n"
            << "LOOKUP_TABLE default\n";
    header5 << "SCALARS p " << STRINGA_REAL << "\n"
            << "LOOKUP_TABLE default\n";
    header6 << "SCALARS Magnitude " << STRINGA_REAL << "\n"
            << "LOOKUP_TABLE default\n";

    MPI_Offset offsetheader1 = header1.str().size(),
               offsetheader2 = header2.str().size(),
               offsetheader3 = header3.str().size(),
               offsetheader4 = header4.str().size(),
               offsetheader5 = header5.str().size(),
               offsetheader6 = header6.str().size(),
               offsetpoints = (3 * sizeof(Real)),
               offsetallpoints = offsetpoints * ((NY + 1) * (NZ + 1)),
               offsetvalue = sizeof(Real),
               offsetallvalue = offsetvalue * ((NY + 1) * (NZ + 1));
    if (rank == 0)
    {
        MPI_File_write_at(fh, 0, header1.str().c_str(), header1.str().size(), MPI_CHAR, MPI_STATUS_IGNORE);
        MPI_File_write_at(fh, offsetheader1 + offsetallpoints, header2.str().c_str(), header2.str().size(), MPI_CHAR, MPI_STATUS_IGNORE);
        MPI_File_write_at(fh, offsetheader1 + offsetallpoints + offsetheader2 + offsetallvalue, header3.str().c_str(), header3.str().size(), MPI_CHAR, MPI_STATUS_IGNORE);
        MPI_File_write_at(fh, offsetheader1 + offsetallpoints + offsetheader2 + 2 * offsetallvalue + offsetheader3, header4.str().c_str(), header4.str().size(), MPI_CHAR, MPI_STATUS_IGNORE);
        MPI_File_write_at(fh, offsetheader1 + offsetallpoints + offsetheader2 + 3 * offsetallvalue + offsetheader3 + offsetheader4, header5.str().c_str(), header5.str().size(), MPI_CHAR, MPI_STATUS_IGNORE);
        MPI_File_write_at(fh, offsetheader1 + offsetallpoints + offsetheader2 + 4 * offsetallvalue + offsetheader3 + offsetheader4 + offsetheader5, header6.str().c_str(), header6.str().size(), MPI_CHAR, MPI_STATUS_IGNORE);
    }
    MPI_Barrier(cart_comm);

    // Sync offset across all processes
    // MPI_Bcast(&offset, 1, MPI_OFFSET, 0, cart_comm);

    //===========================================
    // Data Collection
    //===========================================

    // Find process containing middle slice
    int x_index_x = static_cast<int>(x_middle_x);
    int x_index = static_cast<int>(x_middle);

    int offset_x_x_ = offset_x_x - 1;
    int offset_y_x_ = offset_y_x - 1;
    int offset_x_y_ = offset_x_y - 1;
    int offset_y_y_ = offset_y_y - 1;
    int offset_x_z_ = offset_x_z - 1;
    int offset_y_z_ = offset_y_z - 1;
    int offset_x_p_ = coords[0] * xSize[2] - 1;

    MPI_Barrier(cart_comm);
    if (x_index >= offset_x_x_ + 1 && x_index < dim_x_x + offset_x_x_ + 1)
    {

        int local_x_x = x_index_x - offset_x_x_;
        int local_x_y = x_index - offset_x_y_;
        int local_x_z = x_index - offset_x_z_;
        int local_x_p = x_index - offset_x_p_;

        Real value_x = 0, value_y = 0, value_z = 0, value_p = 0, value_m = 0;
        Real bg_px = 0, bg_py = 0, bg_pz = 0, bg_vx = 0, bg_vy = 0, bg_vz = 0, bg_vp = 0, bg_vm = 0;
        Real point_x = 0, point_y = 0, point_z = 0;
        for (int j = 1; j < newDimY_x - 1; j++)
        {
            for (int k = 0; k < dim_z; k++)
            {

                // Write grid points coordinate
                point_x = SX + LX / 2;
                point_y = static_cast<Real>(SY + (j + offset_y_x_) * DY);
                point_z = static_cast<Real>(SZ + (k)*DZ);

                value_x = grid.u[getx(local_x_x, j, k)];

                if (lby && j == 1)
                {
                    value_y = boundary.boundary_value_v[FRONT]->value(x_index, j + offset_y_y_ - 0.5, k, DT*Nt);
                }
                else if (rby && j == newDimY_x - 2)
                {
                    value_y = boundary.boundary_value_v[BACK]->value(x_index, j + offset_y_y_ - 0.5, k, DT*Nt);
                }
                else
                {
                    value_y = (grid.v[gety(local_x_y, j + resy, k)] + grid.v[gety(local_x_y, j + resy - 1, k)] + 
                                grid.v[gety(local_x_y - 1, j + resy, k)] + grid.v[gety(local_x_y-1, j + resy-1, k)])/4;
                }

                if (k == 0 && lbz)
                {
                    value_z = boundary.boundary_value_w[LOWER]->value(x_index, j + offset_y_z_, k - 0.5, DT*Nt);
                }
                else if (k == dim_z - 1 && rbz)
                {
                    value_z = boundary.boundary_value_w[UPPER]->value(x_index, j + offset_y_z_, k - 0.5, DT*Nt);
                }
                else
                {
                    value_z = (grid.w[getz(local_x_z, j, k)] + grid.w[getz(local_x_z, j, k - 1)] + grid.w[getz(local_x_z-1, j, k)] + grid.w[getz(local_x_z -1, j, k -1)])/4;
                }

                value_p = (grid.p[getHaloP(local_x_p, j, k)] + grid.p[getHaloP(local_x_p -1, j, k)])/2;
                value_m = std::sqrt(value_x * value_x + value_y * value_y + value_z * value_z);

                bg_px = to_big_endian(point_x);
                bg_py = to_big_endian(point_y);
                bg_pz = to_big_endian(point_z);
                bg_vx = to_big_endian(value_x);
                bg_vy = to_big_endian(value_y);
                bg_vz = to_big_endian(value_z);
                bg_vp = to_big_endian(value_p);
                bg_vm = to_big_endian(value_m);

                MPI_File_write_at(fh, offsetheader1 + offsetpoints * ((j + offset_y_x_) * dim_z + k), &bg_px, 1, MPI_REALL, MPI_STATUS_IGNORE);

                MPI_File_write_at(fh, offsetheader1 + offsetpoints * ((j + offset_y_x_) * dim_z + k) + sizeof(Real), &bg_py, 1, MPI_REALL, MPI_STATUS_IGNORE);

                MPI_File_write_at(fh, offsetheader1 + offsetpoints * ((j + offset_y_x_) * dim_z + k) + 2 * sizeof(Real), &bg_pz, 1, MPI_REALL, MPI_STATUS_IGNORE);

                MPI_File_write_at(fh, offsetheader1 + offsetallpoints + offsetheader2 + offsetvalue * ((j + offset_y_x_) * dim_z + k), &bg_vx, 1, MPI_REALL, MPI_STATUS_IGNORE);

                MPI_File_write_at(fh, offsetheader1 + offsetallpoints + offsetheader2 + offsetallvalue + offsetheader3 + offsetvalue * ((j + offset_y_x_) * dim_z + k), &bg_vy, 1, MPI_REALL, MPI_STATUS_IGNORE);

                MPI_File_write_at(fh, offsetheader1 + offsetallpoints + offsetheader2 + 2 * offsetallvalue + offsetheader3 + offsetheader4 + offsetvalue * ((j + offset_y_x_) * dim_z + k), &bg_vz, 1, MPI_REALL, MPI_STATUS_IGNORE);

                MPI_File_write_at(fh, offsetheader1 + offsetallpoints + offsetheader2 + 3 * offsetallvalue + offsetheader3 + offsetheader4 + offsetheader5 + offsetvalue * ((j + offset_y_x_) * dim_z + k), &bg_vp, 1, MPI_REALL, MPI_STATUS_IGNORE);

                MPI_File_write_at(fh, offsetheader1 + offsetallpoints + offsetheader2 + 4 * offsetallvalue + offsetheader3 + offsetheader4 + offsetheader5 + offsetheader6 + offsetvalue * ((j + offset_y_x_) * dim_z + k), &bg_vm, 1, MPI_REALL, MPI_STATUS_IGNORE);
            }
        }
    }
    MPI_Barrier(cart_comm);
    MPI_File_close(&fh);
}

void IcoNS::output_y()
{
    MPI_File fh;
    const Real y_middle = (NY + 1) / 2;
    const Real y_middle_y = NY / 2;
    if (rank == 0)
        std::remove("solution_y.vtk");
    MPI_Barrier(cart_comm);
    MPI_File_open(cart_comm, "solution_y.vtk", MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);

    //===========================================
    // Header Writing (Rank 0 only)
    //===========================================
    std::ostringstream header1, header2, header3, header4, header5, header6;
    header1 << std::fixed << std::setprecision(6);
    header2 << std::fixed << std::setprecision(6);
    header3 << std::fixed << std::setprecision(6);
    header4 << std::fixed << std::setprecision(6);
    header5 << std::fixed << std::setprecision(6);
    header6 << std::fixed << std::setprecision(6);

    // VTK metadata
    header1 << "# vtk DataFile Version 3.0\n"
            << "Solution y\n"
            << "BINARY\n"
            << "DATASET STRUCTURED_GRID\n"
            << "DIMENSIONS " << NX + 1 << " " << 1 << " " << NZ + 1 << "\n"
            << "POINTS " << (NX + 1) * (NZ + 1) << " " << STRINGA_REAL << "\n";

    // Define data format
    header2 << "POINT_DATA " << (NX + 1) * (NZ + 1) << "\n"
            << "SCALARS u " << STRINGA_REAL << "\n"
            << "LOOKUP_TABLE default\n";
    header3 << "SCALARS v " << STRINGA_REAL << "\n"
            << "LOOKUP_TABLE default\n";
    header4 << "SCALARS w " << STRINGA_REAL << "\n"
            << "LOOKUP_TABLE default\n";
    header5 << "SCALARS p " << STRINGA_REAL << "\n"
            << "LOOKUP_TABLE default\n";
    header6 << "SCALARS Magnitude " << STRINGA_REAL << "\n"
            << "LOOKUP_TABLE default\n";

    MPI_Offset offsetheader1 = header1.str().size(),
               offsetheader2 = header2.str().size(),
               offsetheader3 = header3.str().size(),
               offsetheader4 = header4.str().size(),
               offsetheader5 = header5.str().size(),
               offsetheader6 = header6.str().size(),
               offsetpoints = (3 * sizeof(Real)),
               offsetallpoints = offsetpoints * ((NX + 1) * (NZ + 1)),
               offsetvalue = sizeof(Real),
               offsetallvalue = offsetvalue * ((NX + 1) * (NZ + 1));
    if (rank == 0)
    {
        MPI_File_write_at(fh, 0, header1.str().c_str(), header1.str().size(), MPI_CHAR, MPI_STATUS_IGNORE);
        MPI_File_write_at(fh, offsetheader1 + offsetallpoints, header2.str().c_str(), header2.str().size(), MPI_CHAR, MPI_STATUS_IGNORE);
        MPI_File_write_at(fh, offsetheader1 + offsetallpoints + offsetheader2 + offsetallvalue, header3.str().c_str(), header3.str().size(), MPI_CHAR, MPI_STATUS_IGNORE);
        MPI_File_write_at(fh, offsetheader1 + offsetallpoints + offsetheader2 + 2 * offsetallvalue + offsetheader3, header4.str().c_str(), header4.str().size(), MPI_CHAR, MPI_STATUS_IGNORE);
        MPI_File_write_at(fh, offsetheader1 + offsetallpoints + offsetheader2 + 3 * offsetallvalue + offsetheader3 + offsetheader4, header5.str().c_str(), header5.str().size(), MPI_CHAR, MPI_STATUS_IGNORE);
        MPI_File_write_at(fh, offsetheader1 + offsetallpoints + offsetheader2 + 4 * offsetallvalue + offsetheader3 + offsetheader4 + offsetheader5, header6.str().c_str(), header6.str().size(), MPI_CHAR, MPI_STATUS_IGNORE);
    }

    MPI_Barrier(cart_comm);
    MPI_File_sync(fh);
    // Find process containing middle slice
    int y_index = static_cast<int>(y_middle);
    int y_index_y = static_cast<int>(y_middle_y);
    int offset_x_x_ = offset_x_x - 1;
    int offset_y_x_ = offset_y_x - 1;
    int offset_x_y_ = offset_x_y - 1;
    int offset_y_y_ = offset_y_y - 1;
    int offset_x_z_ = offset_x_z - 1;
    int offset_y_z_ = offset_y_z - 1;

    int offset_y_p_ = coords[1] * zSize[1] - 1;


    MPI_Barrier(cart_comm);
    if (y_index >= offset_y_y && y_index < dim_y_y + offset_y_y)
    {
        int local_y_x = y_index - offset_y_x_;
        int local_y_y = y_index_y - offset_y_y_;
        int local_y_z = y_index - offset_y_z_;
        int local_y_p = y_index - offset_y_p_;

        Real value_x = 0, value_y = 0, value_z = 0, value_p = 0, value_m = 0;
        Real bg_px = 0, bg_py = 0, bg_pz = 0, bg_vx = 0, bg_vy = 0, bg_vz = 0, bg_vp = 0, bg_vm = 0;
        Real point_x = 0, point_y = 0, point_z = 0;
        for (int i = 1; i < newDimX_y - 1; i++)
        {
            for (int k = 0; k < dim_z; k++)
            {
                point_x = static_cast<Real>(SX + (i + offset_x_y_) * DX);
                point_y = SY + LY / 2;
                point_z = static_cast<Real>(SZ + (k)*DZ);

                value_y = grid.v[gety(i, local_y_y, k)];

                if (lbx && i == 1)
                {
                    value_x = boundary.boundary_value_u[0]->value(i + offset_x_x_ - 0.5, y_index + offset_y_x_, k, DT*Nt);
                }
                else if (rbx && i == newDimX_y - 2)
                {
                    value_x = boundary.boundary_value_u[1]->value(i + offset_x_x_ - 0.5, y_index + offset_y_x_, k, DT*Nt);
                }
                else
                {
                    value_x = (grid.u[getx(i + resx, local_y_x, k)] + grid.u[getx(i + resx - 1, local_y_x, k)] + grid.u[getx(i + resx, local_y_x - 1, k)] + grid.u[getx(i + resx - 1, local_y_x - 1, k)])/4;
                }

                if (k == 0 && lbz)
                {
                    value_z = boundary.boundary_value_w[4]->value(i + offset_x_z_, y_index, k - 0.5, DT*Nt);
                }
                else if (k == dim_z - 1 && rbz)
                {
                    value_z = boundary.boundary_value_w[5]->value(i + offset_x_z_, y_index, k - 0.5, DT*Nt);
                }
                else
                {
                    value_z = (grid.w[getz(i, local_y_z, k)] + grid.w[getz(i, local_y_z, k-1)] + grid.w[getz(i, local_y_z-1, k)] + grid.w[getz(i, local_y_z-1, k-1)])/4;
                }

                value_p = (grid.p[getHaloP(i, local_y_p, k)] + grid.p[getHaloP(i, local_y_p - 1, k)])/2;

                value_m = std::sqrt(value_x * value_x + value_y * value_y + value_z * value_z);

                bg_px = to_big_endian(point_x);
                bg_py = to_big_endian(point_y);
                bg_pz = to_big_endian(point_z);
                bg_vx = to_big_endian(value_x);
                bg_vy = to_big_endian(value_y);
                bg_vz = to_big_endian(value_z);
                bg_vp = to_big_endian(value_p);
                bg_vm = to_big_endian(value_m);
                MPI_File_write_at(fh, offsetheader1 + offsetpoints * ((i + offset_x_y_) * dim_z + k), &bg_px, 1, MPI_REALL, MPI_STATUS_IGNORE);

                MPI_File_write_at(fh, offsetheader1 + offsetpoints * ((i + offset_x_y_) * dim_z + k) + sizeof(Real), &bg_py, 1, MPI_REALL, MPI_STATUS_IGNORE);

                MPI_File_write_at(fh, offsetheader1 + offsetpoints * ((i + offset_x_y_) * dim_z + k) + 2 * sizeof(Real), &bg_pz, 1, MPI_REALL, MPI_STATUS_IGNORE);

                MPI_File_write_at(fh, offsetheader1 + offsetallpoints + offsetheader2 + offsetvalue * ((i + offset_x_y_) * dim_z + k), &bg_vx, 1, MPI_REALL, MPI_STATUS_IGNORE);

                MPI_File_write_at(fh, offsetheader1 + offsetallpoints + offsetheader2 + offsetallvalue + offsetheader3 + offsetvalue * ((i + offset_x_y_) * dim_z + k), &bg_vy, 1, MPI_REALL, MPI_STATUS_IGNORE);

                MPI_File_write_at(fh, offsetheader1 + offsetallpoints + offsetheader2 + 2 * offsetallvalue + offsetheader3 + offsetheader4 + offsetvalue * ((i + offset_x_y_) * dim_z + k), &bg_vz, 1, MPI_REALL, MPI_STATUS_IGNORE);

                MPI_File_write_at(fh, offsetheader1 + offsetallpoints + offsetheader2 + 3 * offsetallvalue + offsetheader3 + offsetheader4 + offsetheader5 + offsetvalue * ((i + offset_x_y_) * dim_z + k), &bg_vp, 1, MPI_REALL, MPI_STATUS_IGNORE);

                MPI_File_write_at(fh, offsetheader1 + offsetallpoints + offsetheader2 + 4 * offsetallvalue + offsetheader3 + offsetheader4 + offsetheader5 + offsetheader6 + offsetvalue * ((i + offset_x_y_) * dim_z + k), &bg_vm, 1, MPI_REALL, MPI_STATUS_IGNORE);
            }
        }
    }
    MPI_Barrier(cart_comm);
    MPI_File_close(&fh);
}

void IcoNS::output_z()
{
    MPI_File fh;
    const Real z_middle= (NZ + 1) / 2;
    const Real z_middle_z = NZ / 2;
    if (rank == 0)
        std::remove("solution_z.vtk");
    MPI_Barrier(cart_comm);
    MPI_File_open(cart_comm, "solution_z.vtk", MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);



    //===========================================
    // Header Writing (Rank 0 only)
    //===========================================
    std::ostringstream header1, header2, header3, header4, header5, header6;
    header1 << std::fixed << std::setprecision(6);
    header2 << std::fixed << std::setprecision(6);
    header3 << std::fixed << std::setprecision(6);
    header4 << std::fixed << std::setprecision(6);
    header5 << std::fixed << std::setprecision(6);
    header6 << std::fixed << std::setprecision(6);

    // VTK metadata
    header1 << "# vtk DataFile Version 3.0\n"
            << "Solution z\n"
            << "BINARY\n"
            << "DATASET STRUCTURED_GRID\n"
            << "DIMENSIONS " << NX + 1 << " " << NY + 1 << " " << "1" << "\n"
            << "POINTS " << (NY + 1) * (NX + 1) << " " << STRINGA_REAL << "\n";

    // Define data format
    header2 << "POINT_DATA " << (NY + 1) * (NX + 1) << "\n"
            << "SCALARS u " << STRINGA_REAL << "\n"
            << "LOOKUP_TABLE default\n";
    header3 << "SCALARS v " << STRINGA_REAL << "\n"
            << "LOOKUP_TABLE default\n";
    header4 << "SCALARS w " << STRINGA_REAL << "\n"
            << "LOOKUP_TABLE default\n";
    header5 << "SCALARS p " << STRINGA_REAL << "\n"
            << "LOOKUP_TABLE default\n";
    header6 << "SCALARS Magnitude " << STRINGA_REAL << "\n"
            << "LOOKUP_TABLE default\n";

    MPI_Offset offsetheader1 = header1.str().size(),
               offsetheader2 = header2.str().size(),
               offsetheader3 = header3.str().size(),
               offsetheader4 = header4.str().size(),
               offsetheader5 = header5.str().size(),
               offsetheader6 = header6.str().size(),
               offsetpoints = (3 * sizeof(Real)),
               offsetallpoints = offsetpoints * ((NY + 1) * (NX + 1)),
               offsetvalue = sizeof(Real),
               offsetallvalue = offsetvalue * ((NY + 1) * (NX + 1));
    if (rank == 0)
    {
        MPI_File_write_at(fh, 0, header1.str().c_str(), header1.str().size(), MPI_CHAR, MPI_STATUS_IGNORE);
        MPI_File_write_at(fh, offsetheader1 + offsetallpoints, header2.str().c_str(), header2.str().size(), MPI_CHAR, MPI_STATUS_IGNORE);
        MPI_File_write_at(fh, offsetheader1 + offsetallpoints + offsetheader2 + offsetallvalue, header3.str().c_str(), header3.str().size(), MPI_CHAR, MPI_STATUS_IGNORE);
        MPI_File_write_at(fh, offsetheader1 + offsetallpoints + offsetheader2 + 2 * offsetallvalue + offsetheader3, header4.str().c_str(), header4.str().size(), MPI_CHAR, MPI_STATUS_IGNORE);
        MPI_File_write_at(fh, offsetheader1 + offsetallpoints + offsetheader2 + 3 * offsetallvalue + offsetheader3 + offsetheader4, header5.str().c_str(), header5.str().size(), MPI_CHAR, MPI_STATUS_IGNORE);
        MPI_File_write_at(fh, offsetheader1 + offsetallpoints + offsetheader2 + 4 * offsetallvalue + offsetheader3 + offsetheader4 + offsetheader5, header6.str().c_str(), header6.str().size(), MPI_CHAR, MPI_STATUS_IGNORE);
    }

    //===========================================
    // Data Collection
    //===========================================

    // Find process containing middle slice
    int z_index = static_cast<int>(z_middle);
    int z_index_z = static_cast<int>(z_middle_z);
    int offset_x_x_ = offset_x_x - 1;
    int offset_y_x_ = offset_y_x - 1;
    int offset_x_y_ = offset_x_y - 1;
    int offset_y_y_ = offset_y_y - 1;
    int offset_x_z_ = offset_x_z - 1;
    int offset_y_z_ = offset_y_z - 1;

    for (int i = 1; i < newDimX_z - 1; i++)
    {
        for (int j = 1; j < newDimY_z - 1; j++)
        {

            // Write grid points coordinate
            Real point_x = static_cast<Real>(SX + (i + offset_x_z_) * DX),
                 point_y = static_cast<Real>(SY + (j + offset_y_z_) * DY),
                 point_z = SZ + LZ / 2;

            // valuesx[rank*(dim_y_x * dim_z + dim_z) + ((j-1) * dim_z + k)] = grid.v[local_x* newDimY_x * dim_z + j * dim_z + k];
            Real value_x, value_y, value_z, value_p, value_m;
            if (lby && j == 1)
            {
                value_x = (boundary.boundary_value_u[FRONT]->value(i + offset_x_x_ - 0.5, j + offset_y_x_, z_index - 0.5, DT*Nt));
            }
            else if (rby && j == newDimY_z - 2)
            {
                value_x = (boundary.boundary_value_u[BACK]->value(i + offset_x_x_ - 0.5, j + offset_y_x_, z_index - 0.5, DT*Nt));
            }
            else
            {
                value_x = (grid.u[getx(i + resx, j, z_index)] + grid.u[getx(i + resx, j, z_index -1)] + grid.u[getx(i + resx - 1, j, z_index)] + grid.u[getx(i + resx-1, j, z_index-1)])/4;
            }

            if (lbx && i == 1)
            {
                value_y = (boundary.boundary_value_v[LEFT]->value(i + offset_x_y_, j + offset_y_y_ - 0.5, z_index - 0.5, DT*Nt));
            }
            else if (rbx && i == newDimX_z - 2)
            {
                value_y = boundary.boundary_value_v[RIGHT]->value(i + offset_x_y, j + offset_y_y_ - 0.5, z_index - 0.5, DT*Nt);
            }
            else
            {
                value_y = (grid.v[gety(i, j + resy, z_index)] + grid.v[gety(i, j + resy - 1, z_index)] + grid.v[gety(i, j + resy, z_index - 1)] + grid.v[gety(i, j + resy - 1, z_index - 1)])/4;
            }
            // std::cout << i <<" "<< j<< " " << z_index_z << ",   " << grid.w[getz(i, j, z_index_z)] << std::endl;
            // int stop; std::cin >> stop;
            value_z = grid.w[getz(i, j, z_index_z)];

            value_p = (grid.p[getHaloP(i, j, z_index)] + grid.p[getHaloP(i, j, z_index - 1)])/2;

            value_m = std::sqrt(value_x * value_x + value_y * value_y + value_z * value_z);

            Real bg_px = to_big_endian(point_x);
            Real bg_py = to_big_endian(point_y);
            Real bg_pz = to_big_endian(point_z);
            Real bg_vx = to_big_endian(value_x);
            Real bg_vy = to_big_endian(value_y);
            Real bg_vz = to_big_endian(value_z);
            Real bg_vp = to_big_endian(value_p);
            Real bg_vm = to_big_endian(value_m);

            MPI_File_write_at(fh, offsetheader1 + offsetpoints * ((i + offset_x_z_) * dim_z + j + offset_y_z_), &bg_px, 1, MPI_REALL, MPI_STATUS_IGNORE);

            MPI_File_write_at(fh, offsetheader1 + offsetpoints * ((i + offset_x_z_) * dim_z + j + offset_y_z_) + sizeof(Real), &bg_py, 1, MPI_REALL, MPI_STATUS_IGNORE);

            MPI_File_write_at(fh, offsetheader1 + offsetpoints * ((i + offset_x_z_) * dim_z + j + offset_y_z_) + 2 * sizeof(Real), &bg_pz, 1, MPI_REALL, MPI_STATUS_IGNORE);

            MPI_File_write_at(fh, offsetheader1 + offsetallpoints + offsetheader2 + offsetvalue * ((i + offset_x_z_) * dim_z + j + offset_y_z_), &bg_vx, 1, MPI_REALL, MPI_STATUS_IGNORE);

            MPI_File_write_at(fh, offsetheader1 + offsetallpoints + offsetheader2 + offsetallvalue + offsetheader3 + offsetvalue * ((i + offset_x_z_) * dim_z + j + offset_y_z_), &bg_vy, 1, MPI_REALL, MPI_STATUS_IGNORE);

            MPI_File_write_at(fh, offsetheader1 + offsetallpoints + offsetheader2 + 2 * offsetallvalue + offsetheader3 + offsetheader4 + offsetvalue * ((i + offset_x_z_) * dim_z + j + offset_y_z_), &bg_vz, 1, MPI_REALL, MPI_STATUS_IGNORE);

            MPI_File_write_at(fh, offsetheader1 + offsetallpoints + offsetheader2 + 3 * offsetallvalue + offsetheader3 + offsetheader4 + offsetheader5 + offsetvalue * ((i + offset_x_z_) * dim_z + j + offset_y_z_), &bg_vp, 1, MPI_REALL, MPI_STATUS_IGNORE);

            MPI_File_write_at(fh, offsetheader1 + offsetallpoints + offsetheader2 + 4 * offsetallvalue + offsetheader3 + offsetheader4 + offsetheader5 + offsetheader6 + offsetvalue * ((i + offset_x_z_) * dim_z + j + offset_y_z_), &bg_vm, 1, MPI_REALL, MPI_STATUS_IGNORE);
        }
    }

    MPI_File_close(&fh);
}

void IcoNS::output_profile()
{
    const std::string filename = "profile" + std::to_string(testCase) + ".dat";
    const std::string filenametxt = "profile" + std::to_string(testCase) + ".txt";
    if (rank == 0){
        std::remove(filename.c_str());
        std::remove(filenametxt.c_str());
    }
    MPI_File fh;
    MPI_File_open(cart_comm, filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
    MPI_Offset offset = coords[1] * xSize[1] * sizeof(Real) * 7;

    // LINE 1
    if (coords[0] == static_cast<int>(PX / 2))
    {
        const Real xCoord = SX + LX / 2;
        Real yCoord = SY + coords[1] * xSize[1] * DY;
        const Real zCoord = SZ + LZ / 2;
        Real u, v, w, p;
        for (int i = 1; i < xSize[1] + 1; i++)
        {
            MPI_File_write_at(fh, offset, &xCoord, 1, MPI_REALL, MPI_STATUS_IGNORE);
            offset += sizeof(Real);

            MPI_File_write_at(fh, offset, &yCoord, 1, MPI_REALL, MPI_STATUS_IGNORE);
            yCoord += DY;
            offset += sizeof(Real);

            MPI_File_write_at(fh, offset, &zCoord, 1, MPI_REALL, MPI_STATUS_IGNORE);
            offset += sizeof(Real);

            if (PX % 2 == 0)
            {
                u = (grid.u[getx(1, i, (dim_z - 1) / 2)] +
                     grid.u[getx(1, i, (dim_z - 1) / 2 + 1)]) /
                    2;
                if (lby && i == 1)
                {
                    v = boundary.boundary_value_v[FRONT]->value(xCoord, yCoord-DY/2, zCoord, DT*Nt);
                }
                else if (rby && i == xSize[1])
                {
                    v = boundary.boundary_value_v[BACK]->value(xCoord, yCoord-DY/2, zCoord, DT*Nt);
                }
                else
                {
                    v = (grid.v[gety(0, i + resy, (dim_z - 1) / 2)] +
                         grid.v[gety(0, i + resy, (dim_z - 1) / 2 + 1)] +
                         grid.v[gety(0, i + resy - 1, (dim_z - 1) / 2)] +
                         grid.v[gety(0, i + resy - 1, (dim_z - 1) / 2 + 1)] +
                         grid.v[gety(1, i + resy, (dim_z - 1) / 2)] +
                         grid.v[gety(1, i + resy, (dim_z - 1) / 2 + 1)] +
                         grid.v[gety(1, i + resy - 1, (dim_z - 1) / 2)] +
                         grid.v[gety(1, i + resy - 1, (dim_z - 1) / 2 + 1)]) /
                        8;
                }
                w = (grid.w[getz(0, i, (dim_z_z) / 2)] +
                     grid.w[getz(1, i, (dim_z_z) / 2)]) /
                    2;
                p = (grid.p[getHaloP(0, i, (xSize[0] - 1) / 2)] +
                     grid.p[getHaloP(0, i, (xSize[0] - 1) / 2 + 1)] +
                     grid.p[getHaloP(1, i, (xSize[0] - 1) / 2)] +
                     grid.p[getHaloP(1, i, (xSize[0] - 1) / 2 + 1)]) /
                    4;
            }
            else
            { // should distinguish between even or odd number of pressure points on interested processors
              // assuming even for now
              // with 480 cells this is slightly wrong for PX = 7, 11, 13, 21, 23, 27, 29 ...
                u = (grid.u[getx((newDimX_x - 1) / 2, i, (dim_z - 1) / 2)] +
                     grid.u[getx((newDimX_x - 1) / 2, i, (dim_z - 1) / 2 + 1)]) /
                    2;
                if (lby && i == 1)
                { // approximating at the front boundary, could use boundary_value_v[2]->value()
                    v = boundary.boundary_value_v[2]->value(xCoord, yCoord-DY/2, zCoord, DT*Nt);
                }
                else if (rby && i == xSize[1])
                { // approximating at the back boundary, could use boundary_value_v[3]->value()
                    v = boundary.boundary_value_v[3]->value(xCoord, yCoord-DY/2, zCoord, DT*Nt);
                }
                else
                {
                    v = (grid.v[gety((newDimX_y - 1) / 2, i, (dim_z - 1) / 2)] +
                         grid.v[gety((newDimX_y - 1) / 2, i, (dim_z - 1) / 2 + 1)] +
                         grid.v[gety((newDimX_y - 1) / 2, i - 1, (dim_z - 1) / 2)] +
                         grid.v[gety((newDimX_y - 1) / 2, i - 1, (dim_z - 1) / 2 + 1)] +
                         grid.v[gety((newDimX_y - 1) / 2 + 1, i, (dim_z - 1) / 2)] +
                         grid.v[gety((newDimX_y - 1) / 2 + 1, i, (dim_z - 1) / 2 + 1)] +
                         grid.v[gety((newDimX_y - 1) / 2 + 1, i - 1, (dim_z - 1) / 2)] +
                         grid.v[gety((newDimX_y - 1) / 2 + 1, i - 1, (dim_z - 1) / 2 + 1)]) /
                        8;
                }
                w = (grid.w[getz((newDimX_z - 1) / 2, i, (dim_z_z) / 2)] +
                     grid.w[getz((newDimX_z - 1) / 2 + 1, i, (dim_z_z) / 2)]) /
                    2;
                p = (grid.p[getHaloP((xSize[2] - 0) / 2, i, (xSize[0] - 1) / 2)] +
                     grid.p[getHaloP((xSize[2] - 0) / 2, i, (xSize[0] - 1) / 2 + 1)] +
                     grid.p[getHaloP((xSize[2] - 0) / 2 + 1, i, (xSize[0] - 1) / 2)] +
                     grid.p[getHaloP((xSize[2] - 0) / 2 + 1, i, (xSize[0] - 1) / 2 + 1)]) /
                    4;
            }
            MPI_File_write_at(fh, offset, &u, 1, MPI_REALL, MPI_STATUS_IGNORE);
            offset += sizeof(Real);
            MPI_File_write_at(fh, offset, &v, 1, MPI_REALL, MPI_STATUS_IGNORE);
            offset += sizeof(Real);
            MPI_File_write_at(fh, offset, &w, 1, MPI_REALL, MPI_STATUS_IGNORE);
            offset += sizeof(Real);
            MPI_File_write_at(fh, offset, &p, 1, MPI_REALL, MPI_STATUS_IGNORE);
            offset += sizeof(Real);
        }
    }

    // LINE 2
    if (coords[1] == static_cast<int>(PY / 2))
    {
        offset = PY * xSize[1] * sizeof(Real) * 7 + coords[0] * xSize[2] * sizeof(Real) * 7;
        Real xCoord = SX + coords[0] * xSize[2] * DX;
        const Real yCoord = SY + LY / 2;
        const Real zCoord = SZ + LZ / 2;
        Real u, v, w, p;
        for (int i = 1; i < xSize[2] + 1; i++)
        {
            MPI_File_write_at(fh, offset, &xCoord, 1, MPI_REALL, MPI_STATUS_IGNORE);
            xCoord += DX;
            offset += sizeof(Real);

            MPI_File_write_at(fh, offset, &yCoord, 1, MPI_REALL, MPI_STATUS_IGNORE);
            offset += sizeof(Real);

            MPI_File_write_at(fh, offset, &zCoord, 1, MPI_REALL, MPI_STATUS_IGNORE);
            offset += sizeof(Real);

            if (PY % 2 == 0)
            {
                if (lbx && i == 1)
                {
                    u = boundary.boundary_value_u[LEFT]->value(xCoord - DX / 2, yCoord, zCoord, DT*Nt);
                }
                else if (rbx && i == xSize[2])
                {
                    u = boundary.boundary_value_u[RIGHT]->value(xCoord - DX / 2, yCoord, zCoord, DT*Nt);
                }
                else
                {
                    u = (grid.u[getx(i + resx, 0, (dim_z - 1) / 2)] +
                         grid.u[getx(i + resx, 0, (dim_z - 1) / 2 + 1)] +
                         grid.u[getx(i + resx, 1, (dim_z - 1) / 2)] +
                         grid.u[getx(i + resx, 1, (dim_z - 1) / 2 + 1)] +
                         grid.u[getx(i + resx - 1, 0, (dim_z - 1) / 2)] +
                         grid.u[getx(i + resx - 1, 0, (dim_z - 1) / 2 + 1)] +
                         grid.u[getx(i + resx - 1, 1, (dim_z - 1) / 2)] +
                         grid.u[getx(i + resx - 1, 1, (dim_z - 1) / 2 + 1)]) /
                        8;
                }
                v = (grid.v[gety(i, 0 + resy, (dim_z - 1) / 2)] +
                     grid.v[gety(i, 0 + resy, (dim_z - 1) / 2 + 1)]) /
                    2;
                w = (grid.w[getz(i, 0, (dim_z_z) / 2)] +
                     grid.w[getz(i, 1, (dim_z_z) / 2)]) /
                    2;
                p = (grid.p[getHaloP(i, 0, (xSize[0] - 1) / 2)] +
                     grid.p[getHaloP(i, 0, (xSize[0] - 1) / 2 + 1)] +
                     grid.p[getHaloP(i, 1, (xSize[0] - 1) / 2)] +
                     grid.p[getHaloP(i, 1, (xSize[0] - 1) / 2 + 1)]) /
                    4;
            }
            else
            { // same as above
                if (lbx && i == 1)
                {
                    u = boundary.boundary_value_u[0]->value(xCoord - DX / 2, yCoord, zCoord, DT*Nt);
                }
                else if (rbx && i == xSize[2])
                {
                    u = boundary.boundary_value_u[1]->value(xCoord - DX / 2, yCoord, zCoord, DT*Nt);
                }
                else
                {
                    u = (grid.u[getx(i, (newDimY_x - 1) / 2, (dim_z - 1) / 2)] +
                         grid.u[getx(i, (newDimY_x - 1) / 2, (dim_z - 1) / 2 + 1)] +
                         grid.u[getx(i, (newDimY_x - 1) / 2 + 1, (dim_z - 1) / 2)] +
                         grid.u[getx(i, (newDimY_x - 1) / 2 + 1, (dim_z - 1) / 2 + 1)] +
                         grid.u[getx(i - 1, (newDimY_x - 1) / 2, (dim_z - 1) / 2)] +
                         grid.u[getx(i - 1, (newDimY_x - 1) / 2, (dim_z - 1) / 2 + 1)] +
                         grid.u[getx(i - 1, (newDimY_x - 1) / 2 + 1, (dim_z - 1) / 2)] +
                         grid.u[getx(i - 1, (newDimY_x - 1) / 2 + 1, (dim_z - 1) / 2 + 1)]) /
                        8;
                }
                v = (grid.v[gety(i, (newDimY_y - 1) / 2, (dim_z - 1) / 2)] +
                     grid.v[gety(i, (newDimY_y - 1) / 2, (dim_z - 1) / 2 + 1)]) /
                    2;
                w = (grid.w[getz(i, (newDimY_z - 1) / 2, (dim_z_z) / 2)] +
                     grid.w[getz(i, (newDimY_z - 1) / 2 + 1, (dim_z_z) / 2)]) /
                    2;
                p = (grid.p[getHaloP(i, (xSize[1] - 0) / 2, (xSize[0] - 1) / 2)] +
                     grid.p[getHaloP(i, (xSize[1] - 0) / 2, (xSize[0] - 1) / 2 + 1)] +
                     grid.p[getHaloP(i, (xSize[1] - 0) / 2 + 1, (xSize[0] - 1) / 2)] +
                     grid.p[getHaloP(i, (xSize[1] - 0) / 2 + 1, (xSize[0] - 1) / 2 + 1)]) /
                    4;
            }

            MPI_File_write_at(fh, offset, &u, 1, MPI_REALL, MPI_STATUS_IGNORE);
            offset += sizeof(Real);
            MPI_File_write_at(fh, offset, &v, 1, MPI_REALL, MPI_STATUS_IGNORE);
            offset += sizeof(Real);
            MPI_File_write_at(fh, offset, &w, 1, MPI_REALL, MPI_STATUS_IGNORE);
            offset += sizeof(Real);
            MPI_File_write_at(fh, offset, &p, 1, MPI_REALL, MPI_STATUS_IGNORE);
            offset += sizeof(Real);
        }
    }

    // LINE 3
    if (testCase == 2)
    {
        if (coords[0] == static_cast<int>(PX / 2) && coords[1] == static_cast<int>(PY / 2))
        {
            offset = (PX * xSize[2] + PY * xSize[1]) * sizeof(Real) * 7;
            const Real xCoord = SX + LX / 2;
            const Real yCoord = SY + LY / 2;
            Real zCoord = SZ;

            Real u, v, w, p;
            int iX, iY, iZ, iP, jX, jY, jZ, jP;
            if (PX % 2 == 0)
            {
                iX = resx;
                iY = 0;
                iZ = 0;
                iP = 0;
            }
            else
            {
                iX = (newDimX_x - 1) / 2;
                iY = (newDimX_y - 1) / 2;
                iZ = (newDimX_z - 1) / 2;
                iP = (xSize[2]) / 2;
            }
            if (PY % 2 == 0)
            {
                jX = 0;
                jY = resy;
                jZ = 0;
                jP = 0;
            }
            else
            {
                jX = (newDimY_x - 1) / 2;
                jY = (newDimY_y - 1) / 2;
                jZ = (newDimY_z - 1) / 2;
                jP = (xSize[1]) / 2;
            }
            for (int k = 0; k < xSize[0]; k++)
            {
                MPI_File_write_at(fh, offset, &xCoord, 1, MPI_REALL, MPI_STATUS_IGNORE);
                offset += sizeof(Real);

                MPI_File_write_at(fh, offset, &yCoord, 1, MPI_REALL, MPI_STATUS_IGNORE);
                offset += sizeof(Real);

                MPI_File_write_at(fh, offset, &zCoord, 1, MPI_REALL, MPI_STATUS_IGNORE);
                zCoord += DZ;
                offset += sizeof(Real);

                u = (grid.u[getx(iX, jX, k)] +
                     grid.u[getx(iX, jX + 1, k)]) /
                    2;
                v = (grid.v[gety(iY, jY, k)] +
                     grid.v[gety(iY + 1, jY, k)]) /
                    2;
                if (k == 0)
                {
                    w = boundary.boundary_value_w[4]->value(xCoord, yCoord, zCoord-DZ/2, DT*Nt);
                }
                else if (k == xSize[0] - 1)
                {
                    w = boundary.boundary_value_w[5]->value(xCoord, yCoord, zCoord-DZ/2, DT*Nt);
                }
                else
                {
                    w = (grid.w[getz(iZ, jZ, k)] +
                         grid.w[getz(iZ, jZ, k - 1)] +
                         grid.w[getz(iZ, jZ + 1, k)] +
                         grid.w[getz(iZ, jZ + 1, k - 1)] +
                         grid.w[getz(iZ + 1, jZ, k)] +
                         grid.w[getz(iZ + 1, jZ, k - 1)] +
                         grid.w[getz(iZ + 1, jZ + 1, k)] +
                         grid.w[getz(iZ + 1, jZ + 1, k - 1)]) /
                        8;
                }
                p = (grid.p[getHaloP(iP, jP, k)] +
                    grid.p[getHaloP(iP, jP + 1, k)] +
                    grid.p[getHaloP(iP + 1, jP, k)] +
                    grid.p[getHaloP(iP + 1, jP + 1, k)]) /
                    4;
                MPI_File_write_at(fh, offset, &u, 1, MPI_REALL, MPI_STATUS_IGNORE);
                offset += sizeof(Real);
                MPI_File_write_at(fh, offset, &v, 1, MPI_REALL, MPI_STATUS_IGNORE);
                offset += sizeof(Real);
                MPI_File_write_at(fh, offset, &w, 1, MPI_REALL, MPI_STATUS_IGNORE);
                offset += sizeof(Real);
                MPI_File_write_at(fh, offset, &p, 1, MPI_REALL, MPI_STATUS_IGNORE);
                offset += sizeof(Real);
            }
        }
    }
    MPI_Barrier(cart_comm);
    MPI_File_close(&fh);

    if (rank == 0)
    {
        int count = 0;
        std::ifstream input(filename, std::ios::binary);
        std::ofstream output(filenametxt);
        output << "Line 1" << std::endl;
        output << "x y z u v w p" << std::endl;
        Real value;
        while (input.read(reinterpret_cast<char *>(&value), sizeof(Real)))
        {
            count++;
            output << value << " ";
            if (count % 7 == 0)
            {
                output << std::endl;
            }
            if (count == ((NY + 1) * 7))
            {
                output << std::endl;
                output << "Line 2" << std::endl;
                output << "x y z u v w p" << std::endl;
            }
            if (count == ((NX + 1) * 7 + (NY + 1) * 7) && testCase == 2)
            {
                output << std::endl;
                output << "Line 3" << std::endl;
                output << "x y z u v w p" << std::endl;
            }
        }

        input.close();
        output.close();
    }
    MPI_Barrier(cart_comm);
}
