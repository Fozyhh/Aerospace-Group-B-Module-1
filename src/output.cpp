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
    const double x_middle_x = NX / 2;
    const double x_middle = (NX + 1) / 2;

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
            << "POINTS " << (NY + 1) * (NZ + 1) << " double\n";

    // Define data format
    header2 << "POINT_DATA " << (NY + 1) * (NZ + 1) << "\n"
            << "SCALARS u double\n"
            << "LOOKUP_TABLE default\n";
    header3 << "SCALARS v double\n"
            << "LOOKUP_TABLE default\n";
    header4 << "SCALARS w double\n"
            << "LOOKUP_TABLE default\n";
    header5 << "SCALARS p double\n"
            << "LOOKUP_TABLE default\n";
    header6 << "SCALARS Magnitude double\n"
            << "LOOKUP_TABLE default\n";

    MPI_Offset offsetheader1 = header1.str().size(),
               offsetheader2 = header2.str().size(),
               offsetheader3 = header3.str().size(),
               offsetheader4 = header4.str().size(),
               offsetheader5 = header5.str().size(),
               offsetheader6 = header6.str().size(),
               offsetpoints = (3 * sizeof(double)),
               offsetallpoints = offsetpoints * ((NY + 1) * (NZ + 1)),
               offsetvalue = sizeof(double),
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
    // int offset_y_p_ = coords[1] * zSize[1] - 1;

    MPI_Barrier(cart_comm);
    if (x_index >= offset_x_x_ + 1 && x_index < dim_x_x + offset_x_x_ + 1)
    {

        int local_x_x = x_index_x - offset_x_x_;
        int local_x_y = x_index - offset_x_y_;
        int local_x_z = x_index - offset_x_z_;
        int local_x_p = x_index - offset_x_p_;

        Real value_x = 0, value_y = 0, value_z = 0, value_p = 0, value_m = 0;
        double bg_px = 0, bg_py = 0, bg_pz = 0, bg_vx = 0, bg_vy = 0, bg_vz = 0, bg_vp = 0, bg_vm = 0;
        double point_x = 0, point_y = 0, point_z = 0;
        for (int j = 1; j < newDimY_x - 1; j++)
        {
            for (int k = 0; k < dim_z; k++)
            {

                // Write grid points coordinate
                point_x = SX + LX / 2;
                point_y = static_cast<double>(SY + (j + offset_y_x_) * DY);
                point_z = static_cast<double>(SZ + (k)*DZ);

                value_x = grid.u[local_x_x * newDimY_x * dim_z + j * dim_z + k];

                if (lby && j == 1)
                {
                    value_y = (boundary.boundary_value_v[2]->value(x_index, j + offset_y_y_ - 0.5, k, T) + boundary.boundary_value_v[2]->value(x_index + 1, j + offset_y_y_ - 0.5, k, T)) / 2;
                }
                else if (rby && j == newDimY_x - 2)
                {
                    value_y = (boundary.boundary_value_v[3]->value(x_index, j + offset_y_y_ - 0.5, k, T) + boundary.boundary_value_v[3]->value(x_index + 1, j + offset_y_y_ - 0.5, k, T)) / 2;
                }
                else
                {
                    value_y = (grid.v[local_x_y * newDimY_y * dim_z + j * dim_z + k] + grid.v[local_x_y * newDimY_y * dim_z + (j + 1) * dim_z + k] +
                               grid.v[(local_x_y + 1) * newDimY_y * dim_z + j * dim_z + k] + grid.v[(local_x_y + 1) * newDimY_y * dim_z + (j + 1) * dim_z + k]) /
                              4;
                }

                if (k == 0)
                {
                    value_z = boundary.boundary_value_w[4]->value(x_index + 0.5, j + offset_y_z_, k - 0.5, T);
                }
                else if (k == dim_z - 1)
                {
                    value_z = boundary.boundary_value_w[5]->value(x_index + 0.5, j + offset_y_z_, k - 0.5, T);
                }
                else
                {
                    value_z = (grid.w[local_x_z * newDimY_z * dim_z_z + j * dim_z_z + k] + grid.w[local_x_z * newDimY_z * dim_z_z + j * dim_z_z + k + 1] +
                               grid.w[(local_x_z + 1) * newDimY_z * dim_z_z + j * dim_z_z + k] + grid.w[(local_x_z + 1) * newDimY_z * dim_z_z + j * dim_z_z + k + 1]) /
                              4;
                }

                value_p = (halo_p[(local_x_p - resx) * (xSize[1] + 2) * xSize[0] + j * xSize[0] + k] + halo_p[(local_x_p - resx + 1) * (xSize[1] + 2) * xSize[0] + j * zSize[0] + k]) / 2;
                value_m = std::sqrt(value_x * value_x + value_y * value_y + value_z * value_z);

                bg_px = to_big_endian(point_x);
                bg_py = to_big_endian(point_y);
                bg_pz = to_big_endian(point_z);
                bg_vx = to_big_endian(value_x);
                bg_vy = to_big_endian(value_y);
                bg_vz = to_big_endian(value_z);
                bg_vp = to_big_endian(value_p);
                bg_vm = to_big_endian(value_m);

                MPI_File_write_at(fh, offsetheader1 + offsetpoints * ((j + offset_y_x_) * dim_z + k), &bg_px, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);

                MPI_File_write_at(fh, offsetheader1 + offsetpoints * ((j + offset_y_x_) * dim_z + k) + sizeof(double), &bg_py, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);

                MPI_File_write_at(fh, offsetheader1 + offsetpoints * ((j + offset_y_x_) * dim_z + k) + 2 * sizeof(double), &bg_pz, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);

                MPI_File_write_at(fh, offsetheader1 + offsetallpoints + offsetheader2 + offsetvalue * ((j + offset_y_x_) * dim_z + k), &bg_vx, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);

                MPI_File_write_at(fh, offsetheader1 + offsetallpoints + offsetheader2 + offsetallvalue + offsetheader3 + offsetvalue * ((j + offset_y_x_) * dim_z + k), &bg_vy, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);

                MPI_File_write_at(fh, offsetheader1 + offsetallpoints + offsetheader2 + 2 * offsetallvalue + offsetheader3 + offsetheader4 + offsetvalue * ((j + offset_y_x_) * dim_z + k), &bg_vz, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);

                MPI_File_write_at(fh, offsetheader1 + offsetallpoints + offsetheader2 + 3 * offsetallvalue + offsetheader3 + offsetheader4 + offsetheader5 + offsetvalue * ((j + offset_y_x_) * dim_z + k), &bg_vp, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);

                MPI_File_write_at(fh, offsetheader1 + offsetallpoints + offsetheader2 + 4 * offsetallvalue + offsetheader3 + offsetheader4 + offsetheader5 + offsetheader6 + offsetvalue * ((j + offset_y_x_) * dim_z + k), &bg_vm, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);
            }
        }
    }
    MPI_Barrier(cart_comm);
    MPI_File_close(&fh);
}

void IcoNS::output_y()
{
    MPI_File fh;
    const double y_middle = (NY + 1) / 2;
    const double y_middle_y = NY / 2;
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
            << "POINTS " << (NX + 1) * (NZ + 1) << " double\n";

    // Define data format
    header2 << "POINT_DATA " << (NX + 1) * (NZ + 1) << "\n"
            << "SCALARS u double\n"
            << "LOOKUP_TABLE default\n";
    header3 << "SCALARS v double\n"
            << "LOOKUP_TABLE default\n";
    header4 << "SCALARS w double\n"
            << "LOOKUP_TABLE default\n";
    header5 << "SCALARS p double\n"
            << "LOOKUP_TABLE default\n";
    header6 << "SCALARS Magnitude double\n"
            << "LOOKUP_TABLE default\n";

    MPI_Offset offsetheader1 = header1.str().size(),
               offsetheader2 = header2.str().size(),
               offsetheader3 = header3.str().size(),
               offsetheader4 = header4.str().size(),
               offsetheader5 = header5.str().size(),
               offsetheader6 = header6.str().size(),
               offsetpoints = (3 * sizeof(double)),
               offsetallpoints = offsetpoints * ((NX + 1) * (NZ + 1)),
               offsetvalue = sizeof(double),
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
    // int offset_x_p_ = coords[0] * zSize[0] - 1;
    int offset_y_p_ = coords[1] * zSize[1] - 1;

    MPI_Barrier(cart_comm);
    if (y_index >= offset_y_y_ && y_index < dim_y_y + offset_y_y_)
    {
        int local_y_x = y_index - offset_y_x_;
        int local_y_y = y_index_y - offset_y_y_;
        int local_y_z = y_index - offset_y_z_;
        int local_y_p = y_index - offset_y_p_;

        Real value_x = 0, value_y = 0, value_z = 0, value_p = 0, value_m = 0;
        double bg_px = 0, bg_py = 0, bg_pz = 0, bg_vx = 0, bg_vy = 0, bg_vz = 0, bg_vp = 0, bg_vm = 0;
        double point_x = 0, point_y = 0, point_z = 0;
        for (int i = 1; i < newDimX_y - 1; i++)
        {
            for (int k = 0; k < dim_z; k++)
            {
                point_x = static_cast<double>(SX + (i + offset_x_y_) * DX);
                point_y = SY + LY / 2;
                point_z = static_cast<double>(SZ + (k)*DZ);

                value_y = grid.v[i * newDimY_y * dim_z + local_y_y * dim_z + k];

                if (lbx && i == 1)
                {
                    value_x = boundary.boundary_value_u[0]->value(i + offset_x_x_ - 0.5, y_index + offset_y_x_, k, T);
                }
                else if (rbx && i == newDimX_y - 2)
                {
                    value_x = boundary.boundary_value_u[1]->value(i + offset_x_x_ - 0.5, y_index + offset_y_x_, k, T);
                }
                else
                {
                    value_x = grid.u[i * newDimY_x * dim_z + local_y_x * dim_z + k];
                }

                if (k == 0)
                {
                    value_z = boundary.boundary_value_w[4]->value(i + offset_x_z_, y_index + 0.5, k - 0.5, T);
                }
                else if (k == dim_z - 1)
                {
                    value_z = boundary.boundary_value_w[5]->value(i + offset_x_z_, y_index + 0.5, k - 0.5, T);
                }
                else
                {
                    value_z = grid.w[i * newDimY_z * dim_z_z + local_y_z * dim_z_z + k];
                }

                value_p = halo_p[i * (xSize[1] + 2) * xSize[0] + (local_y_p - resy) * xSize[0] + k];

                value_m = std::sqrt(value_x * value_x + value_y * value_y + value_z * value_z);

                bg_px = to_big_endian(point_x);
                bg_py = to_big_endian(point_y);
                bg_pz = to_big_endian(point_z);
                bg_vx = to_big_endian(value_x);
                bg_vy = to_big_endian(value_y);
                bg_vz = to_big_endian(value_z);
                bg_vp = to_big_endian(value_p);
                bg_vm = to_big_endian(value_m);
                MPI_File_write_at(fh, offsetheader1 + offsetpoints * ((i + offset_x_y_) * dim_z + k), &bg_px, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);

                MPI_File_write_at(fh, offsetheader1 + offsetpoints * ((i + offset_x_y_) * dim_z + k) + sizeof(double), &bg_py, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);

                MPI_File_write_at(fh, offsetheader1 + offsetpoints * ((i + offset_x_y_) * dim_z + k) + 2 * sizeof(double), &bg_pz, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);

                MPI_File_write_at(fh, offsetheader1 + offsetallpoints + offsetheader2 + offsetvalue * ((i + offset_x_y_) * dim_z + k), &bg_vx, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);

                MPI_File_write_at(fh, offsetheader1 + offsetallpoints + offsetheader2 + offsetallvalue + offsetheader3 + offsetvalue * ((i + offset_x_y_) * dim_z + k), &bg_vy, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);

                MPI_File_write_at(fh, offsetheader1 + offsetallpoints + offsetheader2 + 2 * offsetallvalue + offsetheader3 + offsetheader4 + offsetvalue * ((i + offset_x_y_) * dim_z + k), &bg_vz, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);

                MPI_File_write_at(fh, offsetheader1 + offsetallpoints + offsetheader2 + 3 * offsetallvalue + offsetheader3 + offsetheader4 + offsetheader5 + offsetvalue * ((i + offset_x_y_) * dim_z + k), &bg_vp, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);

                MPI_File_write_at(fh, offsetheader1 + offsetallpoints + offsetheader2 + 4 * offsetallvalue + offsetheader3 + offsetheader4 + offsetheader5 + offsetheader6 + offsetvalue * ((i + offset_x_y_) * dim_z + k), &bg_vm, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);
            }
        }
    }
    MPI_Barrier(cart_comm);
    MPI_File_close(&fh);
}

void IcoNS::output_z()
{
    MPI_File fh;
    const double z_middle_z = (NZ + 1) / 2;
    const double z_middle = NZ / 2;
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
            << "POINTS " << (NY + 1) * (NX + 1) << " double\n";

    // Define data format
    header2 << "POINT_DATA " << (NY + 1) * (NX + 1) << "\n"
            << "SCALARS u double\n"
            << "LOOKUP_TABLE default\n";
    header3 << "SCALARS v double\n"
            << "LOOKUP_TABLE default\n";
    header4 << "SCALARS w double\n"
            << "LOOKUP_TABLE default\n";
    header5 << "SCALARS p double\n"
            << "LOOKUP_TABLE default\n";
    header6 << "SCALARS Magnitude double\n"
            << "LOOKUP_TABLE default\n";

    MPI_Offset offsetheader1 = header1.str().size(),
               offsetheader2 = header2.str().size(),
               offsetheader3 = header3.str().size(),
               offsetheader4 = header4.str().size(),
               offsetheader5 = header5.str().size(),
               offsetheader6 = header6.str().size(),
               offsetpoints = (3 * sizeof(double)),
               offsetallpoints = offsetpoints * ((NY + 1) * (NX + 1)),
               offsetvalue = sizeof(double),
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
            double point_x = static_cast<double>(SX + (i + offset_x_z_) * DX),
                   point_y = static_cast<double>(SY + (j + offset_y_z_) * DY),
                   point_z = SZ + LZ / 2;

            // valuesx[rank*(dim_y_x * dim_z + dim_z) + ((j-1) * dim_z + k)] = grid.v[local_x* newDimY_x * dim_z + j * dim_z + k];
            Real value_x, value_y, value_z, value_p, value_m;
            if (lby && j == 1)
            {
                value_x = (boundary.boundary_value_u[0]->value(i + offset_x_x_, j + offset_y_x_, z_index, T));
            }
            else if (rby && j == newDimY_z - 2)
            {
                value_x = (boundary.boundary_value_u[1]->value(i + offset_x_x_, j + offset_y_x_, z_index, T));
            }
            else
            {
                value_x = (grid.u[i * newDimY_x * dim_z + j * dim_z + z_index] + grid.u[(i - 1) * newDimY_x * dim_z + j * dim_z + z_index] +
                           grid.u[i * newDimY_x * dim_z + j * dim_z + z_index - 1] + grid.u[(i - 1) * newDimY_x * dim_z + j * dim_z + z_index - 1]) /
                          4;
            }

            if (lbx && i == 1)
            {
                value_y = (boundary.boundary_value_v[0]->value(i + offset_x_y_ + 0.5, j + offset_y_y_ - 0.5, z_index, T));
            }
            else if (rbx && i == newDimX_z - 2)
            {
                value_y = boundary.boundary_value_v[4]->value(i + offset_x_y_ + 0.5, j + offset_y_y_ - 0.5, z_index, T);
                ;
            }
            else
            {
                value_y = (grid.v[i * newDimY_y * dim_z + j * dim_z + z_index] + grid.v[i * newDimY_y * dim_z + j * dim_z + z_index + 1] +
                           grid.v[i * newDimY_y * dim_z + (j + 1) * dim_z + z_index] + grid.v[i * newDimY_y * dim_z + (j + 1) * dim_z + z_index + 1]) /
                          4;
            }

            value_z = grid.w[i * newDimY_z * dim_z_z + j * dim_z_z + z_index_z];

            value_p = (halo_p[i * (xSize[1] + 2) * xSize[0] + j * xSize[0] + z_index] + halo_p[i * (xSize[1] + 2) * xSize[0] + j * xSize[0] + z_index + 1]) / 2;

            value_m = std::sqrt(value_x * value_x + value_y * value_y + value_z * value_z);

            double bg_px = to_big_endian(point_x);
            double bg_py = to_big_endian(point_y);
            double bg_pz = to_big_endian(point_z);
            double bg_vx = to_big_endian(value_x);
            double bg_vy = to_big_endian(value_y);
            double bg_vz = to_big_endian(value_z);
            double bg_vp = to_big_endian(value_p);
            double bg_vm = to_big_endian(value_m);

            MPI_File_write_at(fh, offsetheader1 + offsetpoints * ((i + offset_x_z_) * dim_z + j + offset_y_z_), &bg_px, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);

            MPI_File_write_at(fh, offsetheader1 + offsetpoints * ((i + offset_x_z_) * dim_z + j + offset_y_z_) + sizeof(double), &bg_py, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);

            MPI_File_write_at(fh, offsetheader1 + offsetpoints * ((i + offset_x_z_) * dim_z + j + offset_y_z_) + 2 * sizeof(double), &bg_pz, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);

            MPI_File_write_at(fh, offsetheader1 + offsetallpoints + offsetheader2 + offsetvalue * ((i + offset_x_z_) * dim_z + j + offset_y_z_), &bg_vx, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);

            MPI_File_write_at(fh, offsetheader1 + offsetallpoints + offsetheader2 + offsetallvalue + offsetheader3 + offsetvalue * ((i + offset_x_z_) * dim_z + j + offset_y_z_), &bg_vy, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);

            MPI_File_write_at(fh, offsetheader1 + offsetallpoints + offsetheader2 + 2 * offsetallvalue + offsetheader3 + offsetheader4 + offsetvalue * ((i + offset_x_z_) * dim_z + j + offset_y_z_), &bg_vz, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);

            MPI_File_write_at(fh, offsetheader1 + offsetallpoints + offsetheader2 + 3 * offsetallvalue + offsetheader3 + offsetheader4 + offsetheader5 + offsetvalue * ((i + offset_x_z_) * dim_z + j + offset_y_z_), &bg_vp, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);

            MPI_File_write_at(fh, offsetheader1 + offsetallpoints + offsetheader2 + 4 * offsetallvalue + offsetheader3 + offsetheader4 + offsetheader5 + offsetheader6 + offsetvalue * ((i + offset_x_z_) * dim_z + j + offset_y_z_), &bg_vm, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);
        }
    }

    MPI_File_close(&fh);
}

void IcoNS::output_profile()
{
    const std::string filename = "out" + std::to_string(testCase) + ".dat";
    MPI_File fh;
    MPI_File_open(cart_comm, filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
    MPI_Offset offset = coords[1] * xSize[1] * sizeof(double) * 7;
    // LINE 1
    if (coords[0] == (PX - 1) / 2)
    {
        const double xCoord = SX + LX / 2;
        double yCoord = SY + coords[1] * xSize[1] * DY;
        const double zCoord = SZ + LZ / 2;
        double u, v, w, p;
        for (int i = 1; i < xSize[1] + 1; i++)
        {
            MPI_File_write_at(fh, offset, &xCoord, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);
            offset += sizeof(double);

            MPI_File_write_at(fh, offset, &yCoord, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);
            yCoord += DY;
            offset += sizeof(double);

            MPI_File_write_at(fh, offset, &zCoord, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);
            offset += sizeof(double);

            if (PX % 2 == 0)
            {
                u = (grid.u[getx(newDimX_x - 2, i, (dim_z - 1) / 2)] +
                     grid.u[getx(newDimX_x - 2, i, (dim_z - 1) / 2 + 1)]) /
                    2;
                if (lby && i == 1)
                { // approximating at the front boundary, could use boundary_value_v[2]->value()
                    v = (grid.v[gety(newDimX_y - 2, i, (dim_z - 1) / 2)] +
                         grid.v[gety(newDimX_y - 2, i, (dim_z - 1) / 2 + 1)] +
                         grid.v[gety(newDimX_y - 1, i, (dim_z - 1) / 2)] +
                         grid.v[gety(newDimX_y - 1, i, (dim_z - 1) / 2 + 1)]) /
                        4;
                }
                else if (rby && i == xSize[1])
                { // approximating at the back boundary, could use boundary_value_v[2]->value()
                    v = (grid.v[gety(newDimX_y - 2, i - 1, (dim_z - 1) / 2)] +
                         grid.v[gety(newDimX_y - 2, i - 1, (dim_z - 1) / 2 + 1)] +
                         grid.v[gety(newDimX_y - 1, i - 1, (dim_z - 1) / 2)] +
                         grid.v[gety(newDimX_y - 1, i - 1, (dim_z - 1) / 2 + 1)]) /
                        4;
                }
                else
                {
                    v = (grid.v[gety(newDimX_y - 2, i, (dim_z - 1) / 2)] +
                         grid.v[gety(newDimX_y - 2, i, (dim_z - 1) / 2 + 1)] +
                         grid.v[gety(newDimX_y - 2, i - 1, (dim_z - 1) / 2)] +
                         grid.v[gety(newDimX_y - 2, i - 1, (dim_z - 1) / 2 + 1)] +
                         grid.v[gety(newDimX_y - 1, i, (dim_z - 1) / 2)] +
                         grid.v[gety(newDimX_y - 1, i, (dim_z - 1) / 2 + 1)] +
                         grid.v[gety(newDimX_y - 1, i - 1, (dim_z - 1) / 2)] +
                         grid.v[gety(newDimX_y - 1, i - 1, (dim_z - 1) / 2 + 1)]) /
                        8;
                }
                w = (grid.w[getz(newDimX_z - 2, i, (dim_z_z) / 2)] +
                     grid.w[getz(newDimX_z - 1, i, (dim_z_z) / 2)]) /
                    2;
                p = (grid.p[getp(xSize[2] - 1, i, (xSize[0] - 1) / 2)] +
                     grid.p[getp(xSize[2] - 1, i, (xSize[0] - 1) / 2 + 1)]) /
                    2;
                // approximation, should use halo
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
                    v = (grid.v[gety((newDimX_y - 1) / 2, i, (dim_z - 1) / 2)] +
                         grid.v[gety((newDimX_y - 1) / 2, i, (dim_z - 1) / 2 + 1)] +
                         grid.v[gety((newDimX_y - 1) / 2 + 1, i, (dim_z - 1) / 2)] +
                         grid.v[gety((newDimX_y - 1) / 2 + 1, i, (dim_z - 1) / 2 + 1)]) /
                        4;
                }
                else if (rby && i == xSize[1])
                { // approximating at the back boundary, could use boundary_value_v[3]->value()
                    v = (grid.v[gety((newDimX_y - 1) / 2, i - 1, (dim_z - 1) / 2)] +
                         grid.v[gety((newDimX_y - 1) / 2, i - 1, (dim_z - 1) / 2 + 1)] +
                         grid.v[gety((newDimX_y - 1) / 2 + 1, i - 1, (dim_z - 1) / 2)] +
                         grid.v[gety((newDimX_y - 1) / 2 + 1, i - 1, (dim_z - 1) / 2 + 1)]) /
                        4;
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
                p = (grid.p[getp((xSize[2] - 1) / 2, i, (xSize[0] - 1) / 2)] +
                     grid.p[getp((xSize[2] - 1) / 2, i, (xSize[0] - 1) / 2 + 1)] +
                     grid.p[getp((xSize[2] - 1) / 2 + 1, i, (xSize[0] - 1) / 2)] +
                     grid.p[getp((xSize[2] - 1) / 2 + 1, i, (xSize[0] - 1) / 2 + 1)]) /
                    2;
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
    if (coords[1] == (PY - 1) / 2)
    {
        offset = PY * xSize[1] * sizeof(double) * 7 + coords[0] * xSize[2] * sizeof(double) * 7;
        double xCoord = SX + coords[0] * xSize[2] * DX;
        const double yCoord = SY + LY / 2;
        const double zCoord = SZ + LZ / 2;
        double u, v, w, p;
        for (int i = 1; i < xSize[2] + 1; i++)
        {
            MPI_File_write_at(fh, offset, &xCoord, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);
            xCoord += DX;
            offset += sizeof(double);

            MPI_File_write_at(fh, offset, &yCoord, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);
            offset += sizeof(double);

            MPI_File_write_at(fh, offset, &zCoord, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);
            offset += sizeof(double);

            if (PY % 2 == 0)
            {
                if (lbx && i == 1)
                { // approximating at the left boundary, could use boundary_value_u[0]->value()
                    u = (grid.u[getx(i, newDimY_x - 2, (dim_z - 1) / 2)] +
                         grid.u[getx(i, newDimY_x - 2, (dim_z - 1) / 2 + 1)] +
                         grid.u[getx(i, newDimY_x - 1, (dim_z - 1) / 2)] +
                         grid.u[getx(i, newDimY_x - 1, (dim_z - 1) / 2 + 1)]) /
                        4;
                }
                else if (rbx && i == xSize[2])
                { // approximating at the right boundary, could use boundary_value_u[1]->value()
                    u = (grid.u[getx(i - 1, newDimY_x - 2, (dim_z - 1) / 2)] +
                         grid.u[getx(i - 1, newDimY_x - 2, (dim_z - 1) / 2 + 1)] +
                         grid.u[getx(i - 1, newDimY_x - 1, (dim_z - 1) / 2)] +
                         grid.u[getx(i - 1, newDimY_x - 1, (dim_z - 1) / 2 + 1)]) /
                        4;
                }
                else
                {
                    u = (grid.u[getx(i, newDimY_x - 2, (dim_z - 1) / 2)] +
                         grid.u[getx(i, newDimY_x - 2, (dim_z - 1) / 2 + 1)] +
                         grid.u[getx(i, newDimY_x - 1, (dim_z - 1) / 2)] +
                         grid.u[getx(i, newDimY_x - 1, (dim_z - 1) / 2 + 1)] +
                         grid.u[getx(i - 1, newDimY_x - 2, (dim_z - 1) / 2)] +
                         grid.u[getx(i - 1, newDimY_x - 2, (dim_z - 1) / 2 + 1)] +
                         grid.u[getx(i - 1, newDimY_x - 1, (dim_z - 1) / 2)] +
                         grid.u[getx(i - 1, newDimY_x - 1, (dim_z - 1) / 2 + 1)]) /
                        8;
                }
                v = (grid.v[gety(i, newDimY_y - 2, (dim_z - 1) / 2)] +
                     grid.v[gety(i, newDimY_y - 2, (dim_z - 1) / 2 + 1)]) /
                    2;
                w = (grid.w[getz(i, newDimY_z - 2, (dim_z_z) / 2)] +
                     grid.w[getz(i, newDimY_z - 1, (dim_z_z) / 2)]) /
                    2;
                p = (grid.p[getp(i, xSize[1] - 1, (xSize[0] - 1) / 2)] +
                     grid.p[getp(i, xSize[1] - 1, (xSize[0] - 1) / 2 + 1)]) /
                    2;
                // approximation, should use halo
            }
            else
            { // same as above
                if (lbx && i == 1)
                { // approximating at the left boundary, could use boundary_value_u[0]->value()
                    u = (grid.u[getx(i, (newDimY_x - 1) / 2, (dim_z - 1) / 2)] +
                         grid.u[getx(i, (newDimY_x - 1) / 2, (dim_z - 1) / 2 + 1)] +
                         grid.u[getx(i, (newDimY_x - 1) / 2 + 1, (dim_z - 1) / 2)] +
                         grid.u[getx(i, (newDimY_x - 1) / 2 + 1, (dim_z - 1) / 2 + 1)]) /
                        4;
                }
                else if (rbx && i == xSize[2])
                { // approximating at the right boundary, could use boundary_value_u[1]->value()
                    u = (grid.u[getx(i - 1, (newDimY_x - 1) / 2, (dim_z - 1) / 2)] +
                         grid.u[getx(i - 1, (newDimY_x - 1) / 2, (dim_z - 1) / 2 + 1)] +
                         grid.u[getx(i - 1, (newDimY_x - 1) / 2 + 1, (dim_z - 1) / 2)] +
                         grid.u[getx(i - 1, (newDimY_x - 1) / 2 + 1, (dim_z - 1) / 2 + 1)]) /
                        4;
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
                p = (grid.p[getp(i, (xSize[1] - 1) / 2, (xSize[0] - 1) / 2)] +
                     grid.p[getp(i, (xSize[1] - 1) / 2, (xSize[0] - 1) / 2 + 1)] +
                     grid.p[getp(i, (xSize[1] - 1) / 2 + 1, (xSize[0] - 1) / 2)] +
                     grid.p[getp(i, (xSize[1] - 1) / 2 + 1, (xSize[0] - 1) / 2 + 1)]) /
                    4;
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
    if (testCase == 2)
    {
        if (coords[0] == (PX - 1) / 2 && coords[1] == (PY - 1) / 2)
        {
            offset = (PX * xSize[2] + PY * xSize[1]) * sizeof(double) * 7;
            const double xCoord = SX + LX / 2;
            const double yCoord = SY + LY / 2;
            double zCoord = SZ;

            double u, v, w, p;
            int iX, iY, iZ, iP, jX, jY, jZ, jP;
            if (PX % 2 == 0)
            {
                iX = newDimX_x - 2;
                iY = newDimX_y - 2;
                iZ = newDimX_z - 2;
                iP = xSize[2] - 1;
            }
            else
            {
                iX = (newDimX_x - 1) / 2;
                iY = (newDimX_y - 1) / 2;
                iZ = (newDimX_z - 1) / 2;
                iP = (xSize[2] - 1) / 2;
            }
            if (PY % 2 == 0)
            {
                jX = newDimY_x - 2;
                jY = newDimY_y - 2;
                jZ = newDimY_z - 2;
                jP = xSize[1] - 1;
            }
            else
            {
                jX = (newDimY_x - 1) / 2;
                jY = (newDimY_y - 1) / 2;
                jZ = (newDimY_z - 1) / 2;
                jP = (xSize[1] - 1) / 2;
            }
            for (int k = 0; k < xSize[0]; k++)
            {
                MPI_File_write_at(fh, offset, &xCoord, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);
                offset += sizeof(double);

                MPI_File_write_at(fh, offset, &yCoord, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);
                offset += sizeof(double);

                MPI_File_write_at(fh, offset, &zCoord, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);
                zCoord += DZ;
                offset += sizeof(double);

                u = (grid.u[getx(iX, jX, k)] +
                     grid.u[getx(iX, jX + 1, k)]) /
                    2;
                v = (grid.v[gety(iY, jY, k)] +
                     grid.v[gety(iY + 1, jY, k)]) /
                    2;
                if (k == 0)
                { // approximating at the lower boundary, could use boundary_value_w[4]->value()
                    w = (grid.w[getz(iZ, jZ, k)] +
                         grid.w[getz(iZ, jZ + 1, k)] +
                         grid.w[getz(iZ + 1, jZ, k)] +
                         grid.w[getz(iZ + 1, jZ + 1, k)]) /
                        4;
                }
                else if (k == xSize[0] - 1)
                { // approximating at the upper boundary, could use boundary_value_w[5]->value()
                    w = (grid.w[getz(iZ, jZ, k - 1)] +
                         grid.w[getz(iZ, jZ + 1, k - 1)] +
                         grid.w[getz(iZ + 1, jZ, k - 1)] +
                         grid.w[getz(iZ + 1, jZ + 1, k - 1)]) /
                        4;
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
                if (PX % 2 == 0 && PY % 2 == 0)
                { // no halo, approximate
                    p = grid.p[getp(iP, jP, k)];
                }
                else if (PX % 2 == 0 && PY % 2 == 1)
                { // no halo on x, use y
                    p = (grid.p[getp(iP, jP, k)] +
                         grid.p[getp(iP, jP + 1, k)]) /
                        2;
                }
                else if (PX % 2 == 1 && PY % 2 == 0)
                { // no halo on y, use x
                    p = (grid.p[getp(iP, jP, k)] +
                         grid.p[getp(iP + 1, jP, k)]) /
                        2;
                }
                else
                { // use both
                    p = (grid.p[getp(iP, jP, k)] +
                         grid.p[getp(iP, jP + 1, k)] +
                         grid.p[getp(iP + 1, jP, k)] +
                         grid.p[getp(iP + 1, jP + 1, k)]) /
                        4;
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
    }
    MPI_Barrier(cart_comm);
    MPI_File_close(&fh);

    if (rank == 0)
    {
        int count = 0;
        std::ifstream input(filename, std::ios::binary);
        std::ofstream output("profile" + std::to_string(testCase) + ".dat");
        output << "Line 1" << std::endl;
        output << "x y z u v w p" << std::endl;
        double value;
        while (input.read(reinterpret_cast<char *>(&value), sizeof(double)))
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
        std::remove(filename.c_str());
    }
    MPI_Barrier(cart_comm);
}
