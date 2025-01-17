#include "core.hpp"


Real IcoNS::error_comp_X(const Real t)
{
    Real error = 0.0;
    int offset_x = offset_x_x - 1;
    int offset_y = offset_y_x - 1;
    // first slice (left face)
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
    int offset_x = offset_x_y - 1;
    int offset_y = offset_y_y - 1;
    // first slice (left face)
    if (lbx)
    {
        if (lby)
        {
            error += ((grid.v[(newDimY_y + 1) * dim_z] - exact_solution.value_y(0, 0.5, 0, t)) *
                      (grid.v[(newDimY_y + 1) * dim_z] - exact_solution.value_y(0, 0.5, 0, t)) *
                      DX * DY * DZ / 8);

            for (int k = 1; k < NZ; k++)
            {
                error += ((grid.v[(newDimY_y + 1) * dim_z + k] - exact_solution.value_y(0, 0.5, k, t)) *
                          (grid.v[(newDimY_y + 1) * dim_z + k] - exact_solution.value_y(0, 0.5, k, t)) *
                          DX * DY * DZ / 4);
            }

            error += ((grid.v[(newDimY_y + 1) * dim_z + NZ] - exact_solution.value_y(0, 0.5, NZ, t)) *
                      (grid.v[(newDimY_y + 1) * dim_z + NZ] - exact_solution.value_y(0, 0.5, NZ, t)) *
                      DX * DY * DZ / 8);
        }
        for (int j = 1 + lby; j < newDimY_y - 1 - rby; j++)
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
        if (rby)
        {
            error += ((grid.v[(newDimY_y)*dim_z + (newDimY_y - 2) * (NZ + 1)] - exact_solution.value_y(0, (NY - 0.5), 0, t)) *
                      (grid.v[(newDimY_y)*dim_z + (newDimY_y - 2) * (NZ + 1)] - exact_solution.value_y(0, (NY - 0.5), 0, t)) *
                      DX * DY * DZ / 8);

            for (int k = 1; k < NZ; k++)
            {
                error += ((grid.v[(newDimY_y)*dim_z + (newDimY_y - 2) * (NZ + 1) + k] - exact_solution.value_y(0, (NY - 0.5), k, t)) *
                          (grid.v[(newDimY_y)*dim_z + (newDimY_y - 2) * (NZ + 1) + k] - exact_solution.value_y(0, (NY - 0.5), k, t)) *
                          DX * DY * DZ / 4);
            }

            error += ((grid.v[(newDimY_y)*dim_z + (newDimY_y - 2) * (NZ + 1) + NZ] - exact_solution.value_y(0, (NY - 0.5), NZ, t)) *
                      (grid.v[(newDimY_y)*dim_z + (newDimY_y - 2) * (NZ + 1) + NZ] - exact_solution.value_y(0, (NY - 0.5), NZ, t)) *
                      DX * DY * DZ / 8);
        }
    }

    // middle slices
    {
        for (int i = 1 + lbx; i < newDimX_y - 1 - rbx; i++)
        {
            if (lby)
            {
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
            if (rby)
            {
                error += ((grid.v[i * newDimY_y * (NZ + 1) + (newDimY_y - 2) * (NZ + 1)] - exact_solution.value_y(i + offset_x, (NY - 0.5), 0, t)) *
                          (grid.v[i * newDimY_y * (NZ + 1) + (newDimY_y - 2) * (NZ + 1)] - exact_solution.value_y(i + offset_x, (NY - 0.5), 0, t)) *
                          DX * DY * DZ / 4);

                for (int k = 1; k < NZ; k++)
                {
                    error += ((grid.v[i * newDimY_y * (NZ + 1) + (newDimY_y - 2) * (NZ + 1) + k] - exact_solution.value_y(i + offset_x, (NY - 0.5), k, t)) *
                              (grid.v[i * newDimY_y * (NZ + 1) + (newDimY_y - 2) * (NZ + 1) + k] - exact_solution.value_y(i + offset_x, (NY - 0.5), k, t)) *
                              DX * DY * DZ / 2);
                }

                error += ((grid.v[i * newDimY_y * (NZ + 1) + (newDimY_y - 2) * (NZ + 1) + NZ] - exact_solution.value_y(i + offset_x, (NY - 0.5), NZ, t)) *
                          (grid.v[i * newDimY_y * (NZ + 1) + (newDimY_y - 2) * (NZ + 1) + NZ] - exact_solution.value_y(i + offset_x, (NY - 0.5), NZ, t)) *
                          DX * DY * DZ / 4);
            }
        }
    }

    // last slice (right face)
    if (rbx)
    {
        if (lby)
        {
            error += ((grid.v[(newDimX_y - 2) * newDimY_y * (NZ + 1) + dim_z] - exact_solution.value_y(NX, 0.5, 0, t)) *
                      (grid.v[(newDimX_y - 2) * newDimY_y * (NZ + 1) + dim_z] - exact_solution.value_y(NX, 0.5, 0, t)) *
                      DX * DY * DZ / 8);

            for (int k = 1; k < NZ; k++)
            {
                error += ((grid.v[(newDimX_y - 2) * newDimY_y * (NZ + 1) + dim_z + k] - exact_solution.value_y(NX, 0.5, k, t)) *
                          (grid.v[(newDimX_y - 2) * newDimY_y * (NZ + 1) + dim_z + k] - exact_solution.value_y(NX, 0.5, k, t)) *
                          DX * DY * DZ / 4);
            }

            error += ((grid.v[(newDimX_y - 2) * newDimY_y * (NZ + 1) + dim_z + NZ] - exact_solution.value_y(NX, 0.5, NZ, t)) *
                      (grid.v[(newDimX_y - 2) * newDimY_y * (NZ + 1) + dim_z + NZ] - exact_solution.value_y(NX, 0.5, NZ, t)) *
                      DX * DY * DZ / 8);
        }
        for (int j = 1 + lby; j < newDimY_y - 1 - rby; j++)
        {
            error += ((grid.v[(newDimX_y - 2) * newDimY_y * (NZ + 1) + j * (NZ + 1)] - exact_solution.value_y(NX, j + 0.5 + offset_y, 0, t)) *
                      (grid.v[(newDimX_y - 2) * newDimY_y * (NZ + 1) + j * (NZ + 1)] - exact_solution.value_y(NX, j + 0.5 + offset_y, 0, t)) *
                      DX * DY * DZ / 4);
            for (int k = 1; k < NZ; k++)
            {
                error += ((grid.v[(newDimX_y - 2) * newDimY_y * (NZ + 1) + j * (NZ + 1) + k] - exact_solution.value_y(NX, j + 0.5 + offset_y, k, t)) *
                          (grid.v[(newDimX_y - 2) * newDimY_y * (NZ + 1) + j * (NZ + 1) + k] - exact_solution.value_y(NX, j + 0.5 + offset_y, k, t)) *
                          DX * DY * DZ / 2);
            }
            error += ((grid.v[(newDimX_y - 2) * newDimY_y * (NZ + 1) + j * (NZ + 1) + NZ] - exact_solution.value_y(NX, j + 0.5 + offset_y, NZ, t)) *
                      (grid.v[(newDimX_y - 2) * newDimY_y * (NZ + 1) + j * (NZ + 1) + NZ] - exact_solution.value_y(NX, j + 0.5 + offset_y, NZ, t)) *
                      DX * DY * DZ / 4);
        }
        if (rby)
        {
            error += ((grid.v[(newDimX_y - 2) * newDimY_y * (NZ + 1) + (newDimY_y - 2) * (NZ + 1)] - exact_solution.value_y(NX, (NY - 0.5), 0, t)) *
                      (grid.v[(newDimX_y - 2) * newDimY_y * (NZ + 1) + (newDimY_y - 2) * (NZ + 1)] - exact_solution.value_y(NX, (NY - 0.5), 0, t)) *
                      DX * DY * DZ / 8);

            for (int k = 1; k < NZ; k++)
            {
                error += ((grid.v[(newDimX_y - 2) * newDimY_y * (NZ + 1) + (newDimY_y - 2) * (NZ + 1) + k] - exact_solution.value_y(NX, (NY - 0.5), k, t)) *
                          (grid.v[(newDimX_y - 2) * newDimY_y * (NZ + 1) + (newDimY_y - 2) * (NZ + 1) + k] - exact_solution.value_y(NX, (NY - 0.5), k, t)) *
                          DX * DY * DZ / 4);
            }

            error += ((grid.v[(newDimX_y - 2) * newDimY_y * (NZ + 1) + (newDimY_y - 2) * (NZ + 1) + NZ] - exact_solution.value_y(NX, (NY - 0.5), NZ, t)) *
                      (grid.v[(newDimX_y - 2) * newDimY_y * (NZ + 1) + (newDimY_y - 2) * (NZ + 1) + NZ] - exact_solution.value_y(NX, (NY - 0.5), NZ, t)) *
                      DX * DY * DZ / 8);
        }
    }

    return error;
}

Real IcoNS::error_comp_Z(const Real t)
{
    Real error = 0.0;
    int offset_x = offset_x_z - 1;
    int offset_y = offset_y_z - 1;
    // first slice (left face)
    if (lbx)
    {
        if (lby)
        {
            error += ((grid.w[(newDimY_z + 1) * dim_z_z] - exact_solution.value_z(0, 0, 0.5, t)) *
                      (grid.w[(newDimY_z + 1) * dim_z_z] - exact_solution.value_z(0, 0, 0.5, t)) *
                      DX * DY * DZ / 8);

            for (int k = 1; k < NZ - 1; k++)
            {
                error += ((grid.w[(newDimY_z + 1) * dim_z_z + k] - exact_solution.value_z(0, 0, k + 0.5, t)) *
                          (grid.w[(newDimY_z + 1) * dim_z_z + k] - exact_solution.value_z(0, 0, k + 0.5, t)) *
                          DX * DY * DZ / 4);
            }

            error += ((grid.w[(newDimY_z + 1) * dim_z_z + NZ - 1] - exact_solution.value_z(0, 0, NZ - 0.5, t)) *
                      (grid.w[(newDimY_z + 1) * dim_z_z + NZ - 1] - exact_solution.value_z(0, 0, NZ - 0.5, t)) *
                      DX * DY * DZ / 8);
        }
        for (int j = 1 + lby; j < newDimY_z - 1 - rby; j++)
        {
            error += ((grid.w[(newDimY_z)*dim_z_z + j * NZ] - exact_solution.value_z(0, j + offset_y, 0.5, t)) *
                      (grid.w[(newDimY_z)*dim_z_z + j * NZ] - exact_solution.value_z(0, j + offset_y, 0.5, t)) *
                      DX * DY * DZ / 4);
            for (int k = 1; k < NZ - 1; k++)
            {
                error += ((grid.w[(newDimY_z)*dim_z_z + j * NZ + k] - exact_solution.value_z(0, j + offset_y, k + 0.5, t)) *
                          (grid.w[(newDimY_z)*dim_z_z + j * NZ + k] - exact_solution.value_z(0, j + offset_y, k + 0.5, t)) *
                          DX * DY * DZ / 2);
            }
            error += ((grid.w[(newDimY_z)*dim_z_z + j * NZ + NZ - 1] - exact_solution.value_z(0, j + offset_y, NZ - 0.5, t)) *
                      (grid.w[(newDimY_z)*dim_z_z + j * NZ + NZ - 1] - exact_solution.value_z(0, j + offset_y, NZ - 0.5, t)) *
                      DX * DY * DZ / 4);
        }
        if (rby)
        {
            error += ((grid.w[(newDimY_z)*dim_z_z + (newDimY_z - 2) * dim_z_z] - exact_solution.value_z(0, NY, 0.5, t)) *
                      (grid.w[(newDimY_z)*dim_z_z + (newDimY_z - 2) * dim_z_z] - exact_solution.value_z(0, NY, 0.5, t)) *
                      DX * DY * DZ / 8);

            for (int k = 1; k < NZ - 1; k++)
            {
                error += ((grid.w[(newDimY_z)*dim_z_z + (newDimY_z - 2) * dim_z_z + k] - exact_solution.value_z(0, NY, k + 0.5, t)) *
                          (grid.w[(newDimY_z)*dim_z_z + (newDimY_z - 2) * dim_z_z + k] - exact_solution.value_z(0, NY, k + 0.5, t)) *
                          DX * DY * DZ / 4);
            }

            error += ((grid.w[(newDimY_z)*dim_z_z + (newDimY_z - 2) * dim_z_z + NZ - 1] - exact_solution.value_z(0, NY, NZ - 0.5, t)) *
                      (grid.w[(newDimY_z)*dim_z_z + (newDimY_z - 2) * dim_z_z + NZ - 1] - exact_solution.value_z(0, NY, NZ - 0.5, t)) *
                      DX * DY * DZ / 8);
        }
    }

    // middle slices
    {
        for (int i = 1 + lbx; i < newDimX_z - 1 - rbx; i++)
        {
            if (lby)
            {
                error += ((grid.w[i * newDimY_z * NZ + dim_z_z] - exact_solution.value_z(i + offset_x, 0, 0.5, t)) *
                          (grid.w[i * newDimY_z * NZ + dim_z_z] - exact_solution.value_z(i + offset_x, 0, 0.5, t)) *
                          DX * DY * DZ / 4);

                for (int k = 1; k < NZ - 1; k++)
                {
                    error += ((grid.w[i * newDimY_z * NZ + dim_z_z + k] - exact_solution.value_z(i + offset_x, 0, k + 0.5, t)) *
                              (grid.w[i * newDimY_z * NZ + dim_z_z + k] - exact_solution.value_z(i + offset_x, 0, k + 0.5, t)) *
                              DX * DY * DZ / 2);
                }
                error += ((grid.w[i * newDimY_z * NZ + dim_z_z + NZ - 1] - exact_solution.value_z(i + offset_x, 0, NZ - 0.5, t)) *
                          (grid.w[i * newDimY_z * NZ + dim_z_z + NZ - 1] - exact_solution.value_z(i + offset_x, 0, NZ - 0.5, t)) *
                          DX * DY * DZ / 4);
            }
            for (int j = 1 + lby; j < newDimY_z - 1 - rby; j++)
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
            if (rby)
            {
                error += ((grid.w[i * newDimY_z * NZ + (newDimY_z - 2) * NZ] - exact_solution.value_z(i + offset_x, NY, 0.5, t)) *
                          (grid.w[i * newDimY_z * NZ + (newDimY_z - 2) * NZ] - exact_solution.value_z(i + offset_x, NY, 0.5, t)) *
                          DX * DY * DZ / 4);

                for (int k = 1; k < NZ - 1; k++)
                {
                    error += ((grid.w[i * newDimY_z * NZ + (newDimY_z - 2) * NZ + k] - exact_solution.value_z(i + offset_x, NY, k + 0.5, t)) *
                              (grid.w[i * newDimY_z * NZ + (newDimY_z - 2) * NZ + k] - exact_solution.value_z(i + offset_x, NY, k + 0.5, t)) *
                              DX * DY * DZ / 2);
                }

                error += ((grid.w[i * newDimY_z * NZ + (newDimY_z - 2) * NZ + NZ - 1] - exact_solution.value_z(i + offset_x, NY, NZ - 0.5, t)) *
                          (grid.w[i * newDimY_z * NZ + (newDimY_z - 2) * NZ + NZ - 1] - exact_solution.value_z(i + offset_x, NY, NZ - 0.5, t)) *
                          DX * DY * DZ / 4);
            }
        }
    }

    // last slice (right face)
    if (rbx)
    {
        if (lby)
        {
            error += ((grid.w[(newDimX_z - 2) * (newDimY_z)*NZ + dim_z_z] - exact_solution.value_z(NX, 0, 0.5, t)) *
                      (grid.w[(newDimX_z - 2) * (newDimY_z)*NZ + dim_z_z] - exact_solution.value_z(NX, 0, 0.5, t)) *
                      DX * DY * DZ / 8);

            for (int k = 1; k < NZ - 1; k++)
            {
                error += ((grid.w[(newDimX_z - 2) * (newDimY_z)*NZ + dim_z_z + k] - exact_solution.value_z(NX, 0, k + 0.5, t)) *
                          (grid.w[(newDimX_z - 2) * (newDimY_z)*NZ + dim_z_z + k] - exact_solution.value_z(NX, 0, k + 0.5, t)) *
                          DX * DY * DZ / 4);
            }

            error += ((grid.w[(newDimX_z - 2) * (newDimY_z)*NZ + dim_z_z + NZ - 1] - exact_solution.value_z(NX, 0, NZ - 0.5, t)) *
                      (grid.w[(newDimX_z - 2) * (newDimY_z)*NZ + dim_z_z + NZ - 1] - exact_solution.value_z(NX, 0, NZ - 0.5, t)) *
                      DX * DY * DZ / 8);
        }
        for (int j = 1 + lby; j < newDimY_z - 1 - rby; j++)
        {
            error += ((grid.w[(newDimX_z - 2) * (newDimY_z)*NZ + j * NZ] - exact_solution.value_z(NX, j + offset_y, 0.5, t)) *
                      (grid.w[(newDimX_z - 2) * (newDimY_z)*NZ + j * NZ] - exact_solution.value_z(NX, j + offset_y, 0.5, t)) *
                      DX * DY * DZ / 4);
            for (int k = 1; k < NZ - 1; k++)
            {
                error += ((grid.w[(newDimX_z - 2) * (newDimY_z)*NZ + j * NZ + k] - exact_solution.value_z(NX, j + offset_y, k + 0.5, t)) *
                          (grid.w[(newDimX_z - 2) * (newDimY_z)*NZ + j * NZ + k] - exact_solution.value_z(NX, j + offset_y, k + 0.5, t)) *
                          DX * DY * DZ / 2);
            }
            error += ((grid.w[(newDimX_z - 2) * (newDimY_z)*NZ + j * NZ + NZ - 1] - exact_solution.value_z(NX, j + offset_y, NZ - 0.5, t)) *
                      (grid.w[(newDimX_z - 2) * (newDimY_z)*NZ + j * NZ + NZ - 1] - exact_solution.value_z(NX, j + offset_y, NZ - 0.5, t)) *
                      DX * DY * DZ / 4);
        }
        if (rby)
        {
            error += ((grid.w[(newDimX_z - 2) * (newDimY_z)*NZ + (newDimY_z - 2) * NZ] - exact_solution.value_z(NX, NY, 0.5, t)) *
                      (grid.w[(newDimX_z - 2) * (newDimY_z)*NZ + (newDimY_z - 2) * NZ] - exact_solution.value_z(NX, NY, 0.5, t)) *
                      DX * DY * DZ / 8);

            for (int k = 1; k < NZ - 1; k++)
            {
                error += ((grid.w[(newDimX_z - 2) * (newDimY_z)*NZ + (newDimY_z - 2) * NZ + k] - exact_solution.value_z(NX, NY, k + 0.5, t)) *
                          (grid.w[(newDimX_z - 2) * (newDimY_z)*NZ + (newDimY_z - 2) * NZ + k] - exact_solution.value_z(NX, NY, k + 0.5, t)) *
                          DX * DY * DZ / 4);
            }

            error += ((grid.w[(newDimX_z - 2) * (newDimY_z)*NZ + (newDimY_z - 2) * NZ + NZ - 1] - exact_solution.value_z(NX, NY, NZ - 0.5, t)) *
                      (grid.w[(newDimX_z - 2) * (newDimY_z)*NZ + (newDimY_z - 2) * NZ + NZ - 1] - exact_solution.value_z(NX, NY, NZ - 0.5, t)) *
                      DX * DY * DZ / 8);
        }
    }

    return error;
}

Real IcoNS::error_comp_P(const Real t)
{
    Real error = 0.0;
    // TODO: check as im not sure
    int offset_x = coords[0] * xSize[2] + std::max(0, coords[0] - (PX - (NX + 1) % PX));
    int offset_y = coords[1] * xSize[1] + std::max(0, coords[1] - (PY - (NY + 1) % PY));
    // first slice (left face)
    if (lbx)
    {
        if (lby)
        {
            error += ((grid.p[getp(0, 0, 0)] - exact_solution.value_p(0, 0, 0, t)) *
                      (grid.p[getp(0, 0, 0)] - exact_solution.value_p(0, 0, 0, t)) *
                      DX * DY * DZ / 8);
            for (int k = 1; k < xSize[0] - 1; k++)
            {
                error += ((grid.p[getp(0, 0, k)] - exact_solution.value_p(0, 0, k, t)) *
                          (grid.p[getp(0, 0, k)] - exact_solution.value_p(0, 0, k, t)) *
                          DX * DY * DZ / 4);
            }
            error += ((grid.p[getp(0, 0, xSize[0] - 1)] - exact_solution.value_p(0, 0, NZ, t)) *
                      (grid.p[getp(0, 0, xSize[0] - 1)] - exact_solution.value_p(0, 0, NZ, t)) *
                      DX * DY * DZ / 8);
        }

        for (int j = lby; j < xSize[1] - rby; j++)
        {
            error += ((grid.p[getp(0, j, 0)] - exact_solution.value_p(0, j + offset_y, 0, t)) *
                      (grid.p[getp(0, j, 0)] - exact_solution.value_p(0, j + offset_y, 0, t)) *
                      DX * DY * DZ / 4);
            for (int k = 1; k < xSize[0] - 1; k++)
            {
                error += ((grid.p[getp(0, j, k)] - exact_solution.value_p(0, j + offset_y, k, t)) *
                          (grid.p[getp(0, j, k)] - exact_solution.value_p(0, j + offset_y, k, t)) *
                          DX * DY * DZ / 2);
            }
            error += ((grid.p[getp(0, j, xSize[0] - 1)] - exact_solution.value_p(0, j + offset_y, NZ, t)) *
                      (grid.p[getp(0, j, xSize[0] - 1)] - exact_solution.value_p(0, j + offset_y, NZ, t)) *
                      DX * DY * DZ / 4);
        }
        if (rby)
        {
            error += ((grid.p[getp(0, xSize[1] - 1, 0)] - exact_solution.value_p(0, NY, 0, t)) *
                      (grid.p[getp(0, xSize[1] - 1, 0)] - exact_solution.value_p(0, NY, 0, t)) *
                      DX * DY * DZ / 8);

            for (int k = 1; k < xSize[0] - 1; k++)
            {
                error += ((grid.p[getp(0, xSize[1] - 1, k)] - exact_solution.value_p(0, NY, k, t)) *
                          (grid.p[getp(0, xSize[1] - 1, k)] - exact_solution.value_p(0, NY, k, t)) *
                          DX * DY * DZ / 4);
            }

            error += ((grid.p[getp(0, xSize[1] - 1, xSize[0] - 1)] - exact_solution.value_p(0, NY, NZ, t)) *
                      (grid.p[getp(0, xSize[1] - 1, xSize[0] - 1)] - exact_solution.value_p(0, NY, NZ, t)) *
                      DX * DY * DZ / 8);
        }
    }
    // middle slices
    {
        for (int i = lbx; i < xSize[2] - rbx; i++)
        {
            if (lby)
            {
                error += ((grid.p[getp(i, 0, 0)] - exact_solution.value_p(i + offset_x, 0, 0, t)) *
                          (grid.p[getp(i, 0, 0)] - exact_solution.value_p(i + offset_x, 0, 0, t)) *
                          DX * DY * DZ / 4);

                for (int k = 1; k < xSize[0] - 1; k++)
                {
                    error += ((grid.p[getp(i, 0, k)] - exact_solution.value_p(i + offset_x, 0, k, t)) *
                              (grid.p[getp(i, 0, k)] - exact_solution.value_p(i + offset_x, 0, k, t)) *
                              DX * DY * DZ / 2);
                }
                error += ((grid.p[getp(i, 0, xSize[0] - 1)] - exact_solution.value_p(i + offset_x, 0, NZ, t)) *
                          (grid.p[getp(i, 0, xSize[0] - 1)] - exact_solution.value_p(i + offset_x, 0, NZ, t)) *
                          DX * DY * DZ / 4);
            }
            for (int j = lby; j < xSize[1] - rby; j++)
            {
                error += ((grid.p[getp(i, j, 0)] - exact_solution.value_p(i + offset_x, j + offset_y, 0, t)) *
                          (grid.p[getp(i, j, 0)] - exact_solution.value_p(i + offset_x, j + offset_y, 0, t)) *
                          DX * DY * DZ / 2);

                for (int k = 1; k < xSize[0] - 1; k++)
                {
                    error += ((grid.p[getp(i, j, k)] - exact_solution.value_p(i + offset_x, j + offset_y, k, t)) *
                              (grid.p[getp(i, j, k)] - exact_solution.value_p(i + offset_x, j + offset_y, k, t)) *
                              DX * DY * DZ);
                }

                error += ((grid.p[getp(i, j, xSize[0] - 1)] - exact_solution.value_p(i + offset_x, j + offset_y, NZ, t)) *
                          (grid.p[getp(i, j, xSize[0] - 1)] - exact_solution.value_p(i + offset_x, j + offset_y, NZ, t)) *
                          DX * DY * DZ / 2);
            }
            if (rby)
            {
                error += ((grid.p[getp(i, xSize[1] - 1, 0)] - exact_solution.value_p(i + offset_x, NY, 0, t)) *
                          (grid.p[getp(i, xSize[1] - 1, 0)] - exact_solution.value_p(i + offset_x, NY, 0, t)) *
                          DX * DY * DZ / 4);

                for (int k = 1; k < xSize[0] - 1; k++)
                {
                    error += ((grid.p[getp(i, xSize[1] - 1, k)] - exact_solution.value_p(i + offset_x, NY, k, t)) *
                              (grid.p[getp(i, xSize[1] - 1, k)] - exact_solution.value_p(i + offset_x, NY, k, t)) *
                              DX * DY * DZ / 2);
                }

                error += ((grid.p[getp(i, xSize[1] - 1, xSize[0] - 1)] - exact_solution.value_p(i + offset_x, NY, NZ, t)) *
                          (grid.p[getp(i, xSize[1] - 1, xSize[0] - 1)] - exact_solution.value_p(i + offset_x, NY, NZ, t)) *
                          DX * DY * DZ / 4);
            }
        }
    }
    // last slice (right face)
    if (rbx)
    {
        if (lby)
        {
            error += ((grid.p[getp(xSize[2] - 1, 0, 0)] - exact_solution.value_p(NX, 0, 0, t)) *
                      (grid.p[getp(xSize[2] - 1, 0, 0)] - exact_solution.value_p(NX, 0, 0, t)) *
                      DX * DY * DZ / 8);

            for (int k = 1; k < xSize[0] - 1; k++)
            {
                error += ((grid.p[getp(xSize[2] - 1, 0, k)] - exact_solution.value_p(NX, 0, k, t)) *
                          (grid.p[getp(xSize[2] - 1, 0, k)] - exact_solution.value_p(NX, 0, k, t)) *
                          DX * DY * DZ / 4);
            }

            error += ((grid.p[getp(xSize[2] - 1, 0, xSize[0] - 1)] - exact_solution.value_p(NX, 0, NZ, t)) *
                      (grid.p[getp(xSize[2] - 1, 0, xSize[0] - 1)] - exact_solution.value_p(NX, 0, NZ, t)) *
                      DX * DY * DZ / 8);
        }
        for (int j = lby; j < xSize[1] - rby; j++)
        {
            error += ((grid.p[getp(xSize[2] - 1, j, 0)] - exact_solution.value_p(NX, j + offset_y, 0, t)) *
                      (grid.p[getp(xSize[2] - 1, j, 0)] - exact_solution.value_p(NX, j + offset_y, 0, t)) *
                      DX * DY * DZ / 4);
            for (int k = 1; k < xSize[0] - 1; k++)
            {
                error += ((grid.p[getp(xSize[2] - 1, j, k)] - exact_solution.value_p(NX, j + offset_y, k, t)) *
                          (grid.p[getp(xSize[2] - 1, j, k)] - exact_solution.value_p(NX, j + offset_y, k, t)) *
                          DX * DY * DZ / 2);
            }
            error += ((grid.p[getp(xSize[2] - 1, j, xSize[0] - 1)] - exact_solution.value_p(NX, j + offset_y, NZ, t)) *
                      (grid.p[getp(xSize[2] - 1, j, xSize[0] - 1)] - exact_solution.value_p(NX, j + offset_y, NZ, t)) *
                      DX * DY * DZ / 4);
        }
        if (rby)
        {
            error += ((grid.p[getp(xSize[2] - 1, xSize[1] - 1, 0)] - exact_solution.value_p(NX, NY, 0, t)) *
                      (grid.p[getp(xSize[2] - 1, xSize[1] - 1, 0)] - exact_solution.value_p(NX, NY, 0, t)) *
                      DX * DY * DZ / 8);

            for (int k = 1; k < xSize[0] - 1; k++)
            {
                error += ((grid.p[getp(xSize[2] - 1, xSize[1] - 1, k)] - exact_solution.value_p(NX, NY, k, t)) *
                          (grid.p[getp(xSize[2] - 1, xSize[1] - 1, k)] - exact_solution.value_p(NX, NY, k, t)) *
                          DX * DY * DZ / 4);
            }

            error += ((grid.p[getp(xSize[2] - 1, xSize[1] - 1, xSize[0] - 1)] - exact_solution.value_p(NX, NY, NZ, t)) *
                      (grid.p[getp(xSize[2] - 1, xSize[1] - 1, xSize[0] - 1)] - exact_solution.value_p(NX, NY, NZ, t)) *
                      DX * DY * DZ / 8);
        }
    }

    return error;
}
