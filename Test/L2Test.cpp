# include <iostream>
# include "../includes/core.hpp"
# include "../includes/utils.hpp"
# include "../includes/grid.hpp"
# include "../includes/boundary.hpp"
# define nx 10
# define ny 10
# define nz 10
# define dx 0.1
# define dy 0.1
# define dz 0.1
# define Re 200.0
# define dt 0.1
# define T 1.0
# define lx 1
# define ly 1
# define lz 1
# define t 85.0

int main()
{
    //IcoNS icoNS(lx, ly, lz, nx, ny, nz, dt, T, Re, "input.txt", "output.txt");
    Grid grid(nx, ny, nz);
    ExactSolution exact_solution(dx, dy, dz);
    for (size_t i = 0; i < nx; i++)
    {
        for (size_t j = 0; j < ny+1; j++)
        {
            for (size_t k = 0; k < nz+1; k++)
            {
                grid.u[i*(ny+1)*(nz+1) + j*(nz+1) + k] = exact_solution.value_x(i+0.5, j, k,t);
            }
        }
    }

    for (size_t i = 0; i < nx+1; i++)
    {
        for (size_t j = 0; j < ny; j++)
        {
            for (size_t k = 0; k < nz+1; k++)
            {
                grid.v[i*ny*(nz+1) + j*(nz+1) + k] = exact_solution.value_y(i, j+0.5, k,t);
            }
        }
    }

    for (size_t i = 0; i < nx+1; i++)
    {
        for (size_t j = 0; j < ny+1; j++)
        {
            for (size_t k = 0; k < nz; k++)
            {
                grid.w[i*(ny+1)*nz + j*nz + k] = exact_solution.value_z(i, j, k+0.5,t);
            }
        }
    }

    

    double error = 0.0;
{
    

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
    std::cout << error << std::endl;
}

{
    

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
    std::cout << error << std::endl;
}

{
    
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

    std::cout << error << std::endl;
}
    
}

    

