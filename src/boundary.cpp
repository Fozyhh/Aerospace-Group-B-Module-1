#include "boundary.hpp"

Boundary::Boundary(Grid *grid_, double dx_, double dy_, double dz_)
    : grid(grid_),
      dx(dx_),
      dy(dy_),
      dz(dz_),
      prec(1 / 2)
{
}

// The method, which takes as input only the time step, is called by the program at the start of the time step.
// For each boundary node, it takes the exact value dor each component from the input and updates them.
// It also updates those values that are not directly on a face, but need an approximation.
void Boundary::update_boundary(double t)
{

    // Each face is numbered from 0 to 5 and we treat every face separately

    // LEFT FACE

    size_t face = 0;

    for (size_t k = 0; k < NZ; k++)
    {
        grid->w[k] = boundary_value_w[face]->value(0, 0, k, t);
        grid->v[k] = boundary_value_v[face]->value(0, 0, k, t);
    }
    grid->v[NZ] = boundary_value_v[face]->value(0, 0, NZ, t);

    for (size_t j = 1; j < NY; j++)
    {
        grid->v[j * (NZ + 1)] = boundary_value_v[face]->value(0, j, 0, t);
        grid->w[j * NZ] = boundary_value_w[face]->value(0, j, 0, t);
        for (size_t k = 1; k < NZ; k++)
        {
            grid->v[j * (NZ + 1) + k] = boundary_value_v[face]->value(0, j, k, t);
            grid->w[j * NZ + k] = boundary_value_w[face]->value(0, j, k, t);
            grid->u[j * (NZ + 1) + k] = approximate_boundary_u(0, j, k, t, face);
        }
        grid->v[j * (NZ + 1) + NZ] = boundary_value_v[face]->value(0, j, NZ, t);
    }

    for (size_t k = 0; k < NZ; k++)
    {
        grid->w[NY * NZ + k] = boundary_value_w[face]->value(0, NY, k, t);
    }

    // RIGHT FACE

    face = 1;

    for (size_t k = 0; k < NZ; k++)
    {
        grid->v[NX * NY * (NZ + 1) + k] = boundary_value_v[face]->value(NX, 0, k, t);
        grid->w[NX * (NY + 1) * NZ + k] = boundary_value_w[face]->value(NX, 0, k, t);
    }
    grid->v[NX * NY * (NZ + 1) + NZ] = boundary_value_v[face]->value(NX, 0, NZ, t);

    for (size_t j = 1; j < NY; j++)
    {
        grid->v[NX * NY * (NZ + 1) + j * (NZ + 1)] = boundary_value_v[face]->value(NX, j, 0, t);
        grid->w[NX * (NY + 1) * NZ + j * NZ] = boundary_value_w[face]->value(NX, j, 0, t);
        for (size_t k = 1; k < NZ; k++)
        {
            grid->v[NX * NY * (NZ + 1) + j * (NZ + 1) + k] = boundary_value_v[face]->value(NX, j, k, t);
            grid->w[NX * (NY + 1) * NZ + j * NZ + k] = boundary_value_w[face]->value(NX, j, k, t);
            grid->u[(NX - 1) * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k] = approximate_boundary_u(NX - 1, j, k, t, face);
        }
        grid->v[NX * NY * (NZ + 1) + j * (NZ + 1) + NZ] = boundary_value_v[face]->value(NX, j, NZ, t);
    }

    for (size_t k = 0; k < NZ; k++)
    {
        grid->w[NX * (NY + 1) * NZ + NY * NZ + k] = boundary_value_w[face]->value(NX, NY, k, t);
    }

    // FRONT FACE

    face = 2;

    for (size_t k = 1; k < NZ; k++)
    {
        grid->u[k] = boundary_value_u[face]->value(0, 0, k, t);
    }

    for (size_t i = 1; i < NX; i++)
    {
        grid->w[i * (NY + 1) * NZ] = boundary_value_w[face]->value(i, 0, 0, t);
        for (size_t k = 1; k < NZ; k++)
        {
            grid->u[i * (NY + 1) * (NZ + 1) + k] = boundary_value_u[face]->value(i, 0, k, t);
            grid->w[i * (NY + 1) * NZ + k] = boundary_value_w[face]->value(i, 0, k, t);
            grid->v[i * NY * (NZ + 1) + k] = approximate_boundary_v(i, 0, k, t, face);
        }
    }

    // BACK FACE

    face = 3;

    for (size_t k = 1; k < NZ; k++)
    {
        grid->u[NY * (NZ + 1) + k] = boundary_value_u[face]->value(0, NY, k, t);
    }

    for (size_t i = 1; i < NX; i++)
    {
        grid->w[(NY + 1) * NZ * i + NY * NZ] = boundary_value_w[face]->value(i, NY, 0, t);
        for (size_t k = 1; k < NZ; k++)
        {
            grid->u[NY * (NZ + 1) + (NY + 1) * (NZ + 1) * i + k] = boundary_value_u[face]->value(i, NY, k, t);
            grid->w[NY * NZ + i * (NY + 1) * NZ + k] = boundary_value_w[face]->value(i, NY, k, t);
            grid->v[(NY - 1) * (NZ + 1) + NY * (NZ + 1) * i + k] = approximate_boundary_v(i, NY - 1, k, t, face);
        }
    }

    // LOWER FACE

    face = 4;

    for (size_t j = 0; j < NY + 1; j++)
    {
        grid->u[j * (NZ + 1)] = boundary_value_u[face]->value(0, j, 0, t);
    }

    for (size_t i = 1; i < NX; i++)
    {
        grid->u[i * (NY + 1) * (NZ + 1)] = boundary_value_u[face]->value(i, 0, 0, t);
        grid->v[NY * (NZ + 1) * i] = boundary_value_v[face]->value(i, 0, 0, t);
        for (size_t j = 1; j < NY; j++)
        {
            grid->u[(NY + 1) * (NZ + 1) * i + j * (NZ + 1)] = boundary_value_u[face]->value(i, j, 0, t);
            grid->v[NY * (NZ + 1) * i + j * (NZ + 1)] = boundary_value_v[face]->value(i, j, 0, t);
            grid->w[(NY + 1) * NZ * i + j * NZ] = approximate_boundary_w(i, j, 0, t, face);
        }
        grid->u[(NY + 1) * (NZ + 1) * i + NY * (NZ + 1)] = boundary_value_u[face]->value(i, NY, 0, t);
    }

    // UPPER FACE

    face = 5;

    for (size_t j = 0; j < NY + 1; j++)
    {
        grid->u[j * (NZ + 1) + NZ] = boundary_value_u[face]->value(0, j, NZ, t);
    }

    for (size_t i = 1; i < NX; i++)
    {
        grid->u[i * (NY + 1) * (NZ + 1) + NZ] = boundary_value_u[face]->value(i, 0, NZ, t);
        grid->v[NY * (NZ + 1) * i + NZ] = boundary_value_v[face]->value(i, 0, NZ, t);
        for (size_t j = 1; j < NY; j++)
        {
            grid->u[(NY + 1) * (NZ + 1) * i + j * (NZ + 1) + NZ] = boundary_value_u[face]->value(i, j, NZ, t);
            grid->v[NY * (NZ + 1) * i + j * (NZ + 1) + NZ] = boundary_value_v[face]->value(i, j, NZ, t);
            grid->w[(NY + 1) * NZ * i + j * NZ + (NZ - 1)] = approximate_boundary_w(i, j, NZ - 1, t, face);
        }
        grid->u[(NY + 1) * (NZ + 1) * i + NY * (NZ + 1) + NZ] = boundary_value_u[face]->value(i, NY, NZ, t);
    }
}

// Performs the approximation of the component u that isn't precisely on the boundary.
double Boundary::approximate_boundary_u(size_t x, size_t y, size_t z, double t, size_t face)
{
    double dv = ((boundary_value_v[face]->value(x, y + prec, z, t)) -
                 ((boundary_value_v[face]->value(x, y - prec, z, t)))) /
                (2 * prec);

    double dw = ((boundary_value_w[face]->value(x, y, z + prec, t)) -
                 ((boundary_value_w[face]->value(x, y, z - prec, t)))) /
                (2 * prec);

    return boundary_value_u[face]->value(x, y, z, t) - (dv + dw) * (dx / 2);
}

// Performs the approximation of the component v that isn't precisely on the boundary.
double Boundary::approximate_boundary_v(size_t x, size_t y, size_t z, double t, size_t face)
{
    double du = ((boundary_value_u[face]->value(x + prec, y, z, t)) -
                 ((boundary_value_u[face]->value(x - prec, y, z, t)))) /
                (2 * prec);

    double dw = ((boundary_value_w[face]->value(x, y, z + prec, t)) -
                 ((boundary_value_w[face]->value(x, y, z - prec, t)))) /
                (2 * prec);

    return boundary_value_v[face]->value(x, y, z, t) - (du + dw) * (dy / 2);
}

// Performs the approximation of the component w that isn't precisely on the boundary.
double Boundary::approximate_boundary_w(size_t x, size_t y, size_t z, double t, size_t face)
{
    double du = ((boundary_value_u[face]->value(x + prec, y, z, t)) -
                 ((boundary_value_u[face]->value(x - prec, y, z, t)))) /
                (2 * prec);

    double dv = ((boundary_value_v[face]->value(x, y + prec, z, t)) -
                 ((boundary_value_v[face]->value(x, y - prec, z, t)))) /
                (2 * prec);
    return boundary_value_w[face]->value(x, y, z, t) - (du + dv) * (dz / 2);
}

void Boundary::addFunction(size_t direction, std::shared_ptr<BoundaryFunction> x)
{
    switch (direction)
    {
    case 0:
        boundary_value_u.push_back(x);
        break;
    case 1:
        boundary_value_v.push_back(x);
        break;
    case 2:
        boundary_value_w.push_back(x);
        break;
    }
};