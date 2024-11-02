#include "boundary.hpp"

// The method, which takes as input only the time step, is called by the program at the start of the time step.
// For each boundary node, it takes the exact value dor each component from the input and updates them.
// It also updates those values that are not directly on a face, but need an approximation.
void Boundary::update_boundary(std::vector<double> &Yx, std::vector<double> &Yy, std::vector<double> &Yz, double t)
{

    // Each face is numbered from 0 to 5 and we treat every face separately
    // LEFT FACE

    size_t face = 0;

    for (size_t k = 0; k < NZ; k++)
    {
        Yz[k] = boundary_value_w[face]->value(0, 0, k, t);
        Yy[k] = boundary_value_v[face]->value(0, 0, k, t);
    }
    Yy[NZ] = boundary_value_v[face]->value(0, 0, NZ, t);
    for (size_t j = 1; j < NY; j++)
    {
        Yy[j * (NZ + 1)] = boundary_value_v[face]->value(0, j, 0, t);
        Yz[j * NZ] = boundary_value_w[face]->value(0, j, 0, t);
        for (size_t k = 1; k < NZ; k++)
        {
            Yy[j * (NZ + 1) + k] = boundary_value_v[face]->value(0, j, k, t);
            Yz[j * NZ + k] = boundary_value_w[face]->value(0, j, k, t);
            Yx[j * (NZ + 1) + k] = approximate_boundary_u(0, j, k, t, face, 1);
        }
        Yy[j * (NZ + 1) + NZ] = boundary_value_v[face]->value(0, j, NZ, t);
    }

    for (size_t k = 0; k < NZ; k++)
    {
        Yz[NY * NZ + k] = boundary_value_w[face]->value(0, NY, k, t);
    }

    // RIGHT FACE

    face = 1;

    for (size_t k = 0; k < NZ; k++)
    {
        Yy[NX * NY * (NZ + 1) + k] = boundary_value_v[face]->value(NX, 0, k, t);
        Yz[NX * (NY + 1) * NZ + k] = boundary_value_w[face]->value(NX, 0, k, t);
    }
    Yy[NX * NY * (NZ + 1) + NZ] = boundary_value_v[face]->value(NX, 0, NZ, t);

    for (size_t j = 1; j < NY; j++)
    {
        Yy[NX * NY * (NZ + 1) + j * (NZ + 1)] = boundary_value_v[face]->value(NX, j, 0, t);
        Yz[NX * (NY + 1) * NZ + j * NZ] = boundary_value_w[face]->value(NX, j, 0, t);
        for (size_t k = 1; k < NZ; k++)
        {
            Yy[NX * NY * (NZ + 1) + j * (NZ + 1) + k] = boundary_value_v[face]->value(NX, j, k, t);
            Yz[NX * (NY + 1) * NZ + j * NZ + k] = boundary_value_w[face]->value(NX, j, k, t);
            Yx[(NX - 1) * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k] = approximate_boundary_u(NX, j, k, t, face, -1);
        }
        Yy[NX * NY * (NZ + 1) + j * (NZ + 1) + NZ] = boundary_value_v[face]->value(NX, j, NZ, t);
    }

    for (size_t k = 0; k < NZ; k++)
    {
        Yz[NX * (NY + 1) * NZ + NY * NZ + k] = boundary_value_w[face]->value(NX, NY, k, t);
    }

    // FRONT FACE

    face = 2;

    for (size_t k = 1; k < NZ; k++)
    {
        Yx[k] = boundary_value_u[face]->value(0, 0, k, t);
    }

    for (size_t i = 1; i < NX; i++)
    {
        Yz[i * (NY + 1) * NZ] = boundary_value_w[face]->value(i, 0, 0, t);
        for (size_t k = 1; k < NZ; k++)
        {
            Yx[i * (NY + 1) * (NZ + 1) + k] = boundary_value_u[face]->value(i, 0, k, t);
            Yz[i * (NY + 1) * NZ + k] = boundary_value_w[face]->value(i, 0, k, t);
            Yy[i * NY * (NZ + 1) + k] = approximate_boundary_v(i, 0, k, t, face, 1);
        }
    }

    // BACK FACE

    face = 3;

    for (size_t k = 1; k < NZ; k++)
    {
        Yx[NY * (NZ + 1) + k] = boundary_value_u[face]->value(0, NY, k, t);
    }

    for (size_t i = 1; i < NX; i++)
    {
        Yz[(NY + 1) * NZ * i + NY * NZ] = boundary_value_w[face]->value(i, NY, 0, t);
        for (size_t k = 1; k < NZ; k++)
        {
            Yx[NY * (NZ + 1) + (NY + 1) * (NZ + 1) * i + k] = boundary_value_u[face]->value(i, NY, k, t);
            Yz[NY * NZ + i * (NY + 1) * NZ + k] = boundary_value_w[face]->value(i, NY, k, t);
            Yy[(NY - 1) * (NZ + 1) + NY * (NZ + 1) * i + k] = approximate_boundary_v(i, NY, k, t, face, -1);
        }
    }

    // LOWER FACE

    face = 4;

    for (size_t j = 0; j < NY + 1; j++)
    {
        Yx[j * (NZ + 1)] = boundary_value_u[face]->value(0, j, 0, t);
    }

    for (size_t i = 1; i < NX; i++)
    {
        Yx[i * (NY + 1) * (NZ + 1)] = boundary_value_u[face]->value(i, 0, 0, t);
        Yy[NY * (NZ + 1) * i] = boundary_value_v[face]->value(i, 0, 0, t);
        for (size_t j = 1; j < NY; j++)
        {
            Yx[(NY + 1) * (NZ + 1) * i + j * (NZ + 1)] = boundary_value_u[face]->value(i, j, 0, t);
            Yy[NY * (NZ + 1) * i + j * (NZ + 1)] = boundary_value_v[face]->value(i, j, 0, t);
            Yz[(NY + 1) * NZ * i + j * NZ] = approximate_boundary_w(i, j, 0, t, face, 1);
        }
        Yx[(NY + 1) * (NZ + 1) * i + NY * (NZ + 1)] = boundary_value_u[face]->value(i, NY, 0, t);
    }

    // UPPER FACE

    face = 5;

    for (size_t j = 0; j < NY + 1; j++)
    {
        Yx[j * (NZ + 1) + NZ] = boundary_value_u[face]->value(0, j, NZ, t);
    }

    for (size_t i = 1; i < NX; i++)
    {
        Yx[i * (NY + 1) * (NZ + 1) + NZ] = boundary_value_u[face]->value(i, 0, NZ, t);
        Yy[NY * (NZ + 1) * i + NZ] = boundary_value_v[face]->value(i, 0, NZ, t);
        for (size_t j = 1; j < NY; j++)
        {
            Yx[(NY + 1) * (NZ + 1) * i + j * (NZ + 1) + NZ] = boundary_value_u[face]->value(i, j, NZ, t);
            Yy[NY * (NZ + 1) * i + j * (NZ + 1) + NZ] = boundary_value_v[face]->value(i, j, NZ, t);
            Yz[(NY + 1) * NZ * i + j * NZ + (NZ - 1)] = approximate_boundary_w(i, j, NZ, t, face, -1);
        }
        Yx[(NY + 1) * (NZ + 1) * i + NY * (NZ + 1) + NZ] = boundary_value_u[face]->value(i, NY, NZ, t);
    }
}

void Boundary::update_boundary(std::array<double, NX *(NY + 1) * (NZ + 1)> &Yx, std::array<double, (NX + 1) * NY *(NZ + 1)> &Yy, std::array<double, (NX + 1) * (NY + 1) * NZ> &Yz, double t)
{

    // Each face is numbered from 0 to 5 and we treat every face separately
    // LEFT FACE

    size_t face = 0;

    for (size_t k = 0; k < NZ; k++)
    {
        Yz[k] = boundary_value_w[face]->value(0, 0, k, t);
        Yy[k] = boundary_value_v[face]->value(0, 0, k, t);
    }
    Yy[NZ] = boundary_value_v[face]->value(0, 0, NZ, t);
    for (size_t j = 1; j < NY; j++)
    {
        Yy[j * (NZ + 1)] = boundary_value_v[face]->value(0, j, 0, t);
        Yz[j * NZ] = boundary_value_w[face]->value(0, j, 0, t);
        for (size_t k = 1; k < NZ; k++)
        {
            Yy[j * (NZ + 1) + k] = boundary_value_v[face]->value(0, j, k, t);
            Yz[j * NZ + k] = boundary_value_w[face]->value(0, j, k, t);
            Yx[j * (NZ + 1) + k] = approximate_boundary_u(0, j, k, t, face, 1);
        }
        Yy[j * (NZ + 1) + NZ] = boundary_value_v[face]->value(0, j, NZ, t);
    }

    for (size_t k = 0; k < NZ; k++)
    {
        Yz[NY * NZ + k] = boundary_value_w[face]->value(0, NY, k, t);
    }

    // RIGHT FACE

    face = 1;

    for (size_t k = 0; k < NZ; k++)
    {
        Yy[NX * NY * (NZ + 1) + k] = boundary_value_v[face]->value(NX, 0, k, t);
        Yz[NX * (NY + 1) * NZ + k] = boundary_value_w[face]->value(NX, 0, k, t);
    }
    Yy[NX * NY * (NZ + 1) + NZ] = boundary_value_v[face]->value(NX, 0, NZ, t);

    for (size_t j = 1; j < NY; j++)
    {
        Yy[NX * NY * (NZ + 1) + j * (NZ + 1)] = boundary_value_v[face]->value(NX, j, 0, t);
        Yz[NX * (NY + 1) * NZ + j * NZ] = boundary_value_w[face]->value(NX, j, 0, t);
        for (size_t k = 1; k < NZ; k++)
        {
            Yy[NX * NY * (NZ + 1) + j * (NZ + 1) + k] = boundary_value_v[face]->value(NX, j, k, t);
            Yz[NX * (NY + 1) * NZ + j * NZ + k] = boundary_value_w[face]->value(NX, j, k, t);
            Yx[(NX - 1) * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k] = approximate_boundary_u(NX, j, k, t, face, -1);
        }
        Yy[NX * NY * (NZ + 1) + j * (NZ + 1) + NZ] = boundary_value_v[face]->value(NX, j, NZ, t);
    }

    for (size_t k = 0; k < NZ; k++)
    {
        Yz[NX * (NY + 1) * NZ + NY * NZ + k] = boundary_value_w[face]->value(NX, NY, k, t);
    }

    // FRONT FACE

    face = 2;

    for (size_t k = 1; k < NZ; k++)
    {
        Yx[k] = boundary_value_u[face]->value(0, 0, k, t);
    }

    for (size_t i = 1; i < NX; i++)
    {
        Yz[i * (NY + 1) * NZ] = boundary_value_w[face]->value(i, 0, 0, t);
        for (size_t k = 1; k < NZ; k++)
        {
            Yx[i * (NY + 1) * (NZ + 1) + k] = boundary_value_u[face]->value(i, 0, k, t);
            Yz[i * (NY + 1) * NZ + k] = boundary_value_w[face]->value(i, 0, k, t);
            Yy[i * NY * (NZ + 1) + k] = approximate_boundary_v(i, 0, k, t, face, 1);
        }
    }

    // BACK FACE

    face = 3;

    for (size_t k = 1; k < NZ; k++)
    {
        Yx[NY * (NZ + 1) + k] = boundary_value_u[face]->value(0, NY, k, t);
    }

    for (size_t i = 1; i < NX; i++)
    {
        Yz[(NY + 1) * NZ * i + NY * NZ] = boundary_value_w[face]->value(i, NY, 0, t);
        for (size_t k = 1; k < NZ; k++)
        {
            Yx[NY * (NZ + 1) + (NY + 1) * (NZ + 1) * i + k] = boundary_value_u[face]->value(i, NY, k, t);
            Yz[NY * NZ + i * (NY + 1) * NZ + k] = boundary_value_w[face]->value(i, NY, k, t);
            Yy[(NY - 1) * (NZ + 1) + NY * (NZ + 1) * i + k] = approximate_boundary_v(i, NY, k, t, face, -1);
        }
    }

    // LOWER FACE

    face = 4;

    for (size_t j = 0; j < NY + 1; j++)
    {
        Yx[j * (NZ + 1)] = boundary_value_u[face]->value(0, j, 0, t);
    }

    for (size_t i = 1; i < NX; i++)
    {
        Yx[i * (NY + 1) * (NZ + 1)] = boundary_value_u[face]->value(i, 0, 0, t);
        Yy[NY * (NZ + 1) * i] = boundary_value_v[face]->value(i, 0, 0, t);
        for (size_t j = 1; j < NY; j++)
        {
            Yx[(NY + 1) * (NZ + 1) * i + j * (NZ + 1)] = boundary_value_u[face]->value(i, j, 0, t);
            Yy[NY * (NZ + 1) * i + j * (NZ + 1)] = boundary_value_v[face]->value(i, j, 0, t);
            Yz[(NY + 1) * NZ * i + j * NZ] = approximate_boundary_w(i, j, 0, t, face, 1);
        }
        Yx[(NY + 1) * (NZ + 1) * i + NY * (NZ + 1)] = boundary_value_u[face]->value(i, NY, 0, t);
    }

    // UPPER FACE

    face = 5;

    for (size_t j = 0; j < NY + 1; j++)
    {
        Yx[j * (NZ + 1) + NZ] = boundary_value_u[face]->value(0, j, NZ, t);
    }

    for (size_t i = 1; i < NX; i++)
    {
        Yx[i * (NY + 1) * (NZ + 1) + NZ] = boundary_value_u[face]->value(i, 0, NZ, t);
        Yy[NY * (NZ + 1) * i + NZ] = boundary_value_v[face]->value(i, 0, NZ, t);
        for (size_t j = 1; j < NY; j++)
        {
            Yx[(NY + 1) * (NZ + 1) * i + j * (NZ + 1) + NZ] = boundary_value_u[face]->value(i, j, NZ, t);
            Yy[NY * (NZ + 1) * i + j * (NZ + 1) + NZ] = boundary_value_v[face]->value(i, j, NZ, t);
            Yz[(NY + 1) * NZ * i + j * NZ + (NZ - 1)] = approximate_boundary_w(i, j, NZ, t, face, -1);
        }
        Yx[(NY + 1) * (NZ + 1) * i + NY * (NZ + 1) + NZ] = boundary_value_u[face]->value(i, NY, NZ, t);
    }
}

// Performs the approximation of the component u that isn't precisely on the boundary.
double Boundary::approximate_boundary_u(size_t x, size_t y, size_t z, double t, size_t face, int side)
{

    double dv = (boundary_value_v[face]->value(x, y, z, t) - boundary_value_v[face]->value(x, y - 1.0, z, t)) / DY;
    double dw = (boundary_value_w[face]->value(x, y, z, t) - boundary_value_w[face]->value(x, y, z - 1.0, t)) / DZ;

    return boundary_value_u[face]->value((x - 0.5 /*(DX/2.0)*/), y, z, t) - (dv + dw) * (DX / 2) * side;
}

// Performs the approximation of the component v that isn't precisely on the boundary.
double Boundary::approximate_boundary_v(size_t x, size_t y, size_t z, double t, size_t face, int side)
{
    double du = ((boundary_value_u[face]->value(x, y, z, t)) -
                 ((boundary_value_u[face]->value(x - 1, y, z, t)))) /
                (DX);

    double dw = ((boundary_value_w[face]->value(x, y, z, t)) -
                 ((boundary_value_w[face]->value(x, y, z - 1, t)))) /
                (DZ);

    return boundary_value_v[face]->value(x, y - 0.5, z, t) - (du + dw) * (DY / 2.0) * side;
}

// Performs the approximation of the component w that isn't precisely on the boundary.
double Boundary::approximate_boundary_w(size_t x, size_t y, size_t z, double t, size_t face, int side)
{
    double du = ((boundary_value_u[face]->value(x, y, z, t)) -
                 ((boundary_value_u[face]->value(x - 1.0, y, z, t)))) /
                (DX);

    double dv = ((boundary_value_v[face]->value(x, y, z, t)) -
                 ((boundary_value_v[face]->value(x, y - 1.0, z, t)))) /
                (DY);
    return boundary_value_w[face]->value(x, y, z - 0.5, t) - (du + dv) * (DZ / 2) * side;
}

void Boundary::addFunction(Direction direction, std::shared_ptr<BoundaryFunction> x)
{
    switch (direction)
    {
    case U:
        boundary_value_u.push_back(x);
        break;
    case V:
        boundary_value_v.push_back(x);
        break;
    case W:
        boundary_value_w.push_back(x);
        break;
    }
};
