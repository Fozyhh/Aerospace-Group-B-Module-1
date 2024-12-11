#include "boundary.hpp"

// TODO:
void Boundary::update_boundary(std::vector<Real> &Yx, std::vector<Real> &Yy, std::vector<Real> &Yz, Real t)
{
    // UPDATE EVERY INDEX SUCH THAT YOU SKIP FIRST FACE(GHOST)
    // AND THE FIRST AND LAST ROW OF EACH FACE
    int face, i, j, k;
    // X
    // LEFT FACE

    // face 0 -> lbx left
    // face 1 -> rbx right
    // face 2 -> lby front
    // face 3 -> rby back
    // face 4 5 -> no partition on z

    // Have to adjust indexes in exact values if it is made 
    // by processor with diff coords with global indexes 
    // Have to skip cycles: EITHER by adjusting loops or by hand for lines outside loops
    int offset_x_x = coords[0] * other_dim_x_x;
    int offset_y_x = (PY - 1 - coords[1]) * other_dim_y_x;
    int offset_x_y = coords[0] * other_dim_x_y;
    int offset_y_y = (PY - 1 - coords[1]) * other_dim_y_y;
    int offset_x_z = coords[0] * other_dim_x_z;
    int offset_y_z = (PY - 1 - coords[1]) * other_dim_y_z;

    if (lbx)
    {
        if (lby)
        { // if it's a processor with coords[1]==0 then it has the front face as boundary
            face = 2;
            for (int k = 0; k < NZ + 1; k++)
            {
                Yx[((newDimY_x + 1) * dim_z) + k] = boundary_value_u[face]->value(0, 0, k, t); // FOR X AND Y WE NEED GLOBAL INDICES OR IT DOESNT WORK
            }
        }

        for (int j = 1 + lby; j < newDimY_x - 1 - rby; j++) // Parti da 1(salti ghost) + 1(se devi considerare il boundary) e arrivi a newDimY -1(ghost)-1(se boundary destro)
        {
            face = 4;
            // lby coords = 1 rby coords = 0
            // we move in the opposite direction of the y process count unfortunately
            Yx[((newDimY_x)*dim_z) + j * (NZ + 1)] = boundary_value_u[face]->value(0, j + offset_y_x, 0, t);
            face = 0;
            for (int k = 1; k < NZ; k++)
            {
                Yx[((newDimY_x)*dim_z) + j * (NZ + 1) + k] = approximate_boundary_u(0, j + offset_y_x, k, t, face, 1);
            }
            face = 5;
            Yx[((newDimY_x)*dim_z) + j * (NZ + 1) + NZ] = boundary_value_u[face]->value(0, j + offset_y_x, NZ, t);
        }

        if (rby)
        {
            face = 3;
            for (int k = 0; k < NZ + 1; k++)
            {
                Yx[((newDimY_x)*dim_z) + (newDimY_x - 2) * (NZ + 1) + k] = boundary_value_u[face]->value(0, NY, k, t);
            }
        }
    }
    // CONTINUE LIKE THIS

    // MIDDLE POINTS
    for (int i = 1 + lbx; i < newDimX_x - 1 - rbx; i++)
    {
        if (lby)
        {
            face = 2;
            for (int k = 0; k < NZ + 1; k++)
            {
                Yx[i * newDimY_x * (NZ + 1) + (NZ + 1) + k] = boundary_value_u[face]->value(i + offset_x_x, 0, k, t);
            }
        }

        for (int j = 1 + lby; j < newDimY_x - 1 - rby; j++)
        {
            face = 4;
            Yx[newDimY_x * (NZ + 1) * i + j * (NZ + 1)] = boundary_value_u[face]->value(i + offset_x_x, j + offset_y_x, 0, t);

            face = 5;
            Yx[newDimY_x * (NZ + 1) * i + j * (NZ + 1) + NZ] = boundary_value_u[face]->value(i + offset_x_x, j + offset_y_x, NZ, t);
        }

        if (rby)
        {
            face = 3;
            for (int k = 0; k < NZ + 1; k++)
            {
                Yx[newDimY_x * (NZ + 1) * i + (newDimY_x - 2) * (NZ + 1) + k] = boundary_value_u[face]->value(i + offset_x_x, NY, k, t);
            }
        }
    }

    // RIGHT FACE
    if (rbx)
    {
        if (lby)
        {
            face = 2;
            for (int k = 0; k < NZ + 1; k++)
            {
                Yx[(newDimX_x - 2) * newDimY_x * (NZ + 1) + (NZ + 1) + k] = boundary_value_u[face]->value(NX - 1, 0, k, t);
            }
        }

        for (int j = 1 + lby; j < newDimY_x - 1 - rby; j++)
        {
            face = 4;
            Yx[(newDimX_x - 2) * newDimY_x * (NZ + 1) + j * (NZ + 1)] = boundary_value_u[face]->value(NX - 1, j + offset_y_x, 0, t);

            face = 1;
            for (int k = 1; k < NZ; k++)
            {
                Yx[(newDimX_x - 2) * (newDimY_x) * (NZ + 1) + j * (NZ + 1) + k] = approximate_boundary_u(NX, j + offset_y_x, k, t, face, -1);
            }

            face = 5;
            Yx[(newDimX_x - 2) * (newDimY_x) * (NZ + 1) + j * (NZ + 1) + NZ] = boundary_value_u[face]->value(NX - 1, j + offset_y_x, NZ, t);
        }
        if (rby)
        {
            face = 3;
            for (int k = 0; k < NZ + 1; k++)
            {
                Yx[(newDimX_x - 2) * newDimY_x * (NZ + 1) + (newDimY_x - 2) * (NZ + 1) + k] = boundary_value_u[face]->value(NX - 1, NY, k, t);
            }
        }
    }

    // Y
    // LEFT FACE
    if (lbx)
    {
        face = 0;
        for (int j = 1; j < newDimY_y - 1; j++)
        {
            for (int k = 0; k < NZ + 1; k++)
            {
                Yy[(newDimY_y * dim_z) + j * dim_z + k] = boundary_value_v[face]->value(0, j + offset_y_y, k, t);
            }
        }
    }

    // MIDDLE POINTS
    for (int i = 1 + lbx; i < newDimX_y - 1 - rbx; i++)
    {
        if (lby)
        {
            face = 4;
            Yy[i * newDimY_y * (NZ + 1) + dim_z] = boundary_value_v[face]->value(i + offset_x_y, offset_y_y, 0, t);
            face = 2;
            for (int k = 1; k < NZ; k++)
            {
                Yy[i * newDimY_y * (NZ + 1) + k + dim_z] = approximate_boundary_v(i + offset_x_y, 0, k, t, face, 1);
            }
            face = 5;
            Yy[newDimY_y * (NZ + 1) * i + dim_z + NZ] = boundary_value_v[face]->value(i + offset_x_y, 0 + offset_y_y, NZ, t);
        }
        for (int j = 1 + lby; j < newDimY_y - 1 - rby; j++)
        {
            face = 4;
            Yy[newDimY_y * (NZ + 1) * i + j * (NZ + 1)] = boundary_value_v[face]->value(i + offset_x_y, j + offset_y_y, 0, t);

            face = 5;
            Yy[newDimY_y * (NZ + 1) * i + j * (NZ + 1) + NZ] = boundary_value_v[face]->value(i + offset_x_y, j + offset_y_y, NZ, t);
        }

        if (rby)
        {
            face = 4;
            Yy[newDimY_y * (NZ + 1) * i + (newDimY_y - 2) * (NZ + 1)] = boundary_value_v[face]->value(i + offset_x_y, NY - 1, 0, t);
            face = 3;
            for (int k = 1; k < NZ; k++)
            {
                Yy[newDimY_y * (NZ + 1) * i + (newDimY_y - 2) * (NZ + 1) + k] = approximate_boundary_v(i + offset_x_y, NY, k, t, face, -1);
            }
            face = 5;
            Yy[newDimY_y * (NZ + 1) * i + (newDimY_y - 2) * (NZ + 1) + NZ] = boundary_value_v[face]->value(i + offset_x_y, NY - 1, NZ, t);
        }
    }

    // RIGHT FACE
    if (rbx)
    {
        face = 1;
        for (int j = 1; j < newDimY_y - 1; j++)
        {
            for (int k = 0; k < NZ + 1; k++)
            {
                Yy[(newDimX_y - 2) * newDimY_y * (NZ + 1) + j * (NZ + 1) + k] = boundary_value_v[face]->value(NX, j + offset_y_y, k, t);
            }
        }
    }

    // Z
    // LEFT FACE
    if (lbx)
    {
        face = 0;
        for (int j = 1; j < newDimY_z; j++)
        {
            for (int k = 0; k < NZ; k++)
            {
                Yz[newDimY_z * dim_z_z + j * NZ + k] = boundary_value_w[face]->value(0, j + offset_y_z, k, t);
            }
        }
    }
    // MIDDLE POINTS
    for (int i = 1 + lbx; i < newDimX_z - 1 - rbx; i++)
    {
        if (lby)
        {
            face = 2;
            for (int k = 0; k < NZ; k++)
            {
                Yz[i * (newDimY_z)*NZ + dim_z_z + k] = boundary_value_w[face]->value(i + offset_x_z, 0, k, t);
            }
        }
        for (int j = 1 + lby; j < newDimY_z - 1 - rby; j++)
        {
            face = 4;
            Yz[newDimY_z * NZ * i + j * NZ] = approximate_boundary_w(i + offset_x_z, j + offset_y_z, 0, t, face, 1);

            face = 5;
            Yz[newDimY_z * NZ * i + j * NZ + (NZ - 1)] = approximate_boundary_w(i + offset_x_z, j + offset_y_z, NZ, t, face, -1);
        }
        if (rby)
        {
            face = 3;
            for (int k = 0; k < NZ + 1; k++)
            {
                Yz[i * newDimY_z * NZ + (newDimY_z - 2) * NZ + k] = boundary_value_w[face]->value(i + offset_x_z, NY, k, t);
            }
        }
    }

    // RIGHT FACE
    if (rbx)
    {
        face = 1;
        for (int j = 1; j < newDimY_z - 1; j++)
        {
            for (int k = 0; k < NZ; k++)
            {
                Yz[(newDimX_z - 2) * newDimY_z * NZ + j * NZ + k] = boundary_value_w[face]->value(NX, j + offset_y_z, k, t);
            }
        }
    }
}


Real Boundary::approximate_boundary_u(int x, int y, int z, Real t, int face, int side)
{

    Real dv = (boundary_value_v[face]->value(x, y, z, t) - boundary_value_v[face]->value(x, y - 1.0, z, t)) / DY;
    Real dw = (boundary_value_w[face]->value(x, y, z, t) - boundary_value_w[face]->value(x, y, z - 1.0, t)) / DZ;

    return boundary_value_u[face]->value((x - 0.5 /*(DX/2.0)*/), y, z, t) - (dv + dw) * (DX / 2) * side;
}


Real Boundary::approximate_boundary_v(int x, int y, int z, Real t, int face, int side)
{
    Real du = ((boundary_value_u[face]->value(x, y, z, t)) -
               ((boundary_value_u[face]->value(x - 1, y, z, t)))) /
              (DX);

    Real dw = ((boundary_value_w[face]->value(x, y, z, t)) -
               ((boundary_value_w[face]->value(x, y, z - 1, t)))) /
              (DZ);

    return boundary_value_v[face]->value(x, y - 0.5, z, t) - (du + dw) * (DY / 2.0) * side;
}


Real Boundary::approximate_boundary_w(int x, int y, int z, Real t, int face, int side)
{
    Real du = ((boundary_value_u[face]->value(x, y, z, t)) -
               ((boundary_value_u[face]->value(x - 1.0, y, z, t)))) /
              (DX);

    Real dv = ((boundary_value_v[face]->value(x, y, z, t)) -
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

void Boundary::setBoundaryOffsets(int lbx_, int rbx_, int lby_, int rby_)
{
    lbx = lbx_;
    rbx = rbx_;
    lby = lby_;
    rby = rby_;
}

void Boundary::setCoords(int coords_[2])
{
    coords[0] = coords_[0];
    coords[1] = coords_[1];
}

void Boundary::setOtherDim(int other_dim_x_x_, int other_dim_y_x_, int other_dim_x_y_, int other_dim_y_y_, int other_dim_x_z_, int other_dim_y_z_)
{
    other_dim_x_x = other_dim_x_x_;
    other_dim_y_x = other_dim_y_x_;
    other_dim_x_y = other_dim_x_y_;
    other_dim_y_y = other_dim_y_y_;
    other_dim_x_z = other_dim_x_z_;
    other_dim_y_z = other_dim_y_z_;
}
