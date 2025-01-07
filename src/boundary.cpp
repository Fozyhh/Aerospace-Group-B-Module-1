#include "boundary.hpp"
#include <iostream>

/**
 * @brief The method is called by the program multiple during the time step, in order to update the values of the boundaries at each
 * requested time t, calculating the approximated ones too.
 *
 * @param Yx Boundary x velocities for the Y intermediate function related to the x direction.
 * @param Yy Boundary y velocities for the Y intermediate function related to the y direction.
 * @param Yz Boundary z velocities for the Y intermediate function related to the z direction.
 * @param t Time of the time discretization we are considering.
 */
void Boundary::update_boundary(std::vector<Real> &Yx, std::vector<Real> &Yy, std::vector<Real> &Yz, Real t)
{
    int face;
    // face 0 -> lbx left
    // face 1 -> rbx right
    // face 2 -> lby front
    // face 3 -> rby back
    // face 4 5 -> no partition on z

    // Global offsets
    int offset_x_x = coords[0] * other_dim_x_x -1;
    int offset_y_x = coords[1] * other_dim_y_x -1;
    int offset_x_y = coords[0] * other_dim_x_y -1;
    int offset_y_y = coords[1] * other_dim_y_y -1;
    int offset_x_z = coords[0] * other_dim_x_z -1;
    int offset_y_z = coords[1] * other_dim_y_z -1;

    // X
    // LEFT FACE
    if (lbx)
    {
        if (lby)
        {
            face = 2;
            for (int k = 0; k < dim_z; k++)
            {
                Yx[((newDimY_x + 1) * dim_z) + k] = boundary_value_u[face]->value(0, 0, k, t);
            }
        }

        for (int j = 1 + lby; j < newDimY_x - 1 - rby; j++) 
        {
            face = 4;
            Yx[((newDimY_x)*dim_z) + j * dim_z] = boundary_value_u[face]->value(0, j + offset_y_x, 0, t);
            face = 0;
            for (int k = 1; k < dim_z - 1; k++)
            {
                Yx[((newDimY_x)*dim_z) + j * dim_z + k] = approximate_boundary_u(0, j + offset_y_x, k, t, face, 1);
            }
            face = 5;
            Yx[((newDimY_x)*dim_z) + j * dim_z + dim_z-1] = boundary_value_u[face]->value(0, j + offset_y_x, NZ, t);
        }

        if (rby)
        {
            face = 3;
            for (int k = 0; k < dim_z; k++)
            {
                Yx[((newDimY_x)*dim_z) + (newDimY_x - 2) * dim_z + k] = boundary_value_u[face]->value(0, NY, k, t);
            }
        }
    }

    // MIDDLE POINTS
    for (int i = 1 + lbx; i < newDimX_x - 1 - rbx; i++)
    {
        if (lby)
        {
            face = 2;
            for (int k = 0; k < dim_z; k++)
            {
                Yx[i * newDimY_x * (dim_z) + (dim_z) + k] = boundary_value_u[face]->value(i + offset_x_x, 0, k, t);
            }
        }

        for (int j = 1 + lby; j < newDimY_x - 1 - rby; j++)
        {
            face = 4;
            Yx[newDimY_x * dim_z * i + j * dim_z] = boundary_value_u[face]->value(i + offset_x_x, j + offset_y_x, 0, t);

            face = 5;
            Yx[newDimY_x * dim_z * i + j * dim_z + dim_z-1] = boundary_value_u[face]->value(i + offset_x_x, j + offset_y_x, NZ, t);
        }

        if (rby)
        {
            face = 3;
            for (int k = 0; k < dim_z; k++)
            {
                Yx[newDimY_x * (dim_z) * i + (newDimY_x - 2) * (dim_z) + k] = boundary_value_u[face]->value(i + offset_x_x, NY, k, t);
            }
        }
    }

    // RIGHT FACE
    if (rbx)
    {
        if (lby)
        {
            face = 2;
            for (int k = 0; k < dim_z; k++)
            {
                Yx[(newDimX_x - 2) * newDimY_x * (dim_z) + (dim_z) + k] = boundary_value_u[face]->value(NX - 1, 0, k, t);
            }
        }

        for (int j = 1 + lby; j < newDimY_x - 1 - rby; j++)
        {
            face = 4;
            Yx[(newDimX_x - 2) * newDimY_x * (dim_z) + j * (dim_z)] = boundary_value_u[face]->value(NX - 1, j + offset_y_x, 0, t);

            face = 1;
            for (int k = 1; k < dim_z-1; k++)
            {
                Yx[(newDimX_x - 2) * (newDimY_x) * dim_z + j * dim_z + k] = approximate_boundary_u(NX, j + offset_y_x, k, t, face, -1);
            }

            face = 5;
            Yx[(newDimX_x - 2) * (newDimY_x) * dim_z + j * dim_z + dim_z-1] = boundary_value_u[face]->value(NX - 1, j + offset_y_x, NZ, t);
        }
        if (rby)
        {
            face = 3;
            for (int k = 0; k < dim_z; k++)
            {
                Yx[(newDimX_x - 2) * newDimY_x * (dim_z) + (newDimY_x - 2) * (dim_z) + k] = boundary_value_u[face]->value(NX - 1, NY, k, t);
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
            for (int k = 0; k < dim_z; k++)
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
            Yy[i * newDimY_y * dim_z + dim_z] = boundary_value_v[face]->value(i + offset_x_y, 0, 0, t);
            face = 2;
            for (int k = 1; k < dim_z-1; k++)
            {
                Yy[i * newDimY_y * dim_z + dim_z + k] = approximate_boundary_v(i + offset_x_y, 0, k, t, face, 1);
            }
            face = 5;
            Yy[i * newDimY_y * dim_z + dim_z + dim_z-1] = boundary_value_v[face]->value(i + offset_x_y, 0, NZ, t);
        }
        for (int j = 1 + lby; j < newDimY_y - 1 - rby; j++)
        {
            face = 4;
            Yy[newDimY_y * dim_z * i + j * dim_z] = boundary_value_v[face]->value(i + offset_x_y, j + offset_y_y, 0, t);

            face = 5;
            Yy[newDimY_y * dim_z * i + j * dim_z + dim_z-1] = boundary_value_v[face]->value(i + offset_x_y, j + offset_y_y, NZ, t);
        }

        if (rby)
        {
            face = 4;
            Yy[newDimY_y * dim_z * i + (newDimY_y - 2) * dim_z] = boundary_value_v[face]->value(i + offset_x_y, NY - 1, 0, t);
            face = 3;
            for (int k = 1; k < dim_z -1; k++)
            {
                Yy[newDimY_y * dim_z * i + (newDimY_y - 2) * dim_z + k] = approximate_boundary_v(i + offset_x_y, NY, k, t, face, -1);
            }
            face = 5;
            Yy[newDimY_y * dim_z * i + (newDimY_y - 2) * dim_z + dim_z -1] = boundary_value_v[face]->value(i + offset_x_y, NY - 1, NZ, t);
        }
    }

    // RIGHT FACE
    if (rbx)
    {
        face = 1;
        for (int j = 1; j < newDimY_y - 1; j++)
        {
            for (int k = 0; k < dim_z; k++)
            {
                Yy[(newDimX_y - 2) * newDimY_y * (dim_z) + j * (dim_z) + k] = boundary_value_v[face]->value(NX, j + offset_y_y, k, t);
            }
        }
    }

    // Z
    // LEFT FACE
    if (lbx)
    {
        face = 0;
        for (int j = 1; j < newDimY_z-1 ; j++)
        {
            for (int k = 0; k < dim_z_z; k++)
            {
                Yz[newDimY_z * dim_z_z + j * dim_z_z + k] = boundary_value_w[face]->value(0, j + offset_y_z, k, t);
            }
        }
    }
    // MIDDLE POINTS
    for (int i = 1 + lbx; i < newDimX_z - 1 - rbx; i++)
    {
        if (lby)
        {
            face = 2;
            for (int k = 0; k < dim_z_z; k++)
            {
                Yz[i * (newDimY_z)*dim_z_z + dim_z_z + k] = boundary_value_w[face]->value(i + offset_x_z, 0, k, t);
            }
        }
        for (int j = 1 + lby; j < newDimY_z - 1 - rby; j++)
        {
            face = 4;
            Yz[newDimY_z * dim_z_z * i + j * dim_z_z] = approximate_boundary_w(i + offset_x_z, j + offset_y_z, 0, t, face, 1);

            face = 5;
            Yz[newDimY_z * dim_z_z * i + j * dim_z_z + (dim_z_z - 1)] = approximate_boundary_w(i + offset_x_z, j + offset_y_z, dim_z_z, t, face, -1);
        }
        if (rby)
        {
            face = 3;
            for (int k = 0; k < dim_z_z + 1; k++)
            {
                Yz[i * newDimY_z * dim_z_z + (newDimY_z - 2) * dim_z_z + k] = boundary_value_w[face]->value(i + offset_x_z, NY, k, t);
            }
        }
    }

    // RIGHT FACE
    if (rbx)
    {
        face = 1;
        for (int j = 1; j < newDimY_z - 1; j++)
        {
            for (int k = 0; k < dim_z_z; k++)
            {
                Yz[(newDimX_z - 2) * newDimY_z * dim_z_z + j * dim_z_z + k] = boundary_value_w[face]->value(NX, j + offset_y_z, k, t);
            }
        }
    }
}


Real Boundary::approximate_boundary_u(int x, int y, int z, Real t, int face, int side)
{

    Real dv = (boundary_value_v[face]->value(x, y, z, t) - boundary_value_v[face]->value(x, y - 1.0, z, t)) / DY;
    Real dw = (boundary_value_w[face]->value(x, y, z, t) - boundary_value_w[face]->value(x, y, z - 1.0, t)) / DZ;

    return boundary_value_u[face]->value((x - 0.5), y, z, t) - (dv + dw) * (DX / 2) * side;
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


void Boundary::divergence(std::vector<Real> &Yx, std::vector<Real> &Yy, std::vector<Real> &Yz, double* &Y2_p, Real t, Real c)
{
    //TODO: check the global offsets are the same
    int offset_x_x = coords[0] * other_dim_x_x;
    int offset_y_x = coords[1] * other_dim_y_x;
    int offset_x_y = coords[0] * other_dim_x_y;
    int offset_y_y = coords[1] * other_dim_y_y;
    int offset_x_z = coords[0] * other_dim_x_z;
    int offset_y_z = coords[1] * other_dim_y_z;

    // LEFT FACE
    if(lbx){
        for (int j = lby; j < zSize[1] - rby; j++)
        {
            for(int k = 1; k < zSize[2] - 1; k++)
            {
                Y2_p[getp(0,j,k)] = 120.0 / (c * DT) *
                    ((-8*boundary_value_u[0]->value(-0.5,j + offset_y_x,k,t) + 9*Yx[getx(0,j,k)] - Yx[getx(1,j,k)]) / (3*DX) +
                    (Yy[gety(0,j,k)] - Yy[gety(0,j-1,k)]) / (DY) +
                    (Yz[getz(0,j,k)] - Yz[getz(0,j,k-1)]) / (DZ));
                //Y2_p[getp(0,j,k)] = cos(-0.5*DX)*cos((j+ offset_y_x)*DY)*cos(k*DZ)*(sin(t)-sin(t+64.0/120.0*DT));
            }
        }
    }

    if(rbx){
        // RIGHT FACE
        for (int j = lby; j < zSize[1] - rby; j++)
        {
            for(int k = 1; k < zSize[2] - 1; k++)
            {
                Y2_p[getp(zSize[0]-1,j,k)] = 120.0 / (c * DT) *
                    ((8*boundary_value_u[1]->value(NX-0.5,j + offset_y_x,k,t) - 9*Yx[getx(dim_x_x - 1, j, k)] + Yx[getx(dim_x_x - 2,j,k)]) / (3*DX) +
                    (Yy[gety(dim_x_y-1,j,k)] - Yy[gety(dim_x_y-1,j-1,k)]) / (DY) +
                    (Yz[getz(dim_x_z-1,j,k)] - Yz[getz(dim_x_z-1,j,k-1)]) / (DZ));
                //Y2_p[getp(zSize[0]-1,j,k)] = cos((NX-0.5)*DX)*cos((j+ offset_y_x)*DY)*cos(k*DZ)*(sin(t)-sin(t+64.0/120.0*DT));
            }
        }
    }

    // FRONT FACE
    if(lby){
        for (int i = lbx; i < zSize[0] - rbx; i++)
        {
            for(int k = 1; k < zSize[2] - 1; k++)
            {
                Y2_p[getp(i,0,k)] = 120.0 / (c * DT) *
                    ((Yx[getx(i,0,k)] - Yx[getx(i-1,0,k)]) / (DX) +
                    (-8*boundary_value_v[2]->value(i+ offset_x_y,-0.5,k,t) + 9*Yy[gety(i,0,k)] - Yy[gety(i,1,k)]) / (3*DY) +
                    (Yz[getz(i,0,k)] - Yz[getz(i,0,k-1)]) / (DZ));
                //Y2_p[getp(i,0,k)] = cos((i+ offset_x_y)*DX)*cos(k*DZ)*(sin(t)-sin(t+64.0/120.0*DT));
            }
        }
    }

    // BACK FACE
    if(rby){
        for (int i = lbx; i < zSize[0] - rbx; i++)
        {
            for(int k = 1; k < zSize[2] - 1; k++)
            {
                Y2_p[getp(i,zSize[1] - 1,k)] = 120.0 / (c * DT) *
                    ((Yx[getx(i,dim_y_x-1,k)] - Yx[getx(i-1,dim_y_x-1,k)]) / (DX) +
                    (8*boundary_value_v[3]->value(i + offset_x_y,NY-0.5,k,t) - 9*Yy[gety(i,dim_y_y-1,k)] + Yy[gety(i,dim_y_y-2,k)]) / (3*DY) +
                    (Yz[getz(i,dim_y_z-1,k)] - Yz[getz(i,dim_y_z-1,k-1)]) / (DZ));
                //Y2_p[getp(i,zSize[1] - 1,k)] = cos((i+ offset_x_y)*DX)*cos(k*DZ)*(sin(t)-sin(t+64.0/120.0*DT));
            }
        }
    }

    // LOWER FACE
    for (int i = lbx; i < zSize[0] - rbx; i++)
    {
        for(int j = lby; j < zSize[1] - rby; j++)
        {
            Y2_p[getp(i,j,0)] = 120.0 / (c * DT) *
                ((Yx[getx(i,j,0)] - Yx[getx(i-1,j,0)]) / (DX) +
                (Yy[gety(i,j,0)] - Yy[gety(i,j-1,0)]) / (DY) +
                (-8*boundary_value_w[4]->value(i + offset_x_z,j + offset_y_z,-0.5,t) + 9*Yz[getz(i,j,0)] - Yz[getz(i,j,1)]) / (3*DZ));
            //Y2_p[getp(i,j,0)] = cos((i+ offset_x_z)*DX)*cos((j+ offset_y_z)*DY)*(sin(t)-sin(t+64.0/120.0*DT));
        }
    }
    // UPPER FACE
    for (int i = lbx; i < zSize[0] - rbx; i++)
    {
        for(int j = lby; j < zSize[1] - rby; j++)
        {
            Y2_p[getp(i,j,zSize[2] - 1)] = 120.0 / (c * DT) *
                ((Yx[getx(i,j,dim_z - 1)] - Yx[getx(i-1,j,dim_z - 1)]) / (DX) +
                (Yy[gety(i,j,dim_z-1)] - Yy[gety(i,j-1,dim_z-1)]) / (DY) +
                (8*boundary_value_w[5]->value(i + offset_x_z,j + offset_y_z,NZ-0.5,t) - 9*Yz[getz(i,j,dim_z_z-1)] + Yz[getz(i,j,dim_z_z-2)]) / (3*DZ));
            //Y2_p[getp(i,j,zSize[2] - 1)] = cos((i+ offset_x_z)*DX)*cos((j+ offset_y_z)*DY)*(sin(t)-sin(t+64.0/120.0*DT));
        }
    }
    // 4 X EDGES
    for (int i = lbx; i < zSize[0] - rbx; i++)
    {
        // LOWER FRONT EDGE
        if(lby){
            Y2_p[getp(i,0,0)] = 120.0 / (c * DT) *
                ((Yx[getx(i,0,0)] - Yx[getx(i-1,0,0)]) / (DX) +
                (-8*boundary_value_v[2]->value(i + offset_x_y,-0.5,0,t) + 9*Yy[gety(i,0,0)] - Yy[gety(i,1,0)]) / (3*DY) +
                (-8*boundary_value_w[4]->value(i + offset_x_z,0,-0.5,t) + 9*Yz[getz(i,0,0)] - Yz[getz(i,0,1)]) / (3*DZ));
            //Y2_p[getp(i,0,0)] = cos((i+ offset_x_y)*DX)*cos(DZ)*(sin(t)-sin(t+64.0/120.0*DT));
        // UPPER FRONT EDGE
        Y2_p[getp(i,0,zSize[2]-1)] = 120.0 / (c * DT) *
            ((Yx[getx(i,0,dim_z-1)] - Yx[getx(i-1,0,dim_z-1)]) / (DX) +
            (-8*boundary_value_v[2]->value(i + offset_x_y,-0.5,NZ,t) + 9*Yy[gety(i,0,dim_z-1)] - Yy[gety(i,1,dim_z-1)]) / (3*DY) +
            (8*boundary_value_w[5]->value(i + offset_x_z,0,NZ-0.5,t) - 9*Yz[getz(i,0,dim_z_z-1)] + Yz[getz(i,0,dim_z_z-2)]) / (3*DZ));
        //Y2_p[getp(i,0,zSize[2]-1)] = cos((i+ offset_x_y)*DX)*cos(DZ)*(sin(t)-sin(t+64.0/120.0*DT));
        }
        if(rby){
            // LOWER BACK EDGE
            Y2_p[getp(i,(zSize[1] - 1),0)] = 120.0 / (c * DT) *
                ((Yx[getx(i,dim_y_x - 1,0)] - Yx[getx(i-1,dim_y_x - 1,0)]) / (DX) +
                (8*boundary_value_v[3]->value(i + offset_x_y,NY-0.5,0,t) - 9*Yy[gety(i,dim_y_y - 1,0)] + Yy[gety(i,dim_y_y - 2,0)]) / (3*DY) +
                (-8*boundary_value_w[4]->value(i + offset_x_z,NY,-0.5,t) + 9*Yz[getz(i,dim_y_z - 1,0)] - Yz[getz(i,dim_y_z - 1,1)]) / (3*DZ));
            //Y2_p[getp(i,(zSize[1] - 1),0)] = cos((i+ offset_x_y)*DX)*cos(DZ)*(sin(t)-sin(t+64.0/120.0*DT));
            // UPPER BACK EDGE
            Y2_p[getp(i,zSize[1]-1,zSize[2]-1)] = 120.0 / (c * DT) *
                ((Yx[getx(i,dim_y_x-1,dim_z-1)] - Yx[getx(i-1,dim_y_x-1,dim_z-1)]) / (DX) +
                (8*boundary_value_v[3]->value(i + offset_x_y,NY-0.5,NZ,t) - 9*Yy[gety(i,dim_y_y-1,dim_z-1)] + Yy[gety(i,dim_y_y-2,dim_z-1)]) / (3*DY) +
                (8*boundary_value_w[5]->value(i + offset_x_z,NY,NZ-0.5,t) - 9*Yz[getz(i,dim_y_z-1,dim_z_z-1)] + Yz[getz(i,dim_y_z-1,dim_z_z-2)]) / (3*DZ));
            //Y2_p[getp(i,zSize[1]-1,zSize[2]-1)] = cos((i+ offset_x_y)*DX)*cos(DZ)*(sin(t)-sin(t+64.0/120.0*DT));
        }
    }
    if(lbx){
        // 4 Y EDGES
        for (int j = lby; j < zSize[1] - rby; j++)
        {
            // LOWER LEFT EDGE
            Y2_p[getp(0,j,0)] = 120.0 / (c * DT) *
                ((-8*boundary_value_u[0]->value(-0.5,j + offset_y_x,0,t) + 9*Yx[getx(0,j,0)] - Yx[getx(1,j,0)]) / (3*DX) +
                (Yy[gety(0,j,0)] - Yy[gety(0,j-1,0)]) / (DY) +
                (-8*boundary_value_w[4]->value(0,j + offset_y_z,-0.5,t) + 9*Yz[getz(0,j,0)] - Yz[getz(0,j,1)]) / (3*DZ));
            // UPPER LEFT EDGE
            Y2_p[getp(0,j,zSize[2]-1)] = 120.0 / (c * DT) *
                ((-8*boundary_value_u[0]->value(-0.5,j + offset_y_x,NZ,t) + 9*Yx[getx(0,j,dim_z-1)] - Yx[getx(1,j,dim_z-1)]) / (3*DX) +
                (Yy[gety(0,j,dim_z-1)] - Yy[gety(0,j-1,dim_z-1)]) / (DY) +
                (8*boundary_value_w[5]->value(0,j,NZ-0.5,t) - 9*Yz[getz(0,j,dim_z_z-1)] + Yz[getz(0,j,dim_z_z-2)]) / (3*DZ));
        }// break to exploit locality
    }
    if(rbx){
        for (int j = lby; j < zSize[1] - rby; j++)
        {
            // LOWER RIGHT EDGE
            Y2_p[getp(zSize[0] - 1,j,0)] = 120.0 / (c * DT) *
                ((8*boundary_value_u[1]->value(NX-0.5,j,0,t) - 9*Yx[getx(dim_x_x - 1,j,0)] + Yx[getx(dim_x_x - 2,j,0)]) / (3*DX) +
                (Yy[gety(dim_x_y - 1,j,0)] - Yy[gety(dim_x_y - 1,j-1,0)]) / (DY) +
                (-8*boundary_value_w[4]->value(NX,j,-0.5,t) + 9*Yz[getz(dim_x_z - 1,j,0)] - Yz[getz(dim_x_z - 1,j,1)]) / (3*DZ));
            // UPPER RIGHT EDGE
            Y2_p[getp(zSize[0] - 1,j,zSize[2]-1)] = 120.0 / (c * DT) *
                ((8*boundary_value_u[1]->value(NX-0.5,j,NZ,t) - 9*Yx[getx(dim_x_x - 1,j,dim_z - 1)] + Yx[getx(dim_x_x - 2,j,dim_z - 1)]) / (3*DX) +
                (Yy[gety(dim_x_y - 1,j,dim_z - 1)] - Yy[gety(dim_x_y - 1,j-1,dim_z - 1)]) / (DY) +
                (8*boundary_value_w[5]->value(NX,j,NZ-0.5,t) - 9*Yz[getz(dim_x_z - 1,j,dim_z_z - 1)] + Yz[getz(dim_x_z - 1,j,dim_z - 2)]) / (3*DZ));
        }
    }
    if(lbx && lby){
        // 4 Z EDGES
        for (int k = 1; k < zSize[2] - 1; k++)
        {
            // FRONT LEFT EDGE
            Y2_p[getp(0,0,k)] = 120.0 / (c * DT) *
                ((-8*boundary_value_u[0]->value(-0.5,0,k,t) + 9*Yx[getx(0,0,k)] - Yx[getx(1,0,k)]) / (3*DX) +
                (-8*boundary_value_v[2]->value(0,-0.5,k,t) + 9*Yy[gety(0,0,k)] - Yy[getx(0,1,k)]) / (3*DY) +
                (Yz[getz(0,0,k)] - Yz[getz(0,0,k-1)]) / (DZ));
        }// break to exploit locality
    }
    if(lbx && rby){
        for (int k = 1; k < zSize[2] - 1; k++)
        {
            // BACK LEFT EDGE
            Y2_p[getp(0,zSize[1] - 1, k)] = 120.0 / (c * DT) *
                ((-8*boundary_value_u[0]->value(-0.5,NY,k,t) + 9*Yx[getx(0,dim_y_x - 1,k)] - Yx[getx(1,dim_y_x - 1,k)]) / (3*DX) +
                (8*boundary_value_v[3]->value(0,NY-0.5,k,t) - 9*Yy[gety(0,dim_y_y - 1,k)] + Yy[gety(0,dim_y_y - 2,k)]) / (3*DY) +
                (Yz[getz(0,dim_y_z - 1,k)] - Yz[getz(0,dim_y_z - 1,k-1)]) / (DZ));
        }// break to exploit locality
    }
    if(rbx && lby){
        for (int k = 1; k < zSize[2] - 1; k++)
        {
            // FRONT RIGHT EDGE
            Y2_p[getp(zSize[0] - 1,0,k)] = 120.0 / (c * DT) *
                ((8*boundary_value_u[1]->value(NX-0.5,0,k,t) - 9*Yx[getx(dim_x_x-1,0,k)] + Yx[getx(dim_x_x-2,0,k)]) / (3*DX) +
                (-8*boundary_value_v[2]->value(NX,-0.5,k,t) + 9*Yy[gety(dim_x_y-1,0,k)] - Yy[gety(dim_x_y-1,1,k)]) / (3*DY) +
                (Yz[getz(dim_x_z-1,0,k)] - Yz[getz(dim_x_z-1,0,k-1)]) / (DZ));
        }// break to exploit locality
    }
    if(rbx && rby){
        for (int k = 1; k < zSize[2] - 1; k++)
        {
            // BACK RIGHT EDGE
            Y2_p[getp(zSize[0]-1, zSize[1]-1, k)] = 120.0 / (c * DT) *
                ((8*boundary_value_u[1]->value(NX-0.5,NY,k,t) - 9*Yx[getx(dim_x_x-1,dim_y_x-1,k)] + Yx[getx(dim_x_x-2,dim_y_x-1,k)]) / (3*DX) +
                (8*boundary_value_v[3]->value(NX,NY-0.5,k,t) - 9*Yy[gety(dim_x_y-1,dim_y_y-1,k)] + Yy[gety(dim_x_y-1,dim_y_y-2,k)]) / (3*DY) +
                (Yz[getz(dim_x_z-1,dim_y_z-1,k)] - Yz[getz(dim_x_z-1,dim_y_z-1,k-1)]) / (DZ));
        }
    }
    // 8 VERTICES
    if(lbx && lby){
        // LOWER FRONT LEFT VERTEX (0,0,0)
        Y2_p[getp(0, 0, 0)] = 120.0 / (c * DT) *
            ((-8*boundary_value_u[0]->value(-0.5,0,0,t) + 9*Yx[getx(0,0,0)] - Yx[getx(1,0,0)]) / (3*DX) +
            (-8*boundary_value_v[2]->value(0,-0.5,0,t) + 9*Yy[gety(0,0,0)] - Yy[gety(0,1,0)]) / (3*DY) +
            (-8*boundary_value_w[4]->value(0,0,-0.5,t) + 9*Yz[getz(0,0,0)] - Yz[getz(0,0,1)]) / (3*DZ));
    }
    if(rbx && lby){
        // LOWER FRONT RIGHT VERTEX (1,0,0)
        Y2_p[getp(zSize[0]-1, 0, 0)] = 120.0 / (c * DT) *
            ((8*boundary_value_u[1]->value(NX-0.5,0,0,t) - 9*Yx[getx(dim_x_x-1,0,0)] + Yx[getx(dim_x_x-2,0,0)]) / (3*DX) +
            (-8*boundary_value_v[2]->value(NX,-0.5,0,t) + 9*Yy[gety(dim_x_y-1,0,0)] - Yy[gety(dim_x_y-1,1,0)]) / (3*DY) +
            (-8*boundary_value_w[4]->value(NX,0,-0.5,t) + 9*Yz[getz(dim_x_z-1,0,0)] - Yz[getz(dim_x_z-1,0,1)]) / (3*DZ));
    }
    if(lbx && rby){
        // LOWER BACK LEFT VERTEX (0,1,0)
        Y2_p[getp(0, zSize[1]-1, 0)] = 120.0 / (c * DT) *
            ((-8*boundary_value_u[0]->value(-0.5,NY,0,t) + 9*Yx[getx(0,dim_y_x-1,0)] - Yx[getx(1,dim_y_x-1,0)]) / (3*DX) +
            (8*boundary_value_v[3]->value(0,NY-0.5,0,t) - 9*Yy[gety(0,dim_y_y-1,0)] + Yy[gety(0,dim_y_y-2,0)]) / (3*DY) +
            (-8*boundary_value_w[4]->value(0,NY,-0.5,t) + 9*Yz[getz(0,dim_y_z-1,0)] - Yz[getz(0,dim_y_z-1,1)]) / (3*DZ));
    }
    if(lbx && lby){
        // UPPER FRONT LEFT VERTEX (0,0,1)
        Y2_p[getp(0, 0, zSize[2]-1)] = 120.0 / (c * DT) *
            ((-8*boundary_value_u[0]->value(-0.5,0,NZ,t) + 9*Yx[getx(0,0,dim_z-1)] - Yx[getx(1,0,dim_z-1)]) / (3*DX) +
            (-8*boundary_value_v[2]->value(0,-0.5,NZ,t) + 9*Yy[gety(0,0,dim_z-1)] - Yy[gety(0,1,dim_z-1)]) / (3*DY) +
            (8*boundary_value_w[5]->value(0,0,NZ-0.5,t) - 9*Yz[getz(0,0,dim_z_z-1)] + Yz[getz(0,0,dim_z_z-2)]) / (3*DZ));
    }
    if(lbx && rby){
        // UPPER BACK LEFT VERTEX (0,1,1)
        Y2_p[getp(0, zSize[1]-1, zSize[2]-1)] = 120.0 / (c * DT) *
            ((-8*boundary_value_u[0]->value(-0.5,NY,NZ,t) + 9*Yx[getx(0,dim_y_x-1,dim_z-1)] - Yx[getx(1,dim_y_x-1,dim_z-1)]) / (3*DX) +
            (8*boundary_value_v[3]->value(0,NY-0.5,NZ,t) - 9*Yy[gety(0,dim_y_y-1,dim_z-1)] + Yy[gety(0,dim_y_y-2,dim_z-1)]) / (3*DY) +
            (8*boundary_value_w[5]->value(0,NY,NZ-0.5,t) - 9*Yz[getz(0,dim_y_z-1,dim_z_z-1)] + Yz[getz(0,dim_y_z-1,dim_z_z-2)]) / (3*DZ));
    }
    if(rbx && lby){
        // UPPER FRONT RIGHT VERTEX (1,0,1)
        Y2_p[getp(zSize[0]-1, 0, zSize[2]-1)] = 120.0 / (c * DT) *
            ((8*boundary_value_u[1]->value(NX-0.5,0,NZ,t) - 9*Yx[getx(dim_x_x-1,0,dim_z-1)] + Yx[getx(dim_x_x-2,0,dim_z-1)]) / (3*DX) +
            (-8*boundary_value_v[2]->value(NX,-0.5,NZ,t) + 9*Yy[gety(dim_x_y-1,0,dim_z-1)] - Yy[gety(dim_x_y-1,1,dim_z-1)]) / (3*DY) +
            (8*boundary_value_w[5]->value(NX,0,NZ-0.5,t) - 9*Yz[getz(dim_x_z-1,0,dim_z_z-1)] + Yz[getz(dim_x_z-1,0,dim_z_z-2)]) / (3*DZ));
    }
    if(rbx && rby){
        // LOWER BACK RIGHT VERTEX (1,1,0)
        Y2_p[getp(zSize[0]-1, zSize[1]-1, 0)] = 120.0 / (c * DT) *
            ((8*boundary_value_u[1]->value(NX-0.5,NY,0,t) - 9*Yx[getx(dim_x_x-1,dim_y_x-1,0)] + Yx[getx(dim_x_x-2,dim_y_x-1,0)]) / (3*DX) +
            (8*boundary_value_v[3]->value(NX,NY-0.5,0,t) - 9*Yy[gety(dim_x_y-1,dim_y_y-1,0)] + Yy[gety(dim_x_y-1,dim_y_y-2,0)]) / (3*DY) +
            (-8*boundary_value_w[4]->value(NX,NY,-0.5,t) + 9*Yz[getz(dim_x_z-1,dim_y_z-1,0)] - Yz[getz(dim_x_z-1,dim_y_z-1,1)]) / (3*DZ));
        // UPPER BACK RIGHT VERTEX (1,1,1)
        Y2_p[getp(zSize[0]-1, zSize[1]-1, zSize[2]-1)] = 120.0 / (c * DT) *
            ((8*boundary_value_u[1]->value(NX-0.5,NY,NZ,t) - 9*Yx[getx(dim_x_x-1,dim_y_x-1,dim_z-1)] + Yx[getx(dim_x_x-2,dim_y_x-1,dim_z-1)]) / (3*DX) +
            (8*boundary_value_v[3]->value(NX,NY-0.5,NZ,t) - 9*Yy[gety(dim_x_y-1,dim_y_y-1,dim_z-1)] + Yy[gety(dim_x_y-1,dim_y_y-2,dim_z-1)]) / (3*DY) +
            (8*boundary_value_w[5]->value(NX,NY,NZ-0.5,t) - 9*Yz[getz(dim_x_z-1,dim_y_z-1,dim_z_z-1)] + Yz[getz(dim_x_z-1,dim_y_z-1,dim_z_z-2)]) / (3*DZ));
    }
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
