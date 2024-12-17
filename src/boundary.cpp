#include "boundary.hpp"

/**
 * @brief The method is called by the program multiple during the time step, in order to update the values of the boundaries at each 
 * requested time t, calculating the approximated ones too.
 * 
 * @param Yx Boundary x velocities or the Y size_termediate function related to the x direction.
 * @param Yy Boundary y velocities or the Y size_termediate function related to the y direction.
 * @param Yz Boundary z velocities or the Y size_termediate function related to the z direction.
 * @param t Time of the time discretization we are considering.
*/
void Boundary::update_boundary(std::array<Real, NX *(NY + 1) * (NZ + 1)> &Yx, std::array<Real, (NX + 1) * NY *(NZ + 1)> &Yy, std::array<Real, (NX + 1) * (NY + 1) * NZ> &Yz, std::array<Real, (NX+1) * (NY + 1) * (NZ + 1)> &P, Real t)
{
    int face/* , i, j, k */;
    // X
    // LEFT FACE
    face = 2;
    for (int k=0; k < NZ+1; k++)
    {
        Yx[k] = boundary_value_u[face]->value(0,0,k,t);
    }

    for (int j=1; j < NY; j++)
    {
        face = 4;
        Yx[j*(NZ+1)] = boundary_value_u[face]->value(0,j,0,t);
        face = 0;
        for(int k = 1; k < NZ; k++)
        {
            Yx[j*(NZ+1) + k] = approximate_boundary_u(0,j,k,t,face,1);
        }
        face = 5;
        Yx[j*(NZ+1) + NZ] = boundary_value_u[face]->value(0,j,NZ,t);
    }

    face = 3;
    for (int k=0; k < NZ+1; k++)
    {
        Yx[NY * (NZ+1) + k] = boundary_value_u[face]->value(0,NY,k,t);
    }

    // MIDDLE POINTS
    for (int i=1; i < NX-1; i++)
    {
        face = 2;
        for (int k=0; k < NZ+1; k++)
        {
            Yx[i * (NY+1) * (NZ+1) + k] = boundary_value_u[face]->value(i,0,k,t);
        }
        
        for (int j=1; j < NY; j++)
        {
            face = 4;
            Yx[(NY+1) * (NZ+1) * i + j*(NZ+1)] = boundary_value_u[face]->value(i,j,0,t);

            face = 5;
            Yx[(NY+1) * (NZ+1) * i + j*(NZ+1) + NZ] = boundary_value_u[face]->value(i,j,NZ,t);
        }

        face = 3;
        for (int k=0; k < NZ+1; k++)
        {
            Yx[NY * (NZ + 1) + (NY+1) * (NZ+1) * i + k] = boundary_value_u[face]->value(i,NY,k,t);
        }
    }

    // RIGHT FACE
    face = 2;
    for (int k=0; k < NZ+1; k++)
    {
        Yx[(NX-1) * (NY+1) * (NZ+1) + k] = boundary_value_u[face]->value(NX-1,0,k,t);
    }

    for (int j=1; j < NY; j++)
    {
        face = 4;
        Yx[(NX-1) * (NY+1) * (NZ+1) + j*(NZ+1)] = boundary_value_u[face]->value(NX-1,j,0,t);

        face = 1;
        for(int k = 1; k < NZ; k++)
        {
            Yx[(NX-1)*(NY+1)*(NZ+1) + j*(NZ+1) + k] = approximate_boundary_u(NX,j,k,t,face,-1);
        }

        face = 5;
        Yx[(NX-1) * (NY+1) * (NZ+1) + j*(NZ+1) + NZ] = boundary_value_u[face]->value(NX-1,j,NZ,t);
    }

    face = 3;
    for (int k=0; k < NZ+1; k++)
    {
        Yx[(NX-1) * (NY+1) * (NZ+1) + NY*(NZ+1) + k] = boundary_value_u[face]->value(NX-1,NY,k,t);
    }


    // Y
    // LEFT FACE
    face = 0;
    for (int j = 0; j < NY; j++)
    {
        for(int k = 0; k < NZ+1; k++)
        {
            Yy[j*(NZ+1) + k] = boundary_value_v[face]->value(0,j,k,t);
        }
    }

    // MIDDLE POINTS
    for (int i=1; i < NX; i++)
    {
        face = 4;
        Yy[NY * (NZ+1) * i] = boundary_value_v[face]->value(i,0,0,t); 
        face = 2;
        for (int k=1; k < NZ; k++)
        {
            Yy[i * NY * (NZ+1) + k] = approximate_boundary_v(i,0,k,t,face,1);
        }
        face = 5;
        Yy[NY * (NZ+1) * i + NZ] = boundary_value_v[face]->value(i,0,NZ,t); 

        for (int j=1; j < NY-1; j++)
        {
            face = 4;
            Yy[NY * (NZ+1) * i + j*(NZ+1)] = boundary_value_v[face]->value(i,j,0,t); 

            face = 5;
            Yy[NY * (NZ+1) * i + j*(NZ+1) + NZ] = boundary_value_v[face]->value(i,j,NZ,t); 
        }
        
        face = 4;
        Yy[NY * (NZ+1) * i + (NY-1)*(NZ+1)] = boundary_value_v[face]->value(i,NY-1,0,t); 
        face = 3;
        for(int k = 1; k < NZ; k++)
        {
            Yy[(NY-1) * (NZ+1) + NY * (NZ+1) * i + k] = approximate_boundary_v(i,NY,k,t,face,-1);
        }
        face = 5;
        Yy[NY * (NZ+1) * i + (NY-1)*(NZ+1) + NZ] = boundary_value_v[face]->value(i,NY-1,NZ,t); 
    }

    // RIGHT FACE
    face = 1;
    for (int j = 0; j < NY; j++)
    {
        for(int k = 0; k < NZ+1; k++)
        {
            Yy[NX*NY*(NZ+1) + j*(NZ+1) + k] = boundary_value_v[face]->value(NX,j,k,t);
        }
    }


    // Z
    // LEFT FACE 
    face = 0;
    for (int j = 0; j < NY+1; j++)
    {
        for(int k = 0; k < NZ; k++)
        {
            Yz[j*NZ + k] = boundary_value_w[face]->value(0,j,k,t); 
        }
    }

    // MIDDLE POINTS
    for (int i=1; i < NX; i++)
    {
        face = 2;
        for (int k=0; k < NZ; k++)
        {
            Yz[i * (NY+1) * NZ + k] = boundary_value_w[face]->value(i,0,k,t);
        }

        for(int j = 1; j < NY; j++)
        {
            face = 4;
            Yz[(NY+1) * NZ * i + j*NZ] = approximate_boundary_w(i,j,0,t,face,1);/* boundary_value_w[face]->value(i,j,0,t); */
            
            face = 5;
            Yz[(NY+1) * NZ * i + j*NZ + (NZ - 1)] = approximate_boundary_w(i,j,NZ,t,face,-1); /* boundary_value_w[face]->value(i,j,NZ,t);  */
        } 

        face = 3;
        for(int k = 0; k < NZ+1; k++)
        {
            Yz[NY * NZ + i * (NY+1) * NZ + k] = boundary_value_w[face]->value(i,NY,k,t);
        }
    }

    // RIGHT FACE
    face = 1;
    for (int j = 0; j < NY+1; j++)
    {
        for(int k = 0; k < NZ; k++)
        {
            Yz[NX*(NY+1)*NZ + j*NZ + k] = boundary_value_w[face]->value(NX,j,k,t); 
        }
    }
}

/**
 * @brief Calculate the approximate value of the x velocity in a given posize_t.
 * 
 * @param x,y,z Coordinates of the position in the 3D mesh.
 * @param t Time of the time discretization we are considering.
 * @param face Face of the mesh we are considering. In order: 0 (LEFT), 1 (RIGHT), 2 (FRONT), 3 (BACK), 4 (LOWER), 5 (UPPER)
 * @param side Int used to calculate the correct approximate value. It's 1 for faces 0, 2 and 4, -1 for the others.
 * 
 * @return the approximate value.
*/
Real Boundary::approximate_boundary_u(int x, int y, int z, Real t, int face, int side)
{

    Real dv = (boundary_value_v[face]->value(x, y, z, t) - boundary_value_v[face]->value(x, y - 1.0, z, t)) / DY;
    Real dw = (boundary_value_w[face]->value(x, y, z, t) - boundary_value_w[face]->value(x, y, z - 1.0, t)) / DZ;

    return boundary_value_u[face]->value((x - 0.5 /*(DX/2.0)*/), y, z, t) - (dv + dw) * (DX / 2) * side;
}

/**
 * @brief Calculate the approximate value of the y velocity in a given posize_t.
 * 
 * @param x,y,z Coordinates of the position in the 3D mesh.
 * @param t Time of the time discretization we are considering.
 * @param face Face of the mesh we are considering. In order: 0 (LEFT), 1 (RIGHT), 2 (FRONT), 3 (BACK), 4 (LOWER), 5 (UPPER)
 * @param side Int used to calculate the correct approximate value. It's 1 for faces 0, 2 and 4, -1 for the others.
 * 
 * @return the approximate value.
*/
Real Boundary::approximate_boundary_v(int x, int y, int z, Real t, int face, int side)
{
    Real du = ((boundary_value_u[face]->value(x, y, z, t)) -
                 ((boundary_value_u[face]->value(x - 1.0, y, z, t)))) /
                (DX);

    Real dw = ((boundary_value_w[face]->value(x, y, z, t)) -
                 ((boundary_value_w[face]->value(x, y, z - 1.0, t)))) /
                (DZ);    
    return boundary_value_v[face]->value(x, y - 0.5, z, t) - (du + dw) * (DY / 2.0) * side;
}

/**
 * @brief Calculate the approximate value of the z velocity in a given posize_t.
 * 
 * @param x,y,z Coordinates of the position in the 3D mesh.
 * @param t Time of the time discretization we are considering.
 * @param face Face of the mesh we are considering. In order: 0 (LEFT), 1 (RIGHT), 2 (FRONT), 3 (BACK), 4 (LOWER), 5 (UPPER)
 * @param side Int used to calculate the correct approximate value. It's 1 for faces 0, 2 and 4, -1 for the others.
 * 
 * @return the approximate value.
*/
Real Boundary::approximate_boundary_w(int x, int y, int z, Real t, int face, int side)
{
    Real du = ((boundary_value_u[face]->value(x, y, z, t)) -
                 ((boundary_value_u[face]->value(x - 1.0, y, z, t)))) /
                (DX);

    Real dv = (boundary_value_v[face]->value(x, y, z, t) - boundary_value_v[face]->value(x, y - 1.0, z, t)) / DY;
    return boundary_value_w[face]->value(x, y, z - 0.5, t) - (du + dv) * (DZ / 2.0) * side;
}

/**
 * @brief Add the given function to the selected direction.
 * 
 * @param direction Direction U (length), V (width) or W (height) of the boundary.
 * @param x Function to assign to the boundary
*/
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
