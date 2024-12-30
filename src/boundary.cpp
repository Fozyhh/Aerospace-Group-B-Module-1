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
void Boundary::update_boundary(std::array<Real, NX *(NY + 1) * (NZ + 1)> &Yx, std::array<Real, (NX + 1) * NY *(NZ + 1)> &Yy, std::array<Real, (NX + 1) * (NY + 1) * NZ> &Yz, Real t)
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
        for(int k = 0; k < NZ; k++)
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
 * @brief Calculate the divergence of the velocity at the boundary, used to calculate the input of the poisson solver.
 * 
 * @param Yx  Ux input grid.
 * @param Yy  Uy input grid.
 * @param Yz  Uz input grid.
 * @param Y2_p output grid.
 * @param t Time of the time discretization we are considering.
 * @param c Time constant of the RK intermidiate step (64, 16 or 40).
*/
void Boundary::divergence(std::array<Real, NX * (NY + 1) * (NZ + 1)> &Yx, std::array<Real, (NX + 1) * NY * (NZ + 1)> &Yy, std::array<Real, (NX + 1) * (NY + 1) * NZ> &Yz, std::array<Real, (NX + 1) * (NY + 1) * (NZ + 1)> &Y2_p, Real t, Real c)
{
    // is the denominator 3*DX correct? -> 2*DX ?
    // LEFT FACE
    for (int j = 1; j < NY; j++)
    {
        for(int k = 1; k < NZ; k++)
        {
            Y2_p[j * (NZ + 1) + k] = 120.0 / (c * DT) *
                ((-8*boundary_value_u[0]->value(-0.5,j,k,t) + 9*Yx[j * (NZ + 1) + k] - Yx[1 * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k]) / (3*DX) +
                (Yy[j * (NZ + 1) + k] - Yy[(j-1) * (NZ + 1) + k]) / (DY) +
                (Yz[j * NZ + k] - Yz[j * NZ + k-1]) / (DZ));
        }
    }
    // RIGHT FACE
    for (int j = 1; j < NY; j++)
    {
        for(int k = 1; k < NZ; k++)
        {
            Y2_p[(NX) * (NY + 1) * (NZ + 1) + j * (NZ+1) + k] = 120.0 / (c * DT) *
                ((8*boundary_value_u[1]->value(NX-0.5,j,k,t) - 9*Yx[(NX-1) * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k] + Yx[(NX-2) * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k]) / (3*DX) +
                (Yy[(NX) * NY * (NZ + 1) + j * (NZ + 1) + k] - Yy[(NX) * NY * (NZ + 1) + (j-1) * (NZ + 1) + k]) / (DY) +
                (Yz[(NX) * (NY + 1) * NZ + j * NZ + k] - Yz[(NX) * (NY + 1) * NZ + j * NZ + k-1]) / (DZ));
        }
    }
    // FRONT FACE
    for (int i = 1; i < NX; i++)
    {
        for(int k = 1; k < NZ; k++)
        {
            Y2_p[i * (NY + 1) * (NZ + 1) + k] = 120.0 / (c * DT) *
                ((Yx[i * (NY + 1) * (NZ + 1) + k] - Yx[(i-1) * (NY + 1) * (NZ + 1) + k]) / (DX) +
                (-8*boundary_value_v[2]->value(i,-0.5,k,t) + 9*Yy[i * NY * (NZ + 1) + k] - Yy[i * NY * (NZ + 1) + 1 * (NZ + 1) + k]) / (3*DY) +
                (Yz[i * (NY + 1) * NZ + k] - Yz[i * (NY + 1) * NZ + k-1]) / (DZ));
        }
    }
    // BACK FACE
    for (int i = 1; i < NX; i++)
    {
        for(int k = 1; k < NZ; k++)
        {
            Y2_p[i * (NY + 1) * (NZ + 1) + (NY) * (NZ + 1) + k] = 120.0 / (c * DT) *
                ((Yx[i * (NY + 1) * (NZ + 1) + (NY) * (NZ + 1) + k] - Yx[(i-1) * (NY + 1) * (NZ + 1) + (NY) * (NZ + 1) + k]) / (DX) +
                (8*boundary_value_v[3]->value(i,NY-0.5,k,t) - 9*Yy[i * NY * (NZ + 1) + (NY-1) * (NZ + 1) + k] + Yy[i * NY * (NZ + 1) + (NY-2) * (NZ + 1) + k]) / (3*DY) +
                (Yz[i * (NY + 1) * NZ + (NY) * NZ + k] - Yz[i * (NY + 1) * NZ + (NY) * NZ + k-1]) / (DZ));
        }
    }
    // LOWER FACE
    for (int i = 1; i < NX; i++)
    {
        for(int j = 1; j < NY; j++)
        {
            Y2_p[i * (NY + 1) * (NZ + 1) + j * (NZ + 1)] = 120.0 / (c * DT) *
                ((Yx[i * (NY + 1) * (NZ + 1) + j * (NZ + 1)] - Yx[(i-1) * (NY + 1) * (NZ + 1) + j * (NZ + 1)]) / (DX) +
                (Yy[i * NY * (NZ + 1) + j * (NZ + 1)] - Yy[i * NY * (NZ + 1) + (j-1) * (NZ + 1)]) / (DY) +
                (-8*boundary_value_w[4]->value(i,j,-0.5,t) + 9*Yz[i * (NY + 1) * NZ + j * NZ] - Yz[i * (NY + 1) * NZ + j * NZ + 1]) / (3*DZ));
        }
    }
    // UPPER FACE
    for (int i = 1; i < NX; i++)
    {
        for(int j = 1; j < NY; j++)
        {
            Y2_p[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + (NZ)] = 120.0 / (c * DT) *
                ((Yx[i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + (NZ)] - Yx[(i-1) * (NY + 1) * (NZ + 1) + j * (NZ + 1) + (NZ)]) / (DX) +
                (Yy[i * NY * (NZ + 1) + j * (NZ + 1) + (NZ)] - Yy[i * NY * (NZ + 1) + (j-1) * (NZ + 1) + (NZ)]) / (DY) +
                (8*boundary_value_w[5]->value(i,j,NZ-0.5,t) - 9*Yz[i * (NY + 1) * NZ + j * NZ + (NZ-1)] + Yz[i * (NY + 1) * NZ + j * NZ + (NZ-2)]) / (3*DZ));
        }
    }
    // 4 X EDGES
    for (int i = 1; i < NX; i++)
    {
        // LOWER FRONT EDGE
        Y2_p[i * (NY + 1) * (NZ + 1)] = 120.0 / (c * DT) *
            ((Yx[i * (NY + 1) * (NZ + 1)] - Yx[(i-1) * (NY + 1) * (NZ + 1)]) / (DX) +
            (-8*boundary_value_v[2]->value(i,-0.5,0,t) + 9*Yy[i * NY * (NZ + 1)] - Yy[i * NY * (NZ + 1) + 1 * (NZ + 1)]) / (3*DY) +
            (-8*boundary_value_w[4]->value(i,0,-0.5,t) + 9*Yz[i * (NY + 1) * NZ] - Yz[i * (NY + 1) * NZ + 1]) / (3*DZ));
        // UPPER FRONT EDGE
        Y2_p[i * (NY + 1) * (NZ + 1) + (NZ)] = 120.0 / (c * DT) *
            ((Yx[i * (NY + 1) * (NZ + 1) + (NZ)] - Yx[(i-1) * (NY + 1) * (NZ + 1) + (NZ)]) / (DX) +
            (-8*boundary_value_v[2]->value(i,-0.5,NZ,t) + 9*Yy[i * NY * (NZ + 1) + (NZ)] - Yy[i * NY * (NZ + 1) + 1 * (NZ + 1) + (NZ)]) / (3*DY) +
            (8*boundary_value_w[5]->value(i,0,NZ-0.5,t) - 9*Yz[i * (NY + 1) * NZ + (NZ-1)] + Yz[i * (NY + 1) * NZ + (NZ-2)]) / (3*DZ));
        // LOWER BACK EDGE
        Y2_p[i * (NY + 1) * (NZ + 1) + (NY) * (NZ + 1)] = 120.0 / (c * DT) *
            ((Yx[i * (NY + 1) * (NZ + 1) + (NY) * (NZ + 1)] - Yx[(i-1) * (NY + 1) * (NZ + 1) + (NY) * (NZ + 1)]) / (DX) +
            (8*boundary_value_v[3]->value(i,NY-0.5,0,t) - 9*Yy[i * NY * (NZ + 1) + (NY-1) * (NZ + 1)] + Yy[i * NY * (NZ + 1) + (NY-2) * (NZ + 1)]) / (3*DY) +
            (-8*boundary_value_w[4]->value(i,NY,-0.5,t) + 9*Yz[i * (NY + 1) * NZ + (NY) * NZ] - Yz[i * (NY + 1) * NZ + (NY) * NZ + 1]) / (3*DZ));
        // UPPER BACK EDGE
        Y2_p[i * (NY + 1) * (NZ + 1) + (NY) * (NZ + 1) + (NZ)] = 120.0 / (c * DT) *
            ((Yx[i * (NY + 1) * (NZ + 1) + (NY) * (NZ + 1) + (NZ)] - Yx[(i-1) * (NY + 1) * (NZ + 1) + (NY) * (NZ + 1) + (NZ)]) / (DX) +
            (8*boundary_value_v[3]->value(i,NY-0.5,NZ,t) - 9*Yy[i * NY * (NZ + 1) + (NY-1) * (NZ + 1) + (NZ)] + Yy[i * NY * (NZ + 1) + (NY-2) * (NZ + 1) + (NZ)]) / (3*DY) +
            (8*boundary_value_w[5]->value(i,NY,NZ-0.5,t) - 9*Yz[i * (NY + 1) * NZ + (NY) * NZ + (NZ-1)] + Yz[i * (NY + 1) * NZ + (NY) * NZ + (NZ-2)]) / (3*DZ));
    }
    // 4 Y EDGES
    for (int j = 1; j < NY; j++)
    {
        // LOWER LEFT EDGE
        Y2_p[j * (NZ + 1)] = 120.0 / (c * DT) *
            ((-8*boundary_value_u[0]->value(-0.5,j,0,t) + 9*Yx[j * (NZ + 1)] - Yx[1 * (NY + 1) * (NZ + 1) + j * (NZ + 1)]) / (3*DX) +
            (Yy[j * (NZ + 1)] - Yy[(j-1) * (NZ + 1)]) / (DY) +
            (-8*boundary_value_w[4]->value(0,j,-0.5,t) + 9*Yz[j * NZ] - Yz[j * NZ + 1]) / (3*DZ));
        // UPPER LEFT EDGE
        Y2_p[j * (NZ + 1) + (NZ)] = 120.0 / (c * DT) *
            ((-8*boundary_value_u[0]->value(-0.5,j,NZ,t) + 9*Yx[j * (NZ + 1) + (NZ)] - Yx[1 * (NY + 1) * (NZ + 1) + j * (NZ + 1) + (NZ)]) / (3*DX) +
            (Yy[j * (NZ + 1) + (NZ)] - Yy[(j-1) * (NZ + 1) + (NZ)]) / (DY) +
            (8*boundary_value_w[5]->value(0,j,NZ-0.5,t) - 9*Yz[j * NZ + (NZ-1)] + Yz[j * NZ + (NZ-2)]) / (3*DZ));
    }// break to exploit locality
    for (int j = 1; j < NY; j++)
    {
        // LOWER RIGHT EDGE
        Y2_p[(NX) * (NY + 1) * (NZ + 1) + j * (NZ + 1)] = 120.0 / (c * DT) *
            ((8*boundary_value_u[1]->value(NX-0.5,j,0,t) - 9*Yx[(NX-1) * (NY + 1) * (NZ + 1) + j * (NZ + 1)] + Yx[(NX-2) * (NY + 1) * (NZ + 1) + j * (NZ + 1)]) / (3*DX) +
            (Yy[(NX) * NY * (NZ + 1) + j * (NZ + 1)] - Yy[(NX) * NY * (NZ + 1) + (j-1) * (NZ + 1)]) / (DY) +
            (-8*boundary_value_w[4]->value(NX,j,-0.5,t) + 9*Yz[(NX) * (NY + 1) * NZ + j * NZ] - Yz[(NX) * (NY + 1) * NZ + j * NZ + 1]) / (3*DZ));
        // UPPER RIGHT EDGE
        Y2_p[(NX) * (NY + 1) * (NZ + 1) + j * (NZ + 1) + (NZ)] = 120.0 / (c * DT) *
            ((8*boundary_value_u[1]->value(NX-0.5,j,NZ,t) - 9*Yx[(NX-1) * (NY + 1) * (NZ + 1) + j * (NZ + 1) + (NZ)] + Yx[(NX-2) * (NY + 1) * (NZ + 1) + j * (NZ + 1) + (NZ)]) / (3*DX) +
            (Yy[(NX) * NY * (NZ + 1) + j * (NZ + 1) + (NZ)] - Yy[(NX) * NY * (NZ + 1) + (j-1) * (NZ + 1) + (NZ)]) / (DY) +
            (8*boundary_value_w[5]->value(NX,j,NZ-0.5,t) - 9*Yz[(NX) * (NY + 1) * NZ + j * NZ + (NZ-1)] + Yz[(NX) * (NY + 1) * NZ + j * NZ + (NZ-2)]) / (3*DZ));
    }
    // 4 Z EDGES
    for (int k = 1; k < NZ; k++)
    {
        // FRONT LEFT EDGE
        Y2_p[k] = 120.0 / (c * DT) *
            ((-8*boundary_value_u[0]->value(-0.5,0,k,t) + 9*Yx[k] - Yx[1 * (NY + 1) * (NZ + 1) + k]) / (3*DX) +
            (-8*boundary_value_v[2]->value(0,-0.5,k,t) + 9*Yy[k] - Yy[1 * (NZ + 1) + k]) / (3*DY) +
            (Yz[k] - Yz[k-1]) / (DZ));
    }// break to exploit locality
    for (int k = 1; k < NZ; k++)
    {
        // BACK LEFT EDGE
        Y2_p[(NY) * (NZ + 1) + k] = 120.0 / (c * DT) *
            ((-8*boundary_value_u[0]->value(-0.5,NY,k,t) + 9*Yx[(NY) * (NZ + 1) + k] - Yx[1 * (NY + 1) * (NZ + 1) + (NY) * (NZ + 1) + k]) / (3*DX) +
            (8*boundary_value_v[3]->value(0,NY-0.5,k,t) - 9*Yy[(NY-1) * (NZ + 1) + k] + Yy[(NY-2) * (NZ + 1) + k]) / (3*DY) +
            (Yz[(NY) * NZ + k] - Yz[(NY) * NZ + k-1]) / (DZ));
    }// break to exploit locality
    for (int k = 1; k < NZ; k++)
    {
        // FRONT RIGHT EDGE
        Y2_p[(NX) * (NY + 1) * (NZ + 1) + k] = 120.0 / (c * DT) *
            ((8*boundary_value_u[1]->value(NX-0.5,0,k,t) - 9*Yx[(NX-1) * (NY + 1) * (NZ + 1) + k] + Yx[(NX-2) * (NY + 1) * (NZ + 1) + k]) / (3*DX) +
            (-8*boundary_value_v[2]->value(NX,-0.5,k,t) + 9*Yy[(NX) * NY * (NZ + 1) + k] - Yy[(NX) * NY * (NZ + 1) + 1 * (NZ + 1) + k]) / (3*DY) +
            (Yz[(NX) * (NY + 1) * NZ + k] - Yz[(NX) * (NY + 1) * NZ + k-1]) / (DZ));
    }// break to exploit locality
    for (int k = 1; k < NZ; k++)
    {
        // BACK RIGHT EDGE
        Y2_p[(NX) * (NY + 1) * (NZ + 1) + (NY) * (NZ + 1) + k] = 120.0 / (c * DT) *
            ((8*boundary_value_u[1]->value(NX-0.5,NY,k,t) - 9*Yx[(NX-1) * (NY + 1) * (NZ + 1) + (NY) * (NZ + 1) + k] + Yx[(NX-2) * (NY + 1) * (NZ + 1) + (NY) * (NZ + 1) + k]) / (3*DX) +
            (8*boundary_value_v[3]->value(NX,NY-0.5,k,t) - 9*Yy[(NX) * NY * (NZ + 1) + (NY-1) * (NZ + 1) + k] + Yy[(NX) * NY * (NZ + 1) + (NY-2) * (NZ + 1) + k]) / (3*DY) +
            (Yz[(NX) * (NY + 1) * NZ + (NY) * NZ + k] - Yz[(NX) * (NY + 1) * NZ + (NY) * NZ + k-1]) / (DZ));
    }
    // 8 VERTICES
    // LOWER FRONT LEFT VERTEX (0,0,0)
    Y2_p[0] = 120.0 / (c * DT) *
        ((-8*boundary_value_u[0]->value(-0.5,0,0,t) + 9*Yx[0] - Yx[1 * (NY + 1) * (NZ + 1)]) / (3*DX) +
        (-8*boundary_value_v[2]->value(0,-0.5,0,t) + 9*Yy[0] - Yy[1 * (NZ + 1)]) / (3*DY) +
        (-8*boundary_value_w[4]->value(0,0,-0.5,t) + 9*Yz[0] - Yz[1]) / (3*DZ));
    // LOWER FRONT RIGHT VERTEX (1,0,0)
    Y2_p[(NX) * (NY + 1) * (NZ + 1)] = 120.0 / (c * DT) *
        ((8*boundary_value_u[1]->value(NX-0.5,0,0,t) - 9*Yx[(NX-1) * (NY + 1) * (NZ + 1)] + Yx[(NX-2) * (NY + 1) * (NZ + 1)]) / (3*DX) +
        (-8*boundary_value_v[2]->value(NX,-0.5,0,t) + 9*Yy[(NX) * NY * (NZ + 1)] - Yy[(NX) * NY * (NZ + 1) + 1 * (NZ + 1)]) / (3*DY) +
        (-8*boundary_value_w[4]->value(NX,0,-0.5,t) + 9*Yz[(NX) * (NY + 1) * NZ] - Yz[(NX) * (NY + 1) * NZ + 1]) / (3*DZ));
    // LOWER BACK LEFT VERTEX (0,1,0)
    Y2_p[(NY) * (NZ + 1)] = 120.0 / (c * DT) *
        ((-8*boundary_value_u[0]->value(-0.5,NY,0,t) + 9*Yx[(NY) * (NZ + 1)] - Yx[1 * (NY + 1) * (NZ + 1) + (NY) * (NZ + 1)]) / (3*DX) +
        (8*boundary_value_v[3]->value(0,NY-0.5,0,t) - 9*Yy[(NY-1) * (NZ + 1)] + Yy[(NY-2) * (NZ + 1)]) / (3*DY) +
        (-8*boundary_value_w[4]->value(0,NY,-0.5,t) + 9*Yz[(NY) * NZ] - Yz[(NY) * NZ + 1]) / (3*DZ));
    // UPPER FRONT LEFT VERTEX (0,0,1)
    Y2_p[(NZ)] = 120.0 / (c * DT) *
        ((-8*boundary_value_u[0]->value(-0.5,0,NZ,t) + 9*Yx[(NZ)] - Yx[1 * (NY + 1) * (NZ + 1) + (NZ)]) / (3*DX) +
        (-8*boundary_value_v[2]->value(0,-0.5,NZ,t) + 9*Yy[(NZ)] - Yy[1 * (NZ + 1) + (NZ)]) / (3*DY) +
        (8*boundary_value_w[5]->value(0,0,NZ-0.5,t) - 9*Yz[(NZ-1)] + Yz[(NZ-2)]) / (3*DZ));
    // UPPER BACK LEFT VERTEX (0,1,1)
    Y2_p[(NY) * (NZ + 1) + (NZ)] = 120.0 / (c * DT) *
        ((-8*boundary_value_u[0]->value(-0.5,NY,NZ,t) + 9*Yx[(NY) * (NZ + 1) + (NZ)] - Yx[1 * (NY + 1) * (NZ + 1) + (NY) * (NZ + 1) + (NZ)]) / (3*DX) +
        (8*boundary_value_v[3]->value(0,NY-0.5,NZ,t) - 9*Yy[(NY-1) * (NZ + 1) + (NZ)] + Yy[(NY-2) * (NZ + 1) + (NZ)]) / (3*DY) +
        (8*boundary_value_w[5]->value(0,NY,NZ-0.5,t) - 9*Yz[(NY) * NZ + (NZ-1)] + Yz[(NY) * NZ + (NZ-2)]) / (3*DZ));
    // UPPER FRONT RIGHT VERTEX (1,0,1)
    Y2_p[(NX) * (NY + 1) * (NZ + 1) + (NZ)] = 120.0 / (c * DT) *
        ((8*boundary_value_u[1]->value(NX-0.5,0,NZ,t) - 9*Yx[(NX-1) * (NY + 1) * (NZ + 1) + (NZ)] + Yx[(NX-2) * (NY + 1) * (NZ + 1) + (NZ)]) / (3*DX) +
        (-8*boundary_value_v[2]->value(NX,-0.5,NZ,t) + 9*Yy[(NX) * NY * (NZ + 1) + (NZ)] - Yy[(NX) * NY * (NZ + 1) + 1 * (NZ + 1) + (NZ)]) / (3*DY) +
        (8*boundary_value_w[5]->value(NX,0,NZ-0.5,t) - 9*Yz[(NX) * (NY + 1) * NZ + (NZ-1)] + Yz[(NX) * (NY + 1) * NZ + (NZ-2)]) / (3*DZ));
    // LOWER BACK RIGHT VERTEX (1,1,0)
    Y2_p[(NX) * (NY + 1) * (NZ + 1) + (NY) * (NZ + 1)] = 120.0 / (c * DT) *
        ((8*boundary_value_u[1]->value(NX-0.5,NY,0,t) - 9*Yx[(NX-1) * (NY + 1) * (NZ + 1) + (NY) * (NZ + 1)] + Yx[(NX-2) * (NY + 1) * (NZ + 1) + (NY) * (NZ + 1)]) / (3*DX) +
        (8*boundary_value_v[3]->value(NX,NY-0.5,0,t) - 9*Yy[(NX) * NY * (NZ + 1) + (NY-1) * (NZ + 1)] + Yy[(NX) * NY * (NZ + 1) + (NY-2) * (NZ + 1)]) / (3*DY) +
        (-8*boundary_value_w[4]->value(NX,NY,-0.5,t) + 9*Yz[(NX) * (NY + 1) * NZ + (NY) * NZ] - Yz[(NX) * (NY + 1) * NZ + (NY) * NZ + 1]) / (3*DZ));
    // UPPER BACK RIGHT VERTEX (1,1,1)
    Y2_p[(NX) * (NY + 1) * (NZ + 1) + (NY) * (NZ + 1) + (NZ)] = 120.0 / (c * DT) *
        ((8*boundary_value_u[1]->value(NX-0.5,NY,NZ,t) - 9*Yx[(NX-1) * (NY + 1) * (NZ + 1) + (NY) * (NZ + 1) + (NZ)] + Yx[(NX-2) * (NY + 1) * (NZ + 1) + (NY) * (NZ + 1) + (NZ)]) / (3*DX) +
        (8*boundary_value_v[3]->value(NX,NY-0.5,NZ,t) - 9*Yy[(NX) * NY * (NZ + 1) + (NY-1) * (NZ + 1) + (NZ)] + Yy[(NX) * NY * (NZ + 1) + (NY-2) * (NZ + 1) + (NZ)]) / (3*DY) +
        (8*boundary_value_w[5]->value(NX,NY,NZ-0.5,t) - 9*Yz[(NX) * (NY + 1) * NZ + (NY) * NZ + (NZ-1)] + Yz[(NX) * (NY + 1) * NZ + (NY) * NZ + (NZ-2)]) / (3*DZ));
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
