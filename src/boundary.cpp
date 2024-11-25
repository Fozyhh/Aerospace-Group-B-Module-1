#include "boundary.hpp"

// The method, which takes as input only the time step, is called by the program at the start of the time step.
// For each boundary node, it takes the exact value dor each component from the input and updates them.
// It also updates those values that are not directly on a face, but need an approximation.
void Boundary::update_boundary(std::vector<Real> &Yx, std::vector<Real> &Yy, std::vector<Real> &Yz, Real t)
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


/**
 * @brief The method is called by the program multiple during the time step, in order to update the values of the boundaries at each 
 * requested time t, calculating the approximated ones too.
 * 
 * @param Yx Boundary x velocities or the Y intermediate function related to the x direction.
 * @param Yy Boundary y velocities or the Y intermediate function related to the y direction.
 * @param Yz Boundary z velocities or the Y intermediate function related to the z direction.
 * @param t Time of the time discretization we are considering.
*/
void Boundary::update_boundary(std::array<Real, NX *(NY + 1) * (NZ + 1)> &Yx, std::array<Real, (NX + 1) * NY *(NZ + 1)> &Yy, std::array<Real, (NX + 1) * (NY + 1) * NZ> &Yz, Real t)
{
/*
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
*/

    size_t face, i, j, k;
    // X
    // LEFT FACE
    face = 2;
    for (size_t k=0; k < NZ+1; k++)
    {
        Yx[k] = boundary_value_u[face]->value(0,0,k,t);
    }

//<<<<<<< v2/array1D
    for (size_t j=1; j < NY; j++)
    {
        face = 4;
        Yx[j*(NZ+1)] = boundary_value_u[face]->value(0,j,0,t);
        face = 0;
        for(size_t k = 1; k < NZ; k++)
        {
            Yx[j*(NZ+1) + k] = approximate_boundary_u(0,j,k,t,face,1);
        }
        face = 5;
        Yx[j*(NZ+1) + NZ] = boundary_value_u[face]->value(0,j,NZ,t);
    }

    face = 3;
    for (size_t k=0; k < NZ+1; k++)
    {
        Yx[NY * (NZ+1) + k] = boundary_value_u[face]->value(0,NY,k,t);
    }

    // MIDDLE POINTS
    for (size_t i=1; i < NX-1; i++)
    {
        face = 2;
        for (size_t k=0; k < NZ+1; k++)
        {
            Yx[i * (NY+1) * (NZ+1) + k] = boundary_value_u[face]->value(i,0,k,t);
        }
        
        for (size_t j=1; j < NY; j++)
        {
            face = 4;
            Yx[(NY+1) * (NZ+1) * i + j*(NZ+1)] = boundary_value_u[face]->value(i,j,0,t);

            face = 5;
            Yx[(NY+1) * (NZ+1) * i + j*(NZ+1) + NZ] = boundary_value_u[face]->value(i,j,NZ,t);
        }

        face = 3;
        for (size_t k=0; k < NZ+1; k++)
        {
            Yx[NY * (NZ + 1) + (NY+1) * (NZ+1) * i + k] = boundary_value_u[face]->value(i,NY,k,t);
        }
    }

    // RIGHT FACE
    face = 2;
    for (size_t k=0; k < NZ+1; k++)
    {
        Yx[(NX-1) * (NY+1) * (NZ+1) + k] = boundary_value_u[face]->value(NX-1,0,k,t);
    }
//=======
// Performs the approximation of the component u that isn't precisely on the boundary.
Real Boundary::approximate_boundary_u(size_t x, size_t y, size_t z, Real t,size_t face,int side) {
    // return side == 1 ? boundary_value_u[face]->value(x,y,z,t) : boundary_value_u[face]->value(x- 1.0,y,z,t);

    Real dv = (boundary_value_v[face]->value(x, y , z, t) - boundary_value_v[face]->value(x, y - 1.0, z, t) ) / (dy*1);
    Real dw = (boundary_value_w[face]->value(x, y , z ,t) - boundary_value_w[face]->value(x, y , z - 1.0 , t)) / dz;
    

    // Real dv = (boundary_value_v[face]->value(x, y-(0.5-precision) , z, t) - boundary_value_v[face]->value(x, y - (0.5+precision), z, t) ) / (dy*2*precision);
    // Real dw = (boundary_value_w[face]->value(x, y , z-(0.5-precision) ,t) - boundary_value_w[face]->value(x, y , z - (0.5+precision) , t)) / (dz*2*precision);

    return  boundary_value_u[face]->value((x-0.5), y, z, t) - (dv + dw) * (dx/2)*side;
    
    
}
//>>>>>>> main

    for (size_t j=1; j < NY; j++)
    {
        face = 4;
        Yx[(NX-1) * (NY+1) * (NZ+1) + j*(NZ+1)] = boundary_value_u[face]->value(NX-1,j,0,t);

//<<<<<<< v2/array1D
        face = 1;
        for(size_t k = 1; k < NZ; k++)
        {
            Yx[(NX-1)*(NY+1)*(NZ+1) + j*(NZ+1) + k] = approximate_boundary_u(NX,j,k,t,face,-1);
        }

        face = 5;
        Yx[(NX-1) * (NY+1) * (NZ+1) + j*(NZ+1) + NZ] = boundary_value_u[face]->value(NX-1,j,NZ,t);
    }

    face = 3;
    for (size_t k=0; k < NZ+1; k++)
    {
        Yx[(NX-1) * (NY+1) * (NZ+1) + NY*(NZ+1) + k] = boundary_value_u[face]->value(NX-1,NY,k,t);
    }


    // Y
    // LEFT FACE
    face = 0;
    for (size_t j = 0; j < NY; j++)
    {
        for(size_t k = 0; k < NZ+1; k++)
        {
            Yy[j*(NZ+1) + k] = boundary_value_v[face]->value(0,j,k,t);
        }
    }

    // MIDDLE POINTS
    for (size_t i=1; i < NX; i++)
    {
        face = 4;
        Yy[NY * (NZ+1) * i] = boundary_value_v[face]->value(i,0,0,t); 
        face = 2;
        for (size_t k=1; k < NZ; k++)
        {
            Yy[i * NY * (NZ+1) + k] = approximate_boundary_v(i,0,k,t,face,1);
        }
        face = 5;
        Yy[NY * (NZ+1) * i + NZ] = boundary_value_v[face]->value(i,0,NZ,t); 

        for (size_t j=1; j < NY-1; j++)
        {
            face = 4;
            Yy[NY * (NZ+1) * i + j*(NZ+1)] = boundary_value_v[face]->value(i,j,0,t); 

            face = 5;
            Yy[NY * (NZ+1) * i + j*(NZ+1) + NZ] = boundary_value_v[face]->value(i,j,NZ,t); 
        }
        
        face = 4;
        Yy[NY * (NZ+1) * i + (NY-1)*(NZ+1)] = boundary_value_v[face]->value(i,NY-1,0,t); 
        face = 3;
        for(size_t k = 1; k < NZ; k++)
        {
            Yy[(NY-1) * (NZ+1) + NY * (NZ+1) * i + k] = approximate_boundary_v(i,NY,k,t,face,-1);
        }
        face = 5;
        Yy[NY * (NZ+1) * i + (NY-1)*(NZ+1) + NZ] = boundary_value_v[face]->value(i,NY-1,NZ,t); 
    }

    // RIGHT FACE
    face = 1;
    for (size_t j = 0; j < NY; j++)
    {
        for(size_t k = 0; k < NZ+1; k++)
        {
            Yy[NX*NY*(NZ+1) + j*(NZ+1) + k] = boundary_value_v[face]->value(NX,j,k,t);
        }
    }


    // Z
    // LEFT FACE 
    face = 0;
    for (size_t j = 0; j < NY+1; j++)
    {
        for(size_t k = 0; k < NZ; k++)
        {
            Yz[j*NZ + k] = boundary_value_w[face]->value(0,j,k,t); 
        }
    }

    // MIDDLE POINTS
    for (size_t i=1; i < NX; i++)
    {
        face = 2;
        for (size_t k=0; k < NZ; k++)
        {
            Yz[i * (NY+1) * NZ + k] = boundary_value_w[face]->value(i,0,k,t);
        }

        for(size_t j = 1; j < NY; j++)
        {
            face = 4;
            Yz[(NY+1) * NZ * i + j*NZ] = approximate_boundary_w(i,j,0,t,face,1);
            
            face = 5;
            Yz[(NY+1) * NZ * i + j*NZ + (NZ - 1)] = approximate_boundary_w(i,j,NZ,t,face,-1);
        }

        face = 3;
        for(size_t k = 0; k < NZ+1; k++)
        {
            Yz[NY * NZ + i * (NY+1) * NZ + k] = boundary_value_w[face]->value(i,NY,k,t);
        }
    }

    // RIGHT FACE
    face = 1;
    for (size_t j = 0; j < NY+1; j++)
    {
        for(size_t k = 0; k < NZ; k++)
        {
            Yz[NX*(NY+1)*NZ + j*NZ + k] = boundary_value_w[face]->value(NX,j,k,t); 
        }
    }
//=======
// Performs the approximation of the component v that isn't precisely on the boundary.
Real Boundary::approximate_boundary_v(size_t x, size_t y, size_t z, Real t,size_t face, int side) {
    // return side == 1 ? boundary_value_v[face]->value(x,y,z,t) : boundary_value_v[face]->value(x,y-1.0,z,t);
    Real du = ((boundary_value_u[face]->value(x , y, z,t)) -
                ((boundary_value_u[face]->value(x - 1.0, y, z,t)))) / (dx);
    
    Real dw =  ((boundary_value_w[face]->value(x, y, z,t)) -
                 ((boundary_value_w[face]->value(x, y, z - 1.0,t)))) / (dz);

    // Real du = ((boundary_value_u[face]->value(x - (0.5 - precision), y , z,t)) -
    //             ((boundary_value_u[face]->value(x - (0.5 + precision), y, z,t)))) / (dx*2*precision);
    
    // Real dw =  ((boundary_value_w[face]->value(x, y, z - (0.5 - precision),t)) -
    //              ((boundary_value_w[face]->value(x, y, z - (0.5 + precision),t)))) / (dz*2*precision);

    return  boundary_value_v[face]->value(x, y-0.5, z, t) - (du + dw) * (dy/2.0)*side;
    
//>>>>>>> main
}

/**
 * @brief Calculate the approximate value of the x velocity in a given point.
 * 
 * @param x,y,z Coordinates of the position in the 3D mesh.
 * @param t Time of the time discretization we are considering.
 * @param face Face of the mesh we are considering. In order: 0 (LEFT), 1 (RIGHT), 2 (FRONT), 3 (BACK), 4 (LOWER), 5 (UPPER)
 * @param side Int used to calculate the correct approximate value. It's 1 for faces 0, 2 and 4, -1 for the others.
 * 
 * @return the approximate value.
*/
Real Boundary::approximate_boundary_u(size_t x, size_t y, size_t z, Real t, size_t face, int side)
{

    Real dv = (boundary_value_v[face]->value(x, y, z, t) - boundary_value_v[face]->value(x, y - 1.0, z, t)) / DY;
    Real dw = (boundary_value_w[face]->value(x, y, z, t) - boundary_value_w[face]->value(x, y, z - 1.0, t)) / DZ;

    return boundary_value_u[face]->value((x - 0.5 /*(DX/2.0)*/), y, z, t) - (dv + dw) * (DX / 2) * side;
}

/**
 * @brief Calculate the approximate value of the y velocity in a given point.
 * 
 * @param x,y,z Coordinates of the position in the 3D mesh.
 * @param t Time of the time discretization we are considering.
 * @param face Face of the mesh we are considering. In order: 0 (LEFT), 1 (RIGHT), 2 (FRONT), 3 (BACK), 4 (LOWER), 5 (UPPER)
 * @param side Int used to calculate the correct approximate value. It's 1 for faces 0, 2 and 4, -1 for the others.
 * 
 * @return the approximate value.
*/
Real Boundary::approximate_boundary_v(size_t x, size_t y, size_t z, Real t, size_t face, int side)
{
    Real du = ((boundary_value_u[face]->value(x, y, z, t)) -
                 ((boundary_value_u[face]->value(x - 1, y, z, t)))) /
                (DX);

    Real dw = ((boundary_value_w[face]->value(x, y, z, t)) -
                 ((boundary_value_w[face]->value(x, y, z - 1, t)))) /
                (DZ);

    return boundary_value_v[face]->value(x, y - 0.5, z, t) - (du + dw) * (DY / 2.0) * side;
}

//<<<<<<< v2/array1D
/**
 * @brief Calculate the approximate value of the z velocity in a given point.
 * 
 * @param x,y,z Coordinates of the position in the 3D mesh.
 * @param t Time of the time discretization we are considering.
 * @param face Face of the mesh we are considering. In order: 0 (LEFT), 1 (RIGHT), 2 (FRONT), 3 (BACK), 4 (LOWER), 5 (UPPER)
 * @param side Int used to calculate the correct approximate value. It's 1 for faces 0, 2 and 4, -1 for the others.
 * 
 * @return the approximate value.
*/
Real Boundary::approximate_boundary_w(size_t x, size_t y, size_t z, Real t, size_t face, int side)
{
    Real du = ((boundary_value_u[face]->value(x, y, z, t)) -
                 ((boundary_value_u[face]->value(x - 1.0, y, z, t)))) /
                (DX);

    Real dv = ((boundary_value_v[face]->value(x, y, z, t)) -
                 ((boundary_value_v[face]->value(x, y - 1.0, z, t)))) /
                (DY);
    return boundary_value_w[face]->value(x, y, z - 0.5, t) - (du + dv) * (DZ / 2) * side;
//=======
// Performs the approximation of the component w that isn't precisely on the boundary.
Real Boundary::approximate_boundary_w(size_t x, size_t y, size_t z, Real t,size_t face,int side) {
    // return side == 1 ? boundary_value_w[face]->value(x,y,z,t) : boundary_value_w[face]->value(x,y,z - 1.0,t);
    Real du = ((boundary_value_u[face]->value(x, y, z,t)) -
                (boundary_value_u[face]->value(x - 1.0, y, z,t))) / dx;
    
    Real dv = ((boundary_value_v[face]->value(x, y , z,t)) -
                (boundary_value_v[face]->value(x, y - 1.0, z,t))) / dy;
    // Real du = ((boundary_value_u[face]->value(x - (0.5- precision), y, z,t)) -
    //             ((boundary_value_u[face]->value(x - (0.5 + precision), y, z,t)))) / (dx*2*precision);
    
    // Real dv = ((boundary_value_v[face]->value(x, y - (0.5 - precision), z,t)) -
    //             ((boundary_value_v[face]->value(x, y - (0.5 + precision), z,t)))) / (dy*2*precision);
    return  boundary_value_w[face]->value(x, y, z-0.5, t) - (du + dv) * (dz/2.0)*side;
    
//>>>>>>> main
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
