#include "boundary.hpp"

Boundary::Boundary(Grid& grid_,double dx_, double dy_,double dz_): 
grid(grid_)
,nx(grid.nx)
,ny(grid.ny)
,nz(grid.nz)
,dx(dx_)
,dy(dy_)
,dz(dz_)
,prec(1.0/2.0)
{}

// The method, which takes as input only the time step, is called by the program at the start of the time step.
// For each boundary node, it takes the exact value dor each component from the input and updates them.
// It also updates those values that are not directly on a face, but need an approximation.
void Boundary::update_boundary(double t){

    // Each face is numbered from 0 to 5 and we treat every face separately


    // LEFT FACE

    size_t face=0;

    for (size_t k=0; k < nz; k++)
    {
        grid.w[k] = boundary_value_w[face]->value(0,0,k,t);
        grid.v[k] = boundary_value_v[face]->value(0,0,k,t);
    }
    grid.v[nz] = boundary_value_v[face]->value(0,0,nz,t);
    for (size_t j = 1; j < ny; j++)
    {
        grid.v[j*(nz+1)] = boundary_value_v[face]->value(0,j,0,t);
        grid.w[j*nz] = boundary_value_w[face]->value(0,j,0,t); 
        for(size_t k = 1; k < nz; k++)
        {
            grid.v[j*(nz+1) + k] = boundary_value_v[face]->value(0,j,k,t);
            grid.w[j*nz + k] = boundary_value_w[face]->value(0,j,k,t); 
            grid.u[j*(nz+1) + k] = approximate_boundary_u(0,j,k,t,face);
        }
        grid.v[j*(nz+1) + nz] = boundary_value_v[face]->value(0,j,nz,t);
    }
    
    for (size_t k = 0; k < nz; k++)
    {
        grid.w[ny*nz + k] = boundary_value_w[face]->value(0,ny,k,t);
    }
    


    // RIGHT FACE

    face=1;

    for (size_t k=0; k < nz; k++)
    {
        grid.v[nx*ny*(nz+1) + k] = boundary_value_v[face]->value(nx,0,k,t); 
        grid.w[nx*(ny+1)*nz + k] = boundary_value_w[face]->value(nx,0,k,t); 
    }
    grid.v[nx*ny*(nz+1) + nz] = boundary_value_v[face]->value(nx,0,nz,t);

    for (size_t j = 1; j < ny; j++)
    {
        grid.v[nx*ny*(nz+1) + j*(nz+1)] = boundary_value_v[face]->value(nx,j,0,t);
        grid.w[nx*(ny+1)*nz + j*nz] = boundary_value_w[face]->value(nx,j,0,t); 
        for(size_t k = 1; k < nz; k++)
        {
            grid.v[nx*ny*(nz+1) + j*(nz+1) + k] = boundary_value_v[face]->value(nx,j,k,t);
            grid.w[nx*(ny+1)*nz + j*nz + k] = boundary_value_w[face]->value(nx,j,k,t); 
            grid.u[(nx-1)*(ny+1)*(nz+1) + j*(nz+1) + k] = approximate_boundary_u(nx-1,j,k,t,face);
        }
        grid.v[nx*ny*(nz+1) + j*(nz+1) + nz] = boundary_value_v[face]->value(nx,j,nz,t);
    }

    for (size_t k = 0; k < nz; k++)
    {
        grid.w[nx*(ny+1)*nz + ny*nz + k] = boundary_value_w[face]->value(nx,ny,k,t);
    }



    // FRONT FACE

    face=2;

    for (size_t k=1; k < nz; k++)
    {
        grid.u[k] = boundary_value_u[face]->value(0,0,k,t);
    }
    
    for (size_t i = 1; i < nx; i++)
    {
        grid.w[i * (ny+1) * nz] = boundary_value_w[face]->value(i,0,0,t);
        for(size_t k = 1; k < nz; k++)
        {
            grid.u[i * (ny+1) * (nz+1) + k] = boundary_value_u[face]->value(i,0,k,t);
            grid.w[i * (ny+1) * nz + k] = boundary_value_w[face]->value(i,0,k,t);
            grid.v[i * ny * (nz+1) + k] = approximate_boundary_v(i,0,k,t,face);
        }
    }

    

    // BACK FACE

    face=3;

    for (size_t k=1; k < nz; k++)
    {
        grid.u[ny * (nz+1) + k] = boundary_value_u[face]->value(0,ny,k,t);
    }
    
    for (size_t i = 1; i < nx; i++)
    {
        grid.w[(ny+1) * nz * i + ny * nz] = boundary_value_w[face]->value(i,ny,0,t);
        for(size_t k = 1; k < nz; k++)
        {
            grid.u[ny * (nz + 1) + (ny+1) * (nz+1) * i + k] = boundary_value_u[face]->value(i,ny,k,t);
            grid.w[ny * nz + i * (ny+1) * nz + k] = boundary_value_w[face]->value(i,ny,k,t);
            grid.v[(ny-1) * (nz+1) + ny * (nz+1) * i + k] = approximate_boundary_v(i,ny-1,k,t,face);
        }
    }



    // LOWER FACE

    face=4;

    for (size_t j=0; j < ny+1; j++)
    {
        grid.u[j*(nz+1)] = boundary_value_u[face]->value(0,j,0,t);
    }
    
    for (size_t i = 1; i < nx; i++)
    {
        grid.u[i*(ny+1)*(nz+1)] = boundary_value_u[face]->value(i,0,0,t);
        grid.v[ny * (nz+1) * i] = boundary_value_v[face]->value(i,0,0,t); 
        for(size_t j = 1; j < ny; j++)
        {   
            grid.u[(ny+1) * (nz+1) * i + j*(nz+1)] = boundary_value_u[face]->value(i,j,0,t);
            grid.v[ny * (nz+1) * i + j*(nz+1)] = boundary_value_v[face]->value(i,j,0,t); 
            grid.w[(ny+1) * nz * i + j*nz] = approximate_boundary_w(i,j,0,t,face);
        }
        grid.u[(ny+1) * (nz+1) * i + ny*(nz+1)] = boundary_value_u[face]->value(i,ny,0,t);
    }



    // UPPER FACE

    face=5;

    for (size_t j=0; j < ny+1; j++)
    {
        grid.u[j*(nz+1) + nz] = boundary_value_u[face]->value(0,j,nz,t);
    }
    
    for (size_t i = 1; i < nx; i++)
    {
        grid.u[i * (ny+1) * (nz+1) + nz] = boundary_value_u[face]->value(i,0,nz,t);
        grid.v[ny * (nz+1) * i + nz] = boundary_value_v[face]->value(i,0,nz,t); 
        for (size_t j = 1; j < ny; j++)
        {
            grid.u[(ny+1) * (nz+1) * i + j*(nz+1) + nz] = boundary_value_u[face]->value(i,j,nz,t);
            grid.v[ny * (nz+1) * i + j*(nz+1) + nz] = boundary_value_v[face]->value(i,j,nz,t); 
            grid.w[(ny+1) * nz * i + j*nz + (nz - 1)] = approximate_boundary_w(i,j,nz-1,t,face);
        }
        grid.u[(ny+1) * (nz+1) * i + ny*(nz+1) + nz] = boundary_value_u[face]->value(i,ny,nz,t);
    }


}   


// Performs the approximation of the component u that isn't precisely on the boundary.
double Boundary::approximate_boundary_u(size_t x, size_t y, size_t z, double t,size_t face) {
    double dv = ((boundary_value_v[face]->value(x, y , z,t)) -
                ((boundary_value_v[face]->value(x, y - 1, z,t)))) / (2*prec);

    double dw =  ((boundary_value_w[face]->value(x, y, z ,t)) -
                 ((boundary_value_w[face]->value(x, y, z - 1,t)))) / (2*prec);
    
    return  boundary_value_u[face]->value(x-0.5, y, z, t) - (dv + dw) * (dx/2);
}


// Performs the approximation of the component v that isn't precisely on the boundary.
double Boundary::approximate_boundary_v(size_t x, size_t y, size_t z, double t,size_t face) {
    double du = ((boundary_value_u[face]->value(x , y, z,t)) -
                ((boundary_value_u[face]->value(x - 1, y, z,t)))) / (2*prec);
    
    double dw =  ((boundary_value_w[face]->value(x, y, z,t)) -
                 ((boundary_value_w[face]->value(x, y, z - 1,t)))) / (2*prec);

    return  boundary_value_v[face]->value(x, y-0.5, z, t) - (du + dw) * (dy/2);
}


// Performs the approximation of the component w that isn't precisely on the boundary.
double Boundary::approximate_boundary_w(size_t x, size_t y, size_t z, double t,size_t face) {
    double du = ((boundary_value_u[face]->value(x, y, z,t)) -
                ((boundary_value_u[face]->value(x - 1, y, z,t)))) / (2*prec);
    
    double dv = ((boundary_value_v[face]->value(x, y, z,t)) -
                ((boundary_value_v[face]->value(x, y - 1, z,t)))) / (2*prec);
    return  boundary_value_w[face]->value(x, y, z-0.5, t) - (du + dv) * (dz/2);
}

void Boundary::addFunction(size_t direction, std::shared_ptr<BoundaryFunction> x){
      switch(direction){
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