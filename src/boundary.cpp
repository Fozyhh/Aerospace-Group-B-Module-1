#include "boundary.hpp"

Boundary::Boundary(int nx, int ny, int nz,double dx_, double dy_,double dz_): 
nx(nx)
,ny(ny)
,nz(nz)
,dx(dx_)
,dy(dy_)
,dz(dz_)
{}

// The method, which takes as input only the time step, is called by the program at the start of the time step.
// For each boundary node, it takes the exact value dor each component from the input and updates them.
// It also updates those values that are not directly on a face, but need an approximation.
void Boundary::update_boundary(std::vector<double> Yx,std::vector<double> Yy,std::vector<double> Yz, double t){

    // Each face is numbered from 0 to 5 and we treat every face separately
    // LEFT FACE

    size_t face=0;

    for (size_t k=0; k < nz; k++)
    {
        Yz[k] = boundary_value_w[face]->value(0,0,k,t);
        Yy[k] = boundary_value_v[face]->value(0,0,k,t);
    }
    Yy[nz] = boundary_value_v[face]->value(0,0,nz,t);
    for (size_t j = 1; j < ny; j++)
    {
        Yy[j*(nz+1)] = boundary_value_v[face]->value(0,j,0,t);
        Yz[j*nz] = boundary_value_w[face]->value(0,j,0,t); 
        for(size_t k = 1; k < nz; k++)
        {
            Yy[j*(nz+1) + k] = boundary_value_v[face]->value(0,j,k,t);
            Yz[j*nz + k] = boundary_value_w[face]->value(0,j,k,t); 
            Yx[j*(nz+1) + k] = approximate_boundary_u(0,j,k,t,face,1);
        }
        Yy[j*(nz+1) + nz] = boundary_value_v[face]->value(0,j,nz,t);
    }
    
    for (size_t k = 0; k < nz; k++)
    {
        Yz[ny*nz + k] = boundary_value_w[face]->value(0,ny,k,t);
    }
    


    // RIGHT FACE

    face=1;

    for (size_t k=0; k < nz; k++)
    {
        Yy[nx*ny*(nz+1) + k] = boundary_value_v[face]->value(nx,0,k,t); 
        Yz[nx*(ny+1)*nz + k] = boundary_value_w[face]->value(nx,0,k,t); 
    }
    Yy[nx*ny*(nz+1) + nz] = boundary_value_v[face]->value(nx,0,nz,t);

    for (size_t j = 1; j < ny; j++)
    {
        Yy[nx*ny*(nz+1) + j*(nz+1)] = boundary_value_v[face]->value(nx,j,0,t);
        Yz[nx*(ny+1)*nz + j*nz] = boundary_value_w[face]->value(nx,j,0,t); 
        for(size_t k = 1; k < nz; k++)
        {
            Yy[nx*ny*(nz+1) + j*(nz+1) + k] = boundary_value_v[face]->value(nx,j,k,t);
            Yz[nx*(ny+1)*nz + j*nz + k] = boundary_value_w[face]->value(nx,j,k,t); 
            Yx[(nx-1)*(ny+1)*(nz+1) + j*(nz+1) + k] = approximate_boundary_u(nx,j,k,t,face,-1);
        }
        Yy[nx*ny*(nz+1) + j*(nz+1) + nz] = boundary_value_v[face]->value(nx,j,nz,t);
    }

    for (size_t k = 0; k < nz; k++)
    {
        Yz[nx*(ny+1)*nz + ny*nz + k] = boundary_value_w[face]->value(nx,ny,k,t);
    }



    // FRONT FACE

    face=2;

    for (size_t k=1; k < nz; k++)
    {
        Yx[k] = boundary_value_u[face]->value(0,0,k,t);
    }
    
    for (size_t i = 1; i < nx; i++)
    {
        Yz[i * (ny+1) * nz] = boundary_value_w[face]->value(i,0,0,t);
        for(size_t k = 1; k < nz; k++)
        {
            Yx[i * (ny+1) * (nz+1) + k] = boundary_value_u[face]->value(i,0,k,t);
            Yz[i * (ny+1) * nz + k] = boundary_value_w[face]->value(i,0,k,t);
            Yy[i * ny * (nz+1) + k] = approximate_boundary_v(i,0,k,t,face,1);
        }
    }

    

    // BACK FACE

    face=3;

    for (size_t k=1; k < nz; k++)
    {
        Yx[ny * (nz+1) + k] = boundary_value_u[face]->value(0,ny,k,t);
    }
    
    for (size_t i = 1; i < nx; i++)
    {
        Yz[(ny+1) * nz * i + ny * nz] = boundary_value_w[face]->value(i,ny,0,t);
        for(size_t k = 1; k < nz; k++)
        {
            Yx[ny * (nz + 1) + (ny+1) * (nz+1) * i + k] = boundary_value_u[face]->value(i,ny,k,t);
            Yz[ny * nz + i * (ny+1) * nz + k] = boundary_value_w[face]->value(i,ny,k,t);
            Yy[(ny-1) * (nz+1) + ny * (nz+1) * i + k] = approximate_boundary_v(i,ny,k,t,face,-1);
        }
    }



    // LOWER FACE

    face=4;

    for (size_t j=0; j < ny+1; j++)
    {
        Yx[j*(nz+1)] = boundary_value_u[face]->value(0,j,0,t);
    }
    
    for (size_t i = 1; i < nx; i++)
    {
        Yx[i*(ny+1)*(nz+1)] = boundary_value_u[face]->value(i,0,0,t);
        Yy[ny * (nz+1) * i] = boundary_value_v[face]->value(i,0,0,t); 
        for(size_t j = 1; j < ny; j++)
        {   
            Yx[(ny+1) * (nz+1) * i + j*(nz+1)] = boundary_value_u[face]->value(i,j,0,t);
            Yy[ny * (nz+1) * i + j*(nz+1)] = boundary_value_v[face]->value(i,j,0,t); 
            Yz[(ny+1) * nz * i + j*nz] = approximate_boundary_w(i,j,0,t,face,1);
        }
        Yx[(ny+1) * (nz+1) * i + ny*(nz+1)] = boundary_value_u[face]->value(i,ny,0,t);
    }



    // UPPER FACE

    face=5;

    for (size_t j=0; j < ny+1; j++)
    {
        Yx[j*(nz+1) + nz] = boundary_value_u[face]->value(0,j,nz,t);
    }
    
    for (size_t i = 1; i < nx; i++)
    {
        Yx[i * (ny+1) * (nz+1) + nz] = boundary_value_u[face]->value(i,0,nz,t);
        Yy[ny * (nz+1) * i + nz] = boundary_value_v[face]->value(i,0,nz,t); 
        for (size_t j = 1; j < ny; j++)
        {
            Yx[(ny+1) * (nz+1) * i + j*(nz+1) + nz] = boundary_value_u[face]->value(i,j,nz,t);
            Yy[ny * (nz+1) * i + j*(nz+1) + nz] = boundary_value_v[face]->value(i,j,nz,t); 
            Yz[(ny+1) * nz * i + j*nz + (nz - 1)] = approximate_boundary_w(i,j,nz,t,face,-1);
        }
        Yx[(ny+1) * (nz+1) * i + ny*(nz+1) + nz] = boundary_value_u[face]->value(i,ny,nz,t);
    }
}   


// Performs the approximation of the component u that isn't precisely on the boundary.
double Boundary::approximate_boundary_u(size_t x, size_t y, size_t z, double t,size_t face,int side) {
    
    double dv = (boundary_value_v[face]->value(x, y , z, t) - boundary_value_v[face]->value(x, y - (dy), z, t) ) / (2*(dy/2.0));

    double dw = (boundary_value_w[face]->value(x, y , z ,t) - boundary_value_w[face]->value(x, y , z - (dz) , t)) / (2*(dz/2.0));

    return  boundary_value_u[face]->value((x-(dx/2.0)), y, z, t) - (dv + dw) * (dx/2)*side;
}


// Performs the approximation of the component v that isn't precisely on the boundary.
double Boundary::approximate_boundary_v(size_t x, size_t y, size_t z, double t,size_t face, int side) {
    double du = ((boundary_value_u[face]->value(x , y, z,t)) -
                ((boundary_value_u[face]->value(x - (dx), y, z,t)))) / (dx);
    
    double dw =  ((boundary_value_w[face]->value(x, y, z,t)) -
                 ((boundary_value_w[face]->value(x, y, z - (dz),t)))) / (dz);

    return  boundary_value_v[face]->value(x, y-(dy/2.0), z, t) - (du + dw) * (dy/2.0)*side;
}


// Performs the approximation of the component w that isn't precisely on the boundary.
double Boundary::approximate_boundary_w(size_t x, size_t y, size_t z, double t,size_t face,int side) {
    double du = ((boundary_value_u[face]->value(x, y, z,t)) -
                ((boundary_value_u[face]->value(x - (dx), y, z,t)))) / (dx);
    
    double dv = ((boundary_value_v[face]->value(x, y, z,t)) -
                ((boundary_value_v[face]->value(x, y - (dy), z,t)))) / (dy);
    return  boundary_value_w[face]->value(x, y, z-(dz/2), t) - (du + dv) * (dz/2)*side;
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
