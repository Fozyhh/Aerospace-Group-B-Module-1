#include "boundary.hpp"

Boundary::Boundary(Grid* grid_,double dx_, double dy_,double dz_): 
grid(grid_)
,nx(grid->nx)
,ny(grid->ny)
,nz(grid->nz)
,dx(dx_)
,dy(dy_)
,dz(dz_)
,prec(1/2)
{}

void Boundary::update_boundary(double t){


    std::cout << "Check 0"<< std::endl;
    size_t face=0; // left face

    for (size_t k=0; k < nz; k++)
    {
        grid->w[k] = boundary_value_w[face]->value(0,0,k,t); //0 0 0.5
        grid->v[k] = boundary_value_v[face]->value(0,0,k,t); //0 0.5 0
    }
    grid->v[nz] = boundary_value_v[face]->value(0,0,nz,t);

    for (size_t j = 1; j < ny; j++)
    {
        grid->v[j*(nz+1)] = boundary_value_v[face]->value(0,j,0,t);
        grid->w[j*nz] = boundary_value_w[face]->value(0,j,0,t); 
        for(size_t k = 1; k < nz; k++)
        {
            grid->v[j*(nz+1) + k] = boundary_value_v[face]->value(0,j,k,t);
            grid->w[j*nz + k] = boundary_value_w[face]->value(0,j,k,t); 
            grid->u[j*(nz+1) + k] = approximate_boundary_u(0,j,k,t,face);
        }
        grid->v[j*(nz+1) + nz] = boundary_value_v[face]->value(0,j,nz,t);
    }

    for (size_t k = 0; k < nz; k++)
    {
        grid->w[ny*nz + k] = boundary_value_w[face]->value(0,ny,k,t);
    }
    


    std::cout << "Check 1"<< std::endl;
    face=1; // right face

    for (size_t k=0; k < nz; k++)
    {
        grid->v[nx*ny*(nz+1) + k] = boundary_value_v[face]->value(nx,0,k,t); 
        grid->w[nx*(ny+1)*nz + k] = boundary_value_w[face]->value(nx,0,k,t); 
    }
    grid->v[nx*ny*(nz+1) + nz] = boundary_value_v[face]->value(nx,0,nz,t);

    for (size_t j = 1; j < ny; j++)
    {
        grid->v[nx*ny*(nz+1) + j*(nz+1)] = boundary_value_v[face]->value(nx,j,0,t);
        grid->w[nx*(ny+1)*nz + j*nz] = boundary_value_w[face]->value(nx,j,0,t); 
        for(size_t k = 1; k < nz; k++)
        {
            grid->v[nx*ny*(nz+1) + j*(nz+1) + k] = boundary_value_v[face]->value(nx,j,k,t);
            grid->w[nx*(ny+1)*nz + j*nz + k] = boundary_value_w[face]->value(nx,j,k,t); 
            grid->u[(nx-1)*(ny+1)*(nz+1) + j*(nz+1) + k] = approximate_boundary_u(nx-1,j,k,t,face);
        }
        grid->v[nx*ny*(nz+1) + j*(nz+1) + nz] = boundary_value_v[face]->value(nx,j,nz,t);
    }

    for (size_t k = 0; k < nz; k++)
    {
        grid->w[nx*(ny+1)*nz + ny*nz + k] = boundary_value_w[face]->value(nx,ny,k,t);
    }



    std::cout << "Check 2"<< std::endl;
    face=2; // front face

    for (size_t k=1; k < nz; k++)
    {
        grid->u[k] = boundary_value_u[face]->value(0,0,k,t);
    }
    
    for (size_t i = 1; i < nx; i++)
    {
        grid->w[i * (ny+1) * nz] = boundary_value_w[face]->value(i,0,0,t);
        for(size_t k = 1; k < nz; k++)
        {
            grid->u[i * (ny+1) * (nz+1) + k] = boundary_value_u[face]->value(i,0,k,t);
            grid->w[i * (ny+1) * nz + k] = boundary_value_w[face]->value(i,0,k,t);
            grid->v[i * ny * (nz+1) + k] = approximate_boundary_v(i,0,k,t,face);
        }
    }

    

    std::cout << "Check 3"<< std::endl;
    face=3; //back face

    for (size_t k=1; k < nz; k++)
    {
        grid->u[ny * (nz+1) + k] = boundary_value_u[face]->value(0,ny,k,t);
    }
    
    for (size_t i = 1; i < nx; i++)
    {
        grid->w[(ny+1) * nz * i + ny * nz] = boundary_value_w[face]->value(i,ny,0,t);
        for(size_t k = 1; k < nz; k++)
        {
            grid->u[ny * (nz + 1) + (ny+1) * (nz+1) * i + k] = boundary_value_u[face]->value(i,ny,k,t);
            grid->w[ny * nz + i * (ny+1) * nz + k] = boundary_value_w[face]->value(i,ny,k,t);
            grid->v[(ny-1) * (nz+1) + ny * (nz+1) * i + k] = approximate_boundary_v(i,ny-1,k,t,face);
        }
    }



    face=4; //lower face

    for (size_t j=0; j < ny+1; j++)
    {
        grid->u[j*(nz+1)] = boundary_value_u[face]->value(0,j,0,t);
    }
    
    for (size_t i = 1; i < nx; i++)
    {
        grid->u[i*(ny+1)*(nz+1)] = boundary_value_u[face]->value(i,0,0,t);
        grid->v[ny * (nz+1) * i] = boundary_value_v[face]->value(i,0,0,t); 
        for(size_t j = 1; j < ny; j++)
        {   
            grid->u[(ny+1) * (nz+1) * i + j*(nz+1)] = boundary_value_u[face]->value(i,j,0,t);
            grid->v[ny * (nz+1) * i + j*(nz+1)] = boundary_value_v[face]->value(i,j,0,t); 
            grid->w[(ny+1) * nz * i + j*nz] = approximate_boundary_w(i,j,0,t,face);
        }
        grid->u[(ny+1) * (nz+1) * i + ny*(nz+1)] = boundary_value_u[face]->value(i,ny,0,t);
    }



    face=5; //upper face

    for (size_t j=0; j < ny+1; j++)
    {
        grid->u[j*(nz+1) + nz] = boundary_value_u[face]->value(0,j,nz,t);
    }
    
    for (size_t i = 1; i < nx; i++)
    {
        grid->u[i * (ny+1) * (nz+1) + nz] = boundary_value_u[face]->value(i,0,nz,t);
        grid->v[ny * (nz+1) * i + nz] = boundary_value_v[face]->value(i,0,nz,t); 
        for (size_t j = 1; j < ny; j++)
        {
            grid->u[(ny+1) * (nz+1) * i + j*(nz+1) + nz] = boundary_value_u[face]->value(i,j,nz,t);
            grid->v[ny * (nz+1) * i + j*(nz+1) + nz] = boundary_value_v[face]->value(i,j,nz,t); 
            grid->w[(ny+1) * nz * i + j*nz + (nz - 1)] = approximate_boundary_w(i,j,nz-1,t,face);
        }
        grid->u[(ny+1) * (nz+1) * i + ny*(nz+1) + nz] = boundary_value_u[face]->value(i,ny,nz,t);
    }

}   



double Boundary::approximate_boundary_u(size_t x, size_t y, size_t z, double t,size_t face) {
    double dv = ((boundary_value_v[face]->value(x, y + prec, z,t)) -
                ((boundary_value_v[face]->value(x, y - prec, z,t)))) / (2*prec);
    
    double dw =  ((boundary_value_w[face]->value(x, y, z + prec ,t)) -
                 ((boundary_value_w[face]->value(x, y, z - prec,t)))) / (2*prec);

    return  boundary_value_u[face]->value(x, y, z, t) - (dv + dw) * (dx/2);
}

double Boundary::approximate_boundary_v(size_t x, size_t y, size_t z, double t,size_t face) {
    double du = ((boundary_value_u[face]->value(x + prec, y, z,t)) -
                ((boundary_value_u[face]->value(x- prec, y, z,t)))) / (2*prec);
    
    double dw =  ((boundary_value_w[face]->value(x, y, z + prec,t)) -
                 ((boundary_value_w[face]->value(x, y, z - prec,t)))) / (2*prec);

    return  boundary_value_v[face]->value(x, y, z, t) - (du + dw) * (dy/2);
}

double Boundary::approximate_boundary_w(size_t x, size_t y, size_t z, double t,size_t face) {
    double du = ((boundary_value_u[face]->value(x + prec, y, z,t)) -
                ((boundary_value_u[face]->value(x - prec, y, z,t)))) / (2*prec);
    
    double dv = ((boundary_value_v[face]->value(x, y + prec, z,t)) -
                ((boundary_value_v[face]->value(x, y - prec, z,t)))) / (2*prec);
    return  boundary_value_w[face]->value(x, y, z, t) - (du + dv) * (dz/2);
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