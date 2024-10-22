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

    size_t face=0; // left face
    
    // TODO: check
    for (size_t k=0; k < nz; k++)
    {
        grid->w[k] = boundary_value_w[face]->value(0,0,k,t); //0 0 0.5
        grid->v[k] = boundary_value_v[face]->value(0,0,k,t); //0 0.5 0
    }
    std::cout << "Check 0"<< std::endl;
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
    //Missing: front face, back face, upper face, lower face
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
    grid->u[0] =  boundary_value_u[face]->value(0,0,0,t);
    for (size_t i=1; i < nx; i++)
    {
        grid->u[(ny+1) * (nz+1) * i] = boundary_value_u[face]->value(i,0,0,t);
        grid->v[(nz+1) * ny * i] = boundary_value_v[face]->value(i,0,0,t);
    }
    std::cout << "Check 2,3"<< std::endl;
    for (size_t j = 1; j < ny; j++)
    {
        grid->u[j * (nz+1)] = boundary_value_u[face]->value(0,j,0,t);
        for(size_t i = 1; i < nx; i++)
        {   
            std::cout << "nx " << nx << " ny " << ny << " nz " << nz <<std::endl;
            std::cout << "Check 2,5, indexes : "<< j << "-" << i << std::endl;
            grid->u[(ny+1) * (nz+1) * i + j*(nz+1)] = boundary_value_u[face]->value(i,j,0,t);
            std::cout << "Check 2,55"<< std::endl;
            grid->v[ny * (nz+1) * i + j*(nz+1)] = boundary_value_v[face]->value(i,j,0,t); 
            std::cout << "Check 2,6: " << (ny+1) * nz * i + j*nz << std::endl;
            std::cout << "Size: " << grid->w.size() << std::endl;
            grid->w[(ny+1) * nz * i + j*nz] = approximate_boundary_w(i,j,0,t,face);
        }
        grid->u[(ny+1) * (nz+1) * (nx) + j*(nz+1)] = boundary_value_u[face]->value(nx,j,0,t);
    }
    std::cout << "Check 2,8"<< std::endl;
    for (size_t i = 0; i < nx; i++)
    {
        grid->u[(ny+1) * (nz+1) * i + ny*(nz+1)] = boundary_value_u[face]->value(i,ny,0,t);
    }
    std::cout << "Check 3"<< std::endl;
    face=3; //back face
    grid->u[nz] =  boundary_value_u[face]->value(0,0,nz,t);
    for (size_t i=1; i < nx; i++)
    {
        grid->u[(ny+1) * (nz+1) * i + nz] = boundary_value_u[face]->value(i,0,nz,t);
        grid->v[(nz+1) * ny * i + nz] = boundary_value_v[face]->value(i,0,nz,t);
    }
    
    for (size_t j = 1; j < ny; j++)
    {
        grid->u[j * (nz+1) + nz] = boundary_value_u[face]->value(0,j,nz,t);// THIS IS NOT THE BACK FACE BUT THE UPPER FACE, STILL OK
        for(size_t i = 1; i < nx; i++)
        {
            grid->u[(ny+1) * (nz+1) * i + j*(nz+1) + nz] = boundary_value_u[face]->value(i,j,nz,t);
            grid->v[ny * (nz+1) * i + j*(nz+1) + nz] = boundary_value_v[face]->value(i,j,nz,t); 
            grid->w[(ny+1) * nz * i + j*nz + (nz - 1)] = approximate_boundary_w(i,j,nz-1,t,face);
        }
        grid->u[(ny+1) * (nz+1) * (nx) + j*(nz+1) + nz] = boundary_value_u[face]->value(nx,j,nz,t);
    }

    for (size_t i = 0; i < nx; i++)
    {
        grid->u[(ny+1) * (nz+1) * i + ny*(nz+1) + nz] = boundary_value_u[face]->value(i,ny,nz,t);
    }

    face=4; //lower face

    for (size_t k=1; k < nz; k++)
    {
        grid->u[k] = boundary_value_u[face]->value(0,0,k,t);
    }
    
    for (size_t i = 1; i < nx; i++)
    {
        grid->w[i*(ny+1)*nz] = boundary_value_w[face]->value(i,0,0,t);
        for(size_t k = 1; k < nz; k++)
        {
            grid->u[(ny+1) * (nz+1) * i + k] = boundary_value_u[face]->value(i,0,k,t);
            grid->v[ny * (nz+1) * i + k] = boundary_value_v[face]->value(i,0,k,t); 
            grid->w[(ny+1) * nz * i + k] = approximate_boundary_w(i,0,k,t,face);
        }
    }

    face=5;
    for (size_t k=1; k < nz; k++)
    {
        grid->u[(ny-1) * (nz+1) + k] = boundary_value_u[face]->value(0,ny,k,t);
    }
    
    for (size_t i = 1; i < nx; i++)
    {
        grid->w[(ny+1) * nz * i + ny * nz] = boundary_value_w[face]->value(i,ny,0,t);
        for(size_t k = 1; k < nz; k++)
        {
            grid->u[ny * (nz + 1) + (ny+1) * (nz+1) * i + k] = boundary_value_u[face]->value(i,ny,k,t);
            grid->v[(ny-1) * (nz+1) + ny * (nz+1) * i + k] = boundary_value_v[face]->value(i,ny,k,t); 
            grid->w[(ny-1)*nz + (ny+1) * nz * i + k] = approximate_boundary_w(i,ny-1,k,t,face);
        }
    }

}   



double Boundary::approximate_boundary_u(size_t x, size_t y, size_t z, double t,size_t face) {
    double dv = ((boundary_value_v[face]->value(x * dx, (y + prec) * dy, z,t)) -
                ((boundary_value_v[face]->value(x * dx, (y - prec) * dy, z,t)))) / (2*prec);
    
    double dw =  ((boundary_value_w[face]->value(x * dx, y * dy, (z + prec) * dz ,t)) -
                 ((boundary_value_w[face]->value(x * dx, y * dy, (z - prec) * dz,t)))) / (2*prec);

    return  boundary_value_u[face]->value(x * dx, y * dy, z * dz,t) - (dv + dw) * (dx/2);
}

double Boundary::approximate_boundary_v(size_t x, size_t y, size_t z, double t,size_t face) {
    double du = ((boundary_value_u[face]->value((x + prec) * dx, y * dy, z,t)) -
                ((boundary_value_u[face]->value((x- prec) * dx, y * dy, z,t)))) / (2*prec);
    
    double dw =  ((boundary_value_w[face]->value(x * dx, y * dy, (z + prec) * dz,t)) -
                 ((boundary_value_w[face]->value(x * dx, y * dy, (z - prec) * dz,t)))) / (2*prec);

    return  boundary_value_v[face]->value(x * dx, y * dy, z * dz,t) - (du + dw) * (dy/2);
}

double Boundary::approximate_boundary_w(size_t x, size_t y, size_t z, double t,size_t face) {
    double du = ((boundary_value_u[face]->value((x + prec) * dx, y * dy, z,t)) -
                ((boundary_value_u[face]->value((x - prec) * dx, y * dy, z,t)))) / (2*prec);
    
    double dv = ((boundary_value_v[face]->value(x * dx, (y + prec) * dy, z,t)) -
                ((boundary_value_v[face]->value(x * dx, (y - prec) * dy, z,t)))) / (2*prec);
    return  boundary_value_w[face]->value(x * dx, y * dy, z * dz,t) - (du + dv) * (dz/2);
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





// given x(i) pointing to the right, y(j) pointing into the screen, z(k) pointing up
for (k){
    for(j){
        for(i){
            [i+j*(nx+1)+k*(nx+1)*(ny+1)]
        }
    }
}