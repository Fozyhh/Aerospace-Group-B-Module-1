#include "boundary.hpp"

boundary::boundary(Grid& grid_,double dx_, double dy_,double dz_): 
grid(grid_)
,nx(grid.nx)
,ny(grid.ny)
,nz(grid.nz)
,dx(dx_)
,dy(dy_)
,dz(dz_)
,prec(1/2)
{}
/*
void boundary::update_boundary(double t){

    size_t face=0;*/ // left face
    /*for (size_t k = 0; k < nz-1; k++)
    {
        grid.u[k] = boundary_value_u[face].value(0,0,k,t);
        grid.w[k] = boundary_value_w[face].value(0,0,k,t); 
    }
    grid.u[nz-1] = boundary_value_u[face].value(0,0,nz-1,t);
    for (size_t j = 1; j < ny-1; j++)
    {
        grid.u[j*nz] = boundary_value_u[face].value(0,j,0,t);
        grid.v[j*nz] = boundary_value_v[face].value(0,j,0,t);
        grid.w[j*nz] = boundary_value_w[face].value(0,j,0,t); 
        for(size_t k = 1; k < nz - 1; k++){
            grid.u[j*nz + k] = approximate_boundary_u(0,j,k,t,face);
            grid.v[j*nz + k] = boundary_value_v[face].value(0,j,k,t);
            grid.w[j*nz + k] = boundary_value_w[face].value(0,j,k,t); 
        }
        grid.u[j*nz + nz-1] = boundary_value_u[face].value(0,j,nz-1,t);
        grid.v[j*nz + nz-1] = boundary_value_v[face].value(0,j,nz-1,t);
    }
    for (size_t k = 0; k < nz-1; k++)
    {
        grid.u[(ny-1)*nz + k] = boundary_value_u[face].value(0,(ny-1),k,t);
        grid.w[(ny-1)*nz + k] = boundary_value_w[face].value(0,(ny-1),k,t); 
    }
        grid.u[(nz-1*nz) + nz-1] = boundary_value_u[face].value(0,(nz-1),nz-1,t);*/

    // TODO: check
/*
    for (size_t k=0; k < nz; k++)
    {
        grid.v[k] = boundary_value_v[face].value(0,0,k,t); 
        grid.w[k] = boundary_value_w[face].value(0,0,k,t);
    }
    grid.v[nz] = boundary_value_v[face].value(0,0,nz,t);
    for (size_t j = 1; j < ny; j++)
    {
        grid.v[j*(nz+1)] = boundary_value_v[face].value(0,j,0,t);
        grid.w[j*nz] = boundary_value_w[face].value(0,j,0,t); 
        for(size_t k = 1; k < nz; k++)
        {
            grid.v[j*(nz+1) + k] = boundary_value_v[face].value(0,j,k,t);
            grid.w[j*nz + k] = boundary_value_w[face].value(0,j,k,t); 
            grid.u[j*(nz+1) + k] = approximate_boundary_u(0,j,k,t,face);
        }
        grid.v[j*(nz+1) + nz] = boundary_value_v[face].value(0,j,nz,t);
    }

    for (size_t k = 0; k < nz; k++)
    {
        grid.w[(ny+1)*nz + k] = boundary_value_w[face].value(0,ny,k,t);
    }
    
    
    //Missing: front face, back face, upper face, lower face
    face=1; // right face
    for (size_t k=0; k < nz; k++)
    {
        grid.v[nx*ny*(nz+1) + k] = boundary_value_v[face].value(nx,0,k,t); 
        grid.w[nx*(ny+1)*nz + k] = boundary_value_w[face].value(nx,0,k,t); 
    }
    grid.v[nx*ny*(nz+1) + nz] = boundary_value_v[face].value(nx,0,nz,t);
    for (size_t j = 1; j < ny; j++)
    {
        grid.v[nx*ny*(nz+1) + j*(nz+1)] = boundary_value_v[face].value(nx,j,0,t);
        grid.w[nx*(ny+1)*nz + j*nz] = boundary_value_w[face].value(nx,j,0,t); 
        for(size_t k = 1; k < nz; k++)
        {
            grid.v[nx*ny*(nz+1) + j*(nz+1) + k] = boundary_value_v[face].value(nx,j,k,t);
            grid.w[nx*(ny+1)*nz + j*nz + k] = boundary_value_w[face].value(nx,j,k,t); 
            grid.u[(nx-1)*(ny+1)*(nz+1) + j*(nz+1) + k] = approximate_boundary_u(nx-1,j,k,t,face);
        }
        grid.v[nx*ny*(nz+1) + j*(nz+1) + nz] = boundary_value_v[face].value(nx,j,nz,t);
    }

    for (size_t k = 0; k < nz; k++)
    {
        grid.w[nx*(ny+1)*nz + (ny+1)*nz + k] = boundary_value_w[face].value(nx,ny,k,t);
    }
    
    face=2; // front face
    for (size_t k = 1; k < nx-1; k++)
    {
        grid.w[nx*(ny+1)*nz + (ny+1)*nz + k] = boundary_value_w[face].value(nx,ny,k,t);
    }

    face=3;

    face=4;

    face=5;

}
*/

void boundary::update_boundary(double t){

    size_t face=0; // left face

    // TODO: check
    for (size_t k=0; k < nz; k++)
    {
        grid.w[k] = boundary_value_w[face].value(0,0,k,t); //0 0 0.5
        grid.v[k] = boundary_value_v[face].value(0,0,k,t); //0 0.5 0
    }
    grid.v[nz] = boundary_value_v[face].value(0,0,nz,t);
    for (size_t j = 1; j < ny; j++)
    {
        grid.v[j*(nz+1)] = boundary_value_v[face].value(0,j,0,t);
        grid.w[j*nz] = boundary_value_w[face].value(0,j,0,t); 
        for(size_t k = 1; k < nz; k++)
        {
            grid.v[j*(nz+1) + k] = boundary_value_v[face].value(0,j,k,t);
            grid.w[j*nz + k] = boundary_value_w[face].value(0,j,k,t); 
            grid.u[j*(nz+1) + k] = approximate_boundary_u(0,j,k,t,face);
        }
        grid.v[j*(nz+1) + nz] = boundary_value_v[face].value(0,j,nz,t);
    }

    for (size_t k = 0; k < nz; k++)
    {
        grid.w[(ny+1)*nz + k] = boundary_value_w[face].value(0,ny,k,t);
    }
    
    
    //Missing: front face, back face, upper face, lower face
    face=1; // right face
    for (size_t k=0; k < nz; k++)
    {
        grid.v[nx*ny*(nz+1) + k] = boundary_value_v[face].value(nx,0,k,t); 
        grid.w[nx*(ny+1)*nz + k] = boundary_value_w[face].value(nx,0,k,t); 
    }
    grid.v[nx*ny*(nz+1) + nz] = boundary_value_v[face].value(nx,0,nz,t);
    for (size_t j = 1; j < ny; j++)
    {
        grid.v[nx*ny*(nz+1) + j*(nz+1)] = boundary_value_v[face].value(nx,j,0,t);
        grid.w[nx*(ny+1)*nz + j*nz] = boundary_value_w[face].value(nx,j,0,t); 
        for(size_t k = 1; k < nz; k++)
        {
            grid.v[nx*ny*(nz+1) + j*(nz+1) + k] = boundary_value_v[face].value(nx,j,k,t);
            grid.w[nx*(ny+1)*nz + j*nz + k] = boundary_value_w[face].value(nx,j,k,t); 
            grid.u[(nx-1)*(ny+1)*(nz+1) + j*(nz+1) + k] = approximate_boundary_u(nx-1,j,k,t,face);
        }
        grid.v[nx*ny*(nz+1) + j*(nz+1) + nz] = boundary_value_v[face].value(nx,j,nz,t);
    }

    for (size_t k = 0; k < nz; k++)
    {
        grid.w[nx*(ny+1)*nz + (ny+1)*nz + k] = boundary_value_w[face].value(nx,ny,k,t);
    }
    
    face=2; // front face
    grid.u[0] =  boundary_value_u[face].value(0,0,0,t);
    for (size_t k=1; k < nx; k++)
    {
        grid.u[(ny+1) * (nz+1) * k] = boundary_value_v[face].value(k,0,0,t);
        grid.v[(nz+1) * ny * k] = boundary_value_u[face].value(k,0,0,t);
    }
    
    for (size_t j = 1; j < ny; j++)
    {
        grid.u[j * (nz+1)] = boundary_value_v[face].value(0,j,0,t);
        for(size_t k = 1; k < nx; k++)
        {
            grid.u[(ny+1) * (nz+1) * k + j*(nz+1)] = boundary_value_u[face].value(k,j,0,t);
            grid.v[ny * (nz+1) * k + j*(nz+1)] = boundary_value_v[face].value(k,j,0,t); 
            grid.w[(ny+1) * nz * k + j*nz] = approximate_boundary_w(k,j,0,t,face);
        }
        grid.u[(ny+1) * (nz+1) * (nx) + j*(nz+1)] = boundary_value_u[face].value(nx,j,0,t);
    }

    for (size_t k = 0; k < nx; k++)
    {
        grid.u[(ny+1) * (nz+1) * k + ny*(nz+1)] = boundary_value_w[face].value(k,ny,0,t);
    }

    face=3; //back face
    grid.u[nz] =  boundary_value_u[face].value(0,0,nz,t);
    for (size_t k=1; k < nx; k++)
    {
        grid.u[(ny+1) * (nz+1) * k + nz] = boundary_value_v[face].value(k,0,nz,t);
        grid.v[(nz+1) * ny * k + nz] = boundary_value_u[face].value(k,0,nz,t);
    }
    
    for (size_t j = 1; j < ny; j++)
    {
        grid.u[j * (nz+1) + nz] = boundary_value_v[face].value(0,j,nz,t);
        for(size_t k = 1; k < nx; k++)
        {
            grid.u[(ny+1) * (nz+1) * k + j*(nz+1) + nz] = boundary_value_u[face].value(k,j,nz,t);
            grid.v[ny * (nz+1) * k + j*(nz+1) + nz] = boundary_value_v[face].value(k,j,nz,t); 
            grid.w[(ny+1) * nz * k + j*nz + (nz - 1)] = approximate_boundary_w(k,j,nz,t,face);
        }
        grid.u[(ny+1) * (nz+1) * (nx) + j*(nz+1) + nz] = boundary_value_u[face].value(nx,j,nz,t);
    }

    for (size_t k = 0; k < nx; k++)
    {
        grid.u[(ny+1) * (nz+1) * k + ny*(nz+1) + nz] = boundary_value_w[face].value(k,ny,nz,t);
    }

    face=4; //lower face

    for (size_t k=1; k < nz; k++)
    {
        grid.u[k] = boundary_value_u[face].value(0,0,k,t);
    }
    
    for (size_t j = 1; j < nx; j++)
    {
        grid.w[(ny+1) * nz + j] = boundary_value_w[face].value(j,0,0,t);
        for(size_t k = 1; k < nz; k++)
        {
            grid.u[(ny+1) * (nz+1) * j + k] = boundary_value_u[face].value(j,0,k,t);
            grid.v[ny * (nz+1) * j + k] = boundary_value_v[face].value(j,0,k,t); 
            grid.w[(ny+1) * nz * j + k] = approximate_boundary_w(j,0,k,t,face);
        }
    }

    face=5;
    for (size_t k=1; k < nz; k++)
    {
        grid.u[ny * (nz+1) + k] = boundary_value_u[face].value(0,ny,k,t);
    }
    
    for (size_t j = 1; j < nx; j++)
    {
        grid.w[ny * nz + (ny+1) * nz + j] = boundary_value_w[face].value(j,ny,0,t);
        for(size_t k = 1; k < nz; k++)
        {
            grid.u[ny * (nz + 1) + (ny+1) * (nz+1) * j + k] = boundary_value_u[face].value(j,ny,k,t);
            grid.v[(ny-1) * nz + ny * (nz+1) * j + k] = boundary_value_v[face].value(j,ny,k,t); 
            grid.w[ny*nz + (ny+1) * nz * j + k] = approximate_boundary_w(j,ny,k,t,face);
        }
    }

}   



double boundary::approximate_boundary_u(size_t x, size_t y, size_t z, double t,size_t face) {
    double dv = ((boundary_value_v[face].value(x * dx, y * dy + prec, z,t)) -
                ((boundary_value_v[face].value(x * dx, y * dy - prec, z,t)))) / (2*prec);
    
    double dw =  ((boundary_value_w[face].value(x * dx, y * dy, z * dz + prec,t)) -
                 ((boundary_value_w[face].value(x * dx, y * dy, z * dz - prec,t)))) / (2*prec);

    return  boundary_value_u[face].value(x * dx, y * dy, z * dz,t) - (dv + dw) * (dx/2);
}

double boundary::approximate_boundary_v(size_t x, size_t y, size_t z, double t,size_t face) {
    double du = ((boundary_value_u[face].value(x * dx + prec, y * dy, z,t)) -
                ((boundary_value_u[face].value(x * dx - prec, y * dy, z,t)))) / (2*prec);
    
    double dw =  ((boundary_value_w[face].value(x * dx, y * dy, z * dz + prec,t)) -
                 ((boundary_value_w[face].value(x * dx, y * dy, z * dz - prec,t)))) / (2*prec);

    return  boundary_value_v[face].value(x * dx, y * dy, z * dz,t) - (du + dw) * (dy/2);
}

double boundary::approximate_boundary_w(size_t x, size_t y, size_t z, double t,size_t face) {
    double du = ((boundary_value_u[face].value(x * dx + prec, y * dy, z,t)) -
                ((boundary_value_u[face].value(x * dx - prec, y * dy, z,t)))) / (2*prec);
    
    double dv = ((boundary_value_v[face].value(x * dx, y * dy + prec, z,t)) -
                ((boundary_value_v[face].value(x * dx, y * dy - prec, z,t)))) / (2*prec);
    
    return  boundary_value_w[face].value(x * dx, y * dy, z * dz,t) - (du + dv) * (dz/2);
}

void boundary::addFunction(size_t direction, Function x){
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
}

