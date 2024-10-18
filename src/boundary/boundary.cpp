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


void boundary::update_boundary(double t){

    size_t face=0; // left face
    for (size_t k = 0; k < nz-1; k++)
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
        grid.u[(nz-1*nz) + nz-1] = boundary_value_u[face].value(0,(nz-1),nz-1,t);
    
    
    //Missing: right face, front face, back face, upper face, lower face
    face=1;

    face=2;

    face=3;

    face=4;

    face=5;

}   


double boundary::approximate_boundary_u(size_t x, size_t y, size_t z, double t,size_t face) {
    double dv = ((boundary_value_v[face].value(x * dx, y * dy + prec, z,t)) -
                ((boundary_value_v[face].value(x * dx, y * dy - prec, z,t)))) / (2*prec);
    
    double dw =  ((boundary_value_w[face].value(x * dx, y * dy, z * dz + prec,t)) -
                 ((boundary_value_w[face].value(x * dx, y * dy, z * dz - prec,t)))) / (2*prec);

    return  boundary_value_u[face].value(x * dx, y * dy, z * dz,t) - (dv + dw) * (dx/2);
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

