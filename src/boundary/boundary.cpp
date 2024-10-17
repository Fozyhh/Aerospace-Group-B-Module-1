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


void boundary::update_boundary(){
    //left face only
    
    for (size_t k = 0; k < nz-1; k++)
    {
        grid.u[k] = boundary_value_u(0,0,k);
        grid.w[k] = boundary_value_w(0,0,k); 
    }
        grid.u[nz-1] = boundary_value_u(0,0,nz-1);
    for (size_t j = 1; j < ny-1; j++)
    {
        grid.u[j*nz] = boundary_value_u(0,j,0);
        grid.v[j*nz] = boundary_value_v(0,j,0);
        grid.w[j*nz] = boundary_value_w(0,j,0); 
        for(size_t k = 1; k < nz - 1; k++){
            grid.u[j*nz + k] = approximate_boundary_u(0,j,k);
            grid.v[j*nz + k] = boundary_value_v(0,j,k);
            grid.w[j*nz + k] = boundary_value_w(0,j,k); 
        }
        grid.u[j*nz + nz-1] = boundary_value_u(0,j,nz-1);
        grid.v[j*nz + nz-1] = boundary_value_v(0,j,nz-1);
         
    }
    for (size_t k = 0; k < nz-1; k++)
    {
        grid.u[(ny-1)*nz + k] = boundary_value_u(0,(ny-1),k);
        grid.w[(ny-1)*nz + k] = boundary_value_w(0,(ny-1),k); 
    }
        grid.u[(nz-1*nz) + nz-1] = boundary_value_u(0,(nz-1),nz-1);
    
}


double boundary::approximate_boundary_u(size_t x, size_t y, size_t z) {
    double dv = ((boundary_value_v(x * dx, y * dy + prec, z)) -
                ((boundary_value_v(x * dx, y * dy - prec, z)))) / (2*prec);
    
    double dw =  ((boundary_value_w(x * dx, y * dy, z * dz + prec)) -
                 ((boundary_value_w(x * dx, y * dy, z * dz - prec)))) / (2*prec);

    return  boundary_value_u(x * dx, y * dy, z * dz) - (dv + dw) * (dx/2);
}