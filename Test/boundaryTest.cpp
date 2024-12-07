#include "../includes/core.hpp"
#include "../includes/boundary.hpp"
#include "../includes/grid.hpp"
#include "../includes/utils.hpp"

#define nx 5
#define ny 5
#define nz 5

#define dx 1
#define dy 1
#define dz 1


int main(int argc, char const *argv[])
{
    Grid grid(nx,ny,nz);
    ExactSolution sol;

    Boundary b(grid,dx,dy,dz);

    auto u_func= std::make_shared<Dirichlet>([&](double x, double y, double z, double t){
        return std::sin((x + 1.0/2.0) * dx) *std::cos(y * dy)*std::sin(z * dz)*std::sin(t);
    });
    auto v_func= std::make_shared<Dirichlet>([&](double x, double y, double z, double t){
        return std::cos(x * dx) *std::sin((y + 1.0/2.0) * dy)*std::sin(z * dz)*std::sin(t);
    });
    auto w_func= std::make_shared<Dirichlet>([&](double x, double y, double z, double t){
        return 2*std::cos(x * dx) *std::cos(y * dy)*std::cos((z + 1.0/2.0) * dz)*std::sin(t);
    });
    
    for (size_t i = 0; i < 6/*nfaces*/; i++)
    {
        b.addFunction(U,u_func);
        b.addFunction(V,v_func);
        b.addFunction(W,w_func);
    }
    double t =0.5;

    //std::cout << "u0b:" << b.boundary_value_u[0]->value(3.5, 1, 1, t) << sol.value_x(4,1,1,0.5) << std::endl;

    b.update_boundary(grid.u,grid.v,grid.w,t);



    //std::cout <<"-" << sol.value_x(3.5,1,1,0.5) << std::endl << std::endl;
    for (double i = 0; i < nx; i++)
    {
        for (double j = 0; j < ny+1; j++)
        {
            for (double k = 0; k < nz+1; k++)
            {
                std::cout << i << j << k << "(" << grid.u[i*(ny+1)*(nz+1) + j*(nz+1) +k] << ")-(" << sol.value_x((i + 1.0/2.0 ) *dx,j * dy,k *dz,t) << ") ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
    

    


    return 0;
}
