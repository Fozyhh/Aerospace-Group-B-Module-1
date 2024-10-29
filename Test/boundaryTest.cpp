#include "../includes/core.hpp"
#include "../includes/boundary.hpp"
#include "../includes/grid.hpp"
#include "../includes/utils.hpp"

#define nx 10
#define ny 10
#define nz 10

#define dx 0.01
#define dy 0.01
#define dz 0.01


int main(int argc, char const *argv[])
{
    Grid grid(nx,ny,nz);
    ExactSolution sol;

    Boundary b(grid,dx,dy,dz);

    auto u_func= std::make_shared<Dirichlet>([&](double x, double y, double z, double t){
        return std::sin(x + (dx/2.0)) *std::cos(y)*std::sin(z)*std::sin(t);
    });
    auto v_func= std::make_shared<Dirichlet>([&](double x, double y, double z, double t){
        return std::cos(x) *std::sin(y+(dy/2.0))*std::sin(z)*std::sin(t);
    });
    auto w_func= std::make_shared<Dirichlet>([&](double x, double y, double z, double t){
        return 2*std::cos(x) *std::cos(y)*std::cos(z+(dz/2.0))*std::sin(t);
    });
    
    for (size_t i = 0; i < 6/*nfaces*/; i++)
    {
        b.addFunction(0,u_func);
        b.addFunction(1,v_func);
        b.addFunction(2,w_func);
    }
    double t =0.5;

    std::cout << "u0b:" << b.boundary_value_u[0]->value(3.5, 1, 1, t) << sol.value_x(4,1,1,0.5) << std::endl;

    b.update_boundary(t);



    //std::cout <<"-" << sol.value_x(3.5,1,1,0.5) << std::endl << std::endl;
    for (double i = 0; i < nx; i++)
    {
        for (double j = 0; j < ny+1; j++)
        {
            for (double k = 0; k < nz+1; k++)
            {
                std::cout << i << j << k << "(" << grid.u[i*(ny+1)*(nz+1) + j*(nz) +k] << ")-(" << sol.value_x(i+0.5,j,k,t) << ") ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
    

    


    return 0;
}
