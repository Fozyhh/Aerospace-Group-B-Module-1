#include "../includes/core.hpp"
#include "../includes/boundary.hpp"
#include "../includes/grid.hpp"
#include "../includes/utils.hpp"

#define nx 4
#define ny 4
#define nz 4

int main(int argc, char const *argv[])
{
    Grid grid(nx,ny,nz);
    ExactSolution sol;

    Boundary b(grid,1,1,1);

    auto u_func= std::make_shared<Dirichlet>([&](double x, double y, double z, double t){
        return std::sin(x + 0.5) *std::cos(y)*std::sin(z)*std::sin(t);
    });
    auto v_func= std::make_shared<Dirichlet>([&](double x, double y, double z, double t){
        return std::cos(x) *std::sin(y+0.5)*std::sin(z)*std::sin(t);
    });
    auto w_func= std::make_shared<Dirichlet>([&](double x, double y, double z, double t){
        return 2*std::cos(x) *std::cos(y)*std::cos(z+0.5)*std::sin(t);
    });
    
    for (size_t i = 0; i < 6/*nfaces*/; i++)
    {
        b.addFunction(0,u_func);
        b.addFunction(1,v_func);
        b.addFunction(2,w_func);
    }
    double t =0.5;
    /*double x=0,y=1,z=1,prec=1.0/2.0,t=0.5;


    std::cout << b.boundary_value_v[0]->value(x, y - prec, z,t)<< std::endl;
    double dv = ((b.boundary_value_v[0]->value(x, y + prec, z,t)) - //std::cos(x) *std::sin(y)*std::sin(z)*std::sin(t); 1* 0.0087 * 0,017 * sin(0,5)
                ((b.boundary_value_v[0]->value(x, y - prec, z,t)))) / (2*prec);
    std::cout << dv << std::endl;

    double dw =  ((b.boundary_value_w[0]->value(x, y, z + prec ,t)) - //2*std::cos(x) *std::cos(y)*std::cos(z)*std::sin(t);
                 ((b.boundary_value_w[0]->value(x, y, z - prec,t)))) / (2*prec);
    std::cout << dw << std::endl;

    double bv = b.boundary_value_u[0]->value(x, y, z, t) - (dv + dw) * (1.0/2.0); //std::sin(x) *std::cos(y)*std::sin(z)*std::sin(t);

    std::cout<< "(" << bv<< ")-(" << sol.value_x(x,y,z,t) << ") " << std::endl;*/
    b.update_boundary(t);

    std::cout << sol.value_x(3,0,1,0.5) << std::endl << std::endl;
    for (size_t i = 0; i < nx; i++)
    {
        for (size_t j = 0; j < ny+1; j++)
        {
            for (size_t k = 0; k < nz+1; k++)
            {
                std::cout << i << j << k << "(" << grid.u[i*(ny+1)*(nz+1) + j*(nz+1) +k] << ")-(" << sol.value_x((i+0.5),j,k,t) << ") ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
    

    


    return 0;
}
