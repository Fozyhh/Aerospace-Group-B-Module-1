#include <iostream>
#include "core.hpp"

int main(int argc, char *argv[])
{
    const double lx = 1;
    const double ly = 1;
    const double lz = 1;

    if(argc != 6)
    {
        std::cout << "Invalid arguments." << std::endl;
        std::cout << "Usage: nx ny nz dt T" << std::endl;
        return 0;
    }

    const unsigned int nx = std::atoi(argv[1]);
    const unsigned int ny = std::atoi(argv[2]);
    const unsigned int nz = std::atoi(argv[3]);
    const double dt = std::atof(argv[4]);
    const double T = std::atof(argv[5]);

    // errore al tempo 0.1 con n=10 = 0.147016 ?????????
    // errore al tempo 0.1 con n=20 = 0.149791 ????????? errore massimo tra 40 e 50 celle
    // errore al tempo 0.1 con n=50 = 0.150279
    // errore al tempo 0.1 con n=100 = 0.15024
    // una reuduzione nell'ordine di 10^-5 è giusta
    // secondo me la parte 0.1502 è una costannte dovuta a un errore nostro nel metodo numerico
    
    const double Re = 400.0;
    IcoNS icoNS(lx, ly, lz, nx, ny, nz, dt, T, Re, "input.txt", "output.txt");
    icoNS.preprocessing();
    icoNS.solve();

    return 0;
}
