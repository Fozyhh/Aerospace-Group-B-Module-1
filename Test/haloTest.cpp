#include <mpi.h>
#include <iostream>
#include <vector>
#include <fstream>
#include "../dependencies/2Decomp_C/C2Decomp.hpp"

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    C2Decomp *c2d;
    bool periodss[3] = {false, true, true};
    
    // choose pY and pX, recompile and run with the relative number of processes
    c2d = new C2Decomp(4, 4, 4, 2,2, periodss);
    int xSize[3], ySize[3], zSize[3];

    xSize[0] = c2d->xSize[0];
    xSize[1] = c2d->xSize[1];
    xSize[2] = c2d->xSize[2];

    ySize[0] = c2d->ySize[0];
    ySize[1] = c2d->ySize[1];
    ySize[2] = c2d->ySize[2];   

    zSize[0] = c2d->zSize[0];
    zSize[1] = c2d->zSize[1];
    zSize[2] = c2d->zSize[2];

    double* gridp;
    c2d->allocX(gridp);
    for (size_t i = 0; i < xSize[2]; i++)
    {
        for (size_t j = 0; j < xSize[1]; j++)
        {
            for (size_t k = 0; k < xSize[0]; k++)
            {
                gridp[i*xSize[1]*xSize[0] + j*xSize[0] + k] = (i )*xSize[1]*xSize[0] + j*xSize[0] + k; 
            }
            
        }
        
    }
    // std::cout << rank << " " << c2d->decompMain.zsz[0] << " " << c2d->decompMain.zsz[1]<< " " << c2d->decompMain.zsz[2]<< std::endl;
    //std::cout << rank << " " << xSize[0] << " " << xSize[1]<< " " << xSize[2]<< std::endl<<std::endl;
    double* halo_p;
    c2d->updateHalo(gridp,halo_p,1,0);
        if(rank==0){
            std::cout << std::endl;
            
            for(int k = 0; k < xSize[2]; k++){
                
            for (int j = 0; j < xSize[1]; j++){
                for (int i = 0; i < xSize[0]; i++)
                    {
                        std::cout << gridp[k*xSize[1]*xSize[0] + j*xSize[0] + i] << " ";
                    }
                    std::cout << std::endl;
                }
                std::cout << std::endl;
            }
            std::cout << "Halos: "<< std::endl;
            for(int i = 0; i < xSize[2] +2; i++){
                for (int j = 0; j < xSize[1] + 2; j++){
                    for (int k = 0; k < xSize[0]; k++) //DIMENSIONE LUNGA! costruttore di 2deco come nz,ny,nz?
                        
                        {
                            std::cout << halo_p[(i)*(xSize[1]+2)*(xSize[0]) + (j)*xSize[0] + k] << " ";
                        }
                        std::cout << std::endl;
                    }
                    std::cout << std::endl;
            }
            
            //int stop; std::cin >> stop;
        }
        c2d->deallocXYZ(halo_p);
    MPI_Finalize();
        return 0;
}


