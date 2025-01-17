#include <mpi.h>
#include <iostream>
#include <vector>
#include <fstream>
#include "../dependencies/2Decomp_C/C2Decomp.hpp"

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    C2Decomp *c2d;
    bool periodss[3] = {false, false, false};
    
    // choose pY and pX, recompile and run with the relative number of processes
    c2d = new C2Decomp(4, 6, 4, 2/*p in "our" y direction*/, 2/*p in "our" x direction*/, periodss);
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

    std::ofstream outFile("output_rank_" + std::to_string(rank) + ".txt");

    outFile << "from " << rank << ":" << std::endl;
    outFile << "Neighbours: " << c2d->neighbor[0][0] << " " << c2d->neighbor[0][1] << " " << c2d->neighbor[0][2] << " " << c2d->neighbor[0][3]<< " " << c2d->neighbor[0][4] << " " << c2d->neighbor[0][5] << std::endl;
    outFile << "xSize: " << xSize[0] << " " << xSize[1] << " " << xSize[2] << std::endl;
    outFile << "ySize: " << ySize[0] << " " << ySize[1] << " " << ySize[2] << std::endl;
    outFile << "zSize: " << zSize[0] << " " << zSize[1] << " " << zSize[2] << std::endl;

    double* gridZ, *gridY, *gridX;
    c2d->allocZ(gridZ);
    c2d->allocY(gridY);
    c2d->allocX(gridX);

    outFile << "Displaying usual view (X) from " << rank << ":" << std::endl;
    for (int i = 0; i < xSize[2]; i++) {
        for (int j = 0; j < xSize[1]; j++) {
            for (int k = 0; k < xSize[0]; k++) {
                gridX[i * (xSize[1]) * (xSize[0]) + j * (xSize[0]) + k] = 16*rank + i * (xSize[1]) * (xSize[0]) + j * (xSize[0]) + k;
                outFile << gridX[i * (xSize[1]) * (xSize[0]) + j * (xSize[0]) + k] << " ";
            }
            outFile<< std::endl;
        }
        outFile<< std::endl;
    }
    
    outFile << "-----------------------------" << std::endl;

    c2d->transposeX2Y_MajorIndex(gridX, gridY);

    outFile << "Displaying Y view from " << rank << ":" << std::endl;
    for (int k = 0; k < ySize[0]; k++) {
        for (int i = 0; i < ySize[2]; i++) {
            for (int j = 0; j < ySize[1]; j++) {
                outFile << gridY[k * (ySize[1]) * (ySize[2]) + i * (ySize[1]) + j] << " ";
            }
            outFile<< std::endl;
        }
        outFile << std::endl;
    }
    
    // for (int i = 0; i < ySize[0]*ySize[1]*ySize[2]; i++) {
    //     outFile << gridY[i] << std::endl;
    // }

    outFile << "-----------------------------" << std::endl;
    
    c2d->transposeY2Z_MajorIndex(gridY, gridZ);

    outFile << "Displaying Z view from " << rank << ":" << std::endl;
    for (int j = 0; j < zSize[1]; j++) {
        for (int k = 0; k < zSize[0]; k++) {
            for (int i = 0; i < zSize[2]; i++) {
                outFile << gridZ[j * (zSize[0]) * (zSize[2]) + k * (zSize[2]) + i] << " ";
            }
            outFile<< std::endl;
        }
        outFile << std::endl;
    }
    // for (int i = 0; i < zSize[2]; i++) {
    //     for (int j = 0; j < zSize[1]; j++) {
    //         for (int k = 0; k < zSize[0]; k++) {
    //             outFile << gridZ[j * (zSize[0]) * (zSize[2]) + k * (zSize[2]) + i] << std::endl;
    //         }
    //     }
    // }
    // c2d->transposeZ2Y_MajorIndex(gridX, gridY);
    // outFile << "Back to Y from " << rank << ":" << std::endl;
    // for (int i = 0; i < ySize[0]*ySize[1]*ySize[2]; i++) {
    //     outFile << gridY[i] << std::endl;
    // }

    // c2d->transposeY2X_MajorIndex(gridY, gridZ);

    // outFile << "Back to the start from " << rank << ":" << std::endl;
    // for (int i = 0; i < zSize[0]*zSize[1]*zSize[2]; i++) {
    //     outFile << gridZ[i] << std::endl;
    // }

    outFile.close();

    MPI_Finalize();
    return 0;
}