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

    // TRY WITH ONE PROCESSOR ONLY
    c2d = new C2Decomp(4, 4, 4, 1, 1, periodss);
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

    std::cout << "xSize: " << xSize[0] << " " << xSize[1] << " " << xSize[2] << std::endl;
    std::cout << "ySize: " << ySize[0] << " " << ySize[1] << " " << ySize[2] << std::endl;
    std::cout << "zSize: " << zSize[0] << " " << zSize[1] << " " << zSize[2] << std::endl;

    double* gridZ, *gridY, *gridX;
    c2d->allocZ(gridZ);
    c2d->allocY(gridY);
    c2d->allocX(gridX);

    std::cout << "Displaying usual view (Z)" << std::endl;
    for (int i = 0; i < zSize[0]; i++) {
        for (int j = 0; j < zSize[1]; j++) {
            for (int k = 0; k < zSize[2]; k++) {
                gridZ[i * (zSize[1]) * (zSize[2]) + j * (zSize[2]) + k] = i * (zSize[1]) * (zSize[2]) + j * (zSize[2]) + k;
                std::cout << gridZ[i * (zSize[1]) * (zSize[2]) + j * (zSize[2]) + k] << std::endl;
            }
        }
    }
    std::cout << "-----------------------------" << std::endl;

    std::cout << "Displaying Y view" << std::endl;
    c2d->transposeZ2Y_MajorIndex(gridZ, gridY);

    for (int i = 0; i < ySize[0]; i++) {
        for (int j = 0; j < ySize[1]; j++) {
            for (int k = 0; k < ySize[2]; k++) {
                std::cout << gridY[i * (ySize[1]) * (ySize[2]) + j * (ySize[2]) + k] << std::endl;
            }
        }
    }
    std::cout << "-----------------------------" << std::endl;

    std::cout << "Displaying X view" << std::endl;
    c2d->transposeY2X_MajorIndex(gridY, gridX);

    for (int i = 0; i < xSize[0]; i++) {
        for (int j = 0; j < xSize[1]; j++) {
            for (int k = 0; k < xSize[2]; k++) {
                std::cout << gridX[i * (xSize[1]) * (xSize[2]) + j * (xSize[2]) + k] << std::endl;
            }
        }
    }
/*
    // TRY WITH MULTIPLE PROCESSORS
    c2d = new C2Decomp(4, 4, 4, 2, 1, periodss);
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

    outFile << "from 0:" << std::endl;
    outFile << "xSize: " << xSize[0] << " " << xSize[1] << " " << xSize[2] << std::endl;
    outFile << "ySize: " << ySize[0] << " " << ySize[1] << " " << ySize[2] << std::endl;
    outFile << "zSize: " << zSize[0] << " " << zSize[1] << " " << zSize[2] << std::endl;

    double* gridZ, *gridY, *gridX;
    c2d->allocZ(gridZ);
    c2d->allocY(gridY);
    c2d->allocX(gridX);

    if(rank == 0) {
        outFile << "Displaying usual view (Z) from 0" << std::endl;
        for (int i = 0; i < zSize[0]; i++) {
            for (int j = 0; j < zSize[1]; j++) {
                for (int k = 0; k < zSize[2]; k++) {
                    gridZ[i * (zSize[1]) * (zSize[2]) + j * (zSize[2]) + k] = i * (zSize[1]) * (zSize[2]) + j * (zSize[2]) + k;
                    outFile << gridZ[i * (zSize[1]) * (zSize[2]) + j * (zSize[2]) + k] << std::endl;
                }
            }
        }
    }

    if(rank == 1) {
        outFile << "Displaying usual view (Z) from 1" << std::endl;
        for (int i = 0; i < zSize[0]; i++) {
            for (int j = 0; j < zSize[1]; j++) {
                for (int k = 0; k < zSize[2]; k++) {
                    gridZ[i * (zSize[1]) * (zSize[2]) + j * (zSize[2]) + k] = 32 + i * (zSize[1]) * (zSize[2]) + j * (zSize[2]) + k;
                    outFile << gridZ[i * (zSize[1]) * (zSize[2]) + j * (zSize[2]) + k] << std::endl;
                }
            }
        }
    }
    
    outFile << "-----------------------------" << std::endl;

    c2d->transposeZ2Y_MajorIndex(gridZ, gridY);

    if(rank == 0) {
        outFile << "Displaying Y view from 0" << std::endl;
        for (int i = 0; i < ySize[0]; i++) {
            for (int j = 0; j < ySize[1]; j++) {
                for (int k = 0; k < ySize[2]; k++) {
                    outFile << gridY[i * (ySize[1]) * (ySize[2]) + j * (ySize[2]) + k] << std::endl;
                }
            }
        }
    }

    if(rank == 1) {
        outFile << "Displaying Y view from 1" << std::endl;
        for (int i = 0; i < ySize[0]; i++) {
            for (int j = 0; j < ySize[1]; j++) {
                for (int k = 0; k < ySize[2]; k++) {
                    outFile << gridY[i * (ySize[1]) * (ySize[2]) + j * (ySize[2]) + k] << std::endl;
                }
            }
        }
    }

    outFile << "-----------------------------" << std::endl;
    
    c2d->transposeY2X_MajorIndex(gridY, gridX);

    if(rank == 0) {
        outFile << "Displaying X view from 0" << std::endl;
        for (int i = 0; i < xSize[0]; i++) {
            for (int j = 0; j < xSize[1]; j++) {
                for (int k = 0; k < xSize[2]; k++) {
                    outFile << gridX[i * (xSize[1]) * (xSize[2]) + j * (xSize[2]) + k] << std::endl;
                }
            }
        }
    }

    if(rank == 1) {
        outFile << "Displaying X view from 1" << std::endl;
        for (int i = 0; i < xSize[0]; i++) {
            for (int j = 0; j < xSize[1]; j++) {
                for (int k = 0; k < xSize[2]; k++) {
                    outFile << gridX[i * (xSize[1]) * (xSize[2]) + j * (xSize[2]) + k] << std::endl;
                }
            }
        }
    }

    outFile.close();
*/
    MPI_Finalize();
    return 0;
}