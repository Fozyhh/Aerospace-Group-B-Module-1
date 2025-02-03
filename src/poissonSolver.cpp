#include "poissonSolver.hpp"

void NeumannPoissonSolver::solvePoisson(Real *F)
{
    // 1. X-Direction Transforms
    neumann = FFTW_PREFIX(plan_r2r_1d)(xSize[0], nullptr, nullptr, FFTW_REDFT00, FFTW_ESTIMATE);
    for (int i = 0; i < xSize[2] * xSize[1]; i++)
    {
        FFTW_PREFIX(execute_r2r)(neumann, &F[i * xSize[0]], &F[i * xSize[0]]);
    }
    // FFTW_PREFIX(destroy_plan)(neumann);

    c2d->transposeX2Y_MajorIndex(F, py);

    // 2. Y-Direction Transforms
    neumann = FFTW_PREFIX(plan_r2r_1d)(ySize[1], nullptr, nullptr, FFTW_REDFT00, FFTW_ESTIMATE);
    for (int i = 0; i < ySize[0] * ySize[2]; i++)
    {
        FFTW_PREFIX(execute_r2r)(neumann, &py[i * ySize[1]], &py[i * ySize[1]]);
    }
    // FFTW_PREFIX(destroy_plan)(neumann);

    c2d->transposeY2Z_MajorIndex(py, pz);

    // 3. Z-Direction Transforms
    neumann = FFTW_PREFIX(plan_r2r_1d)(zSize[2], nullptr, nullptr, FFTW_REDFT00, FFTW_ESTIMATE);
    for (int i = 0; i < zSize[1] * zSize[0]; i++)
    {
        FFTW_PREFIX(execute_r2r)(neumann, &pz[i * zSize[2]], &pz[i * zSize[2]]);
    }
    // FFTW_PREFIX(destroy_plan)(neumann);

    // Divide by the eigenvalues
    for (int k = 0; k < zSize[0]; k++)
    {
        for (int j = 0; j < zSize[1]; j++)
        {
            for (int i = 0; i < zSize[2]; i++)
            {
                pz[j * (zSize[1]) * (zSize[2]) + k * (zSize[2]) + i] /= (2 / (DX * DX) * (std::cos(i * M_PI / (c2d->nxGlobal - 1)) - 1) +
                                                                         2 / (DY * DY) * (std::cos((j + c2d->coord[0] * zSize[1]) * M_PI / (c2d->nyGlobal - 1)) - 1) +
                                                                         2 / (DZ * DZ) * (std::cos((k + c2d->coord[1] * zSize[0]) * M_PI / (c2d->nzGlobal - 1)) - 1));
            }
        }
    }

    if (c2d->nRank == 0)
    {
        pz[0] = 0.0;
    }

    // 1. Z-Direction Transforms
    neumann = FFTW_PREFIX(plan_r2r_1d)(zSize[2], nullptr, nullptr, FFTW_REDFT00, FFTW_ESTIMATE);
    for (int i = 0; i < zSize[1] * zSize[0]; i++)
    {
        FFTW_PREFIX(execute_r2r)(neumann, &pz[i * zSize[2]], &pz[i * zSize[2]]);
    }
    // FFTW_PREFIX(destroy_plan)(neumann);

    c2d->transposeZ2Y_MajorIndex(pz, py);

    // 2. Y-Direction Transforms
    neumann = FFTW_PREFIX(plan_r2r_1d)(ySize[1], nullptr, nullptr, FFTW_REDFT00, FFTW_ESTIMATE);
    for (int i = 0; i < ySize[0] * ySize[2]; i++)
    {
        FFTW_PREFIX(execute_r2r)(neumann, &py[i * ySize[1]], &py[i * ySize[1]]);
    }
    // FFTW_PREFIX(destroy_plan)(neumann);

    c2d->transposeY2X_MajorIndex(py, F);

    // 3. X-Direction Transforms
    neumann = FFTW_PREFIX(plan_r2r_1d)(xSize[0], nullptr, nullptr, FFTW_REDFT00, FFTW_ESTIMATE);
    for (int i = 0; i < xSize[2] * xSize[1]; i++)
    {
        FFTW_PREFIX(execute_r2r)(neumann, &F[i * xSize[0]], &F[i * xSize[0]]);
    }
    // FFTW_PREFIX(destroy_plan)(neumann);

    // Normalization
    Real normalization_factor1 = 2.0 * (xSize[0] - 1) * 2.0 * (ySize[1] - 1) * 2.0 * (zSize[2] - 1);
    for (int i = 0; i < xSize[0] * xSize[1] * xSize[2]; i++)
    {
        F[i] /= normalization_factor1;
    }
}

void PeriodicPoissonSolver::solvePoisson(Real *F)
{

    // 1. X-Direction Transforms
    neumann = FFTW_PREFIX(plan_r2r_1d)(xSize[0], nullptr, nullptr, FFTW_REDFT00, FFTW_ESTIMATE);
    for (int i = 0; i < xSize[2] * xSize[1]; i++)
    {
        FFTW_PREFIX(execute_r2r)(neumann, &F[i * xSize[0]], &F[i * xSize[0]]);
    }
    // FFTW_PREFIX(destroy_plan)(neumann);

    c2d->transposeX2Y_MajorIndex(F, py);

    // 2. Y-Direction Transforms
    neumann = FFTW_PREFIX(plan_r2r_1d)(ySize[1], nullptr, nullptr, FFTW_REDFT00, FFTW_ESTIMATE);
    for (int i = 0; i < ySize[0] * ySize[2]; i++)
    {
        FFTW_PREFIX(execute_r2r)(neumann, &py[i * ySize[1]], &py[i * ySize[1]]);
    }
    // FFTW_PREFIX(destroy_plan)(neumann);

    c2d->transposeY2Z_MajorIndex(py, pz);

    // Version using the total complex and the helper

    // // 3. Z-Direction Transforms
    // neumann = FFTW_PREFIX(plan_dft_r2c_1d(zSize[2], nullptr, nullptr, FFTW_ESTIMATE);

    // for (int i = 0; i < zSize[1] * zSize[0]; i++)
    // {
    //     FFTW_PREFIX(execute_dft_r2c(neumann, &pz[i * zSize[2]], &helper[i * (zSize[2]/2+1)]);
    // }
    // FFTW_PREFIX(destroy_plan(neumann);

    // // Divide by the eigenvalues
    // for (int k = 0; k < zSize[0]; k++)
    // {
    //     for (int j = 0; j < zSize[1]; j++)
    //     {
    //         for (int i = 0; i < zSize[2]/2+1; i++)
    //         {
    //             helper[j * (zSize[1]) * (zSize[2]/2+1) + k * (zSize[2]/2+1) + i][0] /= (2 / (DX * DX) * (std::cos(2.0*i * M_PI / (c2d->nxGlobal - 1)) - 1) +
    //                                                                      2 / (DY * DY) * (std::cos((j + c2d->coord[0] * zSize[1]) * M_PI / (c2d->nyGlobal - 1)) - 1) +
    //                                                                      2 / (DZ * DZ) * (std::cos((k + c2d->coord[1] * zSize[0]) * M_PI / (c2d->nzGlobal - 1)) - 1));
    //             helper[j * (zSize[1]) * (zSize[2]/2+1) + k * (zSize[2]/2+1) + i][1] /= (2 / (DX * DX) * (std::cos(2.0*i * M_PI / (c2d->nxGlobal - 1)) - 1) +
    //                                                                             2 / (DY * DY) * (std::cos((j + c2d->coord[0] * zSize[1]) * M_PI / (c2d->nyGlobal - 1)) - 1) +
    //                                                                             2 / (DZ * DZ) * (std::cos((k + c2d->coord[1] * zSize[0]) * M_PI / (c2d->nzGlobal - 1)) - 1));
    //         }
    //     }
    // }

    // if (c2d->nRank == 0)
    // {
    //     helper[0][0] = 0.0;
    //     helper[0][1] = 0.0;
    // }

    // // 1. Z-Direction Transforms
    // neumann = FFTW_PREFIX(plan_dft_c2r_1d(zSize[2], nullptr, nullptr, FFTW_ESTIMATE);

    // for (int i = 0; i < zSize[1] * zSize[0]; i++)
    // {
    //     FFTW_PREFIX(execute_dft_c2r(neumann, &helper[i * (zSize[2] / 2 + 1)], &pz[i * zSize[2]]);
    // }
    // FFTW_PREFIX(destroy_plan(neumann);

    // Version using the Half complex

    // 3. Z-Direction Transforms
    // neumann = FFTW_PREFIX(plan_dft_r2c_1d(zSize[2], nullptr, nullptr, FFTW_ESTIMATE);
    neumann = FFTW_PREFIX(plan_r2r_1d)(zSize[2], nullptr, nullptr, FFTW_R2HC, FFTW_ESTIMATE);

    for (int i = 0; i < zSize[1] * zSize[0]; i++)
    {
        // FFTW_PREFIX(execute_dft_r2c(neumann, &pz[i * zSize[2]], &helper[i * (zSize[2]/2+1)]);
        FFTW_PREFIX(execute_r2r)(neumann, &pz[i * zSize[2]], &pz[i * zSize[2]]);
    }
    // FFTW_PREFIX(destroy_plan)(neumann);

    // Divide by the eigenvalues
    for (int k = 0; k < zSize[0]; k++)
    {
        for (int j = 0; j < zSize[1]; j++)
        {
            pz[j * (zSize[1]) * (zSize[2]) + k * (zSize[2]) + 0] /= (2 / (DX * DX) * (std::cos(2.0 * 0 * M_PI / (c2d->nxGlobal - 1)) - 1) +
                                                                     2 / (DY * DY) * (std::cos((j + c2d->coord[0] * zSize[1]) * M_PI / (c2d->nyGlobal - 1)) - 1) +
                                                                     2 / (DZ * DZ) * (std::cos((k + c2d->coord[1] * zSize[0]) * M_PI / (c2d->nzGlobal - 1)) - 1));

            for (int i = 1; i < zSize[2] / 2; i++)

            {

                pz[j * (zSize[1]) * (zSize[2]) + k * (zSize[2]) + i] /= (2 / (DX * DX) * (std::cos(2.0 * i * M_PI / (c2d->nxGlobal - 1)) - 1) +
                                                                         2 / (DY * DY) * (std::cos((j + c2d->coord[0] * zSize[1]) * M_PI / (c2d->nyGlobal - 1)) - 1) +
                                                                         2 / (DZ * DZ) * (std::cos((k + c2d->coord[1] * zSize[0]) * M_PI / (c2d->nzGlobal - 1)) - 1));
            }
            pz[j * (zSize[1]) * (zSize[2]) + k * (zSize[2]) + zSize[2] / 2] /= (2 / (DX * DX) * (std::cos(2.0 * zSize[2] / 2 * M_PI / (c2d->nxGlobal - 1)) - 1) +
                                                                                2 / (DY * DY) * (std::cos((j + c2d->coord[0] * zSize[1]) * M_PI / (c2d->nyGlobal - 1)) - 1) +
                                                                                2 / (DZ * DZ) * (std::cos((k + c2d->coord[1] * zSize[0]) * M_PI / (c2d->nzGlobal - 1)) - 1));

            for (int i = zSize[2] / 2 + 1; i < zSize[2]; i++)
            {
                pz[j * (zSize[1]) * (zSize[2]) + k * (zSize[2]) + (i)] /= (2 / (DX * DX) * (std::cos(2.0 * (zSize[2] - i) * M_PI / (c2d->nxGlobal - 1)) - 1) +
                                                                           2 / (DY * DY) * (std::cos((j + c2d->coord[0] * zSize[1]) * M_PI / (c2d->nyGlobal - 1)) - 1) +
                                                                           2 / (DZ * DZ) * (std::cos((k + c2d->coord[1] * zSize[0]) * M_PI / (c2d->nzGlobal - 1)) - 1));
            }
        }
    }

    if (c2d->nRank == 0)
    {
        pz[0] = 0.0;
    }

    // 1. Z-Direction Transforms
    neumann = FFTW_PREFIX(plan_r2r_1d)(zSize[2], nullptr, nullptr, FFTW_HC2R, FFTW_ESTIMATE);

    for (int i = 0; i < zSize[1] * zSize[0]; i++)
    {
        FFTW_PREFIX(execute_r2r)(neumann, &pz[i * zSize[2]], &pz[i * zSize[2]]);
    }
    // FFTW_PREFIX(destroy_plan)(neumann);

    // Remaining Poisson Solver

    c2d->transposeZ2Y_MajorIndex(pz, py);

    // 2. Y-Direction Transforms
    neumann = FFTW_PREFIX(plan_r2r_1d)(ySize[1], nullptr, nullptr, FFTW_REDFT00, FFTW_ESTIMATE);
    for (int i = 0; i < ySize[0] * ySize[2]; i++)
    {
        FFTW_PREFIX(execute_r2r)(neumann, &py[i * ySize[1]], &py[i * ySize[1]]);
    }
    // FFTW_PREFIX(destroy_plan)(neumann);

    c2d->transposeY2X_MajorIndex(py, F);

    // 3. X-Direction Transforms
    neumann = FFTW_PREFIX(plan_r2r_1d)(xSize[0], nullptr, nullptr, FFTW_REDFT00, FFTW_ESTIMATE);
    for (int i = 0; i < xSize[2] * xSize[1]; i++)
    {
        FFTW_PREFIX(execute_r2r)(neumann, &F[i * xSize[0]], &F[i * xSize[0]]);
    }
    // FFTW_PREFIX(destroy_plan)(neumann);

    // Normalization
    Real normalization_factor1 = 2.0 * (xSize[0] - 1) * 2.0 * (ySize[1] - 1) * 2.0 * (zSize[2] - 1);
    for (int i = 0; i < xSize[0] * xSize[1] * xSize[2]; i++)
    {
        F[i] /= normalization_factor1;
    }
}