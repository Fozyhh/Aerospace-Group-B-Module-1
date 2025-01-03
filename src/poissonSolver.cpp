#include "poissonSolver.hpp"


//TODO: 2Decomp parallellizazione
void PoissonSolver::solveDirichletPoisson(std::vector<Real>& F_dP, fftw_complex *FD)
{
    bool periodicBC[3] = {periodicX, periodicY, periodicZ};
    // C2Decomp *c2d;
    // c2d = new C2Decomp(NX, NY, NZ, 0, 0, periodicBC);
    
    // dP = fft(fft(fft(F)))
    fftw_plan forward = fftw_plan_dft_r2c_3d(NX, NY, NZ, F_dP.data(), FD, FFTW_ESTIMATE);
    fftw_execute(forward);
    
/*  if(periodicX){
        fftw_plan fftw_plan_dft_3d(int n0, int n1, int n2, fftw_complex *in, fftw_complex *out, int sign, unsigned flags);
    }else{
        dP = dct(F);
    }
    rotateGrid X->Y
    if(periodicY){
        dP = fft(dP);
    }else{
        dP = dct(dP);
    }
    rotateGrid Y->Z
    if(periodicZ){
        dP = fft(dP);
    }else{
        dP = dct(dP);
    }
*/

    // gotta change indexing to accomodate X->Y->Z
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            for (int k = 0; k < (NZ)/2+1; k++) {
                FD[i * (NY) * ((NZ)/2+1) + j * ((NZ)/2+1) + k][0] = FD[i * (NY) * ((NZ)/2+1) + j * ((NZ)/2+1) + k][0] /
                    (2/(DX*DX) * (std::cos(2 * i * M_PI / (NX)) - 1) +
                    2/(DY*DY) * (std::cos(2 * j * M_PI / (NY)) - 1) +
                    2/(DZ*DZ) * (std::cos(2 * k * M_PI / (NZ)) - 1));
                FD[i * (NY) * ((NZ)/2+1) + j * ((NZ)/2+1) + k][1] = FD[i * (NY) * ((NZ)/2+1) + j * ((NZ)/2+1) + k][1] /
                    (2/(DX*DX) * (std::cos(2 * i * M_PI / (NX)) - 1) +
                    2/(DY*DY) * (std::cos(2 * j * M_PI / (NY)) - 1) +
                    2/(DZ*DZ) * (std::cos(2 * k * M_PI / (NZ)) - 1));
            }
        }
    }
    FD[0][0] = 0.0;
    FD[0][1] = 0.0;

    // Inverse Fourier transform
    fftw_plan backward = fftw_plan_dft_c2r_3d(NX, NY, NZ, FD, F_dP.data(), FFTW_ESTIMATE);
    fftw_execute(backward);

    double normalization_factor = (NX) * (NY) * (NZ);
    for (int i=0; i < NX; i++){
        for (int j=0; j < NY; j++){
            for (int k=0; k < NZ; k++){
                F_dP[i * (NY) * (NZ) + j * (NZ) + k] = F_dP[i * (NY) * (NZ) + j * (NZ) + k] / normalization_factor;
            }
        }
    }

    fftw_destroy_plan(forward);
    fftw_destroy_plan(backward);

/*
    if(periodicZ){
        dP = fft(dP);
    }else{
        dP = dct(dP);
    }
    rotateGrid Z->Y
    if(periodicY){
        dP = fft(dP);
    }else{
        dP = dct(dP);
    }
    rotateGrid Y->X
    if(periodicX){
        dP = fft(dP);
    }else{
        dP = dct(dP);
    }
*/
}



void PoissonSolver::solveNeumannPoisson(double* F)
{
    double* py;
    double* px;

    c2d->allocY(py);
    c2d->allocX(px);

     

    for (int i = 0; i < zSize[0]; i++) 
    {
        for (int j = 0; j < zSize[1]; j++)
        {
            fftw_plan neumann = fftw_plan_r2r_1d(zSize[2], &F[i * zSize[1] * zSize[2] + j * zSize[2]], &F[i * zSize[1] * zSize[2] + j * zSize[2]], 
                                                 FFTW_REDFT00, FFTW_ESTIMATE);
            fftw_execute(neumann);
            fftw_destroy_plan(neumann);
        }
    }

    c2d->transposeZ2Y_MajorIndex(F, py);
    for (int i = 0; i < ySize[0]; i++) 
    {
        for (int j = 0; j < ySize[2]; j++)
        {
            fftw_plan neumann = fftw_plan_r2r_1d(ySize[1], &py[i * ySize[2] * ySize[1] + j * ySize[1]], &py[i * ySize[2] * ySize[1] + j * ySize[1]], 
                                                 FFTW_REDFT00, FFTW_ESTIMATE);
            fftw_execute(neumann);
            fftw_destroy_plan(neumann);
        }
    }

    c2d->transposeY2X_MajorIndex(py, px);
    for (int i = 0; i < xSize[1]; i++) 
    {
        for (int j = 0; j < xSize[2]; j++)
        {
            fftw_plan neumann = fftw_plan_r2r_1d(xSize[0], &px[i*xSize[2] * xSize[0] + j * xSize[0]], &px[i*xSize[2] * xSize[0] + j * xSize[0]], 
                                                 FFTW_REDFT00, FFTW_ESTIMATE);
            fftw_execute(neumann);
            fftw_destroy_plan(neumann);
        }
    }

    // Divide by the eigenvalues
    
    for (int i = 0; i < xSize[1]; i++) {
        for (int j = 0; j < xSize[2]; j++) {
            for (int k = 0; k < xSize[0]; k++) {
                px[i * (xSize[2]) * (xSize[0]) + j * (xSize[0]) + k] /= (2/(DX*DX) * (std::cos(i * M_PI / (xSize[0]-1)) - 1) +
                                                                        2/(DY*DY) * (std::cos(j * M_PI / (xSize[1]-1)) - 1) +
                                                                        2/(DZ*DZ) * (std::cos(k * M_PI / (xSize[2]-1)) - 1));
            }
        }
    }

    px[0] = 0.0;

    // Inverse Fourier transform
    for (int i = 0; i < xSize[1]; i++) 
    {
        for (int j = 0; j < xSize[2]; j++)
        {
            fftw_plan neumann = fftw_plan_r2r_1d(xSize[0], &px[i*xSize[2] *xSize[0] + j * xSize[0]], &px[i*xSize[2]* xSize[0] + j* xSize[0]], 
                                                 FFTW_REDFT00, FFTW_ESTIMATE);
            fftw_execute(neumann);
            fftw_destroy_plan(neumann);
        }
    }

    c2d->transposeX2Y_MajorIndex(px, py);
    for (int i = 0; i < ySize[0]; i++) 
    {
        for (int j = 0; j < ySize[2]; j++)
        {
            fftw_plan neumann = fftw_plan_r2r_1d(ySize[1], &py[i * ySize[2]* ySize[1] + j* ySize[1]], &py[i * ySize[2]* ySize[1] + j* ySize[1]], 
                                                 FFTW_REDFT00, FFTW_ESTIMATE);
            fftw_execute(neumann);
            fftw_destroy_plan(neumann);
        }
    }

    c2d->transposeY2Z_MajorIndex(py, F);
    for (int i = 0; i < zSize[0]; i++) 
    {
        for (int j = 0; j < zSize[1]; j++)
        {
            fftw_plan neumann = fftw_plan_r2r_1d(zSize[2], &F[i * zSize[1]* zSize[2] + j* zSize[2]], &F[i * zSize[1]* zSize[2] + j* zSize[2]], 
                                                 FFTW_REDFT00, FFTW_ESTIMATE);
            fftw_execute(neumann);
            fftw_destroy_plan(neumann);
        }
    }

    // Normalization
    double normalization_factor1 = 2.0 * (zSize[0]-1) * 2.0 * (zSize[1]-1) * 2.0 * (zSize[2]-1);
    for(int i = 0; i < zSize[0]; i++) {
        for (int j = 0; j < zSize[1]; j++) {
            for (int k = 0; k < zSize[2]; k++) {
                    F[i * (zSize[1]) * (zSize[2]) + j * (zSize[2]) + k] /= normalization_factor1;
            }
        }
    }
    
}