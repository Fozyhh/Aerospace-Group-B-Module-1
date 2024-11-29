#include "poissonSolver.hpp"

void PoissonSolver::solvePoisson(std::array<Real, NX * NY * NZ>& F_dP, fftw_complex *FD)
{
    bool periodicBC[3] = {periodicX, periodicY, periodicZ};
    //C2Decomp *c2d;
    //c2d = new C2Decomp(NX+1, NY+1, NZ+1, 0, 0, periodicBC);
    
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
    // Compute lambda and solve
    double lambda = 0.0;

    // gotta change indexing to accomodate X->Y->Z
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            for (int k = 0; k < NZ/2+1; k++) {
                FD[i * (NY) * (NZ/2+1) + j * (NZ/2+1) + k][0] = FD[i * (NY) * (NZ/2+1) + j * (NZ/2+1) + k][0] /
                    (2 * (std::cos(2 * i * M_PI / (NX)) - 1) +
                    2 * (std::cos(2 * j * M_PI / (NY)) - 1) +
                    2 * (std::cos(2 * k * M_PI / (NZ)) - 1));
                FD[i * (NY) * (NZ/2+1) + j * (NZ/2+1) + k][1] = FD[i * (NY) * (NZ/2+1) + j * (NZ/2+1) + k][1] /
                    (2 * (std::cos(2 * i * M_PI / (NX)) - 1) +
                    2 * (std::cos(2 * j * M_PI / (NY)) - 1) +
                    2 * (std::cos(2 * k * M_PI / (NZ)) - 1));
            }
        }
    }
    FD[0][0] = 0.0;
    FD[0][1] = 0.0;

    // Inverse Fourier transform
    fftw_plan backward = fftw_plan_dft_c2r_3d(NX, NY, NZ, FD, F_dP.data(), FFTW_ESTIMATE);
    fftw_execute(backward);

    // Normalization
    for (size_t i=0; i < NX; i++){
        for (size_t j=0; j < NY; j++){
            for (size_t k=0; k < NZ; k++){
                F_dP[i * (NY) * (NZ) + j * (NZ) + k] = F_dP[i * (NY) * (NZ) + j * (NZ) + k] / (NX*NY*NZ);
            }
        }
    }


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



void PoissonSolver::solveNeumannPoisson(std::array<Real, (NX+1) * (NY+1) * (NZ+1)>& F)
{
    // Divide all the boundaries by sqrt(2.0)

    // LEFT FACE
    for(size_t j=0; j<NY+1; j++){
        for(size_t k=0; k<NZ+1; k++){
            F[j * (NZ+1) + k] /= sqrt(2.0);
        }
    }

    // LOWER FACE
    for(size_t i=0; i<NX+1; i++){
        for(size_t j=0; j<NY+1; j++){
            F[i * (NY+1) * (NZ+1) + j * (NZ+1)] /= sqrt(2.0);
        }
    } 

    // FRONT FACE
    for(size_t i=0; i<NX+1; i++){
        for(size_t k=0; k<NZ+1; k++){
            F[i * (NY+1) * (NZ+1) + k] /= sqrt(2.0);
        }
    } 

    // RIGHT FACE
    for(size_t j=0; j<NY+1; j++){
        for(size_t k=0; k<NZ+1; k++){
            F[NX * (NY+1) * (NZ+1) + j * (NZ+1) + k] /= sqrt(2.0);
        }
    }

    // UPPER FACE
    for(size_t i=0; i<NX+1; i++){
        for(size_t j=0; j<NY+1; j++){
            F[i * (NY+1) * (NZ+1) + j * (NZ+1) + NZ] /= sqrt(2.0);
        }
    } 

    // LOWER FACE
    for(size_t i=0; i<NX+1; i++){
        for(size_t k=0; k<NZ+1; k++){
            F[i * (NY+1) * (NZ+1) + NY * (NZ+1) + k] /= sqrt(2.0);
        }
    } 

    // Perform the fft

    fftw_plan neumann = fftw_plan_r2r_3d(NX+1, NY+1, NZ+1, F.data(), F.data(), FFTW_REDFT00, FFTW_REDFT00, FFTW_REDFT00, FFTW_ESTIMATE);
    fftw_execute(neumann);

    for (size_t i=0; i < NX+1; i++){
        for (size_t j=0; j < NY+1; j++){
            for (size_t k=0; k < NZ+1; k++){

                std::cout << i+1 << ", " << j+1 << ", " << k+1 << ": " << F[i * (NY+1) * (NZ+1) + j * (NZ+1) + k] << std::endl;
            }   
        }
    }

    for (int i = 0; i < NX+1; i++) {
        for (int j = 0; j < NY+1; j++) {
            for (int k = 0; k < NZ+1; k++) {
                F[i * (NY+1) * (NZ+1) + j * (NZ+1) + k] = F[i * (NY+1) * (NZ+1) + j * (NZ+1) + k] /
                    (2 * (std::cos(i * M_PI / (NX)) - 1) +
                    2 * (std::cos(j * M_PI / (NY)) - 1) +
                    2 * (std::cos(k * M_PI / (NZ)) - 1));
            }
        }
    }
    F[0] = 0.0;

    // Normalization
    double normalization_factor = 1.0 / (2.0 * ((NX+1) * (NY+1) * (NZ+1)) - 1.0);
    for(int i = 0; i < NX+1; i++) {
        for (int j = 0; j < NY+1; j++) {
            for (int k = 0; k < NZ+1; k++) {
                    F[i * (NY+1) * (NZ+1) + j * (NZ+1) + k] *= normalization_factor;
            }
        }
    }


    // Multiplicate all the boundaries by sqrt(2.0)

    // LEFT FACE
    for(size_t j=0; j<NY+1; j++){
        for(size_t k=0; k<NZ+1; k++){
            F[j * (NZ+1) + k] *= sqrt(2.0);
        }
    }

    // LOWER FACE
    for(size_t i=0; i<NX+1; i++){
        for(size_t j=0; j<NY+1; j++){
            F[i * (NY+1) * (NZ+1) + j * (NZ+1)] *= sqrt(2.0);
        }
    } 

    // FRONT FACE
    for(size_t i=0; i<NX+1; i++){
        for(size_t k=0; k<NZ+1; k++){
            F[i * (NY+1) * (NZ+1) + k] *= sqrt(2.0);
        }
    } 

    // RIGHT FACE
    for(size_t j=0; j<NY+1; j++){
        for(size_t k=0; k<NZ+1; k++){
            F[NX * (NY+1) * (NZ+1) + j * (NZ+1) + k] *= sqrt(2.0);
        }
    }

    // UPPER FACE
    for(size_t i=0; i<NX+1; i++){
        for(size_t j=0; j<NY+1; j++){
            F[i * (NY+1) * (NZ+1) + j * (NZ+1) + NZ] *= sqrt(2.0);
        }
    } 

    // LOWER FACE
    for(size_t i=0; i<NX+1; i++){
        for(size_t k=0; k<NZ+1; k++){
            F[i * (NY+1) * (NZ+1) + NY * (NZ+1) + k] *= sqrt(2.0);
        }
    } 

    fftw_destroy_plan(neumann);

}