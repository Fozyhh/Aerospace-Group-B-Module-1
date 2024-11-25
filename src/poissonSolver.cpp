#include "poissonSolver.hpp"
#include <cmath>
#include <numbers>
#include <complex>

void PoissonSolver::solvePoisson(std::array<Real, (NX+1)*(NY+1)*(NZ+1)>& F, std::array<Real, (NX+1)*(NY+1)*(NZ+1)>& dP)
{
    
    constexpr double pi = std::numbers::pi_v<double>;
        
    // dP = fft(fft(fft(F))) TBD
    if(periodicX){
        dP = fft(F);
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

    // Compute lambda and solve
    double lambda = 0.0;
    // gotta change indexing to accomodate X->Y->Z
    dP[0] = 0.0;
    for (int i = 1; i < NX+1; i++) {
        for (int j = 0; j < NY+1; j++) {
            for (int k = 0; k < NZ+1; k++) {
                b[i*(NY+1)*(NZ+1)+j*(NZ+1)+k] = b[i*(NY+1)*(NZ+1)+j*(NZ+1)+k] /
                    (2 * (cos(2 * i * pi / NX+1) - 1) +
                    2 * (cos(2 * j * pi / NY+1) - 1) +
                    2 * (cos(2 * k * pi / NZ+1) - 1));
            }
        }
    }

    // Inverse Fourier transform
    // return ifft(ifft(ifft(b))) TBD
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
}
