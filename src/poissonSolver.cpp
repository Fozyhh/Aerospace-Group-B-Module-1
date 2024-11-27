#include "poissonSolver.hpp"

void PoissonSolver::solvePoisson(std::array<Real, NX * NY * NZ>& F_dP, fftw_complex *FD)
{
    
    constexpr double pi = std::numbers::pi_v<double>;

    bool periodicBC[3] = {periodicX, periodicY, periodicZ};
    //C2Decomp *c2d;
    //c2d = new C2Decomp(NX+1, NY+1, NZ+1, 0, 0, periodicBC);
    
    // dP = fft(fft(fft(F))) TBD
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
    FD[0][0] = 0.0;
    FD[0][1] = 0.0;
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            for (int k = 0; k < NZ; k++) {
                if(i!=0 || j!=0 || k!=0){
                    FD[i * (NY) * (NZ) + j * (NZ) + k][0] = FD[i * (NY) * (NZ) + j * (NZ) + k][0] /
                        (2 * (std::cos(2 * i * pi / (NX)) - 1) +
                        2 * (std::cos(2 * j * pi / (NY)) - 1) +
                        2 * (std::cos(2 * k * pi / (NZ)) - 1));
                    FD[i * (NY) * (NZ) + j * (NZ) + k][1] = FD[i * (NY) * (NZ) + j * (NZ) + k][1] /
                        (2 * (std::cos(2 * i * pi / (NX)) - 1) +
                        2 * (std::cos(2 * j * pi / (NY)) - 1) +
                        2 * (std::cos(2 * k * pi / (NZ)) - 1));
                }
            }
        }
    }

    // Inverse Fourier transform
    fftw_plan backward = fftw_plan_dft_c2r_3d(NX, NY, NZ, FD, F_dP.data(), FFTW_ESTIMATE);
    fftw_execute(backward);
/*  // return ifft(ifft(ifft(b))) TBD
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
