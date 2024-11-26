#include "poissonSolver.hpp"

void PoissonSolver::solvePoisson(std::array<Real, (NX+1)*(NY+1)*(NZ+1)>& F_dP, fftw_complex *FD)
{
    
    constexpr double pi = std::numbers::pi_v<double>;
    
    fftw_complex *in, *out;
    bool periodicBC[3] = {periodicX, periodicY, periodicZ};
    C2Decomp *c2d;
    c2d = new C2Decomp(NX+1, NY+1, NZ+1, 0, 0, periodicBC);
    
    // dP = fft(fft(fft(F))) TBD
    fftw_plan fftw_plan_r2c_3d(int nx, int ny, int nz, Real *in, fftw_complex *out, unsigned flags);
    fftw_execute(fftw_plan_r2c_3d(NX+1, NY+1, NZ+1, F_dP.data(), FD, FFTW_MEASURE));

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
    for (int i = 1; i < NX+1; i++) {
        for (int j = 0; j < NY+1; j++) {
            for (int k = 0; k < NZ+1; k++) {
                FD[i*(NY+1)*(NZ+1)+j*(NZ+1)+k][0] = FD[i*(NY+1)*(NZ+1)+j*(NZ+1)+k][0] /
                    (2 * (cos(2 * i * pi / NX+1) - 1) +
                    2 * (cos(2 * j * pi / NY+1) - 1) +
                    2 * (cos(2 * k * pi / NZ+1) - 1));
                FD[i*(NY+1)*(NZ+1)+j*(NZ+1)+k][1] = FD[i*(NY+1)*(NZ+1)+j*(NZ+1)+k][1] /
                    (2 * (cos(2 * i * pi / NX+1) - 1) +
                    2 * (cos(2 * j * pi / NY+1) - 1) +
                    2 * (cos(2 * k * pi / NZ+1) - 1));
            }
        }
    }

    // Inverse Fourier transform
    fftw_plan fftw_plan_c2r_3d(int nx, int ny, int nz, fftw_complex *in, Real *out, unsigned flags);
    fftw_execute(fftw_plan_c2r_3d(NX+1, NY+1, NZ+1, FD, F_dP.data(), FFTW_MEASURE));
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
