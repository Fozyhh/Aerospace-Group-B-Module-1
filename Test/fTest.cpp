#include <iostream>
#include <array>

struct cell{
    float Xu;
    float Yv;
    float Zw;
    float p;
};

int main() {
    //create grid
    const int I = 10, J = 10, K = 10;
    std::array<std::array<std::array<cell, I>, J>, K> grid;
    //fill grid
    for (int k = 0; k < K; k++){
        for (int j = 0; j < J; j++){
            for (int i = 0; i < I; i++){
                grid[k][j][i].Xu = static_cast<float>(i+1) / I;
                grid[k][j][i].Yv = static_cast<float>(j+1) / J;
                grid[k][j][i].Zw = static_cast<float>(k+1) / K;
                grid[k][j][i].p = 0;
            }
        }
    }

    //set parameters, create variables
    const float Re = 1e3, dx = 0.01, dy = 0.01, dz = 0.01;
    std::array<float, 3> f;
    float Xv, Xw, Yu, Yw, Zu, Zv;
    float udx, udy, udz, vdx, vdy, vdz, wdx, wdy, wdz;
    float uddx, uddy, uddz, vddx, vddy, vddz, wddx, wddy, wddz;

    // !!! SKIPPING BOUNDARY !!!
    for(int k = 1; k < K-1; k++){
        for(int j = 1; j < J-1; j++){
            for (int i = 1; i < I-1; i++){        
                // calc field values
                Xv = (grid[k][j][i].Yv + grid[k][j][i+1].Yv + grid[k][j-1][i].Yv + grid[k][j-1][i+1].Yv)/4;
                Xw = (grid[k][j][i].Zw + grid[k][j][i+1].Zw + grid[k-1][j][i].Zw + grid[k-1][j][i+1].Zw)/4;
                Yu = (grid[k][j][i].Xu + grid[k][j+1][i].Xu + grid[k][j][i-1].Xu + grid[k][j+1][i-1].Xu)/4;
                Yw = (grid[k][j][i].Zw + grid[k][j+1][i].Zw + grid[k-1][j][i].Zw + grid[k-1][j+1][i].Zw)/4;
                Zu = (grid[k][j][i].Xu + grid[k+1][j][i].Xu + grid[k][j][i-1].Xu + grid[k+1][j][i-1].Xu)/4;
                Zv = (grid[k][j][i].Yv + grid[k+1][j][i].Yv + grid[k][j-1][i].Yv + grid[k+1][j-1][i].Yv)/4;
                // calc first derivatives
                udx = (grid[k][j][i+1].Xu - grid[k][j][i-1].Xu)/(2*dx);
                udy = (grid[k][j+1][i].Xu - grid[k][j-1][i].Xu)/(2*dy);
                udz = (grid[k+1][j][i].Xu - grid[k-1][j][i].Xu)/(2*dz);
                vdx = (grid[k][j][i+1].Yv - grid[k][j][i-1].Yv)/(2*dx);
                vdy = (grid[k][j+1][i].Yv - grid[k][j-1][i].Yv)/(2*dy);
                vdz = (grid[k+1][j][i].Yv - grid[k-1][j][i].Yv)/(2*dz);
                wdx = (grid[k][j][i+1].Zw - grid[k][j][i-1].Zw)/(2*dx);
                wdy = (grid[k][j+1][i].Zw - grid[k][j-1][i].Zw)/(2*dy);
                wdz = (grid[k+1][j][i].Zw - grid[k-1][j][i].Zw)/(2*dz);
                // calc second derivatives
                uddx = (grid[k][j][i+1].Xu - 2*grid[k][j][i].Xu + grid[k][j][i-1].Xu)/(dx*dx);
                uddy = (grid[k][j+1][i].Xu - 2*grid[k][j][i].Xu + grid[k][j-1][i].Xu)/(dy*dy);
                uddz = (grid[k+1][j][i].Xu - 2*grid[k][j][i].Xu + grid[k-1][j][i].Xu)/(dz*dz);
                vddx = (grid[k][j][i+1].Yv - 2*grid[k][j][i].Yv + grid[k][j][i-1].Yv)/(dx*dx);
                vddy = (grid[k][j+1][i].Yv - 2*grid[k][j][i].Yv + grid[k][j-1][i].Yv)/(dy*dy);
                vddz = (grid[k+1][j][i].Yv - 2*grid[k][j][i].Yv + grid[k-1][j][i].Yv)/(dz*dz);
                wddx = (grid[k][j][i+1].Zw - 2*grid[k][j][i].Zw + grid[k][j][i-1].Zw)/(dx*dx);
                wddy = (grid[k][j+1][i].Zw - 2*grid[k][j][i].Zw + grid[k][j-1][i].Zw)/(dy*dy);
                wddz = (grid[k+1][j][i].Zw - 2*grid[k][j][i].Zw + grid[k-1][j][i].Zw)/(dz*dz);
                // calc f(), convection + diffusion only
                f[0] = - grid[k][j][i].Xu*udx - Xv*udy - Xw*udz + 1/Re*(uddx + uddy + uddz);
                f[1] = - Yu*vdx - grid[k][j][i].Yv*vdy - Yw*vdz + 1/Re*(vddx + vddy + vddz);
                f[2] = - Zu*wdx - Zv*wdy - grid[k][j][i].Zw*wdz + 1/Re*(wddx + wddy + wddz);

                std::cout << "Components of f found at step " << (I*J*k+I*j+i) << " : [" << f[0] << ", " << f[1] << ", " << f[2] << "]" << std::endl;
            }
        }
    }
    return 0;
}