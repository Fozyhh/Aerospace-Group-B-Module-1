#include "core.hpp"
#include "poissonSolver.hpp"


void IcoNS::solve_time_step(Real time)
{
    int ipencil=0;
    PoissonSolver poissonSolver(false,false,false, c2d);

    // 1) pressure point exchange
    double* halo_p;
    c2d->updateHalo(grid.p, halo_p,1,ipencil);
    
    for (int i = 1 + lbx; i < newDimX_x - 1 - rbx; i++)
    {
        for (int j = 1 + lby; j < newDimY_x - 1 - rby; j++)
        {
            for (int k = 1; k < dim_z - 1; k++)
            {
                Y2_x[getx(i, j, k)] = grid.u[getx(i, j, k)] + 
                                                     64.0 / 120.0 * DT * functionF_u(grid.u, grid.v, grid.w, i, j, k, time) -
                                                     64.0 / 120.0 * DT * (halo_p[getHaloP(i+1,j,k)] - 
                                                     halo_p[getHaloP(i,j,k)]) / (DX);
            }
        }
    }


    for (int i = 1 + lbx; i < newDimX_y - 1 - rbx; i++)
    {
        for (int j = 1 + lby; j < newDimY_y - 1 - rby; j++)
        {
            for (int k = 1; k < dim_z - 1; k++)
            {
                Y2_y[gety(i, j, k)] = grid.v[gety(i, j, k)] + 
                                                     64.0 / 120.0 * DT * functionF_v(grid.u, grid.v, grid.w, i, j, k, time) -
                                                     64.0 / 120.0 * DT * (halo_p[getHaloP(i,j+1,k)] - 
                                                     halo_p[getHaloP(i,j,k)]) / (DY);
            }
        }
    }


    for (int i = 1 + lbx; i < newDimX_z - 1 - rbx; i++)
    {
        for (int j = 1 + lby; j < newDimY_z - 1 - rby; j++)
        {
            for (int k = 1; k < dim_z_z - 1; k++)
            {
                Y2_z[getz(i, j, k)] = grid.w[getz(i, j, k)] + 
                                                     64.0 / 120.0 * DT * functionF_w(grid.u, grid.v, grid.w, i, j, k, time) -
                                                     64.0 / 120.0 * DT * (halo_p[getHaloP(i,j,k+1)] - 
                                                     halo_p[getHaloP(i,j,k)]) / (DZ);
            }
        }
    }

    //TODO: boundary handling?
    boundary.update_boundary(Y2_x, Y2_y, Y2_z, time + 64.0 / 120.0 * DT);
    MPI_Barrier(cart_comm);
    exchangeData(Y2_x, newDimX_x, newDimY_x, dim_z, MPI_face_x_x, MPI_face_y_x,0,1);
    exchangeData(Y2_y, newDimX_y, newDimY_y, dim_z, MPI_face_x_y, MPI_face_y_y,1,0);
    exchangeData(Y2_z, newDimX_z, newDimY_z, dim_z_z, MPI_face_x_z, MPI_face_y_z,1,1);
    MPI_Barrier(cart_comm);

    for (int i = 1 + lbx; i < zSize[0] + 1 - rbx; i++)
    {
        for (int j = 1 + lby; j < zSize[1] + 1 - rby; j++)
        {
            for (int k = 1; k < zSize[2] - 1; k++)
            {
                Y2_p[getp(i-1,j-1,k)] = 120.0 / (64.0 * DT) * ((Y2_x[getx(i, j, k)] - Y2_x[getx(i - 1, j, k)]) / (DX) + (Y2_y[gety(i, j, k)] - Y2_y[gety(i, j - 1, k)]) / (DY) + (Y2_z[getz(i, j, k)] - Y2_z[getz(i, j, k - 1)]) / (DZ));
            }
        }
    }

    boundary.divergence(Y2_x, Y2_y, Y2_z, Y2_p, time + 64.0 / 120.0 * DT, 64.0);
    poissonSolver.solveNeumannPoisson(Y2_p);


    // 2) y2_p pressure point exchange
    c2d->deallocXYZ(halo_p);
    c2d->updateHalo(Y2_p, halo_p, 1, ipencil);
    for (int i = 1; i < newDimX_x - 1; i++)
    {
        for (int j = 1; j < newDimY_x - 1; j++)
        {
            for (int k = 0; k < dim_z; k++)
            {
                Y2_x[getx(i, j, k)] = Y2_x[getx(i, j, k)] - 
                                                     64.0 * DT / (120.0) * (halo_p[getHaloP(i + 1, j, k)] - 
                                                     halo_p[getHaloP(i, j, k)]) / (DX);
            }
        }
    }

    for (int i = 1; i < newDimX_y - 1; i++)
    {
        for (int j = 1; j < newDimY_y - 1; j++)
        {
            for (int k = 0; k < dim_z; k++)
            {
                Y2_y[gety(i, j, k)] = Y2_y[gety(i, j, k)] - 
                                                     64.0 * DT / (120.0) * (halo_p[getHaloP(i, j + 1, k)] - 
                                                     halo_p[getHaloP(i, j, k)]) / (DY);
            }
        }
    }

    for (int i = 1; i < newDimX_z - 1; i++)
    {
        for (int j = 1; j < newDimY_z - 1; j++)
        {
            for (int k = 0; k < dim_z_z; k++)
            {
                Y2_z[getz(i, j, k)] = Y2_z[getz(i, j, k)] - 
                                                     64.0 * DT / (120.0) * (halo_p[getHaloP(i, j, k + 1)] - 
                                                     halo_p[getHaloP(i, j, k)]) / (DZ);
            }
        }
    }

    for (int i = 0; i < zSize[0]; i++)
    {
        for (int j = 0; j < zSize[1]; j++)
        {
            for (int k = 0; k < zSize[2]; k++)
            {
                Phi_p[getp(i,j,k)] = Y2_p[getp(i,j,k)] + grid.p[getp(i,j,k)]; // phi^2
            }
        }
    }


    // 3) Phi_p exchange 
    double* halo_phi;
    c2d->updateHalo(Phi_p, halo_phi, 1, ipencil);
    for (int i = 1 + lbx; i < newDimX_x - 1 - rbx; i++)
    {
        for (int j = 1 + lby; j < newDimY_x - 1 - rby; j++)
        {
            for (int k = 1; k < dim_z - 1; k++)

            {
                Y3_x[getx(i, j, k)] = Y2_x[getx(i, j, k)] +
                                                     50.0 / 120.0 * DT * functionF_u(Y2_x, Y2_y, Y2_z, i, j, k, time + 64.0 / 120.0 * DT) -
                                                     34.0 / 120.0 * DT * functionF_u(grid.u, grid.v, grid.w, i, j, k, time) -
                                                     16.0 / 120.0 * DT * (halo_phi[getHaloP(i+1,j,k)] - halo_phi[getHaloP(i,j,k)]) / (DX);
            }
        }
    }

    for (int i = 1 + lbx; i < newDimX_y - 1 - rbx; i++)
    {
        for (int j = 1 + lby; j < newDimY_y - 1 - rby; j++)
        {
            for (int k = 1; k < dim_z-1; k++)
            {
                Y3_y[gety(i, j, k)] = Y2_y[gety(i, j, k)] +
                                                     50.0 / 120.0 * DT * functionF_v(Y2_x, Y2_y, Y2_z, i, j, k, time + 64.0 / 120.0 * DT) -
                                                     34.0 / 120.0 * DT * functionF_v(grid.u, grid.v, grid.w, i, j, k, time) -
                                                     16.0 / 120.0 * DT * (halo_phi[getHaloP(i,j+1,k)] - halo_phi[getHaloP(i,j,k)]) / (DY);
            }
        }
    }

    for (int i = 1 + lbx; i < newDimX_z - 1 - rbx; i++)
    {
        for (int j = 1 + lby; j < newDimY_z - 1 - rby; j++)
        {
            for (int k = 1; k < dim_z_z - 1; k++)

            {
                Y3_z[getz(i, j, k)] = Y2_z[getz(i, j, k)] +
                                                     50.0 / 120.0 * DT * functionF_w(Y2_x, Y2_y, Y2_z, i, j, k, time + 64.0 / 120.0 * DT) -
                                                     34.0 / 120.0 * DT * functionF_w(grid.u, grid.v, grid.w, i, j, k, time) -
                                                     16.0 / 120.0 * DT * (halo_phi[getHaloP(i,j,k+1)] - halo_phi[getHaloP(i,j,k)]) / (DZ);
            }
        }
    }

    boundary.update_boundary(Y3_x, Y3_y, Y3_z, time + 80.0 / 120.0 * DT);
    MPI_Barrier(cart_comm);
    exchangeData(Y3_x, newDimX_x, newDimY_x, dim_z, MPI_face_x_x, MPI_face_y_x,0,1);
    exchangeData(Y3_y, newDimX_y, newDimY_y, dim_z, MPI_face_x_y, MPI_face_y_y,1,0);
    exchangeData(Y3_z, newDimX_z, newDimY_z, dim_z_z, MPI_face_x_z, MPI_face_y_z,1,1);
    MPI_Barrier(cart_comm);

    for (int i = 1 + lbx; i < zSize[0] + 1 - rbx; i++)
    {
        for (int j = 1 + lby; j < zSize[1] + 1 -rby; j++)
        {
            for (int k = 1; k < zSize[2] - 1; k++)
            {
                Y2_p[getp(i-1,j-1,k)] = 120.0 / (16.0 * DT) * ((Y3_x[getx(i, j, k)] - Y3_x[getx(i - 1, j, k)]) / (DX) + (Y3_y[gety(i, j, k)] - Y3_y[gety(i, j - 1, k)]) / (DY) + (Y3_z[getz(i, j, k)] - Y3_z[getz(i, j, k - 1)]) / (DZ));
            }
        }
    }
    boundary.divergence(Y3_x, Y3_y, Y3_z, Y2_p, time + 80.0 / 120.0 * DT, 16.0);

    poissonSolver.solveNeumannPoisson(Y2_p);

    // 3) y2_p exchange
    c2d->deallocXYZ(halo_p);
    c2d->updateHalo(Y2_p, halo_p, 1, ipencil);
    for (int i = 1; i < newDimX_x - 1; i++)
    {
        for (int j = 1; j < newDimY_x - 1; j++)
        {
            for (int k = 0; k < dim_z; k++)
            {
                Y3_x[getx(i, j, k)] = Y3_x[getx(i, j, k)] - 
                                                     16.0 * DT / (120.0) * (halo_p[getHaloP(i + 1, j, k)] - 
                                                     halo_p[getHaloP(i, j, k)]) / (DX);
            }
        }
    }

    for (int i = 1; i < newDimX_y - 1; i++)
    {
        for (int j = 1; j < newDimY_y - 1; j++)
        {
            for (int k = 0; k < dim_z; k++)
            {
                Y3_y[gety(i, j, k)] = Y3_y[gety(i, j, k)] - 
                                                     16.0 * DT / (120.0) * (halo_p[getHaloP(i, j + 1, k)] - 
                                                     halo_p[getHaloP(i, j, k)]) / (DY);
            }
        }
    }

    for (int i = 1; i < newDimX_z - 1; i++)
    {
        for (int j = 1; j < newDimY_z - 1; j++)
        {
            for (int k = 0; k < dim_z_z; k++)
            {
                Y3_z[getz(i, j, k)] = Y3_z[getz(i, j, k)] - 
                                                     16.0 * DT / (120.0) * (halo_p[getHaloP(i, j, k + 1)] - 
                                                     halo_p[getHaloP(i, j, k)]) / (DZ);
            }
        }
    }

    for (int i = 0; i < zSize[0]; i++)
    {
        for (int j = 0; j < zSize[1]; j++)
        {
            for (int k = 0; k < zSize[2]; k++)
            {
                Phi_p[getp(i,j,k)] = Y2_p[getp(i, j, k)] + Phi_p[getp(i,j,k)]; // Phi_p=phi^3
            }
        }
    }

    MPI_Barrier(cart_comm);

    // 4) Phi_p exchange
    c2d->deallocXYZ(halo_phi);
    c2d->updateHalo(Phi_p, halo_phi, 1, ipencil);

    for (int i = 1 + lbx; i < newDimX_x - 1 - rbx; i++)
    {
        for (int j = 1 + lby; j < newDimY_x - 1 - rby; j++)
        {
            for (int k = 1; k < dim_z-1; k++)

            {
                grid.u[getx(i, j, k)] = Y3_x[getx(i, j, k)] +
                                                       90.0 / 120.0 * DT * functionF_u(Y3_x, Y3_y, Y3_z, i, j, k, time + 80.0 / 120.0 * DT) -
                                                       50.0 / 120.0 * DT * functionF_u(Y2_x, Y2_y, Y2_z, i, j, k, time + 64.0 / 120.0 * DT) -
                                                       40.0 / 120.0 * DT * (halo_phi[getHaloP(i + 1,j,k)] - halo_phi[getHaloP(i,j,k)]) / (DX);
            }
        }
    }

    for (int i = 1 + lbx; i < newDimX_y - 1 - rbx; i++)
    {
        for (int j = 1 + lby; j < newDimY_y - 1 - rby; j++)
        {
            for (int k = 1; k < dim_z-1; k++)
            {
                grid.v[gety(i, j, k)] = Y3_y[gety(i, j, k)] +
                                                       90.0 / 120.0 * DT * functionF_v(Y3_x, Y3_y, Y3_z, i, j, k, time + 80.0 / 120.0 * DT) -
                                                       50.0 / 120.0 * DT * functionF_v(Y2_x, Y2_y, Y2_z, i, j, k, time + 64.0 / 120.0 * DT) -
                                                       40.0 / 120.0 * DT * (halo_phi[getHaloP(i,j+1,k)] - halo_phi[getHaloP(i,j,k)]) / (DY);
            }
        }
    }

    for (int i = 1 + lbx; i < newDimX_z - 1 - rbx; i++)
    {
        for (int j = 1 + lby; j < newDimY_z - 1 - rby; j++)
        {
            for (int k = 1; k < dim_z_z - 1; k++)
            {
                grid.w[getz(i, j, k)] = Y3_z[getz(i, j, k)] +
                                                             90.0 / 120.0 * DT * functionF_w(Y3_x, Y3_y, Y3_z, i, j, k, time + 80.0 / 120.0 * DT) -
                                                             50.0 / 120.0 * DT * functionF_w(Y2_x, Y2_y, Y2_z, i, j, k, time + 64.0 / 120.0 * DT) -
                                                             40.0 / 120.0 * DT * (halo_phi[getHaloP(i,j,k+1)] - halo_phi[getHaloP(i,j,k)]) / (DZ);
           }
        }
    }

    boundary.update_boundary(grid.u, grid.v, grid.w, time + DT);
    MPI_Barrier(cart_comm);
    exchangeData(grid.u, newDimX_x, newDimY_x,dim_z,MPI_face_x_x,MPI_face_y_x,0,1);
    exchangeData(grid.v, newDimX_y, newDimY_y,dim_z,MPI_face_x_y,MPI_face_y_y,1,0);
    exchangeData(grid.w, newDimX_z, newDimY_z,dim_z_z,MPI_face_x_z,MPI_face_y_z,1,1);
    MPI_Barrier(cart_comm);

    for (int i = 1 + lbx; i < zSize[0] + 1 -rbx; i++)
    {
        for (int j = 1 + lby; j < zSize[1] + 1 - rby; j++)
        {
            for (int k = 1; k < zSize[2] - 1; k++)
            {
                Y2_p[getp(i-1,j-1,k)] = 120.0 / (40.0 * DT) * ((grid.u[getx(i, j, k)] - grid.u[getx(i - 1, j, k)]) / (DX) + (grid.v[gety(i, j, k)] - grid.v[gety(i, j - 1, k)]) / (DY) + (grid.w[getz(i, j, k)] - grid.w[getz(i, j, k - 1)]) / (DZ));
            }
        }
    }
    boundary.divergence(grid.u, grid.v, grid.w, Y2_p, time + DT, 40.0);

    poissonSolver.solveNeumannPoisson(Y2_p);

    // 5) y2_p exchange
    c2d->deallocXYZ(halo_p);
    c2d->updateHalo(Y2_p, halo_p, 1, ipencil);

    for(int i = 1; i < newDimX_x - 1; i++)
    {
        for(int j = 1; j < newDimY_x - 1; j++)
        {
            for(int k = 0; k < dim_z; k++)
            {
                grid.u[getx(i,j,k)] -= 40.0 * DT / (120.0) * (halo_p[getHaloP(i + 1, j, k)] - 
                                                      halo_p[getHaloP(i, j, k)]) / (DX); 
            }
        }
    }
    for(int i = 1; i < newDimX_y - 1; i++)
    {
        for(int j = 1; j < newDimY_y - 1; j++)
        {
            for(int k = 0; k < dim_z; k++)
            {
                grid.v[gety(i,j,k)] -= 40.0 * DT / (120.0) * (halo_p[getHaloP(i, j + 1, k)] - 
                                                      halo_p[getHaloP(i, j, k)]) / (DY); 
            }
        }
    }
    for(int i = 1; i < newDimX_z - 1; i++)
    {
        for(int j = 1; j < newDimY_z - 1; j++)
        {
            for(int k = 0; k < dim_z_z; k++)
            {
                grid.w[getz(i,j,k)] -= 40.0 * DT / (120.0) * (halo_p[getHaloP(i, j, k + 1)] - 
                                                      halo_p[getHaloP(i, j, k)]) / (DZ); 
            }
        }
    }
    for(int i = 0; i < zSize[0]; i++)
    {
        for(int j = 0; j < zSize[1]; j++)
        {
            for(int k = 0; k < zSize[2]; k++)
            {
                grid.p[getp(i,j,k)] = Y2_p[getp(i,j, k)]  + Phi_p[getp(i,j,k)]; 
            }
        }
    }
    c2d->deallocXYZ(halo_phi);
    c2d->deallocXYZ(halo_p);
}

//TODO: maybe change everything with the indexing function
Real IcoNS::functionF_u(const std::vector<Real> &u, const std::vector<Real> &v, const std::vector<Real> &w, int i, int j, int k, Real t)
{
    int lu = i * newDimY_x * dim_z + j * dim_z + k;
    int lv = i * newDimY_y * dim_z + j * dim_z + k;
    int lw = i * newDimY_z * dim_z_z + j * dim_z_z + k;

    Real value =-(u[lu] * (u[lu + newDimY_x * dim_z] - u[lu - newDimY_x * dim_z]) / (2.0 * DX) +
             (v[lv] + v[lv + newDimY_y * dim_z] + v[lv - dim_z] + v[lv + newDimY_y * dim_z - dim_z]) / 4.0 * (u[lu + dim_z] - u[lu - dim_z]) / (2.0 * DY) +
             (w[lw] + w[lw + newDimY_z * dim_z_z] + w[lw - 1] + w[lw + newDimY_z * dim_z_z - 1]) / 4.0 * (u[lu + 1] - u[lu - 1]) / (2.0 * DZ)) +
           1 / RE * ((u[lu + newDimY_x * dim_z] - 2 * u[lu] + u[lu - newDimY_x * dim_z]) / (DX * DX) + (u[lu + dim_z] - 2 * u[lu] + u[lu - dim_z]) / (DY * DY) + (u[lu + 1] - 2 * u[lu] + u[lu - 1]) / (DZ * DZ));

    if(testCase==0){
        value += functionG_u(i-1 + coords[0] * other_dim_x_x, j-1 + coords[1] * other_dim_y_x, k, t);
    }
    return value;
    // return -(u[lu] * (u[lu + newDimY_x * dim_z] - u[lu - newDimY_x * dim_z]) / (2.0 * DX) +
    //          (v[lv] + v[lv + newDimY_y * dim_z] + v[lv - dim_z] + v[lv + newDimY_y * dim_z - dim_z]) / 4.0 * (u[lu + dim_z] - u[lu - dim_z]) / (2.0 * DY) +
    //          (w[lw] + w[lw + newDimY_z * dim_z_z] + w[lw - 1] + w[lw + newDimY_z * dim_z_z - 1]) / 4.0 * (u[lu + 1] - u[lu - 1]) / (2.0 * DZ)) +
    //        1 / RE * ((u[lu + newDimY_x * dim_z] - 2 * u[lu] + u[lu - newDimY_x * dim_z]) / (DX * DX) + (u[lu + dim_z] - 2 * u[lu] + u[lu - dim_z]) / (DY * DY) + (u[lu + 1] - 2 * u[lu] + u[lu - 1]) / (DZ * DZ)) +
    //        functionG_u(i-1 + coords[0] * other_dim_x_x, j-1 + coords[1] * other_dim_y_x, k, t);
}

Real IcoNS::functionF_v(const std::vector<Real> &u, const std::vector<Real> &v, const std::vector<Real> &w, int i, int j, int k, Real t)
{
    int lu = i * newDimY_x * dim_z + j * dim_z + k;
    int lv = i * newDimY_y * dim_z + j * dim_z + k;
    int lw = i * newDimY_z * dim_z_z + j * dim_z_z + k;
    
    Real value = -((u[lu] + u[lu + dim_z] + u[lu - newDimY_x * dim_z] + u[lu - newDimY_x * dim_z + dim_z]) / 4.0 * ((v[lv + newDimY_y * dim_z] - v[lv - newDimY_y * dim_z]) / (2.0 * DX)) +
             v[lv] * (v[lv + dim_z] - v[lv - dim_z]) / (2.0 * DY) +
             (w[lw] + w[lw - 1] + w[lw + dim_z_z] + w[lw + dim_z_z - 1]) / 4.0 * (v[lv + 1] - v[lv - 1]) / (2.0 * DZ)) +
           (1.0 / RE) * ((v[lv + newDimY_y * dim_z] - 2.0 * v[lv] + v[lv - newDimY_y * dim_z]) / (DX * DX) +
                         (v[lv + dim_z] - 2.0 * v[lv] + v[lv - dim_z]) / (DY * DY) +
                         (v[lv + 1] - 2.0 * v[lv] + v[lv - 1]) / (DZ * DZ));
    if(testCase==0){
        value += functionG_v(i-1 + coords[0] * other_dim_x_y, j-1 + coords[1] * other_dim_y_y, k, t);
    }

    return value;
    // return -((u[lu] + u[lu + dim_z] + u[lu - newDimY_x * dim_z] + u[lu - newDimY_x * dim_z + dim_z]) / 4.0 * ((v[lv + newDimY_y * dim_z] - v[lv - newDimY_y * dim_z]) / (2.0 * DX)) +
    //          v[lv] * (v[lv + dim_z] - v[lv - dim_z]) / (2.0 * DY) +
    //          (w[lw] + w[lw - 1] + w[lw + dim_z_z] + w[lw + dim_z_z - 1]) / 4.0 * (v[lv + 1] - v[lv - 1]) / (2.0 * DZ)) +
    //        (1.0 / RE) * ((v[lv + newDimY_y * dim_z] - 2.0 * v[lv] + v[lv - newDimY_y * dim_z]) / (DX * DX) +
    //                      (v[lv + dim_z] - 2.0 * v[lv] + v[lv - dim_z]) / (DY * DY) +
    //                      (v[lv + 1] - 2.0 * v[lv] + v[lv - 1]) / (DZ * DZ)) +
    //        functionG_v(i-1 + coords[0] * other_dim_x_y, j-1 + coords[1] * other_dim_y_y, k, t);
}

Real IcoNS::functionF_w(const std::vector<Real> &u, const std::vector<Real> &v, const std::vector<Real> &w, int i, int j, int k, Real t)

{
    int lu = i * newDimY_x * dim_z + j * dim_z + k;
    int lv = i * newDimY_y * dim_z + j * dim_z + k;
    int lw = i * newDimY_z * dim_z_z + j * dim_z_z + k;

    Real value = -((u[lu] + u[lu - newDimY_x * dim_z] + u[lu + 1] + u[lu - dim_z * newDimY_x + 1]) / 4.0 * (w[lw + newDimY_z * dim_z_z] - w[lw - newDimY_z * dim_z_z]) / (2.0 * DX) +
                  (v[lv + 1] + v[lv - dim_z + 1] + v[lv] + v[lv - dim_z]) / 4.0 * (w[lw + dim_z_z] - w[lw - dim_z_z]) / (2.0 * DY) +
                  w[lw] * (w[lw + 1] - w[lw - 1]) / (2.0 * DZ)) +
                  (1.0 / RE) * ((w[lw + newDimY_z * dim_z_z] - 2.0 * w[lw] + w[lw - newDimY_z * dim_z_z]) / (DX * DX) +
                         (w[lw + dim_z_z] - 2.0 * w[lw] + w[lw - dim_z_z]) / (DY * DY) +
                         (w[lw + 1] - 2.0 * w[lw] + w[lw - 1]) / (DZ * DZ));
    
    if(testCase==0){
        value += functionG_w(i-1 + coords[0] * other_dim_x_z, j-1 + coords[1] * other_dim_y_z, k, t);
    }
    return 0;
    // return -((u[lu] + u[lu - newDimY_x * dim_z] + u[lu + 1] + u[lu - dim_z * newDimY_x + 1]) / 4.0 * (w[lw + newDimY_z * dim_z_z] - w[lw - newDimY_z * dim_z_z]) / (2.0 * DX) +
    //          (v[lv + 1] + v[lv - dim_z + 1] + v[lv] + v[lv - dim_z]) / 4.0 * (w[lw + dim_z_z] - w[lw - dim_z_z]) / (2.0 * DY) +
    //          w[lw] * (w[lw + 1] - w[lw - 1]) / (2.0 * DZ)) +
    //        (1.0 / RE) * ((w[lw + newDimY_z * dim_z_z] - 2.0 * w[lw] + w[lw - newDimY_z * dim_z_z]) / (DX * DX) +
    //                      (w[lw + dim_z_z] - 2.0 * w[lw] + w[lw - dim_z_z]) / (DY * DY) +
    //                      (w[lw + 1] - 2.0 * w[lw] + w[lw - 1]) / (DZ * DZ)) +
    //        functionG_w(i-1 + coords[0] * other_dim_x_z, j-1 + coords[1] * other_dim_y_z, k, t);

}

Real IcoNS::functionG_u(int i, int j, int k, Real t)
{
    Real x = i * DX + DX / 2;
    Real y = j * DY;
    Real z = k * DZ;

    return std::sin(x) * std::cos(y) * std::sin(z) * std::cos(t) +
           std::sin(x) * std::cos(x) * std::cos(y) * std::cos(y) * std::sin(z) * std::sin(z) * std::sin(t) * std::sin(t) -
           std::sin(x) * std::cos(x) * std::sin(y) * std::sin(y) * std::sin(z) * std::sin(z) * std::sin(t) * std::sin(t) +
           2 * std::sin(x) * std::cos(x) * std::cos(y) * std::cos(y) * std::cos(z) * std::cos(z) * std::sin(t) * std::sin(t) +
           3.0 / RE * std::sin(x) * std::cos(y) * std::sin(z) * std::sin(t) - std::sin(x) * std::cos(y) * std::cos(z) * std::sin(t);
}

Real IcoNS::functionG_v(int i, int j, int k, Real t)
{
    Real x = i * DX;
    Real y = j * DY + DY / 2;
    Real z = k * DZ;

    return std::cos(x) * std::sin(y) * std::sin(z) * std::cos(t) -
           std::sin(x) * std::sin(x) * std::sin(y) * std::cos(y) * std::sin(z) * std::sin(z) * std::sin(t) * std::sin(t) +
           std::cos(x) * std::cos(x) * std::sin(y) * std::cos(y) * std::sin(z) * std::sin(z) * std::sin(t) * std::sin(t) +
           2.0 * std::cos(x) * std::cos(x) * std::sin(y) * std::cos(y) * std::cos(z) * std::cos(z) * std::sin(t) * std::sin(t) +
           3.0 / RE * std::cos(x) * std::sin(y) * std::sin(z) * std::sin(t) - std::cos(x) * std::sin(y) * std::cos(z) * std::sin(t);
}

Real IcoNS::functionG_w(int i, int j, int k, Real t)
{
    Real x = i * DX;
    Real y = j * DY;
    Real z = k * DZ + DZ / 2;

    return 2 * std::cos(x) * std::cos(y) * std::cos(z) * std::cos(t) -
           2 * std::sin(x) * std::sin(x) * std::cos(y) * std::cos(y) * std::sin(z) * std::cos(z) * std::sin(t) * std::sin(t) -
           2 * std::cos(x) * std::cos(x) * std::sin(y) * std::sin(y) * std::sin(z) * std::cos(z) * std::sin(t) * std::sin(t) -
           4.0 * std::cos(x) * std::cos(x) * std::cos(y) * std::cos(y) * std::sin(z) * std::cos(z) * std::sin(t) * std::sin(t) +
           6.0 / RE * std::cos(x) * std::cos(y) * std::cos(z) * std::sin(t) - std::cos(x) * std::cos(y) * std::sin(z) * std::sin(t);
}
