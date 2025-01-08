#include "core.hpp"
#include "poissonSolver.hpp"
#include <fstream>

void IcoNS::solve_time_step(Real time)
{
    int ipencil=0;
    PoissonSolver poissonSolver(false,false,false, c2d);

    // 1) pressure point exchange
    double* halo_p;
    c2d->updateHalo(grid.p, halo_p,1,ipencil);
    
    for (int i = 1; i < newDimX_x - 1; i++)
    {
        for (int j = 1; j < newDimY_x - 1; j++)
        {
            for (int k = 0; k < dim_z ; k++)
            {
                Y2_x[getx(i, j, k)] = grid.u[getx(i, j, k)] + 
                                                     64.0 / 120.0 * DT * functionF_u(grid.u, grid.v, grid.w, i, j, k, time) -
                                                     64.0 / 120.0 * DT * (halo_p[getHaloP(i+1,j,k)] - 
                                                     halo_p[getHaloP(i,j,k)]) / (DX);
            }
        }
    }


    for (int i = 1 ; i < newDimX_y - 1; i++)
    {
        for (int j = 1 ; j < newDimY_y - 1; j++)
        {
            for (int k = 0; k < dim_z ; k++)
            {
                Y2_y[gety(i, j, k)] = grid.v[gety(i, j, k)] + 
                                                     64.0 / 120.0 * DT * functionF_v(grid.u, grid.v, grid.w, i, j, k, time) -
                                                     64.0 / 120.0 * DT * (halo_p[getHaloP(i,j+1,k)] - 
                                                     halo_p[getHaloP(i,j,k)]) / (DY);
            }
        }
    }


    for (int i = 1 ; i < newDimX_z - 1 ; i++)
    {
        for (int j = 1 ; j < newDimY_z - 1 ; j++)
        {
            for (int k = 0; k < dim_z_z ; k++)
            {
                Y2_z[getz(i, j, k)] = grid.w[getz(i, j, k)] + 
                                                     64.0 / 120.0 * DT * functionF_w(grid.u, grid.v, grid.w, i, j, k, time) -
                                                     64.0 / 120.0 * DT * (halo_p[getHaloP(i,j,k+1)] - 
                                                     halo_p[getHaloP(i,j,k)]) / (DZ);
            }
        }
    }

    //TODO: boundary handling?
    // boundary.update_boundary(Y2_x, Y2_y, Y2_z, time + 64.0 / 120.0 * DT);
    // MPI_Barrier(cart_comm);
    exchangeData(Y2_x, newDimX_x, newDimY_x, dim_z, MPI_face_x_x, MPI_face_y_x,0,1);
    exchangeData(Y2_y, newDimX_y, newDimY_y, dim_z, MPI_face_x_y, MPI_face_y_y,1,0);
    exchangeData(Y2_z, newDimX_z, newDimY_z, dim_z_z, MPI_face_x_z, MPI_face_y_z,1,1);

    for (int i = 1 ; i < zSize[0] + 1 ; i++)
    {
        for (int j = 1; j < zSize[1] + 1 ; j++)
        {
            for (int k = 0; k < zSize[2]; k++)
            {
                // Y2_p[getp(i-1,j-1,k)] = 120.0 / (64.0 * DT) * ((Y2_x[getx(i, j, k)] - Y2_x[getx(i - 1, j, k)]) / (DX) + (Y2_y[gety(i, j, k)] - Y2_y[gety(i, j - 1, k)]) / (DY) + (Y2_z[getz(i, j, k)] - Y2_z[getz(i, j, k - 1)]) / (DZ));
                // Y2_p[i * (NY+1) * (NZ+1) + j * (NZ+1) + k] = 120.0 / (64.0 * DT) * ((Y2_x[getx(i, j, k)] - Y2_x[getx(i - 1, j, k)]) / (DX) + (Y2_y[gety(i, j, k)] - Y2_y[gety(i, j - 1, k)]) / (DY) + (Y2_z[getz(i, j, k)] - Y2_z[getz(i, j, k - 1)]) / (DZ));
                Real diffx = 0.0;
                Real diffy = 0.0;
                Real diffz = 0.0;
 
                if (i == 1)
                {
                    // continue;
                    diffx = (-8 * boundary.boundary_value_u[0]->value(i - 1 - 0.5, j - 1, k, time + 64.0 / 120.0 * DT) + 9 * Y2_x[getx(i, j, k)] - Y2_x[getx(i + 1, j, k)]) / (3 * DX);
                    // diffx =  64.0 * DT / (120.0) * std::cos(i * DX) * std::cos(j * DY) * std::cos(k * DZ) * (sin(time) - sin(time + 64.0 * DT / (120.0)));
                }
                else if (i == zSize[0])
                {
                    // continue;
                    diffx = (8 * boundary.boundary_value_u[0]->value(NX - 0.5, j -1, k, time + 64.0 / 120.0 * DT) - 9 * Y2_x[getx(i - 1, j, k)] + Y2_x[getx(i - 2, j, k)]) / (3 * DX);
                    // diffx =  64.0 * DT / (120.0) * std::cos(i * DX) * std::cos(j * DY) * std::cos(k * DZ) * (sin(time) - sin(time + 64.0 * DT / (120.0)));
                }
                else
                {
                    diffx = (Y2_x[getx(i, j, k)] - Y2_x[getx(i - 1, j, k)]) / (DX);
                    // diffx =  64.0 * DT / (120.0) * std::cos(i * DX) * std::cos(j * DY) * std::cos(k * DZ) * (sin(time) - sin(time + 64.0 * DT / (120.0)));
                }
 
                if (j == 1)
                {
                    // continue;
                    diffy = (-8 * boundary.boundary_value_v[0]->value(i - 1, j - 1.5, k, time + 64.0 / 120.0 * DT) + 9 * Y2_y[gety(i, j, k)] - Y2_y[gety(i, j + 1, k)]) / (3 * DY);
                    // diffy =  64.0 * DT / (120.0) * std::cos(i * DX) * std::cos(j * DY) * std::cos(k * DZ) * (sin(time) - sin(time + 64.0 * DT / (120.0)));
                }
                else if (j == zSize[1])
                {
                    // continue;
                    diffy = (8 * boundary.boundary_value_v[0]->value(i - 1, j - 1.5, k, time + 64.0 / 120.0 * DT) - 9 * Y2_y[gety(i, j - 1, k)] + Y2_y[gety(i, j - 2, k)]) / (3 * DY);
                    // diffy =  64.0 * DT / (120.0) * std::cos(i * DX) * std::cos(j * DY) * std::cos(k * DZ) * (sin(time) - sin(time + 64.0 * DT / (120.0)));
                }
                else
                {
                    diffy = (Y2_y[gety(i, j, k)] - Y2_y[gety(i, j - 1, k)]) / (DY);
                    // diffy =  64.0 * DT / (120.0) * std::cos(i * DX) * std::cos(j * DY) * std::cos(k * DZ) * (sin(time) - sin(time + 64.0 * DT / (120.0)));
                }
 
                if (k == 0)
                {
                    // continue;
                    diffz = (-8 * boundary.boundary_value_w[0]->value(i -1, j - 1, k - 0.5, time + 64.0 / 120.0 * DT) + 9 * Y2_z[getz(i, j, k)] - Y2_z[getz(i, j, k + 1)]) / (3 * DZ);
                    // diffz =  64.0 * DT / (120.0) * std::cos(i * DX) * std::cos(j * DY) * std::cos(k * DZ) * (sin(time) - sin(time + 64.0 * DT / (120.0)));
                }
                else if (k == zSize[2] - 1)
                {
                    // continue;
                    diffz = (8 * boundary.boundary_value_w[0]->value(i - 1, j -1, k - 0.5, time + 64.0 / 120.0 * DT) - 9 * Y2_z[getz(i, j, k - 1)] + Y2_z[getz(i, j, k - 2)]) / (3 * DZ);
                    // diffz =  64.0 * DT / (120.0) * std::cos(i * DX) * std::cos(j * DY) * std::cos(k * DZ) * (sin(time) - sin(time + 64.0 * DT / (120.0)));
                }
                else
                {
                    diffz = (Y2_z[getz(i, j, k)] - Y2_z[getz(i, j, k - 1)]) / (DZ);
                    // diffz =  64.0 * DT / (120.0) * std::cos(i * DX) * std::cos(j * DY) * std::cos(k * DZ) * (sin(time) - sin(time + 64.0 * DT / (120.0)));
                }
 
                Y2_p[getp(i-1,j-1,k)] = 120.0 / (64.0 * DT) * (diffx + diffy + diffz);
            }
        }
    }

    //boundary.divergence(Y2_x, Y2_y, Y2_z, Y2_p, time + 64.0 / 120.0 * DT, 64.0);

    /* for (int i = 0; i < zSize[0]; i++)
    {
        for (int j = 0; j < zSize[1]; j++)
        {
            for (int k = 0; k < zSize[2]; k++)
            {
                Y2_p[getp(i, j, k)] = 3*cos(i*DX)*cos(j*DY)*cos(k*DZ)*(sin(time)-sin(time+64.0/120.0*DT));   
            }
        }
    }
 */
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
    pressionCorrection(Phi_p);
    boundary.update_boundary(Y2_x, Y2_y, Y2_z, time + 64.0 / 120.0 * DT);
    MPI_Barrier(cart_comm);
    exchangeData(Y2_x, newDimX_x, newDimY_x, dim_z, MPI_face_x_x, MPI_face_y_x, 0, 1);
    exchangeData(Y2_y, newDimX_y, newDimY_y, dim_z, MPI_face_x_y, MPI_face_y_y, 1, 0);
    exchangeData(Y2_z, newDimX_z, newDimY_z, dim_z_z, MPI_face_x_z, MPI_face_y_z, 1, 1);

    for (int i = 1; i < newDimX_x - 1 ; i++)
    {
        for (int j = 1 ; j < newDimY_x - 1; j++)
        {
            for (int k = 0; k < dim_z ; k++)

            {
                Y3_x[getx(i, j, k)] = Y2_x[getx(i, j, k)] +
                                                     50.0 / 120.0 * DT * functionF_u(Y2_x, Y2_y, Y2_z, i, j, k, time + 64.0 / 120.0 * DT) -
                                                     34.0 / 120.0 * DT * functionF_u(grid.u, grid.v, grid.w, i, j, k, time) -
                                                     16.0 / 120.0 * DT * (halo_phi[getHaloP(i+1,j,k)] - halo_phi[getHaloP(i,j,k)]) / (DX);
            }
        }
    }

    for (int i = 1; i < newDimX_y - 1 ; i++)
    {
        for (int j = 1 ; j < newDimY_y - 1; j++)
        {
            for (int k = 0; k < dim_z; k++)
            {
                Y3_y[gety(i, j, k)] = Y2_y[gety(i, j, k)] +
                                                     50.0 / 120.0 * DT * functionF_v(Y2_x, Y2_y, Y2_z, i, j, k, time + 64.0 / 120.0 * DT) -
                                                     34.0 / 120.0 * DT * functionF_v(grid.u, grid.v, grid.w, i, j, k, time) -
                                                     16.0 / 120.0 * DT * (halo_phi[getHaloP(i,j+1,k)] - halo_phi[getHaloP(i,j,k)]) / (DY);
            }
        }
    }

    for (int i = 1; i < newDimX_z - 1 ; i++)
    {
        for (int j = 1 ; j < newDimY_z - 1; j++)
        {
            for (int k = 0; k < dim_z_z ; k++)

            {
                Y3_z[getz(i, j, k)] = Y2_z[getz(i, j, k)] +
                                                     50.0 / 120.0 * DT * functionF_w(Y2_x, Y2_y, Y2_z, i, j, k, time + 64.0 / 120.0 * DT) -
                                                     34.0 / 120.0 * DT * functionF_w(grid.u, grid.v, grid.w, i, j, k, time) -
                                                     16.0 / 120.0 * DT * (halo_phi[getHaloP(i,j,k+1)] - halo_phi[getHaloP(i,j,k)]) / (DZ);
            }
        }
    }

    //boundary.update_boundary(Y3_x, Y3_y, Y3_z, time + 80.0 / 120.0 * DT);
    MPI_Barrier(cart_comm);
    exchangeData(Y3_x, newDimX_x, newDimY_x, dim_z, MPI_face_x_x, MPI_face_y_x,0,1);
    exchangeData(Y3_y, newDimX_y, newDimY_y, dim_z, MPI_face_x_y, MPI_face_y_y,1,0);
    exchangeData(Y3_z, newDimX_z, newDimY_z, dim_z_z, MPI_face_x_z, MPI_face_y_z,1,1);

    for (int i = 1 ; i < zSize[0] + 1 ; i++)
    {
        for (int j = 1 ; j < zSize[1] + 1; j++)
        {
            for (int k = 0; k < zSize[2]; k++)
            {
                //Y2_p[getp(i-1,j-1,k)] = 120.0 / (16.0 * DT) * ((Y3_x[getx(i, j, k)] - Y3_x[getx(i - 1, j, k)]) / (DX) + (Y3_y[gety(i, j, k)] - Y3_y[gety(i, j - 1, k)]) / (DY) + (Y3_z[getz(i, j, k)] - Y3_z[getz(i, j, k - 1)]) / (DZ));
                Real diffx = 0.0;
                Real diffy = 0.0;
                Real diffz = 0.0;
 
                if (i == 1)
                {
                    // continue;
                    diffx = (-8 * boundary.boundary_value_u[0]->value(i - 1 - 0.5, j - 1, k, time + 64.0 / 120.0 * DT) + 9 * Y3_x[getx(i, j, k)] - Y3_x[getx(i + 1, j, k)]) / (3 * DX);
                    // diffx =  64.0 * DT / (120.0) * std::cos(i * DX) * std::cos(j * DY) * std::cos(k * DZ) * (sin(time) - sin(time + 64.0 * DT / (120.0)));
                }
                else if (i == zSize[0])
                {
                    // continue;
                    diffx = (8 * boundary.boundary_value_u[0]->value(NX - 0.5, j -1, k, time + 64.0 / 120.0 * DT) - 9 * Y3_x[getx(i - 1, j, k)] + Y3_x[getx(i - 2, j, k)]) / (3 * DX);
                    // diffx =  64.0 * DT / (120.0) * std::cos(i * DX) * std::cos(j * DY) * std::cos(k * DZ) * (sin(time) - sin(time + 64.0 * DT / (120.0)));
                }
                else
                {
                    diffx = (Y3_x[getx(i, j, k)] - Y3_x[getx(i - 1, j, k)]) / (DX);
                    // diffx =  64.0 * DT / (120.0) * std::cos(i * DX) * std::cos(j * DY) * std::cos(k * DZ) * (sin(time) - sin(time + 64.0 * DT / (120.0)));
                }
 
                if (j == 1)
                {
                    // continue;
                    diffy = (-8 * boundary.boundary_value_v[0]->value(i - 1, j - 1.5, k, time + 64.0 / 120.0 * DT) + 9 * Y3_y[gety(i, j, k)] - Y3_y[gety(i, j + 1, k)]) / (3 * DY);
                    // diffy =  64.0 * DT / (120.0) * std::cos(i * DX) * std::cos(j * DY) * std::cos(k * DZ) * (sin(time) - sin(time + 64.0 * DT / (120.0)));
                }
                else if (j == zSize[1])
                {
                    // continue;
                    diffy = (8 * boundary.boundary_value_v[0]->value(i - 1, j - 1.5, k, time + 64.0 / 120.0 * DT) - 9 * Y3_y[gety(i, j - 1, k)] + Y3_y[gety(i, j - 2, k)]) / (3 * DY);
                    // diffy =  64.0 * DT / (120.0) * std::cos(i * DX) * std::cos(j * DY) * std::cos(k * DZ) * (sin(time) - sin(time + 64.0 * DT / (120.0)));
                }
                else
                {
                    diffy = (Y3_y[gety(i, j, k)] - Y3_y[gety(i, j - 1, k)]) / (DY);
                    // diffy =  64.0 * DT / (120.0) * std::cos(i * DX) * std::cos(j * DY) * std::cos(k * DZ) * (sin(time) - sin(time + 64.0 * DT / (120.0)));
                }
 
                if (k == 0)
                {
                    // continue;
                    diffz = (-8 * boundary.boundary_value_w[0]->value(i -1, j - 1, k - 0.5, time + 64.0 / 120.0 * DT) + 9 * Y3_z[getz(i, j, k)] - Y3_z[getz(i, j, k + 1)]) / (3 * DZ);
                    // diffz =  64.0 * DT / (120.0) * std::cos(i * DX) * std::cos(j * DY) * std::cos(k * DZ) * (sin(time) - sin(time + 64.0 * DT / (120.0)));
                }
                else if (k == zSize[2] - 1)
                {
                    // continue;
                    diffz = (8 * boundary.boundary_value_w[0]->value(i - 1, j -1, k - 0.5, time + 64.0 / 120.0 * DT) - 9 * Y3_z[getz(i, j, k - 1)] + Y3_z[getz(i, j, k - 2)]) / (3 * DZ);
                    // diffz =  64.0 * DT / (120.0) * std::cos(i * DX) * std::cos(j * DY) * std::cos(k * DZ) * (sin(time) - sin(time + 64.0 * DT / (120.0)));
                }
                else
                {
                    diffz = (Y3_z[getz(i, j, k)] - Y3_z[getz(i, j, k - 1)]) / (DZ);
                    // diffz =  64.0 * DT / (120.0) * std::cos(i * DX) * std::cos(j * DY) * std::cos(k * DZ) * (sin(time) - sin(time + 64.0 * DT / (120.0)));
                }
 
                Y2_p[getp(i-1,j-1,k)] = 120.0 / (16.0 * DT) * (diffx + diffy + diffz);
            }
        }
    }
    //boundary.divergence(Y3_x, Y3_y, Y3_z, Y2_p, time + 80.0 / 120.0 * DT, 16.0);

    /* for (int i = 0; i < zSize[0]; i++)
    {
        for (int j = 0; j < zSize[1]; j++)
        {
            for (int k = 0; k < zSize[2]; k++)
            {
                Y2_p[getp(i, j, k)] = 3*cos(i*DX)*cos(j*DY)*cos(k*DZ)*(sin(time+64.0/120.0*DT)-sin(time+80.0/120.0*DT));   
            }
        }
    } */

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

    boundary.update_boundary(Y3_x, Y3_y, Y3_z, time + 80.0 / 120.0 * DT);
    MPI_Barrier(cart_comm);
    exchangeData(Y3_x, newDimX_x, newDimY_x, dim_z, MPI_face_x_x, MPI_face_y_x, 0, 1);
    exchangeData(Y3_y, newDimX_y, newDimY_y, dim_z, MPI_face_x_y, MPI_face_y_y, 1, 0);
    exchangeData(Y3_z, newDimX_z, newDimY_z, dim_z_z, MPI_face_x_z, MPI_face_y_z, 1, 1);
    pressionCorrection(Phi_p);


    // 4) Phi_p exchange
    c2d->deallocXYZ(halo_phi);
    c2d->updateHalo(Phi_p, halo_phi, 1, ipencil);

    for (int i = 1; i < newDimX_x - 1 ; i++)
    {
        for (int j = 1 ; j < newDimY_x - 1 ; j++)
        {
            for (int k = 0; k < dim_z; k++)

            {
                grid.u[getx(i, j, k)] = Y3_x[getx(i, j, k)] +
                                                       90.0 / 120.0 * DT * functionF_u(Y3_x, Y3_y, Y3_z, i, j, k, time + 80.0 / 120.0 * DT) -
                                                       50.0 / 120.0 * DT * functionF_u(Y2_x, Y2_y, Y2_z, i, j, k, time + 64.0 / 120.0 * DT) -
                                                       40.0 / 120.0 * DT * (halo_phi[getHaloP(i + 1,j,k)] - halo_phi[getHaloP(i,j,k)]) / (DX);
            }
        }
    }

    for (int i = 1 ; i < newDimX_y - 1; i++)
    {
        for (int j = 1 ; j < newDimY_y - 1; j++)
        {
            for (int k = 0; k < dim_z; k++)
            {
                grid.v[gety(i, j, k)] = Y3_y[gety(i, j, k)] +
                                                       90.0 / 120.0 * DT * functionF_v(Y3_x, Y3_y, Y3_z, i, j, k, time + 80.0 / 120.0 * DT) -
                                                       50.0 / 120.0 * DT * functionF_v(Y2_x, Y2_y, Y2_z, i, j, k, time + 64.0 / 120.0 * DT) -
                                                       40.0 / 120.0 * DT * (halo_phi[getHaloP(i,j+1,k)] - halo_phi[getHaloP(i,j,k)]) / (DY);
            }
        }
    }

    for (int i = 1 ; i < newDimX_z - 1; i++)
    {
        for (int j = 1 ; j < newDimY_z - 1; j++)
        {
            for (int k = 0; k < dim_z_z ; k++)
            {
                grid.w[getz(i, j, k)] = Y3_z[getz(i, j, k)] +
                                                             90.0 / 120.0 * DT * functionF_w(Y3_x, Y3_y, Y3_z, i, j, k, time + 80.0 / 120.0 * DT) -
                                                             50.0 / 120.0 * DT * functionF_w(Y2_x, Y2_y, Y2_z, i, j, k, time + 64.0 / 120.0 * DT) -
                                                             40.0 / 120.0 * DT * (halo_phi[getHaloP(i,j,k+1)] - halo_phi[getHaloP(i,j,k)]) / (DZ);
           }
        }
    }

    //boundary.update_boundary(grid.u, grid.v, grid.w, time + DT);
    exchangeData(grid.u, newDimX_x, newDimY_x,dim_z,MPI_face_x_x,MPI_face_y_x,0,1);
    exchangeData(grid.v, newDimX_y, newDimY_y,dim_z,MPI_face_x_y,MPI_face_y_y,1,0);
    exchangeData(grid.w, newDimX_z, newDimY_z,dim_z_z,MPI_face_x_z,MPI_face_y_z,1,1);

    /* for (int i = 1 + lbx; i < zSize[0] + 1 -rbx; i++)
    {
        for (int j = 1 + lby; j < zSize[1] + 1 - rby; j++)
        {
            for (int k = 1; k < zSize[2]; k++)
            {
                Y2_p[getp(i-1,j-1,k)] = 0.0;
            }
        }
    } */

    for (int i = 1 ; i < zSize[0] + 1 ; i++)
    {
        for (int j = 1 ; j < zSize[1] + 1 ; j++)
        {
            for (int k = 0; k < zSize[2] ; k++)
            {

                //Y2_p[getp(i-1,j-1,k)] = 120.0 / (40.0 * DT) * ((grid.u[getx(i, j, k)] - grid.u[getx(i - 1, j, k)]) / (DX) + (grid.v[gety(i, j, k)] - grid.v[gety(i, j - 1, k)]) / (DY) + (grid.w[getz(i, j, k)] - grid.w[getz(i, j, k - 1)]) / (DZ));
                Real diffx = 0.0;
                Real diffy = 0.0;
                Real diffz = 0.0;
 
                if (i == 1)
                {
                    // continue;
                    diffx = (-8 * boundary.boundary_value_u[0]->value(i - 1 - 0.5, j - 1, k, time + 64.0 / 120.0 * DT) + 9 * grid.u[getx(i, j, k)] - grid.u[getx(i + 1, j, k)]) / (3 * DX);
                    // diffx =  64.0 * DT / (120.0) * std::cos(i * DX) * std::cos(j * DY) * std::cos(k * DZ) * (sin(time) - sin(time + 64.0 * DT / (120.0)));
                }
                else if (i == zSize[0])
                {
                    // continue;
                    diffx = (8 * boundary.boundary_value_u[0]->value(NX - 0.5, j -1, k, time + 64.0 / 120.0 * DT) - 9 * grid.u[getx(i - 1, j, k)] + grid.u[getx(i - 2, j, k)]) / (3 * DX);
                    // diffx =  64.0 * DT / (120.0) * std::cos(i * DX) * std::cos(j * DY) * std::cos(k * DZ) * (sin(time) - sin(time + 64.0 * DT / (120.0)));
                }
                else
                {
                    diffx = (grid.u[getx(i, j, k)] - grid.u[getx(i - 1, j, k)]) / (DX);
                    // diffx =  64.0 * DT / (120.0) * std::cos(i * DX) * std::cos(j * DY) * std::cos(k * DZ) * (sin(time) - sin(time + 64.0 * DT / (120.0)));
                }
 
                if (j == 1)
                {
                    // continue;
                    diffy = (-8 * boundary.boundary_value_v[0]->value(i - 1, j - 1.5, k, time + 64.0 / 120.0 * DT) + 9 * grid.v[gety(i, j, k)] - grid.v[gety(i, j + 1, k)]) / (3 * DY);
                    // diffy =  64.0 * DT / (120.0) * std::cos(i * DX) * std::cos(j * DY) * std::cos(k * DZ) * (sin(time) - sin(time + 64.0 * DT / (120.0)));
                }
                else if (j == zSize[1])
                {
                    // continue;
                    diffy = (8 * boundary.boundary_value_v[0]->value(i - 1, j - 1.5, k, time + 64.0 / 120.0 * DT) - 9 * grid.v[gety(i, j - 1, k)] + grid.v[gety(i, j - 2, k)]) / (3 * DY);
                    // diffy =  64.0 * DT / (120.0) * std::cos(i * DX) * std::cos(j * DY) * std::cos(k * DZ) * (sin(time) - sin(time + 64.0 * DT / (120.0)));
                }
                else
                {
                    diffy = (grid.v[gety(i, j, k)] - grid.v[gety(i, j - 1, k)]) / (DY);
                    // diffy =  64.0 * DT / (120.0) * std::cos(i * DX) * std::cos(j * DY) * std::cos(k * DZ) * (sin(time) - sin(time + 64.0 * DT / (120.0)));
                }
 
                if (k == 0)
                {
                    // continue;
                    diffz = (-8 * boundary.boundary_value_w[0]->value(i -1, j - 1, k - 0.5, time + 64.0 / 120.0 * DT) + 9 * grid.w[getz(i, j, k)] - grid.w[getz(i, j, k + 1)]) / (3 * DZ);
                    // diffz =  64.0 * DT / (120.0) * std::cos(i * DX) * std::cos(j * DY) * std::cos(k * DZ) * (sin(time) - sin(time + 64.0 * DT / (120.0)));
                }
                else if (k == zSize[2] - 1)
                {
                    // continue;
                    diffz = (8 * boundary.boundary_value_w[0]->value(i - 1, j -1, k - 0.5, time + 64.0 / 120.0 * DT) - 9 * grid.w[getz(i, j, k - 1)] + grid.w[getz(i, j, k - 2)]) / (3 * DZ);
                    // diffz =  64.0 * DT / (120.0) * std::cos(i * DX) * std::cos(j * DY) * std::cos(k * DZ) * (sin(time) - sin(time + 64.0 * DT / (120.0)));
                }
                else
                {
                    diffz = (grid.w[getz(i, j, k)] - grid.w[getz(i, j, k - 1)]) / (DZ);
                    // diffz =  64.0 * DT / (120.0) * std::cos(i * DX) * std::cos(j * DY) * std::cos(k * DZ) * (sin(time) - sin(time + 64.0 * DT / (120.0)));
                }
 
                Y2_p[getp(i-1,j-1,k)] = 120.0 / (40.0 * DT) * (diffx + diffy + diffz);
            }
        }
    }
    //boundary.divergence(grid.u, grid.v, grid.w, Y2_p, time + DT, 40.0);

    //vtk file for grid.p
    std::ofstream vtkFile("pressure"+std::to_string(time)+".vtk");
    //header
    vtkFile << "# vtk DataFile Version 3.0\n";
    vtkFile << "Navier-Stokes solution\n";
    vtkFile << "ASCII\n";
    vtkFile << "DATASET STRUCTURED_POINTS\n";
    vtkFile << "DIMENSIONS " <<  zSize[0] << " " << zSize[1] << " " << zSize[2] << "\n";
    vtkFile << "ORIGIN 0 0 0\n";
    vtkFile << "SPACING " << DX << " " << DY << " " << DZ << "\n";
    //scalar data   
    vtkFile << "POINT_DATA " << zSize[0]*zSize[1]*zSize[2] << "\n";
    vtkFile << "SCALARS pressure double 1\n";
    vtkFile << "LOOKUP_TABLE default\n";
    for (int i = 0; i < zSize[0]; i++)
    {
        for (int j = 0; j < zSize[1]; j++)
        {
            for (int k = 0; k < zSize[2]; k++)
            {
                vtkFile << Y2_p[getp(i,j,k)] << "\n";
            }
        }
    }
    vtkFile.close();

    /* for (int i = 0; i < zSize[0]; i++)
    {
        for (int j = 0; j < zSize[1]; j++)
        {
            for (int k = 0; k < zSize[2]; k++)
            {
                Y2_p[getp(i, j, k)] = 3*cos(i*DX)*cos(j*DY)*cos(k*DZ)*(sin(time+80.0/120.0*DT)-sin(time+DT));
            }
        }
    } */


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
    boundary.update_boundary(grid.u, grid.v, grid.w, time + DT);
    MPI_Barrier(cart_comm);
    exchangeData(grid.u, newDimX_x, newDimY_x, dim_z, MPI_face_x_x, MPI_face_y_x, 0, 1);
    exchangeData(grid.v, newDimX_y, newDimY_y, dim_z, MPI_face_x_y, MPI_face_y_y, 1, 0);
    exchangeData(grid.w, newDimX_z, newDimY_z, dim_z_z, MPI_face_x_z, MPI_face_y_z, 1, 1);
    pressionCorrection(grid.p);
    c2d->deallocXYZ(halo_phi);
    c2d->deallocXYZ(halo_p);

    

    /* //vtk file for grid.u
    std::ofstream vtkFile2("grid_u"+std::to_string(time)+".vtk");
    //header
    vtkFile2 << "# vtk DataFile Version 3.0\n";
    vtkFile2 << "Navier-Stokes solution\n";
    vtkFile2 << "ASCII\n";
    vtkFile2 << "DATASET STRUCTURED_POINTS\n";
    vtkFile2 << "DIMENSIONS " <<  newDimX_x << " " << newDimY_x << " " << dim_z << "\n";
    vtkFile2 << "ORIGIN 0 0 0\n";
    vtkFile2 << "SPACING " << DX << " " << DY << " " << DZ << "\n";
    //scalar data
    vtkFile2 << "POINT_DATA " << newDimX_x*newDimY_x*dim_z << "\n";
    vtkFile2 << "SCALARS u double 1\n";
    vtkFile2 << "LOOKUP_TABLE default\n";
    for (int i = 0; i < newDimX_x; i++)
    {
        for (int j = 0; j < newDimY_x; j++)
        {
            for (int k = 0; k < dim_z; k++)
            {
                vtkFile2 << grid.u[getx(i,j,k)] << "\n";
            }
        }
    }
    vtkFile2.close();

    //vtk file for grid.v
    std::ofstream vtkFile3("grid_v"+std::to_string(time)+".vtk");
    //header
    vtkFile3 << "# vtk DataFile Version 3.0\n";
    vtkFile3 << "Navier-Stokes solution\n";
    vtkFile3 << "ASCII\n";
    vtkFile3 << "DATASET STRUCTURED_POINTS\n";
    vtkFile3 << "DIMENSIONS " <<  newDimX_y << " " << newDimY_y << " " << dim_z << "\n";
    vtkFile3 << "ORIGIN 0 0 0\n";
    vtkFile3 << "SPACING " << DX << " " << DY << " " << DZ << "\n";
    //scalar data
    vtkFile3 << "POINT_DATA " << newDimX_y*newDimY_y*dim_z << "\n";
    vtkFile3 << "SCALARS v double 1\n";
    vtkFile3 << "LOOKUP_TABLE default\n";
    for (int i = 0; i < newDimX_y; i++)
    {
        for (int j = 0; j < newDimY_y; j++)
        {
            for (int k = 0; k < dim_z; k++)
            {
                vtkFile3 << grid.v[gety(i,j,k)] << "\n";
            }
        }
    }
    vtkFile3.close();

    //vtk file for grid.w
    std::ofstream vtkFile4("grid_w"+std::to_string(time)+".vtk");
    //header
    vtkFile4 << "# vtk DataFile Version 3.0\n";
    vtkFile4 << "Navier-Stokes solution\n";
    vtkFile4 << "ASCII\n";
    vtkFile4 << "DATASET STRUCTURED_POINTS\n";
    vtkFile4 << "DIMENSIONS " <<  newDimX_z << " " << newDimY_z << " " << dim_z_z << "\n";
    vtkFile4 << "ORIGIN 0 0 0\n";
    vtkFile4 << "SPACING " << DX << " " << DY << " " << DZ << "\n";
    //scalar data
    vtkFile4 << "POINT_DATA " << newDimX_z*newDimY_z*dim_z_z << "\n";
    vtkFile4 << "SCALARS w double 1\n";
    vtkFile4 << "LOOKUP_TABLE default\n";
    for (int i = 0; i < newDimX_z; i++)
    {
        for (int j = 0; j < newDimY_z; j++)
        {
            for (int k = 0; k < dim_z_z; k++)
            {
                vtkFile4 << grid.w[getz(i,j,k)] << "\n";
            }
        }
    }
    vtkFile4.close(); */
}

//TODO: maybe change everything with the indexing function
Real IcoNS::functionF_u(const std::vector<Real> &u, const std::vector<Real> &v, const std::vector<Real> &w, int i, int j, int k, Real t)
{
    /*int lu = i * newDimY_x * dim_z + j * dim_z + k;
    int lv = i * newDimY_y * dim_z + j * dim_z + k;
    int lw = i * newDimY_z * dim_z_z + j * dim_z_z + k;

    Real value =-(u[lu] * (u[lu + newDimY_x * dim_z] - u[lu - newDimY_x * dim_z]) / (2.0 * DX) +
             (v[lv] + v[lv + newDimY_y * dim_z] + v[lv - dim_z] + v[lv + newDimY_y * dim_z - dim_z]) / 4.0 * (u[lu + dim_z] - u[lu - dim_z]) / (2.0 * DY) +
             (w[lw] + w[lw + newDimY_z * dim_z_z] + w[lw - 1] + w[lw + newDimY_z * dim_z_z - 1]) / 4.0 * (u[lu + 1] - u[lu - 1]) / (2.0 * DZ)) +
           1 / RE * ((u[lu + newDimY_x * dim_z] - 2 * u[lu] + u[lu - newDimY_x * dim_z]) / (DX * DX) + (u[lu + dim_z] - 2 * u[lu] + u[lu - dim_z]) / (DY * DY) + (u[lu + 1] - 2 * u[lu] + u[lu - 1]) / (DZ * DZ));

    if(testCase==0){
        value += functionG_u(i-1 + coords[0] * other_dim_x_x, j-1 + coords[1] * other_dim_y_x, k, t);
    }
    return value;  */ 
    Real diffx, diffy, diffz, u_value, v_value, w_value, diffsecx, diffsecy, diffsecz;

    u_value = u[getx(i, j, k)];
    //TODO: global offsets
    if (j == 1 || j == newDimY_x - 1)
    {
        v_value = boundary.boundary_value_v[0]->value(i - 1, j - 1.5, k, t);
    }
    else
    {
        v_value = (v[gety(i, j, k)] + v[gety(i + 1, j, k)] + v[gety(i, j - 1, k)] + v[gety(i + 1, j - 1, k)]) / 4.0;
    }
 
    if (k == 0 || k == NZ)
    {
        w_value = boundary.boundary_value_w[0]->value(i - 1, j - 1, k - 0.5, t);
    }
    else
    {
        w_value = (w[getz(i, j, k)] + w[getz(i + 1, j, k)] + w[getz(i, j, k - 1)] + w[getz(i + 1, j, k - 1)]) / 4.0;
    }
 
    if (i == 1)
    {
        diffx = (-4 * boundary.boundary_value_u[0]->value(i - 1 - 0.5, j - 1, k, t) + 3 * u[getx(i, j, k)] + u[getx(i + 1, j, k)]) / (3.0 * DX);
    }
    else if (i == newDimX_x - 1)
    {
        diffx = (4 * boundary.boundary_value_u[0]->value(NX - 0.5, j - 1, k, t) - 3 * u[getx(i, j, k)] - u[getx(i - 1, j, k)]) / (3.0 * DX);
    }
    else
    {
        diffx = (u[getx(i + 1, j, k)] - u[getx(i - 1, j, k)]) / (2.0 * DX);
    }
 
    if (j == 1)
    {
        diffy = (-3 * u[getx(i, j, k)] + 4 * u[getx(i, j + 1, k)] - u[getx(i, j + 2, k)]) / (2 * DY);
    }
    else if (j == newDimY_x -1)
    {
        diffy = (3 * u[getx(i, j, k)] - 4 * u[getx(i, j - 1, k)] + u[getx(i, j - 2, k)]) / (2 * DY);
    }
    else
    {
        diffy = (u[getx(i, j + 1, k)] - u[getx(i, j - 1, k)]) / (2.0 * DY);
    }
 
    if (k == 0)
    {
        diffz = (-3 * u[getx(i, j, k)] + 4 * u[getx(i, j, k + 1)] - u[getx(i, j, k + 2)]) / (2 * DZ);
    }
    else if (k == NZ)
    {
        diffz = (3 * u[getx(i, j, k)] - 4 * u[getx(i, j, k - 1)] + u[getx(i, j, k - 2)]) / (2 * DZ);
    }
    else
    {
        diffz = (u[getx(i, j, k + 1)] - u[getx(i, j, k - 1)]) / (2.0 * DZ);
    }
 
    if (i == 1)
    {
        diffsecx = (16 * boundary.boundary_value_u[0]->value(i - 1.5, j-1, k, t) - 25 * u[getx(i, j, k)] + 10 * u[getx(i + 1, j, k)] - u[getx(i + 2, j, k)]) / (5.0 * DX * DX);
    }
    else if (i == newDimX_x - 1)
    {
        diffsecx = (-16 * boundary.boundary_value_u[0]->value(i + 1.5, j-1, k, t) + 25 * u[getx(i, j, k)] - 10 * u[getx(i - 1, j, k)] + u[getx(i - 2, j, k)]) / (5.0 * DX * DX);
    }
    else
    {
        diffsecx = (u[getx(i + 1, j, k)] - 2 * u[getx(i, j, k)] + u[getx(i - 1, j, k)]) / (DX * DX);
    }
 
    if (j == 1)
    {
        diffsecy = (2 * u[getx(i, j, k)] - 5 * u[getx(i, j + 1, k)] + 4 * u[getx(i, j + 2, k)] - u[getx(i, j + 3, k)]) / (DY * DY);
    }
    else if (j == newDimY_x - 1)
    {
        diffsecy = (-2 * u[getx(i, j, k)] + 5 * u[getx(i, j - 1, k)] - 4 * u[getx(i, j - 2, k)] + u[getx(i, j - 3, k)]) / (DY * DY);
    }
    else
    {
        diffsecy = (u[getx(i, j + 1, k)] - 2 * u[getx(i, j, k)] + u[getx(i, j - 1, k)]) / (DY * DY);
    }
 
    if (k == 0)
    {
        diffsecz = (2 * u[getx(i, j, k)] - 5 * u[getx(i, j, k + 1)] + 4 * u[getx(i, j, k + 2)] - u[getx(i, j, k + 3)]) / (DZ * DZ);
    }
    else if (k == NZ)
    {
        diffsecz = (-2 * u[getx(i, j, k)] + 5 * u[getx(i, j, k - 1)] - 4 * u[getx(i, j, k - 2)] + u[getx(i, j, k - 3)]) / (DZ * DZ);
    }
    else
    {
        diffsecz = (u[getx(i, j, k + 1)] - 2 * u[getx(i, j, k)] + u[getx(i, j, k - 1)]) / (DZ * DZ);
    }
 
    return -1 * (u_value * diffx + v_value * diffy + w_value * diffz) + 1 / RE * (diffsecx + diffsecy + diffsecz) + functionG_u(i-1 + coords[0] * other_dim_x_x, j-1 + coords[1] * other_dim_y_x, k, t);
}


Real IcoNS::functionF_v(const std::vector<Real> &u, const std::vector<Real> &v, const std::vector<Real> &w, int i, int j, int k, Real t)
{
    // int lu = i * newDimY_x * dim_z + j * dim_z + k;
    // int lv = i * newDimY_y * dim_z + j * dim_z + k;
    // int lw = i * newDimY_z * dim_z_z + j * dim_z_z + k;
    
    // Real value = -((u[lu] + u[lu + dim_z] + u[lu - newDimY_x * dim_z] + u[lu - newDimY_x * dim_z + dim_z]) / 4.0 * ((v[lv + newDimY_y * dim_z] - v[lv - newDimY_y * dim_z]) / (2.0 * DX)) +
    //          v[lv] * (v[lv + dim_z] - v[lv - dim_z]) / (2.0 * DY) +
    //          (w[lw] + w[lw - 1] + w[lw + dim_z_z] + w[lw + dim_z_z - 1]) / 4.0 * (v[lv + 1] - v[lv - 1]) / (2.0 * DZ)) +
    //        (1.0 / RE) * ((v[lv + newDimY_y * dim_z] - 2.0 * v[lv] + v[lv - newDimY_y * dim_z]) / (DX * DX) +
    //                      (v[lv + dim_z] - 2.0 * v[lv] + v[lv - dim_z]) / (DY * DY) +
    //                      (v[lv + 1] - 2.0 * v[lv] + v[lv - 1]) / (DZ * DZ));
    // if(testCase==0){
    //     value += functionG_v(i-1 + coords[0] * other_dim_x_y, j-1 + coords[1] * other_dim_y_y, k, t);
    // }
 
    // return value;
    Real diffx, diffy, diffz, u_value, v_value, w_value, diffsecx, diffsecy, diffsecz;
    v_value = v[gety(i, j, k)];
 
    if (i == 1 || i == newDimX_y - 1)
    {
        u_value = boundary.boundary_value_u[0]->value(i - 1.5, j - 1, k, t);
    }
    else
    {
        u_value = (u[getx(i, j, k)] + u[getx(i, j + 1, k)] + u[getx(i - 1, j, k)] + u[getx(i - 1, j + 1, k)]) / 4.0;
    }
 
    if (k == 0 || k == NZ)
    {
        w_value = boundary.boundary_value_w[0]->value(i -1, j -1, k - 0.5, t);
    }
    else
    {
        w_value = (w[getz(i, j, k)] + w[getz(i, j + 1, k)] + w[getz(i, j, k - 1)] + w[getz(i, j + 1, k - 1)]) / 4.0;
    }
 
    if (i == 1)
    {
        diffx = (-3 * v[gety(i, j, k)] + 4 * v[gety(i + 1, j, k)] - v[gety(i + 2, j, k)]) / (2 * DX);
    }
    else if (i == newDimX_y - 1)
    {
        diffx = (3 * v[gety(i, j, k)] - 4 * v[gety(i - 1, j, k)] + v[gety(i - 2, j, k)]) / (2 * DX);
    }
    else
    {
        diffx = (v[gety(i + 1, j, k)] - v[gety(i - 1, j, k)]) / (2.0 * DX);
    }
 
    if (j == 1)
    {
        diffy = (-4 * boundary.boundary_value_v[0]->value(i -1, j - 1.5, k, t) + 3 * v[gety(i, j, k)] + v[gety(i, j + 1, k)]) / (3.0 * DY);
    }
    else if (j == newDimY_y - 1)
    {
        diffy = (4 * boundary.boundary_value_v[0]->value(i -1, NY - 0.5, k, t) - 3 * v[gety(i, j, k)] - v[gety(i, j - 1, k)]) / (3.0 * DY);
    }
    else
    {
        diffy = (v[gety(i, j + 1, k)] - v[gety(i, j - 1, k)]) / (2.0 * DY);
    }
 
    if (k == 0)
    {
        diffz = (-3 * v[gety(i, j, k)] + 4 * v[gety(i, j, k + 1)] - v[gety(i, j, k + 2)]) / (2 * DZ);
    }
    else if (k == NZ)
    {
        diffz = (3 * v[gety(i, j, k)] - 4 * v[gety(i, j, k - 1)] + v[gety(i, j, k - 2)]) / (2 * DZ);
    }
    else
    {
        diffz = (v[gety(i, j, k + 1)] - v[gety(i, j, k - 1)]) / (2.0 * DZ);
    }
 
    if (j == 1)
    {
        diffsecy = (16 * boundary.boundary_value_v[0]->value(i -1, j - 1.5, k, t) - 25 * v[gety(i, j, k)] + 10 * v[gety(i, j + 1, k)] - v[gety(i, j + 2, k)]) / (5.0 * DY * DY);
    }
    else if (j == newDimY_y - 1)
    {
        diffsecy = (-16 * boundary.boundary_value_v[0]->value(i -1, NY - 0.5, k, t) + 25 * v[gety(i, j, k)] - 10 * v[gety(i, j - 1, k)] + v[gety(i, j - 2, k)]) / (5.0 * DY * DY);
    }
    else
    {
        diffsecy = (v[gety(i, j + 1, k)] - 2 * v[gety(i, j, k)] + v[gety(i, j - 1, k)]) / (DY * DY);
    }
 
    if (i == 0)
    {
        diffsecx = (2 * v[gety(i, j, k)] - 5 * v[gety(i + 1, j, k)] + 4 * v[gety(i + 2, j, k)] - v[gety(i + 3, j, k)]) / (DX * DX);
    }
    else if (i == newDimX_y - 1)
    {
        diffsecx = (-2 * v[gety(i, j, k)] + 5 * v[gety(i - 1, j, k)] - 4 * v[gety(i - 2, j, k)] + v[gety(i - 3, j, k)]) / (DX * DX);
    }
    else
    {
        diffsecx = (v[gety(i + 1, j, k)] - 2 * v[gety(i, j, k)] + v[gety(i - 1, j, k)]) / (DX * DX);
    }
 
    if (k == 0)
    {
        diffsecz = (2 * v[gety(i, j, k)] - 5 * v[gety(i, j, k + 1)] + 4 * v[gety(i, j, k + 2)] - v[gety(i, j, k + 3)]) / (DZ * DZ);
    }
    else if (k == NZ)
    {
        diffsecz = (-2 * v[gety(i, j, k)] + 5 * v[gety(i, j, k - 1)] - 4 * v[gety(i, j, k - 2)] + v[gety(i, j, k - 3)]) / (DZ * DZ);
    }
    else
    {
        diffsecz = (v[gety(i, j, k + 1)] - 2 * v[gety(i, j, k)] + v[gety(i, j, k - 1)]) / (DZ * DZ);
    }
 
    return -1 * (u_value * diffx + v_value * diffy + w_value * diffz) + 1 / RE * (diffsecx + diffsecy + diffsecz) + functionG_v(i-1 + coords[0] * other_dim_x_y, j-1 + coords[1] * other_dim_y_y, k, t);
}

Real IcoNS::functionF_w(const std::vector<Real> &u, const std::vector<Real> &v, const std::vector<Real> &w, int i, int j, int k, Real t)

{
    // int lu = i * newDimY_x * dim_z + j * dim_z + k;
    // int lv = i * newDimY_y * dim_z + j * dim_z + k;
    // int lw = i * newDimY_z * dim_z_z + j * dim_z_z + k;

    // Real value = -((u[lu] + u[lu - newDimY_x * dim_z] + u[lu + 1] + u[lu - dim_z * newDimY_x + 1]) / 4.0 * (w[lw + newDimY_z * dim_z_z] - w[lw - newDimY_z * dim_z_z]) / (2.0 * DX) +
    //               (v[lv + 1] + v[lv - dim_z + 1] + v[lv] + v[lv - dim_z]) / 4.0 * (w[lw + dim_z_z] - w[lw - dim_z_z]) / (2.0 * DY) +
    //               w[lw] * (w[lw + 1] - w[lw - 1]) / (2.0 * DZ)) +
    //               (1.0 / RE) * ((w[lw + newDimY_z * dim_z_z] - 2.0 * w[lw] + w[lw - newDimY_z * dim_z_z]) / (DX * DX) +
    //                      (w[lw + dim_z_z] - 2.0 * w[lw] + w[lw - dim_z_z]) / (DY * DY) +
    //                      (w[lw + 1] - 2.0 * w[lw] + w[lw - 1]) / (DZ * DZ));
    
    // if(testCase==0){
    //     value += functionG_w(i-1 + coords[0] * other_dim_x_z, j-1 + coords[1] * other_dim_y_z, k, t);
    // }
    // return 0;

    Real diffx, diffy, diffz, u_value, v_value, w_value, diffsecx, diffsecy, diffsecz;
    w_value = w[getz(i, j, k)];
    if (i == 1 || i == newDimX_z - 1)
    {
        u_value = boundary.boundary_value_u[0]->value(i - 1.5, j - 1, k, t);
    }
    else
    {
        u_value = (u[getx(i, j, k)] + u[getx(i, j, k + 1)] + u[getx(i - 1, j, k)] + u[getx(i - 1, j, k + 1)]) / 4.0;
    }
    if (j == 1 || j == newDimY_z - 1)
    {
        v_value = boundary.boundary_value_v[0]->value(i - 1, j - 1.5, k, t);
    }
    else
    {
        v_value = (v[gety(i, j, k)] + v[gety(i, j - 1, k)] + v[gety(i, j, k + 1)] + v[gety(i, j - 1, k + 1)]) / 4.0;
    }
 
    if (i == 1)
    {
        diffx = (-3 * w[getz(i, j, k)] + 4 * w[getz(i + 1, j, k)] - w[getz(i + 2, j, k)]) / (2 * DX);
    }
    else if (i == newDimX_z - 1)
    {
        diffx = (3 * w[getz(i, j, k)] - 4 * w[getz(i - 1, j, k)] + w[getz(i - 2, j, k)]) / (2 * DX);
    }
    else
    {
        diffx = (w[getz(i + 1, j, k)] - w[getz(i - 1, j, k)]) / (2.0 * DX);
    }
 
    if (j == 1)
    {
        diffy = (-3 * w[getz(i, j, k)] + 4 * w[getz(i, j + 1, k)] - w[getz(i, j + 2, k)]) / (2 * DY);
    }
    else if (j == newDimY_z -1)
    {
        diffy = (3 * w[getz(i, j, k)] - 4 * w[getz(i, j - 1, k)] + w[getz(i, j - 2, k)]) / (2 * DY);
    }
    else
    {
        diffy = (w[getz(i, j + 1, k)] - w[getz(i, j - 1, k)]) / (2.0 * DY);
    }
 
    if (k == 0)
    {
        diffz = (-4 * boundary.boundary_value_w[0]->value(i - 1, j - 1, k - 0.5, t) + 3 * w[getz(i, j, k)] + w[getz(i, j, k + 1)]) / (3.0 * DZ);
    }
    else if (k == NZ - 1)
    {
        diffz = (4 * boundary.boundary_value_w[0]->value(i - 1, j - 1, NZ - 0.5, t) - 3 * w[getz(i, j, k)] - w[getz(i, j, k - 1)]) / (3.0 * DZ);
    }
    else
    {
        diffz = (w[getz(i, j, k + 1)] - w[getz(i, j, k - 1)]) / (2.0 * DZ);
    }
 
    if (k == 0)
    {
        diffsecz = (16 * boundary.boundary_value_w[0]->value(i - 1, j - 1, k - 0.5, t) - 25 * w[getz(i, j, k)] + 10 * w[getz(i, j, k + 1)] - w[getz(i, j, k + 2)]) / (5.0 * DZ * DZ);
    }
    else if (k == NZ - 1)
    {
        diffsecz = (-16 * boundary.boundary_value_w[0]->value(i - 1, j - 1, NZ - 0.5, t) + 25 * w[getz(i, j, k)] - 10 * w[getz(i, j, k - 1)] + w[getz(i, j, k - 2)]) / (5.0 * DZ * DZ);
    }
    else
    {
        diffsecz = (w[getz(i, j, k + 1)] - 2 * w[getz(i, j, k)] + w[getz(i, j, k - 1)]) / (DZ * DZ);
    }
 
    if (i == 1)
    {
        diffsecx = (2 * w[getz(i, j, k)] - 5 * w[getz(i + 1, j, k)] + 4 * w[getz(i + 2, j, k)] - w[getz(i + 3, j, k)]) / (DX * DX);
    }
    else if (i == newDimX_z - 1)
    {
        diffsecx = (-2 * w[getz(i, j, k)] + 5 * w[getz(i - 1, j, k)] - 4 * w[getz(i - 2, j, k)] + w[getz(i - 3, j, k)]) / (DX * DX);
    }
    else
    {
        diffsecx = (w[getz(i + 1, j, k)] - 2 * w[getz(i, j, k)] + w[getz(i - 1, j, k)]) / (DX * DX);
    }
 
    if (j == 1)
    {
        diffsecy = (2 * w[getz(i, j, k)] - 5 * w[getz(i, j + 1, k)] + 4 * w[getz(i, j + 2, k)] - w[getz(i, j + 3, k)]) / (DY * DY);
    }
    else if (j == newDimY_z -1)
    {
        diffsecy = (-2 * w[getz(i, j, k)] + 5 * w[getz(i, j - 1, k)] - 4 * w[getz(i, j - 2, k)] + w[getz(i, j - 3, k)]) / (DY * DY);
    }
    else
    {
        diffsecy = (w[getz(i, j + 1, k)] - 2 * w[getz(i, j, k)] + w[getz(i, j - 1, k)]) / (DY * DY);
    }
 
    return -1 * (u_value * diffx + v_value * diffy + w_value * diffz) + 1 / RE * (diffsecx + diffsecy + diffsecz) + functionG_w(i-1 + coords[0] * other_dim_x_z, j-1 + coords[1] * other_dim_y_z, k, t);

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
// TODO: exact derivative missing, don't know if needed
void IcoNS::pressionCorrection(double *Y2_p)
{
    // is the denominator 3*DX correct? -> 2*DX ?
    // LEFT FACE
    // TODO: how does iteration on p works? bc rn every processor is skipping first element if it is not boundary
    if (lbx)
    {
        for (int j = lby; j < zSize[1] - rby; j++)
        {
            for (int k = 1; k < zSize[2] - 1; k++)
            {
                Y2_p[getp(0, j, k)] = (4 / 3) * Y2_p[getp(1, j, k)] - (1 / 3) * Y2_p[getp(2, j, k)];
            }
        }
    }

    if (rbx)
    {
        // RIGHT FACE
        for (int j = lby; j < zSize[1] - rby; j++)
        {
            for (int k = 1; k < zSize[2] - 1; k++)
            {
                Y2_p[getp(zSize[0] - 1, j, k)] = (4 / 3) * Y2_p[getp(zSize[0] - 2, j, k)] - (1 / 3) * Y2_p[getp(zSize[0] - 3, j, k)];
            }
        }
    }

    // FRONT FACE
    if (lby)
    {
        for (int i = lbx; i < zSize[0] - rbx; i++)
        {
            for (int k = 1; k < zSize[2] - 1; k++)
            {
                Y2_p[getp(i, 0, k)] = (4 / 3) * Y2_p[getp(i, 1, k)] - (1 / 3) * Y2_p[getp(i, 2, k)];
            }
        }
    }

    // BACK FACE
    if (rby)
    {
        for (int i = lbx; i < zSize[0] - rbx; i++)
        {
            for (int k = 1; k < zSize[2] - 1; k++)
            {
                Y2_p[getp(i, zSize[1] - 1, k)] = (4 / 3) * Y2_p[getp(i, zSize[1] - 2, k)] - (1 / 3) * Y2_p[getp(i, zSize[1] - 3, k)];
            }
        }
    }

    // LOWER FACE
    for (int i = lbx; i < zSize[0] - rbx; i++)
    {
        for (int j = lby; j < zSize[1] - rby; j++)
        {
            Y2_p[getp(i, j, 0)] = (4 / 3) * Y2_p[getp(i, j, 1)] - (1 / 3) * Y2_p[getp(i, j, 2)];
        }
    }
    // UPPER FACE
    for (int i = lbx; i < zSize[0] - rbx; i++)
    {
        for (int j = lby; j < zSize[1] - rby; j++)
        {
            Y2_p[getp(i, j, zSize[2] - 1)] = (4 / 3) * Y2_p[getp(i, j, zSize[2] - 2)] - (1 / 3) * Y2_p[getp(i, j, zSize[2] - 3)];
        }
    }
    // 4 X EDGES
    for (int i = lbx; i < zSize[0] - rbx; i++)
    {
        // LOWER FRONT EDGE
        if (lby)
        {
            Y2_p[getp(i, 0, 0)] = (4 / 3) * Y2_p[getp(i, 1, 0)] - (1 / 3) * Y2_p[getp(i, 2, 0)];
            // UPPER FRONT EDGE
            Y2_p[getp(i, 0, zSize[2] - 1)] = (4 / 3) * Y2_p[getp(i, 1, zSize[2] - 1)] - (1 / 3) * Y2_p[getp(i, 2, zSize[2] - 1)];
        }
        if (rby)
        {
            // LOWER BACK EDGE
            Y2_p[getp(i, (zSize[1] - 1), 0)] = (4 / 3) * Y2_p[getp(i, zSize[1] - 2, 0)] - (1 / 3) * Y2_p[getp(i, zSize[1] - 3, 0)];
            // UPPER BACK EDGE
            Y2_p[getp(i, zSize[1] - 1, zSize[2] - 1)] = (4 / 3) * Y2_p[getp(i, zSize[1] - 2, zSize[2] - 1)] - (1 / 3) * Y2_p[getp(i, zSize[1] - 3, zSize[2] - 1)];
        }
    }
    if (lbx)
    {
        // 4 Y EDGES
        for (int j = lby; j < zSize[1] - rby; j++)
        {
            // LOWER LEFT EDGE
            Y2_p[getp(0, j, 0)] = (4 / 3) * Y2_p[getp(1, j, 0)] - (1 / 3) * Y2_p[getp(2, j, 0)];
            // UPPER LEFT EDGE
            Y2_p[getp(0, j, zSize[2] - 1)] = (4 / 3) * Y2_p[getp(1, j, zSize[2] - 1)] - (1 / 3) * Y2_p[getp(2, j, zSize[2] - 1)];
        } // break to exploit locality
    }
    if (rbx)
    {
        for (int j = lby; j < zSize[1] - rby; j++)
        {
            // LOWER RIGHT EDGE
            Y2_p[getp(zSize[0] - 1, j, 0)] = (4 / 3) * Y2_p[getp(zSize[0] - 2, j, 0)] - (1 / 3) * Y2_p[getp(zSize[0] - 3, j, 0)];
            // UPPER RIGHT EDGE
            Y2_p[getp(zSize[0] - 1, j, zSize[2] - 1)] = (4 / 3) * Y2_p[getp(zSize[0] - 2, j, zSize[2] - 1)] - (1 / 3) * Y2_p[getp(zSize[0] - 3, j, zSize[2] - 1)];
        }
    }
    if (lbx && lby)
    {
        // 4 Z EDGES
        for (int k = 1; k < zSize[2] - 1; k++)
        {
            // FRONT LEFT EDGE
            Y2_p[getp(0, 0, k)] = (4 / 3) * Y2_p[getp(1, 0, k)] - (1 / 3) * Y2_p[getp(2, 0, k)];
        } // break to exploit locality
    }
    if (lbx && rby)
    {
        for (int k = 1; k < zSize[2] - 1; k++)
        {
            // BACK LEFT EDGE
            Y2_p[getp(0, zSize[1] - 1, k)] = (4 / 3) * Y2_p[getp(1, zSize[1] - 1, k)] - (1 / 3) * Y2_p[getp(2, zSize[1] - 1, k)];
        } // break to exploit locality
    }
    if (rbx && lby)
    {
        for (int k = 1; k < zSize[2] - 1; k++)
        {
            // FRONT RIGHT EDGE
            Y2_p[getp(zSize[0] - 1, 0, k)] = (4 / 3) * Y2_p[getp(zSize[0] - 2, 0, k)] - (1 / 3) * Y2_p[getp(zSize[0] - 3, 0, k)];
        } // break to exploit locality
    }
    if (rbx && rby)
    {
        for (int k = 1; k < zSize[2] - 1; k++)
        {
            // BACK RIGHT EDGE
            Y2_p[getp(zSize[0] - 1, zSize[1] - 1, k)] = (4 / 3) * Y2_p[getp(zSize[0] - 2, zSize[1] - 1, k)] - (1 / 3) * Y2_p[getp(zSize[0] - 3, zSize[1] - 1, k)];
        }
    }
    // 8 VERTICES
    if (lbx && lby)
    {
        // LOWER FRONT LEFT VERTEX (0,0,0)
        Y2_p[getp(0, 0, 0)] = (4 / 3) * Y2_p[getp(1, 0, 0)] - (1 / 3) * Y2_p[getp(2, 0, 0)];
    }
    if (rbx && lby)
    {
        // LOWER FRONT RIGHT VERTEX (1,0,0)
        Y2_p[getp(zSize[0] - 1, 0, 0)] = (4 / 3) * Y2_p[getp(zSize[0] - 2, 0, 0)] - (1 / 3) * Y2_p[getp(zSize[0] - 3, 0, 0)];
    }
    if (lbx && rby)
    {
        // LOWER BACK LEFT VERTEX (0,1,0)
        Y2_p[getp(0, zSize[1] - 1, 0)] = (4 / 3) * Y2_p[getp(1, zSize[1] - 1, 0)] - (1 / 3) * Y2_p[getp(2, zSize[1] - 1, 0)];
    }
    if (lbx && lby)
    {
        // UPPER FRONT LEFT VERTEX (0,0,1)
        Y2_p[getp(0, 0, zSize[2] - 1)] = (4 / 3) * Y2_p[getp(1, 0, zSize[2] - 1)] - (1 / 3) * Y2_p[getp(2, 0, zSize[2] - 1)];
    }
    if (lbx && rby)
    {
        // UPPER BACK LEFT VERTEX (0,1,1)
        Y2_p[getp(0, zSize[1] - 1, zSize[2] - 1)] = (4 / 3) * Y2_p[getp(1, zSize[1] - 1, zSize[2] - 1)] - (1 / 3) * Y2_p[getp(2, zSize[1] - 1, zSize[2] - 1)];
    }
    if (rbx && lby)
    {
        // UPPER FRONT RIGHT VERTEX (1,0,1)
        Y2_p[getp(zSize[0] - 1, 0, zSize[2] - 1)] = (4 / 3) * Y2_p[getp(zSize[0] - 2, 0, zSize[2] - 1)] - (1 / 3) * Y2_p[getp(zSize[0] - 3, 0, zSize[2] - 1)];
    }
    if (rbx && rby)
    {
        // LOWER BACK RIGHT VERTEX (1,1,0)
        Y2_p[getp(zSize[0] - 1, zSize[1] - 1, 0)] = (4 / 3) * Y2_p[getp(zSize[0] - 2, zSize[1] - 1, 0)] - (1 / 3) * Y2_p[getp(zSize[0] - 3, zSize[1] - 1, 0)];
        // UPPER BACK RIGHT VERTEX (1,1,1)
        Y2_p[getp(zSize[0] - 1, zSize[1] - 1, zSize[2] - 1)] = (4 / 3) * Y2_p[getp(zSize[0] - 2, zSize[1] - 1, zSize[2] - 1)] - (1 / 3) * Y2_p[getp(zSize[0] - 3, zSize[1] - 1, zSize[2] - 1)];
    }

    // compute the average value of the pressure field and subtract it from the pressure field
    Real sum = 0.0;
    for (int i = lbx; i < zSize[0] - rbx; i++)
    {
        for (int j = lby; j < zSize[1] - rby; j++)
        {
            for (int k = 1; k < zSize[2] - 1; k++)
            {
                sum += Y2_p[getp(i, j, k)];
            }
        }
    }
    Real average = sum / (zSize[0] * zSize[1] * zSize[2]);
    for (int i = lbx; i < zSize[0] - rbx; i++)
    {
        for (int j = lby; j < zSize[1] - rby; j++)
        {
            for (int k = 1; k < zSize[2] - 1; k++)
            {
                Y2_p[getp(i, j, k)] -= average;
            }
        }
    }

}