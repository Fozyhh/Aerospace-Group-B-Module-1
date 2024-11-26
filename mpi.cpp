#include <iostream>
#include <array>
#include <mpi.h>
#include <vector>

#define NX 4
#define NY 4
#define NZ 4

#define PX 2
#define PY 2
#define PZ 1

int main(int argc, char **argv)
{

    // std::array<double, NX * (NY+1) * (NZ+1)> grid{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,
    //                                                    16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,
    //                                                  32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,
    //                                                48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,
    //                                              64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,
    //                                             80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,
    //                                           96,97,98,99};
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Create a Cartesian topology (2D)
    MPI_Comm cart_comm;
    int dims[2] = {PX, PY};
    int periods[2] = {0, 0};
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &cart_comm);
    /*
    // x-component
        {
            constexpr int local_x = NX / PX;
            constexpr int local_y = (NY + 1) / PY + 1;
            // constexpr int local_y_last = (NY + 1) / PY + 1;
            // int local_y = ((NY+1) % PY != 0 && rank == size-1) ? local_y_last : local_y;
            constexpr int local_z = NZ + 1;
            std::vector<double> local_grid(local_x * local_y * local_z);

            int coords[2];
            MPI_Cart_coords(cart_comm, rank, 2, coords);
            int coord_x = coords[0];
            int coord_y = coords[1];

            for (int i = 0; i < local_x; i++)
            {
                for (int j = 0; j < local_y; j++)
                {
                    for (int k = 0; k < local_z; k++)
                    {
                        int global_i = coord_x * local_x + i;
                        int global_j = coord_y * local_y + j;
                        local_grid[i * local_y * local_z + j * local_z + k] = grid[global_i * (NZ + 1) * (NY + 1) + global_j * (NZ + 1) + k];
                    }
                }
            }

            // print local grid for each process
            for (int i = 0; i < size; i++)
            {
                if (rank == i)
                {
                    std::cout << "Rank: " << rank << std::endl;
                    std::cout << "Coords: " << coord_x << " " << coord_y << std::endl;

                    for (int i = 0; i < local_x; i++)
                    {
                        for (int j = 0; j < local_y; j++)
                        {
                            for (int k = 0; k < local_z; k++)
                            {
                                std::cout << local_grid[i * local_y * local_z + j * local_z + k] << " ";
                            }
                            std::cout << std::endl;
                        }
                        // std::cout << std::endl;
                    }
                }
                MPI_Barrier(MPI_COMM_WORLD);
            }
        }
    */

    int dim_x = NX / PX;
    int dim_y = (NY + 1) / PY;
    int dim_z = NZ + 1;
    int coords[2];
    MPI_Cart_coords(cart_comm, rank, 2, coords);
    int neighbors[4];

    MPI_Cart_shift(cart_comm, 0, 1, &neighbors[0], &neighbors[2]);

    MPI_Cart_shift(cart_comm, 1, 1, &neighbors[1], &neighbors[3]);

    if (NX % PX != 0 && coords[0] == dim_x - 1)
        dim_x++;
    if ((NY + 1) % PY != 0 && coords[1] == 1/*dim_y - 1*/){
        dim_y++;
    }
        

    int res_x = 0;
    int res_y = 0;

    for (int i = 0; i < 4; i++)
    {
        if (neighbors[i] != -2)
        {
            if (i == 0 || i == 2)
            {
                res_x++;
            }
            else
            {
                res_y++;
            }
        }
    }
    int newDimX = dim_x + res_x;
    int newDimY = dim_y + res_y;
    std::vector<int> grid_loc(newDimX * newDimY * dim_z);

    


    for (int i = 0; i < newDimX; i++)
    {
        for (int j = 0; j < newDimY ; j++)
        {
            for (int k = 0; k < dim_z; k++)
            {

                // ---> z   
                // ^
                // | y,  uscente x 
                int glob_address_x = i +coords[0]*dim_x;
                int glob_address_y =j+ coords[1]*(dim_y-2); // il -1 andrebbe solo nell'ultimo processore dato che ha la dim_y diversa dagli altri(o si fa una var diversa)
                //boundary points segnati nella griglia
                if(glob_address_x == 0 || glob_address_x==NX)
                    grid_loc[i * newDimY * dim_z + j * dim_z + k] += 2;
                if(glob_address_y== 0 || glob_address_y==NY)
                    grid_loc[i * newDimY * dim_z + j * dim_z + k] += 3;
                if(k==0 || k==NZ)
                    grid_loc[i * newDimY * dim_z + j * dim_z + k] += 4;

                
                if((rank == 0 || rank ==1) && i == dim_x) grid_loc[i * newDimY * dim_z + j * dim_z + k] =0;
                if((rank == 2 || rank== 3) && i == 0) grid_loc[i * newDimY * dim_z + j * dim_z + k] =0;
                //Setting ghost points to 0
                if((rank == 0 || rank ==2) && j == dim_y) grid_loc[i * newDimY * dim_z + j * dim_z + k] =0;
                if((rank == 1 || rank ==3) && j == 0) grid_loc[i * newDimY * dim_z + j * dim_z + k] =0;
            }
        }
    }
    //may be dim_x instead of newDimX - 1 but that will come with more processes
    for (int i = 1; i < newDimX - 1; i++)
    {
        for (int j = 1; j < newDimY - 1; j++)
        {
            for (int k = 1; k < dim_z - 1; k++)
            {
                grid_loc[i * newDimY * dim_z + j * dim_z + k] = 1;
                
            }
        }
    }


    // Probabilmente il problema è qui. Considerare che ogni rank definisce 
    // la propria MPI_face (per il momento non è un problema perchè i processi che scambiano dati:
    // 0 e 2, 1 e 3, sono adiacenti e quindi hanno la stessa dimensione di face)
    
    //TODO: riguardare perchè testata con 2 proc e basta
    MPI_Datatype MPI_face_x;
    MPI_Type_vector(1,dim_z* dim_y, res_y, MPI_INT, &MPI_face_x);
    MPI_Type_commit(&MPI_face_x);
    // MPI_Status status1;

    // if (neighbors[0] != -2)
    // {
    //     MPI_Send(&grid_loc[newDimY * dim_z], 1, MPI_face_x, neighbors[0], rank, cart_comm);
    // }

    // if (neighbors[2] != -2)
    // {
    //     MPI_Recv(&grid_loc[(newDimX - 1) * newDimY * dim_z], 1, MPI_face_x, neighbors[2], neighbors[2], cart_comm, &status1);
    // }

    //LO STRIDE INIZIA A CONTARE DALL INIZIO DEL BLOCCO PRECEDENTE!!!!
    MPI_Datatype MPI_face_y;
    MPI_Type_vector(dim_x,dim_z, (dim_y)*dim_z + dim_z, MPI_INT, &MPI_face_y);
    MPI_Type_commit(&MPI_face_y);
    
    //processore 1 invia a 0
    MPI_Status status2;
    if (rank == 1)//(neighbors[1] != -2)
    {
        MPI_Send(&grid_loc[dim_z], 1, MPI_face_y, neighbors[1], rank, cart_comm);
    }

    if (rank == 0) //(neighbors[3] != -2)
    {
        MPI_Recv(&grid_loc[dim_z*(dim_y)], 1, MPI_face_y, neighbors[3], neighbors[3], cart_comm, &status2);
    }

    // 3 invia a 2
    MPI_Status status3;
    if (rank == 3)//(neighbors[1] != -2)
    {
        MPI_Send(&grid_loc[dim_z * newDimY + dim_z], 1, MPI_face_y, neighbors[1], rank, cart_comm);
    }

    if (rank == 2) //(neighbors[3] != -2)
    {
        MPI_Recv(&grid_loc[dim_z*(dim_y) + (newDimY) * dim_z], 1, MPI_face_y, neighbors[3], neighbors[3], cart_comm, &status3);
    }

    MPI_Barrier(cart_comm);
    if (rank == 2){
        std::cout << " rank: " << rank /*rank << "N0: " << neighbors[0] << " N2: " <<neighbors[2]*/ << std::endl;
        for (int i = 0; i < newDimX; i++)
        {
            for (int j = 0; j < newDimY; j++)
            {
                for (int k = 0; k < dim_z; k++)
                {
                    std::cout << grid_loc[i * newDimY * dim_z + j * dim_z + k] << " ";
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }
    }

    MPI_Barrier(cart_comm);

    // if (rank == 1){
    //     std::cout << " rank: " << rank /*rank << "N0: " << neighbors[0] << " N2: " <<neighbors[2]*/ << std::endl;
    //     for (int i = 0; i < newDimX; i++)
    //     {
    //         for (int j = 0; j < newDimY; j++)
    //         {
    //             for (int k = 0; k < dim_z; k++)
    //             {
    //                 std::cout << grid_loc[i * newDimY * dim_z + j * dim_z + k] << " ";
    //             }
    //             std::cout << std::endl;
    //         }
    //         std::cout << std::endl;
    //     }
    // }

    MPI_Barrier(cart_comm);

    if (rank ==0)
    {
        std::cout << "Rank: " << rank << std::endl;
        std::cout << "Dim_x: " << dim_x << " Dim_y: " << dim_y << std::endl;
        std::cout << "Res_x: " << res_x << " Res_y: " << res_y << std::endl;
        std::cout << "NewDimX: " << newDimX << " NewDimY: " << newDimY << std::endl;
        std::cout << "upper neighbour: " << neighbors[1] << " lower: " << neighbors[3] << std::endl;
        std::cout << "left neighbour: " << neighbors[0] << " right: " << neighbors[2] << std::endl;
        std::cout << "Coords: " << coords[0] << "," << coords[1] << std::endl<< std::endl;

    }
    MPI_Barrier(cart_comm);
    if (rank == 1)
    {
        std::cout << "Rank: " << rank << std::endl;
        std::cout << "Dim_x: " << dim_x << " Dim_y: " << dim_y << std::endl;
        std::cout << "Res_x: " << res_x << " Res_y: " << res_y << std::endl;
        std::cout << "NewDimX: " << newDimX << " NewDimY: " << newDimY << std::endl;
        std::cout << "upper neighbour: " << neighbors[1] << " lower: " << neighbors[3] << std::endl;
        std::cout << "left neighbour: " << neighbors[0] << " right: " << neighbors[2] << std::endl;
        std::cout << "Coords: " << coords[0] << "," << coords[1] << std::endl;
    }   

    MPI_Finalize();
    return 0;
}
