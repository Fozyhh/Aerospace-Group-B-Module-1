#include <iostream>
#include <array>
#include <mpi.h>
#include <vector>

#define NX 4
#define NY 4
#define NZ 4

#define PX 2
#define PY 1
#define PZ 1

/**
 * 2D Cartesian Topology with periodic boundaries, following 2decomp partitioning.
 * Only the general case is neeeded to handle all the communications. 
 */


// Check if the process is in the middle of the topology.
bool isMiddleProcess(int coords[2])
{
    return coords[0] != 0 && coords[0] != PX - 1 && coords[1] != 0 && coords[1] != PY - 1;
}

// Check if the process is on edge of the topology.
bool isEdgeProcess(int coords [2])
{
    return coords[0] == 0 || coords[0] == PX - 1 || coords[1] == 0 || coords[1] == PY - 1;
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Create a Cartesian topology (2D)
    MPI_Comm cart_comm;
    int dims[2] = {PX, PY};
    int periods[2] = {1, 0};
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &cart_comm);

    int dim_x = (NX+1) / PX;
    int dim_y = (NY) / PY;
    int dim_z = NZ + 1;
    int coords[2];
    MPI_Cart_coords(cart_comm, rank, 2, coords);
    int neighbors[4];

    MPI_Cart_shift(cart_comm, 0, 1, &neighbors[0], &neighbors[2]);

    MPI_Cart_shift(cart_comm, 1, 1, &neighbors[1], &neighbors[3]);

    if (NX % PX != 0 && coords[0] == dim_x - 1)
        dim_x++;
    if ((NY + 1) % PY != 0 && coords[1] == 1 /*dim_y - 1*/)
    {
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
    

    
    // may be dim_x instead of newDimX - 1 but that will come with more processes
    for (int i = 1; i < newDimX -1 ; i++)
    {
        for (int j = 1; j < newDimY - 1; j++)
        {
            for (int k = 1; k < dim_z - 1; k++)
            {
                grid_loc[i * newDimY * dim_z + j * dim_z + k] = 1;
            }
        }
    }

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
                // if(glob_address_x == 1 || glob_address_x==NX)
                //     grid_loc[i * newDimY * dim_z + j * dim_z + k] += 2;
                if(glob_address_x == 1 || glob_address_x==NX)
                    grid_loc[i * newDimY * dim_z + j * dim_z + k] += 2;

                if(glob_address_y== 1 || glob_address_y==NY)
                    grid_loc[i * newDimY * dim_z + j * dim_z + k] += 3;
                if(k==0 || k==NZ)
                    grid_loc[i * newDimY * dim_z + j * dim_z + k] += 4;

                if(i==0 || i == dim_x+1) grid_loc[i * newDimY * dim_z + j * dim_z + k] =0; //only x
                // if(j==0 || j==dim_y+1) grid_loc[i * newDimY * dim_z + j * dim_z + k] =0;//only y
                //if((rank == 0 || rank ==1) && i == dim_x) grid_loc[i * newDimY * dim_z + j * dim_z + k] =0;
                // if((rank == 0 ) && i == dim_x) grid_loc[i * newDimY * dim_z + j * dim_z + k] =0;
                // //if((rank == 0 || rank ==1) && i == dim_x) grid_loc[i * newDimY * dim_z + j * dim_z + k] =0;
                // //if((rank == 2 || rank== 3) && i == 0) grid_loc[i * newDimY * dim_z + j * dim_z + k] =0;
                // //Setting ghost points to 0
                // //if((rank == 0 || rank ==2) && j == dim_y) grid_loc[i * newDimY * dim_z + j * dim_z + k] =0;
                // if((rank == 1 || rank ==3) && i/*j*/ == 0) grid_loc[i * newDimY * dim_z + j * dim_z + k] =0;
            }
        }
    }



    

    // Probabilmente il problema è qui. Considerare che ogni rank definisce
    // la propria MPI_face (per il momento non è un problema perchè i processi che scambiano dati:
    // 0 e 2, 1 e 3, sono adiacenti e quindi hanno la stessa dimensione di face)

    // LO STRIDE INIZIA A CONTARE DALL INIZIO DEL BLOCCO PRECEDENTE!!!!
    MPI_Datatype MPI_face_x;
    MPI_Type_vector(dim_x, dim_z, (dim_y)*dim_z + 2*dim_z, MPI_INT, &MPI_face_x);
    MPI_Type_commit(&MPI_face_x);

   

    

    // TODO: riguardare perchè testata con 2 proc e basta
    MPI_Datatype MPI_face_y;
    MPI_Type_vector(1, dim_z * newDimY, 0/*2 *dim_z for multiple faces communication if needed*/, MPI_INT, &MPI_face_y);
    MPI_Type_commit(&MPI_face_y);
    MPI_Status status1;

 

    MPI_Status status11;

   
    MPI_Barrier(cart_comm);
// CORNERS, make ranks general
//{
    /****** TOP LEFT CORNER ******/
    // 0
    // |
    // v
    // 1
    MPI_Status status4;
    if (rank == 0) //(neighbors[1] != -2)
    {
        MPI_Send(&grid_loc[dim_z], 1, MPI_face_x, neighbors[3], rank, cart_comm);
    }

    if (rank == 1) //(neighbors[3] != -2)
    {
        MPI_Recv(&grid_loc[0], 1, MPI_face_x, neighbors[1], neighbors[1], cart_comm, &status4);
    }

    /****** BOTTOM LEFT CORNER ******/
    // 0
    // ^
    // |
    // 1
    MPI_Status status2;
    if (rank == 1) //(neighbors[1] != -2)
    {
        MPI_Send(&grid_loc[dim_z], 1, MPI_face_x, neighbors[1], rank, cart_comm);
    }

    if (rank == 0) //(neighbors[3] != -2)
    {
        MPI_Recv(&grid_loc[dim_z * (dim_y)], 1, MPI_face_x, neighbors[3], neighbors[3], cart_comm, &status2);
    }

    /****** TOP RIGHT CORNER ******/
    // 2
    // |
    // v
    // 3
    MPI_Status status5;
    if (rank == 2) //(neighbors[1] != -2)
    {
        MPI_Send(&grid_loc[dim_z * newDimY + dim_z], 1, MPI_face_x, neighbors[3], rank, cart_comm);
    }

    if (rank == 3) //(neighbors[3] != -2)
    {
        MPI_Recv(&grid_loc[dim_z * newDimY], 1, MPI_face_x, neighbors[1], neighbors[1], cart_comm, &status5);
    }
    
     /****** BOTTOM RIGHT CORNER ******/
    2
    ^
    |
    3


    General
    MPI_Status status3;
    if (rank == 3) //(neighbors[1] != -2)
    {
        MPI_Send(&grid_loc[dim_z * newDimY + dim_z], 1, MPI_face_x, neighbors[1], rank, cart_comm);
    }

    if (rank == 2) //(neighbors[3] != -2)
    {
        MPI_Recv(&grid_loc[dim_z * (dim_y) + (newDimY)*dim_z], 1, MPI_face_x, neighbors[3], neighbors[3], cart_comm, &status3);
    }

    /**** LEFT CORNERS TO THE PROCESSES ON THEIR RIGHT ****/
    // 0 -> 2 and 1 -> 3
    if (neighbors[2] != -2)
    {
        MPI_Send(&grid_loc[(dim_x - 1) * newDimY * dim_z], 1, MPI_face_y, neighbors[2], rank, cart_comm);
    }

    if (neighbors[0] != -2)
    {
        MPI_Recv(&grid_loc[0], 1, MPI_face_y, neighbors[0], neighbors[0], cart_comm, &status1);
    }

    /**** RIGHT CORNERS TO THE PROCESSES ON THEIR LEFT ****/
    // 0 <- 2 and 1 <- 3
    if (neighbors[0] != -2)
    {
        MPI_Send(&grid_loc[newDimY * dim_z], 1, MPI_face_y, neighbors[0], rank, cart_comm);
    }

    if (neighbors[2] != -2)
    {
        MPI_Recv(&grid_loc[(dim_x)*newDimY * dim_z], 1, MPI_face_y, neighbors[2], neighbors[2], cart_comm, &status1);
    }

}
if(false){
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
} 

/*************************** GENERAL CASE ***************************/
if (true)//isMiddleProcess(coords))
{
    // TO THE LEFT
    // (x, y) -> (x-1, y)
    MPI_Status status6;
    if (neighbors[0] != -2)
    {
        // MPI_Send(&grid_loc[newDimY * dim_z /*+ dim_z*/], 1, MPI_face_y, neighbors[0], rank, cart_comm); // caso generale
        MPI_Send(&grid_loc[(2*newDimY) * dim_z /*+ dim_z*/], 1, MPI_face_y, neighbors[0], rank, cart_comm); //limite sinistro
    }
    if (neighbors[2] != -2)
    {
        MPI_Recv(&grid_loc[(dim_x+1)*newDimY*dim_z /*+ dim_z*/], 1, MPI_face_y, neighbors[2], neighbors[2], cart_comm, &status6);
    }

    // TO THE RIGHT
    // (x, y) -> (x+1, y)
    MPI_Status status7;
    if (neighbors[2] != -2)
    {
        MPI_Send(&grid_loc[(dim_x) * newDimY * dim_z/*+ dim_z */], 1, MPI_face_y, neighbors[2], rank, cart_comm);//caso generale
        // MPI_Send(&grid_loc[(dim_x-1) * newDimY * dim_z/*+ dim_z*/], 1, MPI_face_y, neighbors[2], rank, cart_comm);// limite destro     
    }
    if(neighbors[0] != -2){
        MPI_Recv(&grid_loc[0/*+dim_z*/], 1, MPI_face_y, neighbors[0], neighbors[0], cart_comm, &status7);
    }

    // TO THE TOP
    // (x, y) -> (x, y+1)
    MPI_Status status8;
    if (neighbors[1] != -2)
    {
        // First row of the local grid: skip the first face of ghos points, skip the first column of ghost points.
        MPI_Send(&grid_loc[/*newDimY * dim_z +*/ dim_z], 1, MPI_face_x, neighbors[1], rank, cart_comm);//general case
        // MPI_Send(&grid_loc[/*newDimY * dim_z +*/ 2*dim_z], 1, MPI_face_x, neighbors[1], rank, cart_comm);//limite sopra
    }
    if (neighbors[3] != -2)
    {
        MPI_Recv(&grid_loc[/*newDimY*dim_z +*/ (newDimY-1)*dim_z], 1, MPI_face_x, neighbors[3], neighbors[3], cart_comm, &status8);
    }

    // TO THE BOTTOM
    // (x, y) -> (x, y-1)
    MPI_Status status9;
    if (neighbors[3] != -2)
    {
        MPI_Send(&grid_loc[/*dim_z*newDimY +*/ dim_z*(newDimY-2)/*dim_z * (2 * newDimY - 1)*/], 1, MPI_face_x, neighbors[3], rank, cart_comm);//caso generale
        // MPI_Send(&grid_loc[/*dim_z*newDimY +*/ dim_z*(newDimY-3)/*dim_z * (2 * newDimY - 1)*/], 1, MPI_face_x, neighbors[3], rank, cart_comm);//limite sotto
    }
    if (neighbors[1] != -2)
    {
        MPI_Recv(&grid_loc[/*dim_z * newDimY+*/0], 1, MPI_face_x, neighbors[1], neighbors[1], cart_comm, &status9);
    }
}

if (rank == 1){
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


/*************************** EDGES ***************************/
if (isEdgeProcess(coords))
{

}

    MPI_Finalize();
    return 0;
}
