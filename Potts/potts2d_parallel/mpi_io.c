#include"mpi_io.h"

void save_grid(int g[m][n], int s, int e){
        MPI_Datatype subarray;
        int ndims=2;
        int sizes[2] = {m,n};
        int subsizes[2] = {(e-s+1),n};
        int starts[2] = {s,0};
        MPI_Type_create_subarray(ndims, sizes, subsizes, starts, MPI_ORDER_C, MPI_INT, &subarray);
        MPI_Type_commit(&subarray);

        MPI_File fh;
        MPI_File_open(MPI_COMM_WORLD, "configuration.txt", MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
        MPI_File_set_view(fh, /*s*n*/0, MPI_INT, subarray, "native", MPI_INFO_NULL);
        MPI_File_write_all(fh, &g[s][0]/*[s][0]*/, (e-s+1)*n, /*subarray*/ MPI_INT, MPI_STATUS_IGNORE);
        MPI_File_close(&fh);
        MPI_Type_free(&subarray);
}

void save_grid2d(int g[m][n], int s, int e, int s2, int e2){
        MPI_Datatype subarray;
        int ndims=2;
        int sizes[2] = {m,n};
        int subsizes[2] = {(e-s+1),(e2-s2+1)};
        int starts[2] = {s,s2};
        MPI_Type_create_subarray(ndims, sizes, subsizes, starts, MPI_ORDER_C, MPI_INT, &subarray);
        MPI_Type_commit(&subarray);

        MPI_File fh;
        MPI_File_open(MPI_COMM_WORLD, "configuration.txt", MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
        MPI_File_set_view(fh, /*s*n*/0, MPI_INT, subarray, "native", MPI_INFO_NULL);
        MPI_File_write_all(fh, &g[s][s2], /*1*/(e-s+1)*n, /*subarray*/ MPI_INT, MPI_STATUS_IGNORE);
        MPI_File_close(&fh);
        MPI_Type_free(&subarray);
}

void load_grid(int g[m][n], int s, int e){
        MPI_Datatype subarray;
        int ndims=2;
        int sizes[2] = {m,n};
        int subsizes[2] = {(e-s+1),n};
        int starts[2] = {s,0};
        MPI_Type_create_subarray(ndims, sizes, subsizes, starts, MPI_ORDER_C, MPI_INT, &subarray);
        MPI_Type_commit(&subarray);

        MPI_File fh;
        MPI_File_open(MPI_COMM_WORLD, "configuration.txt", MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
        MPI_File_set_view(fh, 0, subarray, subarray, "native", MPI_INFO_NULL);
        MPI_File_read_all(fh, &g[0][0], 1, subarray, MPI_STATUS_IGNORE);
        MPI_File_close(&fh);
        MPI_Type_free(&subarray);



//      Alternative code to load full grid instead of just relevant part

//      MPI_File fh1;
//      MPI_File_open(MPI_COMM_WORLD, "configuration", MPI_MODE_RDONLY, MPI_INFO_NULL, &fh1);
//      MPI_File_set_view(fh1, 0, MPI_INT, MPI_INT, "native", MPI_INFO_NULL);
//      MPI_File_read_all(fh1, &g[0][0], m*n, MPI_INT, MPI_STATUS_IGNORE);
//      MPI_File_close(&fh1);

}

void load_grid2d(int g[m][n], int s, int e, int s2, int e2){
        MPI_Datatype subarray;
        //printf("%d %d %d %d\n", s, e, s2, e2);
        int ndims=2;
        int sizes[2] = {m,n};
        int subsizes[2] = {(e-s+1),(e2-s2+1)};
        int starts[2] = {s,s2};
        MPI_Type_create_subarray(ndims, sizes, subsizes, starts, MPI_ORDER_C, MPI_INT, &subarray);
        MPI_Type_commit(&subarray);

        MPI_File fh;
        MPI_File_open(MPI_COMM_WORLD, "configuration.txt", MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
        MPI_File_set_view(fh, /*((s*n)+(e2-s2+1))*sizeof(int)*/0, /*MPI_INT*/subarray, /*MPI_INT*/subarray, "native", MPI_INFO_NULL);
        MPI_File_read_all(fh, &g[s][s2], (e-s+1)*n/*(e2-s2+1)*sizeof(int)*/, MPI_INT, MPI_STATUS_IGNORE);
        //MPI_File_read_all(fh, &g[s][s2], (e-s+1)*(e2-s2+1)*sizeof(int), /*subarray*/MPI_INT, MPI_STATUS_IGNORE);
        MPI_File_close(&fh);
        MPI_Type_free(&subarray);

//      Alternative code to load full grid instead of just relevant part

//      MPI_File fh1;
//      MPI_File_open(MPI_COMM_WORLD, "configuration.txt", MPI_MODE_RDONLY, MPI_INFO_NULL, &fh1);
//      MPI_File_set_view(fh1, 0, MPI_INT, MPI_INT, "native", MPI_INFO_NULL);
//      MPI_File_read_all(fh1, &g[0][0], m*n, MPI_INT, MPI_STATUS_IGNORE);
//      MPI_File_close(&fh1);
}

