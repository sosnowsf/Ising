#include"mpi_io.h"

void save_grid(int g[x][y][z], int s, int e){
        MPI_Datatype subarray;
        int ndims=3;
        int sizes[3] = {x,y,z};
        int subsizes[3] = {(e-s+1),y,z};
        int starts[3] = {s,0,0};
        MPI_Type_create_subarray(ndims, sizes, subsizes, starts, MPI_ORDER_C, MPI_INT, &subarray);
        MPI_Type_commit(&subarray);

        MPI_File fh;
        MPI_File_open(MPI_COMM_WORLD, "configuration.txt", MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
        MPI_File_set_view(fh, 0, MPI_INT, subarray, "native", MPI_INFO_NULL);
        MPI_File_write_all(fh, &g[s][0][0], (e-s+1)*y*z, MPI_INT, MPI_STATUS_IGNORE);
        MPI_File_close(&fh);
        MPI_Type_free(&subarray);
}

void save_grid2d(int g[x][y][z], int s, int e, int s2, int e2){
        MPI_Datatype subarray;
        int ndims=3;
        int sizes[3] = {x,y,z};
        int subsizes[3] = {(e-s+1),(e2-s2+1),z};
        int starts[3] = {s,s2,0};
        MPI_Type_create_subarray(ndims, sizes, subsizes, starts, MPI_ORDER_C, MPI_INT, &subarray);
        MPI_Type_commit(&subarray);

        MPI_File fh;
        MPI_File_open(MPI_COMM_WORLD, "configuration.txt", MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
        MPI_File_set_view(fh, 0, MPI_INT, subarray, "native", MPI_INFO_NULL);
        MPI_File_write_all(fh, &g[s][s2][0], (e-s+1)*y*z, MPI_INT, MPI_STATUS_IGNORE);
        MPI_File_close(&fh);
        MPI_Type_free(&subarray);
}

void save_grid3d(int g[x][y][z], int s, int e, int s2, int e2, int s3, int e3){
        MPI_Datatype subarray;
        int ndims=3;
        int sizes[3] = {x,y,z};
        int subsizes[3] = {(e-s+1),(e2-s2+1),(e3-s3+1)};
        int starts[3] = {s,s2,s3};
        MPI_Type_create_subarray(ndims, sizes, subsizes, starts, MPI_ORDER_C, MPI_INT, &subarray);
        MPI_Type_commit(&subarray);

        MPI_File fh;
        MPI_File_open(MPI_COMM_WORLD, "configuration.txt", MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
        MPI_File_set_view(fh, 0, MPI_INT, subarray, "native", MPI_INFO_NULL);
        //MPI_File_write_all(fh, &g[s][s2][s3], x*y*z, MPI_INT, MPI_STATUS_IGNORE);
	MPI_File_write_all(fh, &g[0][0][0], 1, subarray, MPI_STATUS_IGNORE);
        MPI_File_close(&fh);
        MPI_Type_free(&subarray);
}


void load_grid(int g[x][y][z], int s, int e){
        MPI_Datatype subarray;
        int ndims=3;
        int sizes[3] = {x,y,z};
        int subsizes[3] = {(e-s+1),y,z};
        int starts[3] = {s,0,0};
        MPI_Type_create_subarray(ndims, sizes, subsizes, starts, MPI_ORDER_C, MPI_INT, &subarray);
        MPI_Type_commit(&subarray);

        MPI_File fh;
        MPI_File_open(MPI_COMM_WORLD, "configuration.txt", MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
        MPI_File_set_view(fh, 0, subarray, subarray, "native", MPI_INFO_NULL);
        MPI_File_read_all(fh, &g[0][0][0], 1, subarray, MPI_STATUS_IGNORE);
        MPI_File_close(&fh);
        MPI_Type_free(&subarray);



//      Alternative code to load full grid instead of just relevant part

//      MPI_File fh1;
//      MPI_File_open(MPI_COMM_WORLD, "configuration", MPI_MODE_RDONLY, MPI_INFO_NULL, &fh1);
//      MPI_File_set_view(fh1, 0, MPI_INT, MPI_INT, "native", MPI_INFO_NULL);
//      MPI_File_read_all(fh1, &g[0][0], m*n, MPI_INT, MPI_STATUS_IGNORE);
//      MPI_File_close(&fh1);

}

void load_grid2d(int g[x][y][z], int s, int e, int s2, int e2){
        MPI_Datatype subarray;
        int ndims=3;
        int sizes[3] = {x,y,z};
        int subsizes[3] = {(e-s+1),(e2-s2+1),z};
        int starts[3] = {s,s2,0};
        MPI_Type_create_subarray(ndims, sizes, subsizes, starts, MPI_ORDER_C, MPI_INT, &subarray);
        MPI_Type_commit(&subarray);

        MPI_File fh;
        MPI_File_open(MPI_COMM_WORLD, "configuration.txt", MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
        MPI_File_set_view(fh, 0, subarray, subarray, "native", MPI_INFO_NULL);
        MPI_File_read_all(fh, &g[s][s2][0], (e-s+1)*y*z, MPI_INT, MPI_STATUS_IGNORE);
	//MPI_File_read_all(fh, &g[0][0][0], 1/*(e2-s2+1)*sizeof(int)*/, subarray, MPI_STATUS_IGNORE);
        //MPI_File_read_all(fh, &g[s][s2], (e-s+1)*(e2-s2+1)*sizeof(int), /*subarray*/MPI_INT, MPI_STATUS_IGNORE);
        MPI_File_close(&fh);
        MPI_Type_free(&subarray);
}

void load_grid3d(int g[x][y][z], int s, int e, int s2, int e2, int s3, int e3){
        MPI_Datatype subarray;
        int ndims=3;
        int sizes[3] = {x,y,z};
        int subsizes[3] = {(e-s+1),(e2-s2+1),(e3-s3+1)};
        int starts[3] = {s,s2,s3};
        MPI_Type_create_subarray(ndims, sizes, subsizes, starts, MPI_ORDER_C, MPI_INT, &subarray);
        MPI_Type_commit(&subarray);

        MPI_File fh;
        MPI_File_open(MPI_COMM_WORLD, "configuration.txt", MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
        MPI_File_set_view(fh, 0, subarray, subarray, "native", MPI_INFO_NULL);
        //MPI_File_read_all(fh, &g[s][s2][s3], x*y*z, MPI_INT, MPI_STATUS_IGNORE);
        MPI_File_read_all(fh, &g[0][0][0], 1, subarray, MPI_STATUS_IGNORE);
        //MPI_File_read_all(fh, &g[s][s2], (e-s+1)*(e2-s2+1)*sizeof(int), /*subarray*/MPI_INT, MPI_STATUS_IGNORE);
        MPI_File_close(&fh);
        MPI_Type_free(&subarray);
}

