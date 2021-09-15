#include"exchange.h"

void exchange(int g[x][y][z], int s, int e, int nbrtop, int nbrbot, MPI_Comm comm){
        MPI_Request reqs[4];
        MPI_Datatype slice;
        MPI_Type_vector(y*z, 1, 1, MPI_INT, &slice);
        MPI_Type_commit(&slice);

        MPI_Irecv(&g[modulo(s-1,x)][0][0], 1, slice, nbrtop, 0, comm, &reqs[0]);
        MPI_Irecv(&g[modulo(e+1,x)][0][0], 1, slice, nbrbot, 0, comm, &reqs[1]);
        MPI_Isend(&g[modulo(e,x)][0][0], 1, slice, nbrbot, 0, comm, &reqs[2]);
        MPI_Isend(&g[modulo(s,x)][0][0], 1, slice, nbrtop, 0, comm, &reqs[3]);


//      MPI_Irecv(&g[modulo(s-1,x)][0][0], 100, MPI_INT, nbrtop, 0, comm, &reqs[0]);
//      MPI_Irecv(&g[modulo(e+1,x)][0][0], 100, MPI_INT, nbrbot, 0, comm, &reqs[1]);
//      MPI_Isend(&g[modulo(e,x)][0][0], 100, MPI_INT, nbrbot, 0, comm, &reqs[2]);
//      MPI_Isend(&g[modulo(s,x)][0][0], 100, MPI_INT, nbrtop, 0, comm, &reqs[3]);

        MPI_Waitall(4, reqs, MPI_STATUSES_IGNORE);
        MPI_Type_free(&slice);
}

//exchange along rows and cols
void exchange2d(int g[x][y][z], int s, int e, int s2, int e2, int nbrtop, int nbrbot, int nbrleft, int nbrright, MPI_Comm comm){
        MPI_Request reqs[8];
        MPI_Datatype slice1;
        MPI_Type_vector((e-s), z, y*z, MPI_INT, &slice1);
        MPI_Type_commit(&slice1);

        //Send to left and right nbrs
        MPI_Irecv(&g[s][modulo(e2+1,x)][0], 1, slice1, nbrright, 0, comm, &reqs[0]);
        MPI_Irecv(&g[s][modulo(s2-1,x)][0], 1, slice1, nbrleft, 1, comm, &reqs[1]);
        MPI_Isend(&g[s][modulo(s2,x)][0], 1, slice1, nbrleft, 0, comm, &reqs[2]);
        MPI_Isend(&g[s][modulo(e2,x)][0], 1, slice1, nbrright, 1, comm, &reqs[3]);

        //Send to top and bot nbrs
        MPI_Datatype slice2;
        MPI_Type_vector((e2-s2+1)*z, 1, 1, MPI_INT, &slice2);
        MPI_Type_commit(&slice2);

        MPI_Irecv(&g[modulo(e+1,x)][s2][0], 1, slice2, nbrtop, 2, comm, &reqs[4]);
        MPI_Irecv(&g[modulo(s-1,x)][s2][0], 1, slice2, nbrbot, 3, comm, &reqs[5]);
        MPI_Isend(&g[modulo(s,x)][s2][0], 1, slice2, nbrbot, 2, comm, &reqs[6]);
        MPI_Isend(&g[modulo(e,x)][s2][0], 1, slice2, nbrtop, 3, comm, &reqs[7]);

        MPI_Waitall(8, reqs, MPI_STATUSES_IGNORE);
        MPI_Type_free(&slice1);
        MPI_Type_free(&slice2);
}

//exchange along all three dimensions
void exchange3d(int g[x][y][z], int s, int e, int s2, int e2, int s3, int e3, int nbrtop, int nbrbot, int nbrleft, int nbrright, int nbrfront, int nbrback, MPI_Comm comm){
        MPI_Request reqs[12];
        MPI_Datatype slice1;
        //MPI_Type_vector((e-s)*z, 1, z, MPI_INT, &slice1);
        MPI_Type_vector((e-s+1), (e3-s3+1), y*z, MPI_INT, &slice1);
        MPI_Type_commit(&slice1);

        //Send to left and right nbrs
        MPI_Irecv(&g[s][modulo(e2+1,x)][s3], 1, slice1, nbrright, 0, comm, &reqs[0]);
        MPI_Irecv(&g[s][modulo(s2-1,x)][s3], 1, slice1, nbrleft, 1, comm, &reqs[1]);
        MPI_Isend(&g[s][modulo(s2,x)][s3], 1, slice1, nbrleft, 0, comm, &reqs[2]);
        MPI_Isend(&g[s][modulo(e2,x)][s3], 1, slice1, nbrright, 1, comm, &reqs[3]);

        //Send to top and bot nbrs
        MPI_Datatype slice2;
        //MPI_Type_vector((e2-s2+1)*(e3-s3-1)*z, 1, 1, MPI_INT, &slice2);
        MPI_Type_vector((e2-s2+1), (e3-s3+1), z, MPI_INT, &slice2);
        MPI_Type_commit(&slice2);

        MPI_Irecv(&g[modulo(e+1,x)][s2][s3], 1, slice2, nbrtop, 2, comm, &reqs[4]);
        MPI_Irecv(&g[modulo(s-1,x)][s2][s3], 1, slice2, nbrbot, 3, comm, &reqs[5]);
        MPI_Isend(&g[modulo(s,x)][s2][s3], 1, slice2, nbrbot, 2, comm, &reqs[6]);
        MPI_Isend(&g[modulo(e,x)][s2][s3], 1, slice2, nbrtop, 3, comm, &reqs[7]);

        //Send front and back nbrs
        MPI_Datatype slice3_1;
        MPI_Datatype slice3_2;
        MPI_Datatype slice3_3;
        MPI_Datatype slice3_4;
        int ndims=3;
        int sizes[3] = {x,y,z};
        int subsizes[3] = {(e-s-1),(e2-s2-1), 1};//{(e-s+1), (e2-s2+1), (e3-s3+1)};
        int starts1[3] = {(s+1),(s2+1),s3};//{(s+1), (s2+1), s3};
        int starts2[3] = {(s+1), (s2+1), e3};
        int starts3[3] = {(s+1), (s2+1), (modulo(e3+1,x))};
        int starts4[3] = {(s+1), (s2+1), (modulo(s3-1,x))};
        MPI_Type_create_subarray(ndims, sizes, subsizes, starts1, MPI_ORDER_C, MPI_INT, &slice3_1);
        MPI_Type_create_subarray(ndims, sizes, subsizes, starts2, MPI_ORDER_C, MPI_INT, &slice3_2);
        MPI_Type_create_subarray(ndims, sizes, subsizes, starts3, MPI_ORDER_C, MPI_INT, &slice3_3);
        MPI_Type_create_subarray(ndims, sizes, subsizes, starts4, MPI_ORDER_C, MPI_INT, &slice3_4);
        MPI_Type_commit(&slice3_1);
        MPI_Type_commit(&slice3_2);
        MPI_Type_commit(&slice3_3);
        MPI_Type_commit(&slice3_4);

        MPI_Irecv(&g[0][0][0], 1, slice3_3, nbrfront, 4, comm, &reqs[8]);
        MPI_Irecv(&g[0][0][0], 1, slice3_4, nbrback, 5, comm, &reqs[9]);
        MPI_Isend(&g[0][0][0], 1, slice3_1, nbrback, 4, comm, &reqs[10]);
        MPI_Isend(&g[0][0][0], 1, slice3_2, nbrfront, 5, comm, &reqs[11]);

        MPI_Waitall(12, reqs, MPI_STATUSES_IGNORE);
        MPI_Type_free(&slice1);
        MPI_Type_free(&slice2);
        MPI_Type_free(&slice3_1);
        MPI_Type_free(&slice3_2);
        MPI_Type_free(&slice3_3);
        MPI_Type_free(&slice3_4);
}

