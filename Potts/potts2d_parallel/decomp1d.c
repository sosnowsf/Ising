#include"decomp1d.h"

/*Function to divide the vector s.th. they differ in size by at most one */
int decomp1d(int n, int p, int myid, int *s, int *e){
        int r = n%p; /*Find division remainder*/
        int d = (n-r)/p;
        if(myid<r){
                *s = myid*(d)+myid;
                *e = (myid+1)*(d)+myid;
        }
        if(myid>=r){
                *s = myid*d+r;
                *e = (myid+1)*d+r-1;
        }
	/*if(myid==0){
		*s=1;
	}*/
return 0;
}

