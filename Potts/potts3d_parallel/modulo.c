#include"modulo.h"

int modulo(int a, int b){
        int r = a%b;
        if(r<0) r+=b;
return r;
}

