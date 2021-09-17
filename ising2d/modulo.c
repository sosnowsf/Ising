#include"modulo.h"

//Modulo function for handling -1(mod)n
int modulo(int a, int b){
        int r = a%b;
        if(r<0) r+=b;
return r;
}
