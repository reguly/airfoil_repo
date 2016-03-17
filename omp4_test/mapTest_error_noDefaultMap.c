#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

int main(){
  int N = 10000000;
  int *input, *res,*tmp;
  res   = malloc(N*sizeof(int));
  input = malloc(N*sizeof(int));
  tmp   = malloc(N*sizeof(int));
/* //Igy lehet helyesen target data env-et letrehozni. (minden hasznalt valtozora map..) 
  #pragma omp target data map(to: input[:N]) map(from: res[:N]) map(alloc:tmp[:N])
  {
    printf("Default device is %d\nNumber of devices: %d\r\n",omp_get_default_device(),omp_get_num_devices());
    #pragma omp target
   {
 #pragma omp parallel for
    for(int i = 0; i < N; ++i){
      tmp[i]=i;
    }
}
    #pragma omp target
{
    #pragma omp parallel for
     for(int i = 0; i < N; ++i){
      res[i] = 2*tmp[i];
    }
}    
  }

    #pragma omp target
    #pragma omp parallel for
     for(int i = 0; i < N; ++i){
      res[i] = 2*tmp[i];
    } 
*/
///Ez a resz nem mukodik nincs alapertelmezett map (hivatalosan tofromnak kene lennie).
  for(int i = 0; i < N;++i){
    res[i] = i;
  }
  #pragma omp target map(to:N)
#pragma omp parallel for
  for(int i = 0; i < N; ++i){
    res[i]*=i;
  }
printf("%d\n",res[2]);
  return 0;
}
