#include<stdio.h>


int main(){

  int a = -1;

  if ( a == -1 ){

    printf("a=%d, line=%d\n", a, __LINE__);
    a++;
  }
  else if ( a == 0 ){

    printf("a=%d, line=%d\n", a, __LINE__);
    a++;
  }
  else{

    printf("a=%d, line=%d\n", a, __LINE__);

  }


  return 0;
}
