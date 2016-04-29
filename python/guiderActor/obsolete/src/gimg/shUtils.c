

#include <stdio.h>
#include <stdlib.h>
#include "shLegacy.h"



void shError(char *fmt, ...){
  printf("%s \n", fmt);
}

void shDebug(int level, char *fmt, ...){
   if (level > 0) {
      printf("%d: %s \n", level, fmt);
   }
}

void *shMalloc(size_t size){
   return malloc(size);
}

void shFree(void *ptr){
  free(ptr);
}

