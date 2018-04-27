#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define UNLIMIT

struct myStringStruct {
  char qstring[128];
};

int compare(const void *elem1, const void *elem2)
{
  int result;
  
  result = strcmp((*((struct myStringStruct *)elem1)).qstring, (*((struct myStringStruct *)elem2)).qstring);

  return (result < 0) ? 1 : ((result == 0) ? 0 : -1);
}


int
main(int argc, char *argv[]) {
 
  struct myStringStruct *array = (struct myStringStruct*)malloc(1*sizeof(struct myStringStruct));
  if (array == NULL) {
    printf("malloc failed.\n");
    exit(1);
  }

  FILE *fp;
  int i,count=0;
  
  if (argc<2) {
    fprintf(stderr,"Usage: qsort_small <file>\n");
    exit(-1);
  }
  else {
    fp = fopen(argv[1],"r");
    while(fscanf(fp, "%s", &array[count].qstring) == 1) {
      count++;
      struct myStringStruct *reallocated_array = realloc(array, (count+1)*sizeof(struct myStringStruct));
      if (reallocated_array == NULL) {
        printf("malloc failed.\n");
        exit(1);
      }
      array = reallocated_array;
    }
  }
  printf("\nSorting %d elements.\n\n",count);
  qsort(array,count,sizeof(struct myStringStruct),compare);
  
  for(i=0;i<count;i++)
    printf("%s\n", array[i].qstring);
  
  free(array);

  return 0;
}
