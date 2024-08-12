#include"main.h"
#include <time.h>

int main(int argc, char const* argv[]){
    srand(time(NULL) + clock());
    int threadNumber = 1 ;
    int size_sbox = 3 ;
    int rows = 1 ;

  for (int i = 0; i < argc; i++) {
    if (!strcmp(argv[i], "-n")) size_sbox = atoi(argv[i + 1]);
    if (!strcmp(argv[i], "-t")) threadNumber = atoi(argv[i + 1]);
    if (!strcmp(argv[i], "-r")) rows = atoi(argv[i + 1]);
  }
   
	lat_sbox(size_sbox, threadNumber, rows);
	
	return 0;
}


