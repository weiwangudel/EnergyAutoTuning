 #include <stdio.h>
 #include <stdlib.h>
 
 int main()
 {
     int i,j,temp[10];
     int a, b;
     for (i=0; i<=9; i++)     // only this loop considered SCoP by PolyOpt
      temp[i] = i*2;
 #pragma scop
     for(i=0;i<9;i++) {
       a = temp[i];
       b = temp[i+1];
       for (j=a; j<b; j++)    // only this loop considerred SCoP by PolyOpt, not the whole loop 
         temp[i] = b;
     }
 #pragma endscop
     for (i=0; i<9; i++)
	printf("temp[i]:%d\n", temp[i]);
     return 1;
 }

