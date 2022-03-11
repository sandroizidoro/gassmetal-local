#pragma once

int partition(Atomo **tempop, float a[], int l, int r) {
   float pivot, t;
   int i, j;
   Atomo z;
   pivot = a[l];
   i = l; j = r+1;
   while(1) {
	do ++i; while( a[i] <= pivot && i <= r );
	do --j; while( a[j] > pivot );
	if( i >= j ) break;
	t = a[i]; a[i] = a[j]; a[j] = t;
	for(int k = 0; k<tam_sitiof; k++){
	   z = tempop[i][k];
	   tempop[i][k] = tempop[j][k];
	   tempop[j][k] = z;
	}
   }
   t = a[l]; a[l] = a[j]; a[j] = t;
   for(int k = 0; k<tam_sitiof; k++){
	   z = tempop[l][k];
	   tempop[l][k] = tempop[j][k];
	   tempop[j][k] = z;
   }
   return j;
}

void quickSort(Atomo **tempop, float a[], int l, int r){
   int j;
   if( l < r ){
	 j = partition( tempop, a, l, r);
	 quickSort(tempop, a, l, j-1);
	 quickSort(tempop, a, j+1, r);
   }
}

