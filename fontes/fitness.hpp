#pragma once

void fitness(Atomo **p, float fitnessvet[], float distanciasitio[], float distanciasref[], int tam_dist){
   float fit = 0.0, aux;
   int pos;
   for(int k=0; k < tpopulacao; k++){
	  pos = 0;
	  for(int i=0; i<tam_sitiof; i++)
		 for(int j=i+1; j<tam_sitiof; j++){
			aux = sqrt(pow(p[k][i].x - p[k][j].x, 2)  + pow(p[k][i].y - p[k][j].y, 2)  + pow(p[k][i].z - p[k][j].z, 2));
			distanciasitio[pos++] = aux;
		 }
	  for(int i=0;i<tam_dist;i++)
		 fit = fit + fabs(distanciasitio[i]-distanciasref[i]);
	  fitnessvet [k]= fit;
	  fit = 0.0;

	  for(int j=0; j<tam_sitiof; j++)
		  for(int i=j+1; i < tam_sitiof; i++)
			 //if (p[k][j].atomo_ID == p[k][i].atomo_ID){
		 if ((p[k][j].atomo_ID == p[k][i].atomo_ID) && (p[k][j].cadeia == p[k][i].cadeia)){
				  fitnessvet[k] = fitnessvet[k] + 100.0;
				  i = tam_sitiof;
				  j = tam_sitiof;
			  }
   }
}
