#pragma once

void mutation (Atomo **pop, int *numero_aminos, Atomo **repositorio){
   char find[5];
   int pos, j;
   for(int i=0; i < tpopulacao; i++){
	  	if(rand()%100<=txmutacao){
			for(j=0; j<tam_sitiof; j++){
		   		if(rand()%100<=50){
			 		strcpy(find, pop[i][j].amino);
			 		pos = posicaonalista(teste, find);
			 		pop[i][j] = repositorio[pos][rand() % numero_aminos[pos]];
		   		}
			}
	  	}	
   	}
}
