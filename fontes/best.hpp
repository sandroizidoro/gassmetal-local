#pragma once

int melhordageracao(float fitnessvet[]){
	float melhor = fitnessvet[0];
	int posicao = 0;
	for (int i = 1; i<tpopulacao; i++)
		 if (melhor > fitnessvet[i]){
			melhor = fitnessvet[i];
			posicao = i;
		 }
	return posicao;
}

// FOUND THE BEST OF GENERATION
int melhordageracaof(float fitnessvet[]){
	  float melhor = fitnessvet[0];
	  int posicao = 0;
	  for (int i = 1; i<ngeracoes; i++)
		   if (melhor > fitnessvet[i]){
			  melhor = fitnessvet[i];
			  posicao = i;
		   }
	   return posicao;
}
// WORST OF GENERATION
int piordageracao(float fitnessvet[]){
	float melhor = fitnessvet[0];
	int posicao = 0;
	for (int i = 1; i<tpopulacao; i++)
		 if (melhor < fitnessvet[i]){
			melhor = fitnessvet[i];
			posicao = i;
		 }
	return posicao;
}
