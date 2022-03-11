#pragma once

int funcaotorneio(float fitnessvet[], int escolhidos_torneio[]){
	 int escolhido = escolhidos_torneio[0];
	 float menor = fitnessvet[escolhido];
	 for(int i=1;i<torneio;i++)
		 if(fitnessvet[escolhidos_torneio[i]]< menor){
			 menor = fitnessvet[escolhidos_torneio[i]];
			 escolhido = escolhidos_torneio[i];
		 }
return escolhido;
}