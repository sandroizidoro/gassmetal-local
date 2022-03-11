#pragma once

float fitnessRef (Atomo *refe, float distanciasref[]){
   float fit = 0.0, aux;
   int pos = 0;
   for(int i=0;i<tam_sitiof;i++)
	  for(int j=i+1;j<tam_sitiof;j++){
		 aux = sqrt(pow(refe[i].x - refe[j].x, 2) + pow(refe[i].y - refe[j].y,2)  +  pow(refe[i].z - refe[j].z,2));
		 distanciasref[pos++] = aux;
		 fit = fit + aux;
	 }
   return fit;
}
