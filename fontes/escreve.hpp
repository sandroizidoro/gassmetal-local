#pragma once

void EscreveAtomo(Atomo atomo){
cout << atomo.amino << "  " << atomo.atomo_ID << "  " << atomo.x << "  " << atomo.y << "  " << atomo.z << "  " << atomo.cadeia << "  -  ";
}
inline int MostraQuadrante(Atomo amino, float *centroide){
	if(amino.x > centroide[0] && amino.y > centroide[1]){ //primeiro quadrante
		return 1;
	}
	else if(amino.x > centroide[0] && amino.y < centroide[1]){
		return 2;
	}
	else if(amino.x < centroide[0] && amino.y < centroide[1]){
		return 3;
	}
	else{
		return 4;
	}   
}
