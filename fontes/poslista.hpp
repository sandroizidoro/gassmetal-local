#pragma once

int posicaonalista (char teste[20][4], char find[5]){
	  for(int i=0; i<20; i++)
		 if (!strcmp(find, teste[i]))
			return i;
	return -1;
}
