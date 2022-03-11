bool AbrirArquivo(ifstream &ifs, char nomeArquivo[]){
	ifs.open(nomeArquivo, ios::binary);
	return ifs.is_open();
}
