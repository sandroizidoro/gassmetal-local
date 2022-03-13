#python rungass.py reference_pdb  template_site target_pdb reference_atom
#python rungass.py reference_pdb template_site mutations target_pdb reference_atom


import os
import re
from re import template
import subprocess
import sys
import requests
import json
from os.path import exists

# Verify template format and return template size as string
def verify_1x1_template(template):
	listaAA = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
	i = 0
	p = ""
	res = ""
	AAok = False
	seq = False
	while(i < len(template)):
		p = p + template[i]
		if(template[i] == ','):
			res = p.replace(',','')
			p = ''
			if res in listaAA and AAok == False:
				AAok = True
				print(res)
			elif res.isdigit() and AAok == True and seq == False:
				seq = True
				#print(res)
			else:
				#print(template[i])
				#print(i)
				#print("nenhum")
				return "-1"
		if(template[i] == ';'):
			cad = p.replace(';','')
			p = ''
			AAok = False
			seq = False
			if cad.isalpha() and len(cad) == 1:
				print(cad)
			else:
				return "-1"
				#print("nenhum")
				
		i = i + 1
	tamanho_sitio = template.count(';')
	return str(tamanho_sitio)

def verify_1x1_mutation(mutation):
	listaAA = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
	i = 0
	p = ''
	res = ''
	achou = False
	while(i < len(mutation)):
		p = p + mutation[i]
		if(mutation[i] == ','):
			res = p.replace(',','').replace(';','')
			p = ''
			if res in listaAA and achou == False:
				achou = True
				#print(res)
			else:
				#print("errado")
				return False
				#return False
		if( mutation[i] == ';'):
			res = p.replace(',','').replace(';','')
			p = ''
			if res in listaAA and achou == True:
				achou = False
				#print(res)
			else:
				#print("errado")
				return False
		i = i + 1
	return True

def DownloadPdB(code, path_to_save):
	if(re.match("^\d\w\w\w$", code)):
		cmd = "wget  http://www.rcsb.org/pdb/files/" + code + ".pdb -O " + path_to_save 
	else:
		page = requests.get("https://www.alphafold.ebi.ac.uk/api/prediction/"+code.upper())
		print(page)
		try:
			jsonParsed = json.loads(page.text)[0]
		except:
			return 0
		cmd ="wget "+jsonParsed['pdbUrl']+" -O " + path_to_save 
	os.system(cmd)
	cmd = "grep \"requested file is not available\" " + path_to_save  + " | wc -l"
	error = subprocess.check_output(cmd, shell=True)
	if int(error) == 0:
		return 1
	else:
		return 0

def CalculaCentroide(nome_pdb):
    arq_pdb = open(nome_pdb, "r")
    linhas = arq_pdb.readlines()
    arq_pdb.close()

    maior_XYZ = [float(linhas[4][30:38]), float(
        linhas[4][38:46]), float(linhas[4][46:54])]
    menor_XYZ = [float(linhas[4][30:38]), float(
        linhas[4][38:46]), float(linhas[4][46:54])]

    for l in linhas[5:]:
        # busca a menor e a maior coordenada X
        if float(l[30:38]) > maior_XYZ[0]:
            maior_XYZ[0] = float(l[30:38])
        elif float(l[30:38]) < menor_XYZ[0]:
            menor_XYZ[0] = float(l[30:38])

        # busca a menor e a maior coordenada Y
        if float(l[38:46]) > maior_XYZ[1]:
            maior_XYZ[1] = float(l[38:46])
        elif float(l[38:46]) < menor_XYZ[1]:
            menor_XYZ[1] = float(l[38:46])

        # busca a menor e a maior coordenada Z
        if float(l[46:54]) > maior_XYZ[2]:
            maior_XYZ[2] = float(l[46:54])
        elif float(l[46:54]) < menor_XYZ[2]:
            menor_XYZ[2] = float(l[46:54])

    media_X = round((maior_XYZ[0]+menor_XYZ[0])/2, 3)
    media_Y = round((maior_XYZ[1]+menor_XYZ[1])/2, 3)
    media_Z = round((maior_XYZ[2]+menor_XYZ[2])/2, 3)

    centroide = []
    centroide.append("{:.3f}".format(media_X))
    centroide.append("{:.3f}".format(media_Y))
    centroide.append("{:.3f}".format(media_Z))

    return centroide

def prepareRunFolder(run_folder):
	cmd="rm "+run_folder+"*.txt && rm "+run_folder+"*.dat&& rm "+run_folder+"*.pdb"
	os.system(cmd)

def prepareReferenceTxt(reference, template_site, run_folder, atom, mutations):
	if template_site[-1] != ';':
		template_site = template_site + ';'
		template_size = verify_1x1_template(template_site)
		template_site = template_site[:-1]
		if template_size == "-1":
			print("Template out of format!!")
			quit()

			# Check mutation format
	if mutations != "":
		if mutations[-1] != ';' :
			mutations = mutations + ';'
			if verify_1x1_mutation(mutations) == False:
				print("Mutations out of format!!")
				quit()

	if(re.match("^\d\w\w\w$",reference, re.I) or re.match("[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$", reference, re.I)):
		DownloadPdB(reference, run_folder+"refe.pdb")
	else:
		cmd = "cp "+reference+" "+run_folder+"refe.pdb"
		os.system(cmd)
	pdb2txt(run_folder+"refe.pdb",run_folder+"refe_.txt",atom,template_site)
	tamanho_sitio = template_site.count(';') 
	arquivo = open(run_folder + "templates/Templates.txt", "w")
	#alterar o nome da proteina .dat no template
	arquivo.write("Proteina\nrefe_.dat\nTamanhodoSitio\n")
	arquivo.write("%s\n" % tamanho_sitio)
	arquivo.write("Residuos\n")

	for i in range(0, len(template_site)):
		if template_site[i] == ',' or template_site[i] == ';':
			arquivo.write('\n')
		else:
			arquivo.write("%s" % template_site[i])
	arquivo.close()
	
	#gera o arquivo da matriz de substituicao
	#excluir o arquivo padrao que esta na pasta de scripts
	#132l,ASP,SER
	
	if mutations != "":
		arquivo = open(run_folder + "templates/SubstitutionMatrix.txt", "w")
		arquivo.write("refe,")
		for i in range(0, len(mutations)-1):
			if mutations[i] == ";":
				arquivo.write("\nrefe,")
			else:
				arquivo.write("%s" % mutations[i])
		arquivo.close()

def prepareTargetTxt(target, run_folder, atom):
	if(re.match("^\d\w\w\w$",target, re.I) or re.match("[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$", target, re.I)):
		DownloadPdB(target, run_folder+"targ.pdb")
	else:
		cmd = "cp "+target+" "+run_folder+"targ.pdb"
		os.system(cmd)
	pdb2txt(run_folder+"targ.pdb",run_folder+"targ_.txt", atom)

def pdb2txt(pdb,path_to_save, atom, template=""):
	dict_lha = {"ALA":"CB", "ARG":"CZ", "ASN":"CG", "ASP":"CG", "CYS":"SG", "GLN":"CD", "GLU":"CD", "GLY":"CA", "HIS":"CE1", "ILE":"CD1", "LEU":"CD1", "LYS":"NZ", "MET":"CE", "PHE":"CZ", "PRO":"CG", "SER":"OG", "THR":"CG2", "TRP":"CH2", "TYR":"OH", "VAL":"CG1"}

	reference_atom = atom

	



	arq_pdb = open(pdb, 'r')
	pdb_lines = arq_pdb.readlines()


	arq_txt = open(path_to_save, 'w')
    
	nome_proteina = 'NULL'
	ec_proteina = 'NULL'
	unp_proteina = 'NULL'
	reso_proteina = 'NULL'

	for l_pdb in pdb_lines:
		if l_pdb.find('HEADER') != -1:
			nome_proteina = l_pdb[62:66]
		if l_pdb.find('EC:') != -1 and l_pdb.find('COMPND') != -1:
			ec_proteina = l_pdb[15:].strip().replace(';','')
		if l_pdb.find('UNP') != -1 and l_pdb.find('DBREF') != -1:
			unp_proteina = l_pdb[33:39]
		if l_pdb.find('RESOLUTION.') != -1 and l_pdb.find('REMARK') != -1:
			reso_proteina = l_pdb[26:31].strip()
   
	regex = re.match("^\d\w\w\w$",nome_proteina, re.I)
	if not regex:
		nome_proteina = 'NULL'
	arq_txt.write(nome_proteina)
	arq_txt.write('\n')
	arq_txt.write(ec_proteina)
	arq_txt.write('\n')
	arq_txt.write(unp_proteina)
	arq_txt.write('\n')
	arq_txt.write(reso_proteina)
	arq_txt.write('\n')
	arq_pdb.close()

    # loops que percorrem todos os atomos de CA e verifica se existem dados duplicados
	for i in range(0,len(pdb_lines)):
		encontrou = False
		if pdb_lines[i][0:4] == 'ATOM':
			res = pdb_lines[i][17:20]
           #print(res)
			if reference_atom == 'CA':
				ra = 'CA'
			else:
				ra = dict_lha[res]
			if pdb_lines[i][12:16].strip() == ra:
				for j in range(i+1,len(pdb_lines)):
					if pdb_lines[j][0:4] == 'ATOM' and pdb_lines[j][12:16].strip() == ra:
						if pdb_lines[i][21:22] == pdb_lines[j][21:22] and pdb_lines[i][22:26].strip() == pdb_lines[j][22:26].strip():
							encontrou = True  
				if encontrou == False:
					arq_txt.write("%s" % pdb_lines[i])
       
	arq_txt.close()


def generateDats(run_folder):
	cmd="cd "+run_folder+" && ./filtroC_Dados"
	os.system(cmd)

def moveFiles(run_folder):
	cmd ="mv "+run_folder+"refe_.dat "+run_folder+"templates/refe_.dat"
	os.system(cmd)
	cmd ="mv "+run_folder+"targ_.dat "+run_folder+"target/protein.dat"
	os.system(cmd)
	cmd ="cp scripts/GA_Conf.txt "+run_folder+"templates/GA_Conf.txt"
	os.system(cmd)
	cmd ="cp scripts/Substitution_matrix.txt "+run_folder+"templates/SubstitutionMatrix.txt"
	os.system(cmd)

def setup(run_folder):
	if not exists(run_folder):
		os.system("mkdir "+run_folder)
	if not exists(run_folder+"target/"):
		os.system("mkdir "+run_folder+"target/")
	if not exists(run_folder+"templates/"):
		os.system("mkdir "+run_folder+"templates/")
	if not exists(run_folder+"filtroC_Dados"):
		os.system("g++ filtroC_Dados.cpp -o "+run_folder+"filtroC_Dados")
	if not exists(run_folder+"LISTAENZIMAS.TXT"):
		os.system("echo -e 'refe_.txt\ntarg_.txt' > "+run_folder+"LISTAENZIMAS.TXT")
	if not exists("thread_GASS"):
		os.system("g++ thread_GASS.cpp -o thread_GASS")
	
	
	
if __name__=="__main__":
	run_folder="run/"
	mutations=""
	if len(sys.argv) == 5:
		reference = sys.argv[1]
		template_site = sys.argv[2]
		target = sys.argv[3]
		atom = sys.argv[4]
	
	elif len(sys.argv) == 6:
		reference = sys.argv[1]
		template_site = sys.argv[2]
		mutations = sys.argv[3]
		target = sys.argv[4]
		atom = sys.argv[5]
	
	else:
		print("Error: Wrong number of arguments")
		quit()
	
	setup(run_folder)
	prepareRunFolder(run_folder)
	prepareReferenceTxt(reference, template_site, run_folder, atom, mutations)
	prepareTargetTxt(target, run_folder, atom)
	generateDats(run_folder)
	centroide = CalculaCentroide(run_folder+"targ_.txt")
	moveFiles(run_folder)
	
	process = subprocess.Popen(['./thread_GASS',run_folder+'target/',run_folder+'templates/',centroide[0], centroide[1], centroide[2]])
	process.communicate()
	cmd="mv "+run_folder+"target/ActiveSitesFound.txt ActiveSitesFound.txt"
	os.system(cmd)
