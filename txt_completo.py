'''
Script que gera um arquivo .txt completo a partir de um arquivo .pdb.
Salva as informacoes de nome, ec number, uniprot e resolucao.
Todos os atomos de CA de todos os residuos da proteina sao salvos.
Ao final do script, atualiza o arquivo 'LISTAENZIMAS.TXT' com os nomes .txt que serao usados pelo ./filtroC_dados
'''

import os
import sys


list_pdbs = []
if len(sys.argv) == 1:
    caminho = ""
    arq_dir = os.listdir()

elif len(sys.argv) == 2:
    caminho = sys.argv[1]
    arq_dir = os.listdir(caminho)

elif len(sys.argv) == 3:
    caminho = sys.argv[1]
    arq_dir = os.listdir(caminho)
    list_pdbs.append(sys.argv[2])

else:
    print("Argumentos invalidos")
    exit()

if len(sys.argv) != 3:
    for a in arq_dir:
        if '.pdb' in a and len(a) == 8:
            list_pdbs.append(a)

if not list_pdbs:
    print("Nao existem arquivos PDBs no diretorio")
    exit()
cont = 1
for nome_pdb in list_pdbs:
    print("{}: {}".format(cont, nome_pdb))
    cont = cont + 1
    arq_pdb = open(caminho + nome_pdb, 'r')
    linhas_pdb = arq_pdb.readlines()
    linhas_pdb = [linha.replace('\n', '') for linha in linhas_pdb]
    nome_txt = nome_pdb[:4] + '+.txt'
    arq_txt = open(caminho + nome_txt, 'w')
    ec_proteina = 'NULL'
    unp_proteina = 'NULL'
    for l_pdb in linhas_pdb:
        if l_pdb.find('HEADER') != -1:
            nome_proteina = l_pdb.strip().split()[-1]
        if l_pdb.find('EC:') != -1 and l_pdb.find('COMPND') != -1:
            ec_proteina = l_pdb.strip().split()[-1].replace(';', '')
        if l_pdb.find('UNP') != -1 and l_pdb.find('DBREF') != -1:
            unp_proteina = l_pdb.strip().split()[6]
        if l_pdb.find('RESOLUTION.') != -1 and l_pdb.find('REMARK') != -1:
            reso_proteina = l_pdb.strip().split()[3]
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
    for i in range(0,len(linhas_pdb)):
        encontrou = False
        aux_i = linhas_pdb[i].strip().split()
        if aux_i[0] == 'ATOM' and aux_i[2] == 'CA':
            for j in range(i+1,len(linhas_pdb)):
                aux_j = linhas_pdb[j].strip().split() 
                if aux_j[0] == 'ATOM' and aux_j[2] == 'CA' :
                    if aux_i[4] == aux_j[4] and aux_i[5] == aux_j[5]:
                        encontrou = True  
            if encontrou == False:
                if(len(aux_i[3]) > 3):
                    linhas_pdb[i] = linhas_pdb[i][:16] + " " + linhas_pdb[i][17:]
                arq_txt.write("%s\n" % linhas_pdb[i])
        
    arq_txt.close()

if len(sys.argv) != 1:
    arq_dir = os.listdir(caminho)
else:
    arq_dir = os.listdir()
list_txt = []
for a in arq_dir:
    if ('_.txt' in a or '+.txt' in a) and len(a) == 9:
        list_txt.append(a)
arq_LISTAENZIMAS = open(caminho + 'LISTAENZIMAS.TXT','w')
for txt in list_txt:
    arq_LISTAENZIMAS.write("{}\n" .format(txt))
arq_LISTAENZIMAS.close()