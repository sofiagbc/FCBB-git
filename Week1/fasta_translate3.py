import sys
    
table={'ATT':'I', 'ATC':'I', 'ATA':'I', 'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L', 'TTA':'L', 'TTG':'L', 'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V', 'TTT':'F', 'TTC':'F', 'ATG':'M', 'TGT':'C', 'TGC':'C', 'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A', 'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G', 'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P', 'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T', 'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S', 'AGT':'S', 'AGC':'S', 'TAT':'Y', 'TAC':'Y', 'TGG':'W', 'CAA':'Q', 'CAG':'Q', 'AAT':'N', 'AAC':'N', 'CAT':'H', 'CAC':'H', 'GAA':'E', 'GAG':'E', 'GAT':'D', 'GAC':'D', 'AAA':'K', 'AAG':'K', 'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R', 'AGA':'R', 'AGG':'R', 'TAA':'*', 'TAG':'*', 'TGA':'*'}    


genes = []
descriptions =[]

gene = "" #string of the complete gene
line = sys.stdin.readline()
line = line.replace('\n', '')
line = line.replace('>', '')
descriptions.append(line)
for line in sys.stdin : #make it general
	if (line[0] == '>') :
		gene = gene.replace('\n', '')
		line = line.replace('\n', '')
		line = line.replace('>', '')
		descriptions.append(line)
		genes.append(gene)
		gene = "" #string of the complete gene

	else:
		gene += line
		gene = gene.replace('\n', '')

genes.append(gene)

# genes translation
for i in range(0,len(genes)) :
	print(f"\n>{descriptions[i]}")
	aminoacid = ''

	for j in range(0, len(genes[i])-3, 3):
		if (j+3) > len(genes[i]):
			continue
		else:
			codon = genes[i][j:j+3]
			aminoacid += table[codon]
			
	print(aminoacid, '\n')
