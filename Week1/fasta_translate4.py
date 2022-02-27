import sys

a=open('codon_table_hard.txt','r')

# dictionnary implementation
text=a.readlines()
table={}

for line in text:
    [name,aa,last]=line.split('\t')
    last=last.replace('\n','')
    codon=last.split(',')
    for i in codon:
        table[i]=aa

a.close()

genes = []
descriptions =[]

# genes and descriptions storage
gene = "" #string of the complete gene
line = sys.stdin.readline()
line = line.replace('\n', '')
line = line.replace('>', '')
descriptions.append(line)
for line in sys.stdin :
    if (line[0] == '>') :
        gene = gene.replace('\n', '')
        line = line.replace('\n', '')
        line = line.replace('>', '')
        descriptions.append(line)
        genes.append(gene)
        gene = ""

    else:
        gene += line
        gene = gene.replace('\n', '')

genes.append(gene)

# genes translation
for i in range(0,len(genes)) :
    print(f"\n>{descriptions[i]}")
    aminoacid = ''
    missing = []

    for j in range(0, len(genes[i]), 3):
        if (j+3) > len(genes[i]):
            continue
        else:
            codon = genes[i][j:j+3]
            if(codon in table) :
               aminoacid += table[codon]
            else :
               sys.stderr.write(f"Invalid: {codon}\n")
				
    print(aminoacid, '\n')  
