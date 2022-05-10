import argparse

my_parser = argparse.ArgumentParser(description='Counting of bases')

# add arguments
my_parser.add_argument('-f')
my_parser.add_argument('-o')
files = my_parser.parse_args()

# dictionary implementation
training_set=open(files.f,'r')
text=training_set.readlines()

total=len(text);
count={'A': [0,0,0,0,0,0], 'C': [0,0,0,0,0,0], 'G':[0,0,0,0,0,0], 'T':[0,0,0,0,0,0]};
for line in text:
	line=line.replace('\n','')
	pos=0;
	for letter in line:
		count[letter][pos]=count[letter][pos]+1;
		pos=pos+1;
training_set.close()

output_file=open(files.o,'w')
output_file.write('Position\t1\t2\t3\t4\t5\t6\n');
for letter in count :
	output_file.write('%s' % (letter));
	for num in range(0,len(count[letter])) :
		output_file.write('\t%d/%d' %(count[letter][num],total));
	output_file.write('\n');

output_file.close()
