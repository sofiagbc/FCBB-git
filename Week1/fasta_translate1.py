import sys


line = sys.stdin.readline()

text = sys.stdin.readlines()
val = ''
for i in range(0,len(text)):
    val += text[i]

val = val.replace('\n', '')
line = line.replace('\n', '')
line = line.replace('>', '')
sequence = {line : val }


print(sequence)
