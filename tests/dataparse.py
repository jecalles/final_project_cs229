from Bio import SeqIO
from context import src
from src import dataparse

# get test data
filename ='../res/GCF_000146045.2_R64_genomic.gbff'
outfile='../res/out.pickle'
records = [record for record in SeqIO.parse(filename, 'genbank')]
chr3 = records[2]
gene = chr3.features[42]

# tests
print('testing: dataparse.__init__() without arguments')
x = dataparse()
# print('testing: dataparse.__init__() with infile specified')
# x = dataparse(infile=filename)
# x = dataparse(filename)
# print('testing: dataparse.__init__() with parsing')
# print('get_genes=True')
# x = dataparse(filename, get_genes=True)
# print('get_promoters=True')
# x = dataparse(filename, get_promoters=True)
# print('parse=True')
# x = dataparse(filename, parse=True)
# print('testing: dataparse.write')
x.write(outfile)
print('testing: dataparse.from_pickle()')
x = dataparse.from_pickle(outfile)
