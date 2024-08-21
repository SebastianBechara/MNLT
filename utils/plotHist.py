#! /usr/bin/env python3

import sys
import pandas as pd
import matplotlib.pyplot as plt

def main():
	IRSFile = sys.argv[1]
	if not IRSFile.lower().endswith('tsv'):
		print('Provided file is not a TSV file.')
		sys.exit()

	kds = sys.argv[2].split(',')

	df = pd.read_table(IRSFile, header=0, sep='\t')
	cols = df.columns

	for kd in kds:
		if kd not in cols:
			print('Specified KD not in File')
			sys.exit()

	cols = ['forestgreen','dodgerblue','red']
	stuff = [ df[kd] for kd in kds]
	labels = kds

	plt.figure(figsize=(9,6))
	plt.hist(stuff, bins=52, color=cols, label=labels)
	plt.legend()
	plt.savefig(f'./IRSHist{"_".join(kds).replace("/","_")}.svg', bbox_inches='tight')

if __name__ == '__main__':
	print(len(sys.argv))
	if len(sys.argv) == 1:
		print('USAGE: plotHist.py <Infile> <KDs>\n\
    Infile: The TSV file containing the IRSs\n\
    KDs: Kds you want in the plot. Specify them them as a comma separated string.\n\
         For better readability I would suggest to limit each plot to 3 KDs\n\n\
It will save the plot in the current working directory')
		sys.exit()
	if len(sys.argv) > 3:
		print('Too many arumgnet!\nUsage: plotHist.py <infile> <KDs>')
		sys.exit()

	main()
