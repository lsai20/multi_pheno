# make identity matrix for specified number of individuals

import sys

def writeIdentity(outfName, N):
	with open(outfName, 'w') as outf:
		for i in range(N):
			lineL = ['0' for j in range(N)]
			lineL[i] = '1'
			outf.write('\t'.join(lineL) + '\n')

	return

def usage():
	print('Usage: python makeIdentityMatrix.py [output filename] [num indivs]')
	return


outfName = sys.argv[1]
N = sys.argv[2]
writeIdentity(outfName, int(N))
