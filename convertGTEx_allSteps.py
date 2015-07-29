# Process GTEx data into plink format

# TODO should not have two genos per emma. 
# also redundant with the version in grat_proj, except paths/input/output changed

# This script will:
#   tfam - remove header for genotype files and convert to .tfam for each chunk
# 	map and tped - write only 1000 snps per file
#   map - place chr and snp columns in right order
#   tped - add chr, snp, 0, bp-coordinate info from .map
# 	tped - split each genotype into two values per individual, 
#			in plink 0/1/2 (missing/var1/var2) format
#	map and tped - write only filtered snps to file (see "Filtering")
#   make snp_list_xxxx for each file
#   optional - useEMMA option will output EMMA format instead of tped, and no tfam

# Assumes:
#	input geno file is rsID geno1 geno2 .... (no map cols)

# Filtering
#   tped - remove snps with MAF < 0.05
# 	tped - keep imputed values (converted to two val per indiv)
# 	tped - keep snps with missing values

import sys

def makeTfam(headerLn, outBase, outPath = "."):
    headerL = headerLn.rstrip().split()
    # header looks like:
    # Id      GTEX-N7MS       GTEX-N7MT       GTEX-NFK9 ...
    with open(outPath + '/' + outBase + '.tfam', 'w') as tfamf:
        for indivID in headerL[1:]:
        # TODO are these assumptions actually true?
        # assume indivs unrelated, sex unknown, pheno unknown
            tfamf.write(indivID + '\t' + indivID + '\t0\t0\t0\t0\n')

	return


def maf(gtexGenoL):
	'''given gtex genotypes, find maf. treat imputed values as fractional towards maf'''
	# TODO would it be better to treat imputed as nearest whole value
	# TODO handle missing value??????
	countA = sum([float(g) for g in gtexGenoL]) 	# count A alleles (including fracs)
	freqA = countA/(2*len(gtexGenoL))
	return min(freqA, (1-freqA))


def recodeGeno(gtexGeno):
# TODO some imputed values are probably so close to an int that they were rounded
# leading to a handful of 1/2's in a snp that's mostly floats

	# recode genotypes
	# gtex |   plink  | explanation
	# -----------------------------
	# 	2  =	2 2    = homo
	# 	1  = 	1 2    = het (should be '1 2')
	# 	0  = 	1 1    = homo
	#  -9  = 	0 0    = missing

	# imputed genotype conversion examples
	#  1.5 =  1.75 + 1.75 = in-between het/homo 
	#  0.5 = 1.25  1.25 = in-between het/homo

	if gtexGeno == '-9': # not sure if snps.txt even has missing vals or if they are imputed
		return '0 0' 

	if gtexGeno == '2':
		return '2 2'
	if gtexGeno == '1':
		return '1 2'
	if gtexGeno == '0':
		return '1 1'

 	# split imputed value btwn two alleles, add one to convert to 1/2
	gtexVal = float(gtexGeno)
	plinkVal = 1.0*gtexVal/2 + 1 
	return '%f %f' % (plinkVal, plinkVal)


# TODO add tfam from header
def fixFormat(genoToFix, mapToFix, inPath, outBase, outPath, outPath_snp_list, useEMMA = False): #, nonMissingThresh = 1):
	lineNo = -1
	chunkNo = 0
	CHUNKNO = str(chunkNo).zfill(4) # pad to 4 digit w leading zeros

	inGenof = open(inPath + "/" + genoToFix, 'r') #  .snps.txt files
	inMapf = open(inPath + "/" + mapToFix, 'r')

	# output files will be per chunk

	outTpedExtension = ".tped"
	if useEMMA:
		outTpedExtension = ".emma"
	outTpedf = open(outPath + "/" + outBase + "_" + CHUNKNO + outTpedExtension, 'w')
	#outTpedf = open(outPath + "/" + outBase + "_" + CHUNKNO + ".tped", 'w')
	outMapf = open(outPath + "/" + outBase + "_" + CHUNKNO + ".map", 'w')
	outListf = open(outPath_snp_list + "/snp_list_" + CHUNKNO, 'w')

	# remove header (list of indivs) and make into .tfam
	headerLn = inGenof.readline()

	# make .tped and filter .map
	for genoLine in inGenof:
		gtexGenoL = genoLine.rstrip().split()
		geno_rsid = gtexGenoL[0]
		gtexGenoL = gtexGenoL[1:]   # genotypes only


		# .snps.map has entries like "rs79010578 1 736289 9"
		# snp_no i seems to indicate ith snp in file 
		# (plink doesn't use it, and filter out some lines anyway)
		rsid, chro, bp_coord, snp_no = inMapf.readline().split()

		if geno_rsid != rsid:
			print("ERROR: rsID mismatch between input geno and map file\n %s: %s\n %s: %d")
			# another possibility is that map cols already in order. either way, quit
			# print("ERROR: columns in %d are already in correct order" % mapToFilter)
			inGenof.close()
			inMapf.close()
			outTpedf.close()
			outMapf.close()
			return

		# continue has to be after reading line from both snps.txt and snps.map

		# filter out rare variants
		if maf(gtexGenoL) < MAF_CUTOFF:
			continue

		recodedGenoL = [recodeGeno(g) for g in gtexGenoL]

		# if too many indivs have missing vals, don't write snp to file
		#if not (all(g != '0 0' for g in recodedGenoL)):
		#if sum([g != '0 0' for g in recodedGenoL]) < nonMissingThresh:
			#print("All values missing for snp %s\n" % rsid)
			#continue

		lineNo += 1 	# else we're adding a line

		# if first line in new chunk, close current file and open new one before writing
		#if lineNo % 1000 == 0 and chunkNo != 0: # if chunkNo == 0, everything just init before for loop
		if lineNo == 1000:
			outTpedf.close()
			outMapf.close()
			outListf.close()

			lineNo = 0
			chunkNo += 1
			CHUNKNO = str(chunkNo).zfill(4) # pad to 4 digit w leading zeros

			outTpedExtension = ".tped"
			if useEMMA:
				outTpedExtension = ".emma"
			outTpedf = open(outPath + "/" + outBase + "_" + CHUNKNO + outTpedExtension, 'w')
			#outTpedf = open(outPath + "/" + outBase + "_" + CHUNKNO + ".tped", 'w')
			outMapf = open(outPath + "/" + outBase + "_" + CHUNKNO + ".map", 'w')
			outListf = open(outPath_snp_list + "/snp_list_" + CHUNKNO, 'w')


		# write output
		mapStr = " ".join([chro, rsid, "0", bp_coord])
		recodedGenos = " ".join(recodedGenoL)
		#print(mapStr + " " + recodedGenos + "\n")

		if useEMMA:
			outTpedf.write(recodedGenos + "\n") # no map cols in .emma
		else:
			outTpedf.write(mapStr + " " + recodedGenos + "\n") # add map cols to tped
		
		outMapf.write(mapStr + "\n")
		outListf.write(rsid + " " + PRIOR + "\n")

		# also make a Tfam for this chunk
		if not useEMMA:
			makeTfam(headerLn, outBase + "_" + CHUNKNO, outPath)



	inGenof.close()
	inMapf.close()
	
	outTpedf.close()
	outMapf.close()
	outListf.close()

	return



#nonMissingThresh = 96 # min number of non-missing integer-geno indivs allowed

MAF_CUTOFF = 0.005 #0.05
PRIOR = "1e-6" 
# TODO should prior be same as alpha? don't think so since grat example runs used 1e-6 prior and other alphas

genoToFix = "Lung1k.snps.txt"
mapToFilter = "Lung1k.snps.map"

#inPath = "/u/home/f/fhormoz/project/data/GTEx/jaehoon/genotypes"
#inPath = "/u/scratch/l/lgai/GTEx_data/Lung_geno_unsplit"
inPath = "/u/home/l/lgai/multi_pheno/GTEx_data"
#outBase = "Lung"
#outPath = "/u/scratch/l/lgai/GTEx_data/Lung_input_plink" # put plink formatted split data here
#outPath_snp_list = "/u/scratch/l/lgai/GTEx_data/Lung_input_2fat" # snp prior list
outBase = "Lung1k"
outPath = "/u/home/l/lgai/multi_pheno/GTEx_data_plink"
outPath_snp_list = outPath

if len(sys.argv) == 1: # if run with no parameters, use tped
	fixFormat(genoToFix, mapToFilter, inPath, outBase, outPath, outPath_snp_list) #, nonMissingThresh = nonMissingThresh)

elif len(sys.argv) == 2 and sys.argv[1] == "useEMMA":
    outPath = "/u/home/l/lgai/multi_pheno/GTEx_data_emma"
    outPath_snp_list = outPath
    #outPath = "/u/scratch/l/lgai/GTEx_data/Lung_input_emma" # put plink formatted split data here
    #outPath_snp_list = "/u/scratch/l/lgai/GTEx_data/Lung_input_2fat_rotated" # snp prior list
    fixFormat(genoToFix, mapToFilter, inPath, outBase, outPath, outPath_snp_list, useEMMA=True) #, nonMissingThresh = nonMissingThresh)

else:
	print("ERROR: parameter(s) not recognized")


