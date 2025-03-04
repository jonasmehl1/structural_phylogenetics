#!/usr/bin/env python

#Script used to obtain the BRH of a pair of reciprocal blast results

# import argparse

#Load blast results into memory after filtering based on evalue and overlap thresholds
def filter_blast(BlastFile, threshold, evalue_thr=1e-05, overlap_thr=0.2):
	hits = {}
	stats = {}
	# if evalue_thr > 0.1:
	# 	print("WARNING: High e-value threshold!")
	for line in open(BlastFile):
		line = line.strip()
		dades = line.split("\t")
		if "#" not in line:
			# try:
			# 	overlap = (float(dades[7]) - float(dades[6])) / float(len(seqs[dades[0]]))
			# except:
			# 	overlap = 0.0
			# if overlap > overlap_thr and float(dades[10]) < evalue_thr:
			if float(dades[10]) < float(evalue_thr):
				if dades[0] not in hits:
					hits[dades[0]] = []
					# stats[dades[0]] = {}
				if dades[1] not in hits[dades[0]] and len(hits[dades[0]]) < threshold:
					hits[dades[0]].append(dades[1])
					# stats[dades[0]][dades[1]] = {}
					# stats[dades[0]][dades[1]]["overlap"] = overlap
					# stats[dades[0]][dades[1]]["eval"] = float(dades[10])
					# stats[dades[0]][dades[1]]["ident"] = float(dades[2])
	return hits

#Get best reciprocal hits
def get_BRH(hits1,hits2):
	BBH = set([])
	taken = set([])
	for c1 in hits1:
		h1 = hits1[c1][0]
		if h1 in hits2:
			if hits2[h1][0] == c1:
				# if c1 > h1:
				pair = c1 +"\t"+ h1
				# else:
				# 	pair = h1 + "\t"+c1
				BBH.add(pair)
				taken.add(h1)
				taken.add(c1)
	return BBH,taken


if __name__ == "__main__":
	# parser = argparse.ArgumentParser(description="Obtains pairs of BRH")
	# parser.add_argument('-h1', dest='h1',required=True,
	# 					help=('blast where seed species is species 1.'),
	# 					nargs='*')
	# parser.add_argument('-h2', dest='h2',required=True,
	# 					help=('blast where seed species is species 2.'),
	# 					nargs='*')
	# parser.add_argument("-e", dest="evalue",
	# 					help="minimum evalue", type=float, default=1e-3)
	# parser.add_argument("-o", dest="out", required=True,
	# 					help="outfileName")

	# args = parser.parse_args()
	with open(snakemake.output[0], "w") as outfile:
		for idx, blast1 in enumerate(snakemake.input["a"]):
			hits1 = filter_blast(snakemake.input["a"][idx], 1, evalue_thr=snakemake.params[0])
			hits2 = filter_blast(snakemake.input["b"][idx], 1, evalue_thr=snakemake.params[0])
			pairs,taken = get_BRH(hits1, hits2)
			outfile.write('\n'.join(pairs)+'\n')
