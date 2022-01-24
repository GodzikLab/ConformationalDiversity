#!/usr/bin/env python

import argparse
import os
from datetime import datetime

def readFileToDF(file, exclPDBs):

	import pandas as pd
	
	df = pd.read_csv(file, sep='\t', names = ['st1', 'st2', 'cov', 'rmsd'])
	st1 = set([i for i in df.loc[:,'st1']])
	st2 = set([i for i in df.loc[:, 'st2']])

	pdbs = list(st1 | st2)
	df2 = pd.DataFrame(0.0, index = pdbs, columns = pdbs)
	# print(df2)

	for i in range(df.shape[0]):
		row = df.iloc[i]
		pdb1 = row['st1']
		pdb2 = row['st2']
		cov = row['cov']
		rmsd = row['rmsd']

		df2.loc[pdb1][pdb2] = rmsd
		df2.loc[pdb2][pdb1] = rmsd
	

	if exclPDBs:
		exclPDBs = list(set(exclPDBs) & set(pdbs))
		df2.drop(labels = exclPDBs, axis = 0, inplace = True)
		df2.drop(labels = exclPDBs, axis = 1, inplace = True)
	

	return df2
	

def getCluster(df, simThresh):

	pdb_list = df.columns.tolist()
	similarPairs = []
	similarPairs_dict = {}
	neighborCount = {}
	for pdb1 in pdb_list:
		neighborCount[pdb1] = 0
		similarPairs_dict[pdb1] = []
		for pdb2 in pdb_list:
			if df.loc[pdb1][pdb2] < simThresh:
				similarPairs.append((pdb1, pdb2))
				similarPairs_dict[pdb1].append(pdb2)
				neighborCount[pdb1] += 1

	maxn = max(neighborCount.values())
	for k, v in neighborCount.items():
		if v == maxn:
			clusterCenter = k
			break

	clusterMembers = similarPairs_dict[clusterCenter]
	df_modified = (df.drop(labels = clusterMembers, axis = 0)).drop(labels = clusterMembers, axis = 1)
	clusterMembers.remove(clusterCenter)

	return clusterCenter, clusterMembers, df_modified


if __name__ == '__main__':


	parser = argparse.ArgumentParser(prog='Subclusters the PDBFlex clusters based on RMSD. Try:\n./createSubclusters.py --help for help.')
	parser.add_argument('--covrmsDir', type=str, required=True, help='Directory containing the all-by-all covrms files.')
	parser.add_argument('--similarityThreshold', type=float, required=True, help='Provide an RMSD threshold below which structures would be considered similar.')
	parser.add_argument('--outputFile', type=str, required=True, help='Path to output file.')
	parser.add_argument('--logFile', type=str, required=True, help='Path to log file.')
	parser.add_argument('--exclude', type=str, required=False, help='Path to file containing list of PDB chains to exclude.')
	args = parser.parse_args()

	covrmsDir = args.covrmsDir
	similarityThreshold = args.similarityThreshold
	outputFile = args.outputFile
	logFile = args.logFile
	excludeFile = args.exclude
	if excludeFile:
		excludePDBs = []
		e = open(excludeFile, 'r')
		for line in e:
			excludePDBs.append(line.strip())
	else:
		excludePDBs = None

	f = open(outputFile, 'w')
	l = open(logFile, 'w')
	f.write('CLUSTER_REP,SUBCLUSTER_REP,SUBCLUSTER_MEMBERS\n')
	for filename in os.listdir(covrmsDir):
		if filename.endswith('.covrms'):
			continue
		else:
			try:
				l.write('Working on {} @ {}\n'.format(filename, datetime.now()))
				df_rmsd = readFileToDF(os.path.join(covrmsDir, filename), excludePDBs)
				while df_rmsd.shape[0] > 0:
					clusterCenterPDB, clusterMembersPDB, df_rmsd = getCluster(df_rmsd, similarityThreshold)
					f.write('{},{},"{}"\n'.format(filename,clusterCenterPDB,clusterMembersPDB))
			except Exception as e:
				print(e)
				l.write('Exception occurred for {}!\n'.format(filename) )

	l.close()
	f.close()
	

