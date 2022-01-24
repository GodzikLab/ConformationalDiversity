#!/usr/bin/env python

class PDBChainError(Exception):
	pass
class EmptyAlignmentFileException(Exception):
	pass
class NoAlignmentException(Exception):
	pass
class PDBSequenceError(Exception):
	pass
class AlignmentCorrectionError(Exception):
	pass
class ResidueEquivalenceCorrectionError(Exception):
	pass
class CalphaNotFoundError(Exception):
	pass
class UnequalMatrixShapesError(Exception):
	pass

def parseAlignmentFile(alignFile, verbose):

    if verbose:
        print("Parsing alignment file: ", alignFile)

    import re

    f = open(alignFile, 'r')
    lines = []
    for line in f:
        lines.append(line)
    if len(lines) == 0:
        raise EmptyAlignmentFileException("EmptyAlignmentFileException: {} is empty!\n".format(alignFile))
    if re.search('NO ALIGNMENT', lines[0]):
        raise NoAlignmentException("NoAlignmentException: No alignment in: {}\n".format(alignFile))
    pdbChain1_alignOffset = int(lines[3].strip().split()[0])
    pdbChain1_alignSeq = lines[3].strip().split()[1]
    pdbChain2_alignOffset = int(lines[4].strip().split()[0])
    pdbChain2_alignSeq = lines[4].strip().split()[1]

    seq1_len = 0
    seq1_alignedLen = 0
    seq2_len = 0
    seq2_alignedLen = 0
    for i in range(len(pdbChain1_alignSeq)):
        s1 = pdbChain1_alignSeq[i]
        s2 = pdbChain2_alignSeq[i]
        if s1 != '-' and s2 != '-':
            seq1_len += 1
            seq2_len += 1
            seq1_alignedLen += 1
            seq2_alignedLen += 1
        elif s1 != '-' and s2 == '-':
            seq1_len += 1
        elif s1 == '-' and s2 != '-':
            seq2_len += 1
        else:
            continue
    seq1_cov = float(seq1_alignedLen) / seq1_len * 100
    seq2_cov = float(seq2_alignedLen) / seq2_len * 100

    return pdbChain1_alignOffset, pdbChain1_alignSeq, seq1_cov, pdbChain2_alignOffset, pdbChain2_alignSeq, seq2_cov

def getPDBStructure(pdbChain, pdbDir, verbose):

    import os
    import Bio.PDB
    from Bio.PDB.PDBParser import PDBParser
    from Bio.PDB.PDBExceptions import PDBConstructionException

    if verbose:
        print("Getting PDB chain structure for: ", pdbChain)

    pdbID = pdbChain[:4]
    chainName = pdbChain[4:]

    structure =  PDBParser().get_structure("struc", os.path.join(pdbDir, pdbID+'.pdb'))
    chain = None
    for c in structure.get_chains():
        if c.get_id() == chainName:
            chain = c

    if chain == None:
        raise PDBChainError("PDBChainError: {} chain is None!\n".format(pdbChain))

    return chain

def compareSequenceToPDB(alignSeq, alignOff, alignLen, pdbChStruc, verbose):

    import Bio.PDB

    if verbose:
        print("Comparing alignment sequence to PDB sequence.")

    aaDict = {"GLY": "G", "ALA": "A", "VAL": "V", "LEU": "L", "ILE": "I", "PHE": "F", "TRP": "W", "PRO": "P", "LYS": "K", "ARG": "R", "HIS": "H", "ASP": "D", "GLU": "E", "SER": "S", "THR": "T", "TYR": "Y", "MET": "M", "CYS": "C", "ASN": "N", "GLN": "Q", "MSE": "M"}

    alignSeq_woGapsList = []
    for c in alignSeq:
        if c != '-':
            alignSeq_woGapsList.append(c)

    alignSeq_woGaps = ''.join(alignSeq_woGapsList)

    pdbResidues_pre = []
    for r in pdbChStruc.get_residues():
        if r.id[0] == ' ':
            pdbResidues_pre.append(r)
        else:
            if r.get_resname() == 'MSE':
                pdbResidues_pre.append(r)
    pdbResidues = pdbResidues_pre[alignOff - 1: alignOff - 1 + alignLen]

    pdbSeqList = []
    for r in pdbResidues:
        if r.get_resname() in aaDict.keys():
            pdbSeqList.append(aaDict[r.get_resname()])
    pdbSeq = ''.join(pdbSeqList)
    returnValue = alignSeq_woGaps == pdbSeq

    return returnValue

def getResidueEquivalence(pdbCh1, alignOff1, alignSeq1, pdbChStruc1, pdbCh2, alignOff2, alignSeq2, pdbChStruc2, verbose):

    if verbose:
        print("Getting residue equivalence for {} and {}".format(pdbCh1, pdbCh2))

    seq1_alignedLen = 0
    for c in alignSeq1:
        if c != '-':
            seq1_alignedLen += 1
    seq2_alignedLen = 0
    for c in alignSeq2:
        if c != '-':
            seq2_alignedLen += 1

    seqs_same = compareSequenceToPDB(alignSeq1, alignOff1, seq1_alignedLen, pdbChStruc1, verbose)
    if not seqs_same:
        raise PDBSequenceError("PDBSequenceError: PDB sequence != aligned sequence for {}\n".format(pdbCh1))

    seqs_same = compareSequenceToPDB(alignSeq2, alignOff2, seq2_alignedLen, pdbChStruc2, verbose)
    if not seqs_same:
        raise PDBSequenceError("PDBSequenceError: PDB sequence != aligned sequence for {}\n".format(pdbCh2))

    if verbose:
        print('PDBSeq same as aligned seq w/o gaps')

    residues1_pre = []
    for r in pdbChStruc1.get_residues():
        if r.id[0] == ' ':
            residues1_pre.append(r)
        else:
            if r.get_resname() == 'MSE':
                residues1_pre.append(r)
    residues1 = residues1_pre[alignOff1 - 1: alignOff1 - 1 + seq1_alignedLen]

    residues2_pre = []
    for r in pdbChStruc2.get_residues():
        if r.id[0] == ' ':
            residues2_pre.append(r)
        else:
            if r.get_resname() == 'MSE':
                residues2_pre.append(r)
    residues2 = residues2_pre[alignOff2 - 1: alignOff2 - 1 + seq2_alignedLen]

    residueEq1_2 = []
    residueEq2_1 = []
    seq1_pos = 0
    seq2_pos = 0
    for i in range(len(alignSeq1)):
        r1 = alignSeq1[i]
        r2 = alignSeq2[i]
        if r1 != '-':
            res1 = residues1[seq1_pos]
            seq1_pos += 1
            if r2 != '-':
                res2 = residues2[seq2_pos]
                seq2_pos += 1
                residueEq1_2.append((res1, res2)) 
                residueEq2_1.append((res2, res1))
            else:
                res2 = None
                residueEq1_2.append((res1,res2))
        else:
            res1 = None
            if r2 != '-':
                res2 = residues2[seq2_pos]
                seq2_pos += 1
                residueEq2_1.append((res2, res1))
            else:
                res2 = None

    return residueEq1_2, residueEq2_1


def correctResidueEquivalences(residueEq_list1, residueEq_list2):

    commonResidues = set([i[0] for i in residueEq_list1]) & set([i[0] for i in residueEq_list2])
    residueEq_list1_corr = []
    residueEq_list2_corr = []
    for x, y in residueEq_list1:
        if x in commonResidues:
            residueEq_list1_corr.append((x, y))
    for x, y in residueEq_list2:
        if x in commonResidues:
            residueEq_list2_corr.append((x, y))

    if not set(residueEq_list1_corr) <= set(residueEq_list1):
        raise ResidueEquivalenceCorrectionError("ResidueEquivalenceCorrectionError: set(residueEq_list1_corr) is not a subset of set(residueEq_list1)\n")
    if not set(residueEq_list2_corr) <= set(residueEq_list2):
        raise ResidueEquivalenceCorrectionError("ResidueEquivalenceCorrectionError: set(residueEq_list2_corr) is not a subset of set(residueEq_list2)\n")

    return residueEq_list1_corr, residueEq_list2_corr

def getDistanceMatrix(residues):

    '''
    residues = [<Residue GLN het=  resseq=230 icode= >, <Residue LEU het=  resseq=228 icode= >]
    '''
    from scipy.spatial.distance import pdist,squareform

    Calpha_coords_list = []
    for r in residues:
        rID = r.id
        Calpha_coord = []
        for a in r.get_atoms():
            if a.get_id() == 'CA':
                Calpha_coord = a.coord
        if len(Calpha_coord )!= 0:
            Calpha_coords_list.append(Calpha_coord)
        else:
            raise CalphaNotFoundError("CalphaNotFoundError: Calpha not found for residue with ID: {}\n".format(rID))

    dm = squareform(pdist(Calpha_coords_list, metric='euclidean'))

    return dm

def plotMatrixCorrelation(m1_elts, m2_elts, file_name, verbose):

    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot as plt

    if verbose:
        print("Creating matrix correlation figure {}".format(file_name))

    fig=plt.figure(figsize = (5, 5), dpi = 300)
    plt.scatter(m1_elts, m2_elts)
    plt.xlabel("Value from DDM 1")
    plt.ylabel("Value from DDM 2")
    plt.savefig(file_name, dpi = 300, bbox_inches = "tight")
    plt.close(fig)

    return True

def calculateMatrixCorrelation(m1, m2, file_name, verbose):

    from scipy.stats import pearsonr, spearmanr

    if verbose:
        print("Calculating matrix correlation")

    if not m1.shape == m2.shape:
        raise UnequalMatrixShapesError("UnequalMatrixShapesError: m1.shape {} != m2.shape {}\n".format(m1.shape, m2.shape))

    m1_elements = []
    m2_elements = []
    for i in range(m1.shape[0]):
        for j in range(m1.shape[0]):
            m1_elements.append(m1[i, j])
            m2_elements.append(m2[i, j])

    plotMatrixCorrelation(m1_elements, m2_elements, file_name, verbose)

    pearson_r = pearsonr(m1_elements, m2_elements)[0]
    pearson_p = pearsonr(m1_elements, m2_elements)[1]
    spearman_r = spearmanr(m1_elements, m2_elements)[0]
    spearman_p = spearmanr(m1_elements, m2_elements)[1]
    return pearson_r, pearson_p, spearman_r, spearman_p

def MatrixFigure(m, file_name, verbose):

    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot as plt
    from matplotlib.colors import ListedColormap, BoundaryNorm

    if verbose:
        print("Creating matrix figure {} ".format(file_name))

    max_min_values = [np.amax(m), -np.amin(m)]
    abs_max = max(max_min_values)

    fig, ax = plt.subplots(figsize = (5, 5), dpi = 300)
    heatmap = ax.imshow(m, cmap='bwr', interpolation=None, vmin=-abs_max, vmax=abs_max)
    ax.figure.colorbar(heatmap, ax=ax)
    plt.xlabel("Alignment position")
    plt.ylabel("Alignment position")
    plt.savefig(file_name, dpi = 300, bbox_inches = "tight")
    plt.close(fig)

    return True

if __name__ == '__main__':

    import sys
    import os
    import numpy as np
    from scipy.stats import pearsonr
    import pandas as pd
    import ast
    import argparse
    import Bio.PDB
    from Bio.PDB.PDBExceptions import PDBConstructionException

    parser = argparse.ArgumentParser(prog='Compare number of subclusters for homologs. Try:\n./countSubclustersForHomologs.py --help for help.')
    parser.add_argument('--inputFile', type=str, required=True, help='Input file')
    parser.add_argument('--outputFile', type=str, required=True, help='Output file')
    parser.add_argument('--logFile', type=str, required=True, help='Log file')
    parser.add_argument('--alignDir', type=str, required=True, help='Dir containing alignment files')
    parser.add_argument('--pdbDir', type=str, required=True, help='Dir containing PDB files')
    parser.add_argument('--makeFigures', action='store_true', help='Figures will be made for the DMs and DDMs.')
    parser.add_argument('--figDir', type=str, required=False, help='Dir in which to write fig files.')
    parser.add_argument('--verbose', action='store_true', help='Verbose mode on.')

    args = parser.parse_args()

    inputFile = args.inputFile
    outputFile = args.outputFile
    logFile = args.logFile
    alignmentDir = args.alignDir
    pdbFilesDir = args.pdbDir
    makeFigures = args.makeFigures
    figureDir = args.figDir
    verb = args.verbose
    outputFileHandle = open(outputFile, 'w')
    outputFileHandle.write('QUERY_CLUSTER\tSUBJECT_CLUSTER\tCLUSTER_SEQUENCE_ID\tQUERY_SUBCLUSTERS\tSUBJECT_SUBCLUSTERS\tDDM_CORRELATION_PEARSON\tP-VALUE_PEARSON\tDDM_CORRELATION_SPEARMAN\tP-VALUE_SPEARMAN\tFIXED_SEQ_ID\tIG_PAIR\tSEQ_ID_RANGE\n')

    logFileHandle = open(logFile, 'w')

    inputFileHandle = open(inputFile, 'r')
    header = inputFileHandle.readline()
    for line in inputFileHandle:
        qCluster = line.strip().split('\t')[0]
        sCluster = line.strip().split('\t')[1]
        clusterSeqID = line.strip().split('\t')[2]
        qSubclusters = ast.literal_eval(ast.literal_eval(line.strip().split('\t')[3]))
        sSubclusters = ast.literal_eval(ast.literal_eval(line.strip().split('\t')[4]))
        igPair = line.strip().split('\t')[5]
        seqID_range = line.strip().split('\t')[6]
        pdbChain1 = qSubclusters[0]
        pdbChain2 = qSubclusters[1]
        pdbChain3 = sSubclusters[0]
        pdbChain4 = sSubclusters[1]
        try:
            logFileHandle.write("Working on pdbChain1: {}, pdbChain2: {}, pdbChain3: {}, pdbChain4: {}\n".format(pdbChain1, pdbChain2, pdbChain3, pdbChain4))
            chainStruc1 = getPDBStructure(pdbChain1, pdbFilesDir, verb)
            chainStruc2 = getPDBStructure(pdbChain2, pdbFilesDir, verb)
            chainStruc3 = getPDBStructure(pdbChain3, pdbFilesDir, verb)
            chainStruc4 = getPDBStructure(pdbChain4, pdbFilesDir, verb)

            # For subcluster pair of first cluster:
            if verb:
                logFileHandle.write("Working on subcluster pair of first cluster\n")

            alignmentFile12 = os.path.join(alignmentDir, '{}.{}.fixed.mu'.format(pdbChain1, pdbChain2))
            alignOffset1, alignSequence1, alignCoverage1, alignOffset2, alignSequence2, alignCoverage2 = parseAlignmentFile(alignmentFile12, verb)
            if alignCoverage1 < 90 or alignCoverage2 < 90:
                logFileHandle.write("LowCoverageWarning: Coverage less than 90% for subcluster representatives: {} and {}\n".format(pdbChain1, pdbChain2))

            residueEq1_2_list, residueEq2_1_list = getResidueEquivalence(pdbChain1, alignOffset1, alignSequence1, chainStruc1, pdbChain2, alignOffset2, alignSequence2, chainStruc2, verb)

            residueEq1_2_list_corr = []
            for x, y in residueEq1_2_list:
                if y != None:
                    residueEq1_2_list_corr.append((x, y))
            residueEq2_1_list_corr = []
            for x, y in residueEq2_1_list:
                if y != None:
                    residueEq2_1_list_corr.append((x, y))

            if len(residueEq1_2_list_corr) != len(residueEq2_1_list_corr):
                raise AlignmentCorrectionError("AlignmentCorrectionError: Corrected residue equivalences are not equal for: {} and {}\n".format(pdbChain1, pdbChain2))

            for i in range(len(residueEq1_2_list_corr)):
                x1, y1 = residueEq1_2_list_corr[i]
                x2, y2 = residueEq2_1_list_corr[i]
                if x1 != y2 or y1 != x2:
                    raise AlignmentCorrectionError("AlignmentCorrectionError: Corrected residue equivalences are not equal for: {} and {}\n".format(pdbChain1, pdbChain2))
            if verb:
                logFileHandle.write("Alignments corrected successfully for pdbChains {} and {}\n".format(pdbChain1, pdbChain2))

            distanceMatrix1 = getDistanceMatrix([i[0] for i in residueEq1_2_list_corr])
            distanceMatrix2 = getDistanceMatrix([i[0] for i in residueEq2_1_list_corr])
            if not distanceMatrix1.shape == distanceMatrix2.shape:
                raise UnequalMatrixShapesError("UnequalMatrixShapesError: distanceMatrix1.shape {} != distanceMatrix2.shape {}\n".format(distanceMatrix1.shape, distanceMatrix2.shape))
            ddm12 = distanceMatrix1 - distanceMatrix2

            if makeFigures:
                MatrixFigure(ddm12, os.path.join(figureDir, '{}.{}_DDM_orig.png'.format(pdbChain1, pdbChain2)), verb)

            # For subcluster pair of second cluster:
            if verb:
                logFileHandle.write("Working on subcluster pair of second cluster\n")

            alignmentFile34 = os.path.join(alignmentDir, '{}.{}.fixed.mu'.format(pdbChain3, pdbChain4))
            alignOffset3, alignSequence3, alignCoverage3, alignOffset4, alignSequence4, alignCoverage4 = parseAlignmentFile(alignmentFile34, verb)
            if alignCoverage3 < 90 or alignCoverage4 < 90:
                logFileHandle.write("LowCoverageWarning: Coverage less than 90% for subcluster representatives: {} and {}\n".format(pdbChain3, pdbChain4))

            residueEq3_4_list, residueEq4_3_list = getResidueEquivalence(pdbChain3, alignOffset3, alignSequence3, chainStruc3, pdbChain4, alignOffset4, alignSequence4, chainStruc4, verb)

            residueEq3_4_list_corr = []
            for x, y in residueEq3_4_list:
                if y != None:
                    residueEq3_4_list_corr.append((x, y))
            residueEq4_3_list_corr = []
            for x, y in residueEq4_3_list:
                if y != None:
                    residueEq4_3_list_corr.append((x, y))

            if len(residueEq3_4_list_corr) != len(residueEq4_3_list_corr):
                raise AlignmentCorrectionError("AlignmentCorrectionError: Corrected residue equivalences are not equal for: {} and {}\n".format(pdbChain3, pdbChain4))

            for i in range(len(residueEq3_4_list_corr)):
                x1, y1 = residueEq3_4_list_corr[i]
                x2, y2 = residueEq4_3_list_corr[i]
                if x1 != y2 or y1 != x2:
                    raise AlignmentCorrectionError("AlignmentCorrectionError: Corrected residue equivalences are not equal for: {} and {}\n".format(pdbChain3, pdbChain4))
            if verb:
                logFileHandle.write("Alignments corrected successfully for pdbChains {} and {}\n".format(pdbChain3, pdbChain4))

            distanceMatrix3 = getDistanceMatrix([i[0] for i in residueEq3_4_list_corr])
            distanceMatrix4 = getDistanceMatrix([i[0] for i in residueEq4_3_list_corr])
            if not distanceMatrix3.shape == distanceMatrix4.shape:
                raise UnequalMatrixShapesError("UnequalMatrixShapesError: distanceMatrix3.shape {} != distanceMatrix4.shape {}\n".format(distanceMatrix3.shape, distanceMatrix4.shape))
            ddm34 = distanceMatrix3 - distanceMatrix4

            if makeFigures:
                MatrixFigure(ddm34, os.path.join(figureDir, '{}.{}_DDM_orig.png'.format(pdbChain3, pdbChain4)), verb)

            # For homolog alignment:
            if verb:
                logFileHandle.write("Working on homolog alignment\n")

            alignmentFile13 = os.path.join(alignmentDir, '{}.{}.fixed.mu'.format(pdbChain1, pdbChain3))
            alignOffset1_h, alignSequence1_h, alignCoverage1_h, alignOffset3_h, alignSequence3_h, alignCoverage3_h = parseAlignmentFile(alignmentFile13, verb)
            if alignCoverage1_h < 90 or alignCoverage3_h < 90:
                 logFileHandle.write("LowCoverageWarning: Coverage less than 90% for homologous subcluster representatives: {} and {}\n".format(pdbChain1, pdbChain3))

            residueEq1_3_list, residueEq3_1_list = getResidueEquivalence(pdbChain1, alignOffset1_h, alignSequence1_h, chainStruc1, pdbChain3, alignOffset3_h, alignSequence3_h, chainStruc3, verb)

            residueEq1_3_list_corr = []
            for x, y in residueEq1_3_list:
                if y != None:
                    residueEq1_3_list_corr.append((x, y))
            residueEq3_1_list_corr = []
            for x, y in residueEq3_1_list:
                if y != None:
                    residueEq3_1_list_corr.append((x, y))

            if len(residueEq1_3_list_corr) != len(residueEq3_1_list_corr):
                raise AlignmentCorrectionError("AlignmentCorrectionError: Corrected residue equivalences are not equal for: {} and {}\n".format(pdbChain1, pdbChain3))

            for i in range(len(residueEq1_3_list_corr)):
                x1, y1 = residueEq1_3_list_corr[i]
                x2, y2 = residueEq3_1_list_corr[i]
                if x1 != y2 or y1 != x2:
                    raise AlignmentCorrectionError("AlignmentCorrectionError: Corrected residue equivalences are not equal for: {} and {}\n".format(pdbChain1, pdbChain3))
            if verb:
                logFileHandle.write("Alignments corrected successfully for pdbChains {} and {}\n".format(pdbChain1, pdbChain3))


            # Correct residue equivalences for both subcluster pairs based on all alignments
            residueEq1_3_list_corr2, residueEq1_2_list_corr2 = correctResidueEquivalences(residueEq1_3_list_corr, residueEq1_2_list_corr)
            residueEq3_1_list_corr2 = [(y, x) for x, y in residueEq1_3_list_corr2]
            if not set(residueEq3_1_list_corr2) <= set(residueEq3_1_list_corr):
                raise ResidueEquivalenceCorrectionError("ResidueEquivalenceCorrectionError: set(residueEq3_1_list_corr2) is not a subset of set(residueEq3_1_list_corr)\n")
            residueEq2_1_list_corr2 = [(y, x) for x, y in residueEq1_2_list_corr2]
            if not set(residueEq2_1_list_corr2) <= set(residueEq2_1_list_corr):
                raise ResidueEquivalenceCorrectionError("ResidueEquivalenceCorrectionError: set(residueEq2_1_list_corr2) is not a subset of set(residueEq2_1_list_corr)\n")

            residueEq3_1_list_corr3, residueEq3_4_list_corr2 = correctResidueEquivalences(residueEq3_1_list_corr2, residueEq3_4_list_corr)
            residueEq1_3_list_corr3 = [(y, x) for x, y in residueEq3_1_list_corr3]
            if not set(residueEq1_3_list_corr3) <= set(residueEq1_3_list_corr2):
                raise ResidueEquivalenceCorrectionError("ResidueEquivalenceCorrectionError: set(residueEq1_3_list_corr3) is not a subset of set(residueEq1_3_list_corr2)\n")
            residueEq4_3_list_corr2 = [(y, x) for x, y in residueEq3_4_list_corr2]
            if not set(residueEq4_3_list_corr2) <= set(residueEq4_3_list_corr):
                raise ResidueEquivalenceCorrectionError("ResidueEquivalenceCorrectionError: set(residueEq4_3_list_corr2) is not a subset of set(residueEq4_3_list_corr)\n")

            residueEq1_3_list_corr4, residueEq1_2_list_corr3 = correctResidueEquivalences(residueEq1_3_list_corr3, residueEq1_2_list_corr2)
            if not residueEq1_3_list_corr4 == residueEq1_3_list_corr3:
                raise ResidueEquivalenceCorrectionError("ResidueEquivalenceCorrectionError: residueEq1_3_list_corr4 != residueEq1_3_list_corr3\n")
            residueEq3_1_list_corr4 = [(y, x) for x, y in residueEq1_3_list_corr4]
            if not residueEq3_1_list_corr4 == residueEq3_1_list_corr3:
                raise ResidueEquivalenceCorrectionError("ResidueEquivalenceCorrectionError: residueEq3_1_list_corr4 != residueEq3_1_list_corr3\n")
            residueEq2_1_list_corr3 = [(y, x) for x, y in residueEq1_2_list_corr3]
            if not set(residueEq2_1_list_corr3) <= set(residueEq2_1_list_corr2):
                raise ResidueEquivalenceCorrectionError("ResidueEquivalenceCorrectionError: set(residueEq2_1_list_corr3 is not a subset of set(residueEq2_1_list_corr2\n")

            if verb:
                logFileHandle.write("All residue equivalences corrected successfully!\n")

            distanceMatrix1 = getDistanceMatrix([i[0] for i in residueEq1_3_list_corr4])
            distanceMatrix2 = getDistanceMatrix([i[0] for i in residueEq2_1_list_corr3])
            if not distanceMatrix1.shape == distanceMatrix2.shape:
                raise UnequalMatrixShapesError("UnequalMatrixShapesError: distanceMatrix1.shape {} != distanceMatrix2.shape {}\n".format(distanceMatrix1.shape, distanceMatrix2.shape))
            ddm12 = distanceMatrix1 - distanceMatrix2
            distanceMatrix3 = getDistanceMatrix([i[0] for i in residueEq3_4_list_corr2])
            distanceMatrix4 = getDistanceMatrix([i[0] for i in residueEq4_3_list_corr2])
            if not distanceMatrix3.shape == distanceMatrix4.shape:
                raise UnequalMatrixShapesError("UnequalMatrixShapesError: distanceMatrix3.shape {} != distanceMatrix4.shape {}\n".format(distanceMatrix3.shape, distanceMatrix4.shape))
            ddm34 = distanceMatrix3 - distanceMatrix4
            pearson_cor, pearson_pvalue, spearman_cor, spearman_pvalue = calculateMatrixCorrelation(ddm12, ddm34, os.path.join(figureDir, '{}.{}_matrixCorr.png'.format(qCluster, sCluster)), verb)

            # Calculate SeqID from corrected alignments:
            sameResidues = 0
            for x, y in residueEq1_3_list_corr4:
                if x.get_resname() == y.get_resname():
                    sameResidues += 1
            fixedSeqID = float(sameResidues) / len(residueEq1_3_list_corr4) * 100

            outputFileHandle.write('{}\t{}\t{}\t"{}"\t"{}"\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(qCluster, sCluster, clusterSeqID, qSubclusters, sSubclusters, pearson_cor, pearson_pvalue, spearman_cor, spearman_pvalue, fixedSeqID, igPair, seqID_range))

            if makeFigures:
                MatrixFigure(distanceMatrix1, os.path.join(figureDir, '{}_DM_corr.png'.format(pdbChain1)), verb)
                MatrixFigure(distanceMatrix2, os.path.join(figureDir, '{}_DM_corr.png'.format(pdbChain2)), verb)
                MatrixFigure(distanceMatrix3, os.path.join(figureDir, '{}_DM_corr.png'.format(pdbChain3)), verb)
                MatrixFigure(distanceMatrix4, os.path.join(figureDir, '{}_DM_corr.png'.format(pdbChain4)), verb)
                MatrixFigure(ddm12, os.path.join(figureDir, '{}.{}_DDM_corr_accordingTo_{}.{}.png'.format(pdbChain1, pdbChain2, pdbChain3, pdbChain4)), verb)
                MatrixFigure(ddm34, os.path.join(figureDir, '{}.{}_DDM_corr_accordingTo_{}.{}.png'.format(pdbChain3, pdbChain4, pdbChain1, pdbChain2)), verb)

            logFileHandle.write("Finished working on pdbChain1: {}, pdbChain2: {}, pdbChain3: {}, pdbChain4: {}\n".format(pdbChain1, pdbChain2, pdbChain3, pdbChain4))

        except IOError as e:
            logFileHandle.write("IOError: {}\n".format(str(e)))
        except PDBChainError as e:
            logFileHandle.write(str(e))
        except PDBSequenceError as e:
            logFileHandle.write(str(e))
        except AlignmentCorrectionError as e:
            logFileHandle.write(str(e))
        except ResidueEquivalenceCorrectionError as e:
            logFileHandle.write(str(e))
        except CalphaNotFoundError as e:
            logFileHandle.write(str(e))
        except EmptyAlignmentFileException as e:
            logFileHandle.write(str(e))
        except NoAlignmentException as e:
            logFileHandle.write(str(e))
        except PDBConstructionException as e:
            logFileHandle.write("PDBConstructionException: {}\n".format(str(e)))
        except UnequalMatrixShapesError as e:
            logFileHandle.write(str(e))
