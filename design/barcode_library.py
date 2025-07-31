#!/usr/bin/env python
#-*- encoding: utf-8 -*-

import os
import sys
import random
import numpy as np
from Bio import SeqIO
from collections import Counter


"""select balanced barcode from barcode pool
"""


def readBarcodePool(barcodeFile):
    lstBarcode = []
    for barcode in SeqIO.parse(barcodeFile, 'fasta'):
        barcode_name = barcode.name
        barcode_seq = str(barcode.seq)
        lstBarcode.append((barcode_name, barcode_seq))
    return lstBarcode


def calcScore(lstSelect):
    """calculate the score of barcode balance
    """
    length = len(lstSelect[0][1])
    score = 0
    lstPosBase = []
    for i in range(length):
        posBase = [b[1][i] for b in lstSelect]
        count = Counter(posBase)
        freqs = [count.get(base, 0) for base in 'ACGT']
        lstPosBase.append(freqs)
        std = np.std(freqs)
        score += std
    return score, lstPosBase


def selectBalancedBarcode(lstBarcode, numSelect, outFile, numTrail=100000):
    bestScore = float('inf') # smallest score is better
    bestSelect = None
    bestPosBase = None
    for i in range(numTrail):
        if i % 1e4 == 0:
            print('sample {}'.format(i), flush=True)
        lstSelect = random.sample(lstBarcode, numSelect)
        score, lstPosBase = calcScore(lstSelect)
        if score < bestScore:
            bestScore = score
            bestPosBase = lstPosBase
            bestSelect = lstSelect
    
    # write the result to file
    with open(outFile, 'w') as f:
        for i, (name, seq) in enumerate(bestSelect):
            f.write("barcode{:03}\t{}\n".format(i+1, seq))
        f.write('\n')
        f.write('score: {}\n'.format(bestScore))
        f.write('\tA\tC\tG\tT\n')
        for i, posBase in enumerate(bestPosBase):
            f.write('{}\t{}\t{}\t{}\t{}\n'.format(i+1, *bestPosBase[i]))


def main():
    if len(sys.argv) != 4:
        print("Usage: %s <barcodeFile> <selectBarcodeFile> <numSelect>" % sys.argv[0])
        exit(1)
        
    barcodeFile = sys.argv[1]
    selectBarcodeFile = sys.argv[2]
    numSelect = int(sys.argv[3])
    
    lstBarcode = readBarcodePool(barcodeFile)
    selectBalancedBarcode(lstBarcode, numSelect, selectBarcodeFile)
    
    
if __name__ == '__main__':
    main()
