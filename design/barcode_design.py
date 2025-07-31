#!/usr/bin/env python
#-*- encoding=utf-8 -*-

__author__ = 'yishenxing'
__version__ = '1.0'


import re
import random
import logging
import argparse
import Levenshtein
from multiprocessing import Pool, Manager, Lock, Process
from collections import Counter


class BarcodeDesigner:
    """BarcodeDesigner: A class for designing barcode pools for high throughput sequencing.

    This class provides functionality to generate optimized barcode pools for NGS applications, 
    with features to minimize cross-talk and maximize barcode diversity.
    """

    def __init__(self, logger, barcodeLen, polymerLen, minGCRate, maxGCRate, dictRepeats, 
        innerPairLen, outerPairLen, minEditDistance, numCandidateBarcode, numGeneratedBarcode, 
        numProcess, outputSeqFile, outputDistanceFile, outputLogFile):
        """Create a BarcodeDesigner object.
        
        Arguments:
            logger: logging.Logger      the logging handler
            barcodeLen: int             the length of barcode
            polymerLen: int             the maximum length of polymer
            minGCRate: float            the minimum GC rate of barcode
            maxGCRate: float            the maximum GC rate of barcode
            dictRepeats: dict           the maximum repeat unit and its number in barcode. 
                                        {repeatUnit: repeatNum, ...}
            innerPairLen: int           the maximum length of inner pair
            outerPairLen: int           the maximum length of outer pair
            minEditDistance: int        the minimum edit distance between barcodes
            numCandidateBarcode: int    the number of candidate barcode
            numGeneratedBarcode: int    the number of final generated barcode
            numProcess: int             the number of process
            outputSeqFile: str          the output barcode sequence file
            outputDistanceFile: str     the output barcode distance infomation file
            outputLogFile: str          the output log file
            
        >>> barcodeDesigner = BarcodeDesigner(
        ...     logger, barcodeLen, polymerLen, minGCRate, maxGCRate, dictRepeats, innerPairLen, 
        ...     outerPairLen, minEditDistance, numCandidateBarcode, numGeneratedBarcode, 
        ...     numProcess, outputSeqFile, outputDistanceFile, logFile
        ... )
        >>> barcodeDesigner.runDesignTaskParallel()
        """
        self.logger = logger
        self.barcodeLen = barcodeLen
        self.polymerLen = polymerLen
        self.polymerPattern = re.compile(r'(A{{{},}}|C{{{},}}|G{{{},}}|T{{{},}})'.format(
            self.polymerLen, self.polymerLen, self.polymerLen, self.polymerLen
        ))
        self.minGCRate = minGCRate
        self.maxGCRate = maxGCRate
        self.dictRepeats = [_.split(':') for _ in dictRepeats.strip().split(',')]
        self.dictRepeats = {_[0]: int(_[1]) for _ in self.dictRepeats}
        self.innerPairLen = innerPairLen
        self.outerPairLen = outerPairLen
        self.minEditDistance = minEditDistance
        
        self.numCandidateBarcode = numCandidateBarcode
        # lock.  for multiprocessing
        self.lock = Lock()
        # shared list.  for multiprocessing
        self.lstCandidateBarcode = Manager().list()
        # shared dict.  store the sub sequence of barcode that already in candidate list, 
        # for avoiding sub sequence between barcodes
        self.dictSubSequence = Manager().dict()
        # add the sequence that need to be excluded to the sub sequence dict
        self._getExistSeq()

        self.numGeneratedBarcode = numGeneratedBarcode
        self.lstGeneratedBarcode = []
        
        self.numProcess = numProcess
        self.outputSeqFile = outputSeqFile
        self.outputDistanceFile = outputDistanceFile
        self.logFile = open(outputLogFile, 'w')
        
        # record the filter type, for debug
        self.dictErrorType = Manager().dict({
            'GC rate error': 0,
            'polymer error': 0,
            'repeat error': 0, 
            'inner pair error': 0,
            'outer pair error': 0,
            'editDistance error': 0
        })


    def _checkGC(self, seqBarcode):
        """Check the GC rate of barcode."""
        gcCount = seqBarcode.count('G') + seqBarcode.count('C')
        gcRate = float(gcCount) / len(seqBarcode)
        if self.minGCRate <= gcRate and gcRate <= self.maxGCRate:
            return True
        return False


    def _filtEditDistance(self, ):
        """check the edit distance between barcodes. Filter barcodes that do 
        not meet the edit distance requirement."""
        # collect all edit distance to check the edit distance between barcodes.
        lstEditDistanceCollector = []
        # create the edit distance matrix
        lstEditDistance = [[0]*self.numCandidateBarcode for _ in range(self.numCandidateBarcode)]
        # barcode filtered status list
        lstFiltedBarcodeStat = [1]*self.numCandidateBarcode

        # calculate the edit distance. if the edit distance is greater than or equal to the 
        # minimum requirement, mark as 1
        for i in range(self.numCandidateBarcode):
            for j in range(i+1, self.numCandidateBarcode):
                editDistance = Levenshtein.distance(
                    self.lstCandidateBarcode[i], self.lstCandidateBarcode[j]
                )
                lstEditDistanceCollector.append(editDistance)
                if editDistance >= self.minEditDistance:
                    lstEditDistance[i][j] = 1
                    lstEditDistance[j][i] = 1

        # start the loop to filter the barcode
        while True:
            # if the number of filtered barcode is equal to the number of generated barcode, break
            if sum(lstFiltedBarcodeStat) == self.numGeneratedBarcode:
                break
            # traverse all candidate barcodes, find the smallest edit distance, kick it out
            index = -1
            # there will be at least one zero in the matrix, so the minimum meet is 1
            minMeet = sum(lstFiltedBarcodeStat) - 1
            for i in range(self.numCandidateBarcode):
                if lstFiltedBarcodeStat[i] == 0:
                    continue
                if sum(lstEditDistance[i]) < minMeet:
                    minMeet = sum(lstEditDistance[i])
                    index = i
            # if do not find any, it means all barcodes are qualified, break
            if index == -1:
                break
            # set the status to 0 in the status list, and set the corresponding row and column to 
            # 0 in the matrix
            lstFiltedBarcodeStat[index] = 0
            for i in range(self.numCandidateBarcode):
                lstEditDistance[i][index] = 0
                lstEditDistance[index][i] = 0

            with self.lock:
                self.dictErrorType['editDistance error'] += 1

        # add the qualified barcode to the final barcode list
        for i in range(self.numCandidateBarcode):
            if lstFiltedBarcodeStat[i] == 1:
                self.lstGeneratedBarcode.append(self.lstCandidateBarcode[i])

        c = sorted(Counter(lstEditDistanceCollector).items(), key=lambda x:x[0])
        print('edit distance distribution:', c, file=self.logFile)


    def _getRevCmp(self, seqBarcode):
        """get the reverse complementary sequence."""
        dictBasePairs = {'A':'T', 'C':'G', 'T':'A', 'G':'C'}
        return ''.join([dictBasePairs[base] for base in seqBarcode[::-1]])


    def _checkInnerPair(self, seqBarcode):
        """check the inner pair of barcode.
        the minimum length of inner pair is self.innerPairLen, the loop length must be greater 
        than or equal to 3"""
        pairLen = self.innerPairLen
        loopLen = 3

        for i in range(len(seqBarcode)-pairLen*2-loopLen+1):
            seqForward = seqBarcode[i:i+pairLen]
            for j in range(i+pairLen+loopLen, len(seqBarcode)-pairLen+1):
                seqReverse = seqBarcode[j:j+pairLen]
                if seqForward == self._getRevCmp(seqReverse):
                    return False
        return True


    def _getExistSeq(self, ):
        """add the sequence that need to be excluded to the sub sequence dict.
        """
        existFwdSeq = ''
        existFwdSeq = existFwdSeq.upper()
        existRevSeq = self._getRevCmp(existFwdSeq)

        with self.lock:
            for i in range(0, len(existFwdSeq)-self.outerPairLen+1):
                self.dictSubSequence[existFwdSeq[i:i+self.outerPairLen]] = None
                self.dictSubSequence[existRevSeq[i:i+self.outerPairLen]] = None


    def _generateBarcode(self, ):
        """generate a barcode.
        """
        gcRate = random.uniform(self.minGCRate, self.maxGCRate)/2
        atRate = 0.5 - gcRate
        
        seqBarcode = ''
        while True:
            if len(seqBarcode) == self.barcodeLen:
                return seqBarcode
            
            rBase = random.choices(
                ('A', 'T', 'C', 'G'), 
                weights=(atRate, atRate, gcRate, gcRate), 
                k=1
            )[0]

            # exclude polymer
            if len(seqBarcode) >= self.polymerLen - 1:
                if rBase * (self.polymerLen-1) == seqBarcode[-self.polymerLen+1:]:
                    with self.lock:
                        self.dictErrorType['polymer error'] += 1
                    continue
            
            # exclude repeat structure
            flag = True
            for repeatUnit, repeatNum in self.dictRepeats.items():
                repeatLen = len(repeatUnit)*repeatNum
                repeatSeq = repeatUnit*repeatNum
                if len(seqBarcode) >= repeatLen-1:
                    if repeatSeq == seqBarcode[-repeatLen+1:]+rBase:
                        flag = False
                        with self.lock:
                            self.dictErrorType['repeat error'] += 1
                        break

            if not flag:
                continue

            seqBarcode += rBase


    def _writeBarcode(self, ):
        """write the generated barcode and edit distance to output files.
        """
        with open(self.outputSeqFile, 'w') as f:
            for i, seqBarcode in enumerate(self.lstGeneratedBarcode):
                f.write('>barcode_{:03} len={};GC={:.2f}\n'.format(
                    i+1, 
                    len(seqBarcode), 
                    (seqBarcode.count('G') + seqBarcode.count('C')) / len(seqBarcode)
                ))
                f.write(seqBarcode+'\n')

        with open(self.outputDistanceFile, 'w') as f:
            f.write('Barcode1\tBarcode2\tEditDistance\n')
            for i in range(len(self.lstGeneratedBarcode)):
                for j in range(i+1, len(self.lstGeneratedBarcode)):
                    distance = Levenshtein.distance(
                        self.lstGeneratedBarcode[i], self.lstGeneratedBarcode[j]
                    )
                    f.write('barcode_{:03}\tbarcode_{:03}\t{}\n'.format(i+1, j+1, distance))

        if self.outputSeqFile.endswith('.fa'):
            outputCandidateSeqFile = self.outputSeqFile.replace('.fa', '_candidate.fa')
        elif self.outputSeqFile.endswith('.fasta'):
            outputCandidateSeqFile = self.outputSeqFile.replace('.fasta', '_candidate.fasta')
        else:
            outputCandidateSeqFile = self.outputSeqFile + '_candidate.fa'
        with open(outputCandidateSeqFile, 'w') as f:
            for i, seqBarcode in enumerate(self.lstCandidateBarcode):
                f.write('>barcode_{:05}\n'.format(i+1))
                f.write(seqBarcode+'\n')


    def _runProcessTask(self, seed, lst, number, dictSubSeq):
        """run the process task to generate candidate barcode.
        Arguments:
            seed: int       seed for random
            lst: list       list for candidate barcode
            number: int     number of requried candidate barcode
        """
        random.seed(seed)

        self.logger.info('process start - {}'.format(seed))
        while True:
            # check if enough candidate sequences have been generated
            with self.lock:
                if len(lst) == number:
                    self.logger.info('process end - {}'.format(seed))
                    return
            
            seqBarcode = self._generateBarcode()

            if not self._checkGC(seqBarcode):
                with self.lock:
                    self.dictErrorType['GC rate error'] += 1
                continue
            if not self._checkInnerPair(seqBarcode):
                with self.lock:
                    self.dictErrorType['inner pair error'] += 1
                continue
            
            # check if there is no sub sequence repeat
            setSubSeq = set()
            for i in range(self.barcodeLen-self.outerPairLen+1):
                seq = seqBarcode[i:i+self.outerPairLen]
                seqRC = self._getRevCmp(seq)
                setSubSeq.add(seq)
                setSubSeq.add(seqRC)
                if seq in dictSubSeq or seqRC in dictSubSeq:
                    with self.lock:
                        self.dictErrorType['outer pair error'] += 1
                    break
            # if there is no sub sequence repeat, add it to the candidate sequence list
            else:
                with self.lock:
                    if len(lst) >= number:
                        self.logger.info('process end - {}'.format(seed))
                        return
                    for seq in setSubSeq:
                        dictSubSeq[seq] = None
                    lst.append(seqBarcode)


    def runDesignTaskParallel(self, ):
        """run the design task parallelly.
        
        Program flow:
        1. create multiple processes, each process generates some candidate barcodes;
        2. after generating enough candidate sequences, then filter them;
        3. calculate the edit distance between the candidate sequences, and filter out the 
          sequences with large differences.
        """
        lstProcess = []
        for i in range(self.numProcess):
            p = Process(
                target=self._runProcessTask, 
                args=(i, self.lstCandidateBarcode, self.numCandidateBarcode, self.dictSubSequence)
            )
            p.start()
            lstProcess.append(p)
        for p in lstProcess:
            p.join()

        self._filtEditDistance()
        self.logger.info('barcodes generated: {}'.format(len(self.lstGeneratedBarcode)))
        self._writeBarcode()
        print('error type:', self.dictErrorType, file=self.logFile)


def getLogger(name, level=logging.INFO):
    """get logger
    """
    logger = logging.getLogger(name)
    logger.setLevel(level)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    console_handler = logging.StreamHandler()
    console_handler.setLevel(level)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)
    return logger


def getParser():
    """get parser
    """
    parser = argparse.ArgumentParser(description='Barcode designer')
    parser.add_argument('--barcodeLen', type=int, required=True, help='Barcode length')
    parser.add_argument('--polymerLen', type=int, required=True, help='Polymer length')
    parser.add_argument('--minGCRate', type=float, required=True, help='Min GC ratio')
    parser.add_argument('--maxGCRate', type=float, required=True, help='Max GC ratio')
    parser.add_argument('--dictRepeats', type=str, required=True, help='Repeat dict, exp: "GC:3,CG:3,AT:3,TA:3"')
    parser.add_argument('--innerPairLen', type=int, required=True, help='Inner pair length')
    parser.add_argument('--outerPairLen', type=int, required=True, help='Outer pair length')
    parser.add_argument('--minEditDistance', type=int, required=True, help='Min edit distance')
    parser.add_argument('--numCandidateBarcode', type=int, required=True, help='Number of candidate barcode')
    parser.add_argument('--numGeneratedBarcode', type=int, required=True, help='Number of generated barcode')
    parser.add_argument('--numProcess', type=int, default=8, help='Number of process')
    parser.add_argument('--outputSeqFile', type=str, required=True, help='Output barcode sequence file')
    parser.add_argument('--outputDistanceFile', type=str, required=True, help='Output barcode distance info file')
    parser.add_argument('--logFile', type=str, required=True, help='Log file')
    args = parser.parse_args()
    return args


def main():
    logger = getLogger('BarcodeDesigner', logging.INFO)
    args = getParser()

    designer = BarcodeDesigner(
        logger, 
        args.barcodeLen, 
        args.polymerLen, 
        args.minGCRate, 
        args.maxGCRate, 
        args.dictRepeats, 
        args.innerPairLen, 
        args.outerPairLen, 
        args.minEditDistance, 
        args.numCandidateBarcode, 
        args.numGeneratedBarcode, 
        args.numProcess, 
        args.outputSeqFile, 
        args.outputDistanceFile, 
        args.logFile
    )
    designer.runDesignTaskParallel()


if __name__ == '__main__':
    main()