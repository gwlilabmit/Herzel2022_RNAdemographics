### removes linker sequence from Illumina sequencing results ###

import numpy as np
import os, subprocess
import sys

def removeLinker(fqIN,matchThres=80,linker='CTGTAGGCACCATCAAT'):
    inFile = open(fqIN, 'r')
    outFile = open(fqIN[:-4] + '-stripped.fq', 'w')
    dumpFile = open(fqIN[:-4] + '-notpassed.fq', 'w')
    line = inFile.readline()
    processed = 0
    passed = 0
    noinsert = 0
    while line != '':
        seq = inFile.readline()[:-1]
        line3 = inFile.readline()
        line4 = inFile.readline()
        length = findLinker(seq, linker, matchThres)
        if length > 3:
            outFile.write(line)
            outFile.write(seq[:length] + '\n')
            outFile.write(line3)
            outFile.write(line4[:length] + '\n')
            passed += 1
        elif length == -1:
            dumpFile.write(line)
            dumpFile.write(seq + '\n')
            dumpFile.write(line3)
            dumpFile.write(line4)
        else:
            noinsert += 1
        line = inFile.readline()
        processed += 1
    inFile.close()
    outFile.close()
    dumpFile.close()
    print 'Total # reads: '+str(processed)+'. Total # passed: '+str(passed)+'.'
    print 'Percent passed: '+str(passed*100/processed)+'%'
    print 'Percent no insert: '+str(noinsert*100/processed)+'%'
    return fqIN[:-4] + '-stripped.fq'

def findLinker(seq, linker, matchThres=80):
    exactmatch = seq.find(linker)
    if exactmatch + 1:
        return exactmatch
    else:
        s, l = [len(seq), len(linker)]
        maxpercent, maxindex = [0, 0]
        for i in range(s-1): #added -1 on 11/2/11
            percentmatch = sum([seq[i+j]==linker[j] for j in range(min(l, s-i))])*100/min(l, s-i)
            if percentmatch > maxpercent:
                maxpercent = percentmatch
                maxindex = i
        if maxpercent > matchThres:
            return maxindex
        else:
            return -1

def stepAlign(fqIN, refs = ['ecoli/MG-rrna', 'ecoli/MG-chr']):
    steps = len(refs)
    fq = fqIN
    for i in range(steps):
        alignedFile = fq[:fq.find('.')] + refs[i][refs[i].find('-'):] + '.align'
        unalignedFile = fq[:fq.find('.')] + refs[i][refs[i].find('-'):] + '.unalign'
        output=subprocess.Popen('bowtie -v 1 -k 1 --quiet /home/herzel/annotations/'+refs[i]+' '+fq+' '+alignedFile+' --un '+unalignedFile, shell=True, stdout=subprocess.PIPE) # adjusted by LH 09/05/2018, changed v=1 to v=0
        output.wait()
	fq = unalignedFile
	print fq
        

if __name__ == '__main__':
    linker='CTGTAGGCACCATCAAT'
    fqIN = sys.argv[1]
    fqOUT = removeLinker(fqIN, 70, linker)
    fqOUT = fqIN[:-4] + '-stripped.fq'


