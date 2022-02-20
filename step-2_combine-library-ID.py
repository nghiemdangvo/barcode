# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 15:50:51 2019

@author: Du Jiang

tesed operating system: Windows 7, Windows 10

task:
    1. quality check based on PCR handle sequence
    2. demultiplex based on virus library IDs

input:
    1. step 1 results, i.e. txt files named as 'NNNNNN_mmddyy.txt'
    2. sample information file
    3. "library_ID.txt" file containing library IDs, tab-delimited

variables subject to change:
    1. step1_location: directory where the step1 results are stored. default: 'step-1'
    2. sample_info: tab-delimited txt files containing sample information.
                    see "readme.pdf" file on how to prepare it
    3. distance_allowed: Levenshtein edit distance allowed when determining whether two sequences are legitimately the same
    4. step2_output: folder to save output files. default 'step-2'

output:
    1. bin files storing results, which is a dictionary with unique sequences as keys and copy number as items.
        these bin files will be the input for step-3
    2. stat.txt files storing quality report for each sequencing sample

"""

#variables subject to change
step1_location = r'step-1'
sample_info = r'sample info 041519.txt'
distance_allowed = 4
step2_output = r'step-2'

#constant
library_ID_file = 'library_ID.txt'
end_seq = 'AGATCGGAAGAGCTCGT'

#packages
import os
import sys
import time
import pickle
import multiprocessing as mp
from Levenshtein import distance

def main():
    #check if the output folder exists. If not, will create one
    if not os.path.exists(step2_output):
        os.makedirs(step2_output)
    
    #read "library_ID.txt" file
    with open(library_ID_file, 'r') as rf:
        vlines = rf.readlines()
    virus = {} #key:library ID; item:virus library number
    for i in range(len(vlines)):
        virus[vlines[i].split('\t')[1].strip()] = vlines[i].split('\t')[0].strip()
    print ("There are"), len(vlines), "virus libraries;"
    
    #read sample info file
    try:
        rf = open(sample_info, 'r')
    except IOError:
        print ("ERROR: can't find sample info file.")
        sys.exit()
    multiplex = {} #key:multiplex index; item:[targeted libraries]
    sample = {} #key:[mutiplex index, targeted library]; item:sample name
    slines = rf.readlines()
    rf.close()
    print ("There are"), len(slines)-1, "samples to process;"
    for i in range(1, len(slines)):
        sline = slines[i].strip().split('\t')
        if not sline[1] in multiplex.keys():
            multiplex[sline[1]] = []
        multiplex[sline[1]].append(sline[3])
        sample[sline[1], sline[3]] = sline[2]
    print (sample)
    print (multiplex)
    
    #sort virus for better visualization later
    sorted_virus = []
    for key in virus:
        sorted_virus.append(int(virus[key]))
    sorted_virus.sort()
    
    #read step 1 results
    file_numbers = 0
    for fil in os.listdir(step1_location):
        #process each sequencing sample
        if not fil.endswith('.txt'):
            continue
        fname = fil.split('_')[0] #fname is the multiplex index
        if not fname in multiplex:
            print ("NO %s sample info, skipped.") % fname
            continue
        file_numbers += 1
        print ("\nProcessing:"), fname, multiplex[fname]
        
        #record the library IDs that are expected from this sequencing sample
        expected_virus = {} #key:6bp library ID; item:corresponding library number 1-24
        for v in virus.keys():
            if virus[v] in multiplex[fname]:
                expected_virus[v] = virus[v]
        
        #initiate a few dictionaries to store the result
        valid_virus_dic = {} #key:expected library number; item:dict(barcode: copy number)
        all_virus_dic = {} #key:all library number 1-24; item:number of total reads in that lib
        for v in expected_virus:
            valid_virus_dic[str(expected_virus[v])] = {}
        for v in virus:
            all_virus_dic[str(virus[v])] = 0
        
        #--------------------Parallel Computing.
        
        #Create a pool of processes
        pool = mp.Pool(mp.cpu_count())
        
        #Use the pool to apply function Map() to 1,000,000 sized chunks of reads
        results = pool.imap(Map, generateChunk(os.path.join(step1_location, fil)))
        
        #Make sure processes finish before continuing
        pool.close()
        pool.join()
        
        #--------------------End Parallel processing
        
        #further demultiplex reads based on their library IDs
        print ("\nOrganizing barcodes by library ID ...")
        ctall,ctend = Reduce(valid_virus_dic, all_virus_dic, results, expected_virus, virus)
        print ("Done!\n")
        
        #create bin files for each cell source
        for v in valid_virus_dic:
            fp = open(os.path.join(step2_output, sample[fname, v] + '.bin'), 'wb')
            pickle.dump(valid_virus_dic[v], fp)
            fp.close()
        
        #write stat files
        wf = open(os.path.join(step2_output, fname + '_stats.txt'), 'w')
        wf.writelines('library ID\tsample\tnumber of reads\tnumber of unique reads\tnote\n')
        all_virus,valid_virus = 0,0
        for i in range(len(sorted_virus)):
            wf.writelines(str(sorted_virus[i]) +'\t')
            if str(sorted_virus[i]) in multiplex[fname]:
                valid_virus += all_virus_dic[str(sorted_virus[i])]
                all_virus += all_virus_dic[str(sorted_virus[i])]
                wf.writelines(sample[fname,str(sorted_virus[i])] +'\t'+ str(all_virus_dic[str(sorted_virus[i])]) +'\t'+ str(len(valid_virus_dic[str(sorted_virus[i])])) +'\t\n')
            else:
                all_virus += all_virus_dic[str(sorted_virus[i])]
                wf.writelines('no sample\t'+ str(all_virus_dic[str(sorted_virus[i])]) +'\t\t')
                exp_val = int(1.0*ctend/(4**6))
                if all_virus_dic[str(sorted_virus[i])] > exp_val:
                    wf.writelines('higher than expected value: '+ str(exp_val) + '\n')
                else:
                    wf.writelines('\n')
        
        wf.writelines('\nTotal reads:\t' + str(ctall) + '\n')
        wf.writelines('\nReads with the expected 17bp ending (allowing 2bp misreads):\t' + str(ctend))
        wf.writelines('\n% valid reads based on 17bp ending:\t' + str(100.0*ctend/ctall) + '%\n')
        wf.writelines('\nReads with all virus library IDs:\t' +str(all_virus))
        wf.writelines('\n% reads with virus library IDs:\t' + str(100.0*all_virus/ctall) + '%\n')
        wf.writelines('\nReads with expected virus library IDs:\t' + str(valid_virus))
        wf.writelines('\n% reads with expected virus library IDs:\t' + str(100.0*valid_virus/ctall) + '%\n')
        wf.close()
        
        #print on screen
        print ('\nTotal reads:', ctall)
        print ('Reads with expected 17bp ending (allowing 2bp misreads):', ctend)
        print ('% valid reads based on 17bp ending:', str(100.0*ctend/ctall)+'%')
        print ('Reads with all virus library IDs:', all_virus)
        print ('% reads with virus library IDs:', str(100.0*all_virus/ctall)+'%')
        print ('Reads with expected virus library IDs:', valid_virus)
        print ('% reads with expected virus library IDs:', str(100.0*valid_virus/ctall)+'%')
        print ('\n')
    
    #final
    print ('ALL DONE!')

def generateChunk(path):
    """
    input the step-1 result files containing reads
    iterate and produce chunks of one million reads for parallel processing
    """
    rf = open(path, 'rU')
    eof = False
    while True:
        results = []
        for iii in range(1000000):
            line = rf.readline().strip()
            if line:
                results.append(line)
            else:
                eof = True
                break
        yield results
        if eof:
            break
    rf.close()

def Map(reads):
    """
    input a chunk of reads
    retain the reads whose last 17bp matches predesigned sequence in a temporary file
    """
    #create a temporary file
    pid = str(mp.current_process().pid)
    filename = 'temporary_' + pid + '_' + str(time.time()) + '.txt'
    wf = open(filename, 'w')
    #process reads
    count = 0 #count the number of reads from the input chunk
    for read in reads:
        count += 1
        #only retain reads whose Levenshtein distance between its last 17bp and 
        #predesigned sequence is less than the variable "distance_allowed"
        if distance(end_seq, read[-17:]) <= distance_allowed:
            wf.write(read + '\n')
        if count % 100000 == 0:
            print ("Process %s: %d is done!" %(pid, count))
    wf.close()
    #return the temporary files containing reads whose last 17bp matches predesigned sequence,
    #and the number of reads in the input chunk
    return (filename, count)

def Reduce(d1, d2, L, vlist, v2):
    """
    demultiplex reads based on its first 6bp:
        1. if the 6bp does not exactly match any of the library IDs in library_ID_file, the read will be excluded
        2. if the 6bp matches exactly an expected library ID, the read will be recorded and counted in dictionary d1
        3. if the 6bp matches exactly an unexpected library ID, the read will be excluded, but will be counted in dictionary d2
    """
    count_end = 0 #count the number of reads whose last 17bp matches predesigned sequence
    count_all = 0 #count the number of total reads
    count_p = 0 #count the number of chunks for monitoring progress
    for pair in L:
        filename, count = pair
        count_all += count
        for read in [line.strip() for line in open(filename, 'r')]:
            count_end += 1
            if read[:6] in vlist:
                if read in d1[vlist[read[:6]]]:
                    d1[vlist[read[:6]]][read] += 1
                else:
                    d1[vlist[read[:6]]][read] = 1
            if read[:6] in v2:
                d2[v2[read[:6]]] += 1
        os.remove(filename) #after finishing analysis, delete the temporary file
        count_p += 1
        print ("Done with %d million..." % count_p)
    return (count_all, count_end)

if __name__ == '__main__':
    main()