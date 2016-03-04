#! /usr/bin/env python
################################################################################
'''
Given a list of tag directories, summarizes the mapping statistics using the 
tagInfo and the log file
'''

### header ###
__author__ = "Jenhan Tao"
__license__ = "BSD"
__email__ = "jenhantao@gmail.com"


### imports ###
import sys 
import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt 
import fnmatch
import seaborn as sns 

### functions ###
def readTagInfoFile(filePath):
    with open(filePath) as f:
        data = f.readlines()
    genome = data[1].split()[0].split('=')[1]
    uniquePositions  =  int(data[1].split()[1])
    fragmentLengthEstimate = float(data[2].split('=')[1])
    tagsPerBP = float(data[4].split('=')[1])
    averageTagsPerPosition = float(data[5].split('=')[1])
    averageTagLength = float(data[6].split('=')[1])
    averageFragmentGCcontent = float(data[8].split('=')[1])
    return genome, uniquePositions, fragmentLengthEstimate, tagsPerBP, averageTagsPerPosition, averageTagLength, averageFragmentGCcontent

def readStarLog(filePath):
    with open(filePath) as f:
        data = f.readlines()
    totalReads = int(data[5].split()[5])
    uniquelyMappedReads = int(data[8].split()[5])
    multiMappedReads = int(data[23].split()[8])
    unmappedReads = totalReads - uniquelyMappedReads - multiMappedReads
    return totalReads, uniquelyMappedReads, multiMappedReads, unmappedReads

def readBowtieLog(filePath):
    with open(filePath) as f:
        data = f.readlines()
    totalReads = float(data[0].split()[0])
    unmappedReads = float(data[2].split()[0])
    uniquelyMappedReads = float(data[3].split()[0])
    multiMappedReads = float(data[4].split()[0])
    return totalReads, uniquelyMappedReads, multiMappedReads, unmappedReads

if __name__ == '__main__':
    ### check that there are sufficient number of arguments ###
    if len(sys.argv) < 4:
        print("Usage:")
        print("summarize_logs.py <rna|chip|atac|gro>> <output_directory> <tagDir1> ..<tagDirN>")
        sys.exit(1)

    ### parse arguments ###
    experimentType = sys.argv[1].lower()
    output_directory = sys.argv[2]
    input_directories = sys.argv[3:]
    
    if not os.path.isdir(output_directory):
        os.mkdirs(output_directory)

    # intitialize arrays
    _genome= []
    _uniquePositions= []
    _fragmentLengthEstimate= []
    _tagsPerBP= []
    _averageTagsPerPosition= []
    _averageTagLength= []
    _averageFragmentGCcontent= []
    _totalReads= []
    _uniquelyMappedReads= []
    _multiMappedReads= []
    _unmappedReads= []
    _sampleNames = []
    # for each tag directory
    for td in input_directories:
        ### read tag info file ###
        _sampleNames.append(td.split('/')[-1])
        genome, uniquePositions, fragmentLengthEstimate, tagsPerBP, averageTagsPerPosition, averageTagLength, averageFragmentGCcontent = readTagInfoFile(td + '/tagInfo.txt')
        for f in os.listdir(td):
            if fnmatch.fnmatch(f, '*.log'):
                logFile = f
                break
            
        ### read in mapping log file
        if experimentType == 'rna':
            totalReads, uniquelyMappedReads, multiMappedReads, unmappedReads = readStarLog(td + '/' + logFile)
        elif experimentType == 'chip':
            totalReads, uniquelyMappedReads, multiMappedReads, unmappedReads = readBowtieLog(td + '/' + logFile)
        elif experimentType == 'atac':
            totalReads, uniquelyMappedReads, multiMappedReads, unmappedReads = readBowtieLog(td + '/' + logFile)
        elif experimentType == 'gro':
            totalReads, uniquelyMappedReads, multiMappedReads, unmappedReads = readBowtieLog(td + '/' + logFile)
        else:
            pass
        _genome.append(genome)
        _uniquePositions.append(uniquePositions)
        _fragmentLengthEstimate.append(fragmentLengthEstimate)
        _tagsPerBP.append(tagsPerBP)
        _averageTagsPerPosition.append(averageTagsPerPosition)
        _averageTagLength.append(averageTagLength)
        _averageFragmentGCcontent.append(averageFragmentGCcontent)
        _totalReads.append(totalReads)
        _uniquelyMappedReads.append(uniquelyMappedReads)
        _multiMappedReads.append(multiMappedReads)
        _unmappedReads.append(unmappedReads)

    
    ### create data frame ####
    stat_dict = {'sample':_sampleNames,
    'genome':_genome,
    'uniquePositions':_uniquePositions,
    'fragmentLengthEstimate':_fragmentLengthEstimate,
    'tagsPerBP':_tagsPerBP,
    'clonality':_averageTagsPerPosition,
    'averageTagLength':_averageTagLength,
    'GC Content':_averageFragmentGCcontent,
    'totalReads':_totalReads,
    'uniquelyMappedReads':_uniquelyMappedReads,
    'multiMappedReads':_multiMappedReads,
    'unmappedReads':_unmappedReads
    }
    stat_frame = pd.DataFrame(stat_dict)
    stat_frame = stat_frame[
        [
    'sample',
    'genome',
    'uniquePositions',
    'fragmentLengthEstimate',
    'tagsPerBP',
    'clonality',
    'averageTagLength',
    'GC Content',
    'totalReads',
    'uniquelyMappedReads',
    'multiMappedReads',
    'unmappedReads'
        ]]
    stat_frame['uniquelyMappedFraction'] = stat_frame['uniquelyMappedReads'] / stat_frame['totalReads']
    stat_frame['mappedFraction'] = (stat_frame['uniquelyMappedReads'] + stat_frame['multiMappedReads']) / stat_frame['totalReads']
    stat_frame.to_csv(output_directory + '/mapping_stats.tsv',sep='\t', index=False)

    ### create plots ###
    sns.factorplot(data=stat_frame, y='clonality', kind='box')
    plt.savefig(output_directory + '/clonality_boxplot.pdf')
    plt.close()

    sns.factorplot(data=stat_frame, y='GC Content', kind='box')
    plt.savefig(output_directory + '/GC_Content_boxplot.pdf')
    plt.close()

    sns.factorplot(data=stat_frame, y='clonality', kind='box')
    plt.savefig(output_directory + '/clonality_boxplot.pdf')
    plt.close()

    sns.factorplot(data=stat_frame, y='uniquelyMappedReads', kind='box')
    plt.savefig(output_directory + '/uniquelyMappedReads_boxplot.pdf')
    plt.close()

    sns.factorplot(data=stat_frame, y='uniquelyMappedFraction', kind='box')
    plt.savefig(output_directory + '/uniquelyMappedFraction_boxplot.pdf')
    plt.close()

    sns.factorplot(data=stat_frame, y='mappedFraction', kind='box')
    plt.savefig(output_directory + '/mappedFraction_boxplot.pdf')
    plt.close()

    sns.regplot(data=stat_frame, x='totalReads', y='uniquelyMappedReads')
    plt.savefig(output_directory + '/sequencingDepth.pdf')
    plt.close()
