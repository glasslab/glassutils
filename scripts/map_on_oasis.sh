################################################################################

### SUMMARY OF FUNCTIONS ###
# One stop destination for all your mapping needs. This script, which is meant
# to be executed on the TSCC, will create qsub scripts (and execute them) to
# perform the following:
# 1. Copy raw data (fastq.gz) from the Glassome to the TSCC (without using qsub)
#    into an automatically generated directory on the Oasis file system
# 2. Decompress the files if necessary
# 3. Map the fastq files to a reference genome, producing a *.sam file
#    * using STAR or Bowtie2
# 4. Creates tag directories
# 5. Calculates the PBC coefficient - which involves removing all non-uniquely 
#    mapped reads, duplicate reads, and calculating the pileups
#    * read more about the PBC coefficient here: 
#      http://genome.ucsc.edu/ENCODE/qualityMetrics.html
###

#Copyright (c) 2015, Jenhan Tao 
#All rights reserved. 
#
#Redistribution and use in source and binary forms, with or without 
#modification, are permitted provided that the following conditions are met: 
#
# * Redistributions of source code must retain the above copyright notice, 
#   this list of conditions and the following disclaimer. 
# * Redistributions in binary form must reproduce the above copyright 
#   notice, this list of conditions and the following disclaimer in the 
#   documentation and/or other materials provided with the distribution. 
# * Neither the name of UC San Diego nor the names of its contributors may 
#   be used to endorse or promote products derived from this software 
#   without specific prior written permission. 
#
#THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
#AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
#IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
#ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
#LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
#CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
#SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
#INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
#CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
#ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
#POSSIBILITY OF SUCH DAMAGE. 



################################################################################
