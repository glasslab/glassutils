import sys, os, re, logging, subprocess, glob #GEOparse, 
from pysradb.sraweb import SRAweb
from datetime import datetime
db = SRAweb()

class suppress_output:
    def __init__(self, suppress_stdout=False, suppress_stderr=False):
        self.suppress_stdout = suppress_stdout
        self.suppress_stderr = suppress_stderr
        self._stdout = None
        self._stderr = None

    def __enter__(self):
        devnull = open(os.devnull, "w")
        if self.suppress_stdout:
            self._stdout = sys.stdout
            sys.stdout = devnull

        if self.suppress_stderr:
            self._stderr = sys.stderr
            sys.stderr = devnull

    def __exit__(self, *args):
        if self.suppress_stdout:
            sys.stdout = self._stdout
        if self.suppress_stderr:
            sys.stderr = self._stderr

def MsgHelp():
  print("Download public data from GEO")
  print("Please provide accession number")
  exit()
  
def MsgError(strMsg=""):
  print(strMsg)
  exit()
    
def dwOneGSM(gsm,uID,tryN):
  name_regex = r"[\s\*\?\(\),\.;]"
  strDIR = "%s_%s_%s"%("Supp",gsm.get_accession(),re.sub(name_regex, "_", gsm.metadata["title"][0]))
  if os.path.isdir(strDIR):
    extF = os.listdir(strDIR)
    if any(s.endswith('.gz') for s in extF) and not any(s.endswith('.sra') for s in extF):
      print("\tSkip! *.gz existed in %s"%strDIR)
      return()  
  try:
    with suppress_output(suppress_stdout=True, suppress_stderr=True):
      logging.getLogger("GEOparse").setLevel(logging.CRITICAL)
      gsm.download_SRA('%s@biogen.com'%uID,filetype="fastq",keep_sra=True)
  except Exception as e:
    print(e)
    if tryN<3:
      print("\tTry again! Max 3 times")
      dwOne(gsm,uID,tryN+1)
    else:
      print("Finished %s\n"%gsm.get_accession())
    
def downloadGEO(strGEO):
  srp = db.gse_to_srp(strGEO)
  for i in range(srp.shape[0]):
    print("\n***** ",srp.study_accession[i]," *****")
    downloadPRJNA(srp.study_accession[i])

def downloadPRJNA(strPRJNA,core=4):
  df = db.sra_metadata(strPRJNA,detailed=True)
  df.to_csv("%s.csv"%strPRJNA)
  print("Downloading all sra ...")
  if 'sra_url' in df.columns:
    RMdf = df[df["sra_url"].isna()]
    if(RMdf.shape[0]>0):
      print("\n\n******\nThe downloading links are missing for the following samples:\n\t%s\n"%';\n\t'.join(list(RMdf['experiment_title'])))
      df = df[~df["sra_url"].isna()]
  db.download(df=df,skip_confirmation=True,out_dir=os.getcwd())
  print("Obtain fastqs from sra ...")
  for i in range(df.shape[0]):
    srr=df.run_accession[i]
    srp=df.study_accession[i]
    srx=df.experiment_accession[i]
    strD=os.path.join(srp,srx)
    strTime = datetime.now().strftime("%Y %B %d %H:%M:%S")
    print("\t%s - %s:%s"%(strTime,srx,srr))
    
    if len(glob.glob("%s/*fastq.gz"%strD))>0:
      print("\t\tSKIP: fastq.gz exists")
      continue
    if not os.path.exists(os.path.join(strD,"%s.sra"%srr)):
      raise Exception("\t\tERROR: missing sra file in %s"%strD)
    cmd= "cd %s;fasterq-dump --threads %d --split-3 %s.sra;pigz -p %d *fastq"%(strD,core,srr,core)
    #cmd="cd %s/%s;fastq-dump --gzip --skip-technical --readids --split-3 %s.sra"%(srp,srx,srr)
    cmdR=subprocess.run(cmd,shell=True,check=True,stdout=subprocess.PIPE)
  print("Rename sample folder")
  name_regex = r"[\s\*\?\(\),\.;:\/]"
  for i in range(df.shape[0]):
    srp=df.study_accession[i]
    srx=df.experiment_accession[i]
    srxTitle=srx+"_"+re.sub(name_regex,"_",df.experiment_title[i])    
    strPath=os.path.join(srp,srx)
    if os.path.exists(strPath):
      os.rename(strPath,os.path.join(srp,srxTitle))
  print("Rename study folder")
  for i in range(df.shape[0]):
    srp=df.study_accession[i]
    srpTitle=srp+"_"+re.sub(name_regex,"_",df.study_title[i])    
    if os.path.exists(srp):
      os.rename(srp,srpTitle)

def main():
  if len(sys.argv)<2:
    MsgHelp()
  
  prjID = sys.argv[1]
  uID = sys.argv[2]
  print("Processing %s"%prjID)
  if prjID.startswith("GSE"):
    downloadGEO(prjID)
  elif prjID.startswith("PRJNA") or prjID.startswith("SRP"):
    downloadPRJNA(prjID)
  else:
    MsgError("Currently only support GSE/PRJNA accessions.")
  
  print("Thanks for using downloadSRA!")
  

if __name__ == "__main__":
  main()
