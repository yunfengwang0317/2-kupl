#/bin/bash
from Bio import SeqIO
import gzip,sys,os,subprocess,time
import pandas as pd
import numpy as np
import itertools
from tempfile import *
from tqdm import *
import subprocess,sys, getopt,os,time,datetime
from multiprocessing import Process,Pool,Manager,Value
#set default parameters
mutant1,mutant2,wild1,wild2='','','',''
outdir='mutcaller_result'
datatype='WES'
method='hamming'
distance=1
lowdepth=20
cutoff=0.05
threads=20
support=5
nb_kmers=1
genome='None'
rmLQB='cutadapt'
help=''
#output USAGE information
def usage():
    print ('''
    USAGE
    MANDATORY PARAMETERS:
    --m1 mutant1\tThe compressed fq format file of mutant type library\n
    --m2 mutant2\tThe compressed fq format file of the second mutant type library in case the data is paired-end\n
    --w1 wild1\tThe compressed fq format file of wild type library\n
    --w2 wild2\tThe compressed fq format file of the second wild type library in case the data is paired-end\n
    -o outdir\tThe directory of output results (default mutcaller_result)

    OPTIONAL PARAMETERS:
    --type datatype\tThe data type of input, RNAseq, WES or genome (default WES)\n.
    -m method\tThe algorithm to calculate the distance (edit or hamming), which is a main parameter for bbduk (default hamming)\n
    -d distance\tThe max distance between kmer and reads (default 1)\n
    -l lowdepth\tThe contigs whose numbers of corresponding reads retrieved by bbduk are lower than this parameter will be removed (default 20)\n
    -c cutoff\tCutoff of relative coverage and contigs with AF larger than this value will be regarded as variants(default 0.05)\n
    -t thread\tThe number of threads used in searching pair kmers and variant calling process(default 20)\n
    -f filter\tMethod used to remove low quanlity bases (default cutadapt)\n
    --sp support\tKmers whose counts are larger than this value and absent in the control or reference will be used as case specific kmers (default 5)\n
    --nbk nb_kmers\tContigs with nb_kmers larger than this value will be used in the downstream analysis (default 1)\n
    --genome genome\tIf genome.fa provided, the contigs in vcf file will be replaced with corresponding chrom and position on genome (default None)\n
    --help \tHelp information
     ''')
if len(sys.argv)<2:
    usage()
    sys.exit(1)
try:
    opts, args = getopt.getopt(sys.argv[1:], "ho:m:d:l:c:t:f:",['m1=','m2=','w1=','w2=','type=','sp=','nbk=','genome=','help'])
    for op,value in opts:
        if op=='--m1':
            rawmutant1=value
        elif op=='--m2':
            rawmutant2=value
        elif op=='--w1':
            rawwild1=value
        elif op=='--w2':
            rawwild2=value
        elif op=='-o':
            outdir=value
        elif op=='--type':
            datatype=value
        elif op=='-m':
            method=value
        elif op=='-d':
            distance=value
        elif op=='-l':
            lowdepth=value
        elif op=='-c':
            cutoff=value
        elif op=='-t':
            threads=value
        elif op=='-f':
            rmLQB=value
        elif op=='--sp':
            support=int(value)
        elif op=='--nbk':
            nb_kmers=value
        elif op=='--genome':
            genome=value
        elif op=='-h' or '--help':
            usage()
            sys.exit(1)
except getopt.GetoptError:
    usage()
    sys.exit(1)

casename=rawmutant1.split('/')[-1].split('_')[0]
controlname=rawwild1.split('/')[-1].split('_')[0]
PAIRED,SINGLE=False,False
if rawmutant2!='':
    PAIRED=True
elif rawmutant2=='':
    SINGLE=True

def trimLBQ(fq1,fq2):
    if os.path.exists(outdir)is False:os.system('mkdir %s'%outdir)
    if os.path.exists(outdir+'/'+fq1.split('/')[-1]) and os.path.getsize(outdir+'/'+fq1.split('/')[-1])>50000000 or os.path.exists(outdir+'/variant_result/SNV_alignments.vcf'):return True
    subprocess.call(r'''cutadapt --cores=0 -q 10,10 -m 31 --pair-filter=any -o %s/%s -p %s/%s %s %s'''%(outdir,fq1.split('/')[-1],outdir,fq2.split('/')[-1],fq1,fq2),shell=True)

def OverrideN(fqfile):
    if os.path.exists(outdir)is False:os.system('mkdir %s'%outdir)
    if int(os.popen('ls %s/*fa.gz|wc -c'%outdir).readline().strip())>0 or os.path.exists(outdir+'/variant_result/SNV_alignments.vcf'):return True
    fi=NamedTemporaryFile(delete=True)
    seed=fi.name.split('/')[-1]
    export=[]
    with open(fqfile, "r") as handle:
        for record in tqdm(list(SeqIO.parse(handle,'fastq'))):
            qual=np.array(record.letter_annotations["phred_quality"])
            if min(qual)>10:
                export.append('>%s\n'%record.id)
                export.append(str(record.seq)+'\n')
                continue
            reads=record.seq.tomutable()
            for i in np.where(qual<10)[0]:reads[int(i)]='N'
            export.append('>%s'%record.id+'\n')
            export.append(''.join(reads)+'\n')
    out=open(fqfile+'.fa','w')
    for i in export:out.write(i)
    out.close()

def generate_kmers():
    fpath='%s/case_specific_kmers/%s.txt.gz'%(outdir,casename)
    if os.path.isfile(fpath) and os.path.getsize(fpath) > 1000:return True
    os.system('mkdir -p %s/case_specific_kmers'%outdir)
    jfdir=os.popen('which jellyfish').readline().strip()
    jfdir='jellyfish'
    if PAIRED:
        controlsh="%s count -m 31 -s 100000 -t 10 -o %s/case_specific_kmers/%s.jf -F 2 -C <(gunzip -c %s) <(gunzip -c %s)\n%s dump -c %s/case_specific_kmers/%s.jf | sort -k 1 -S 3000M --parallel 10|pigz -p 10 -c > %s/case_specific_kmers/%s.txt.gz\n"%(jfdir,outdir,controlname,wild1,wild2,jfdir,outdir,controlname,outdir,controlname)
        casesh="%s count -m 31 -s 100000 -t 10 -o %s/case_specific_kmers/%s.jf -F 2 -C -L %s <(gunzip -c %s) <(gunzip -c %s)\n%s dump -c %s/case_specific_kmers/%s.jf | sort -k 1 -S 3000M --parallel 10|pigz -p 10 -c > %s/case_specific_kmers/%s.txt.gz\n"%(jfdir,outdir,casename,support,mutant1,mutant2,jfdir,outdir,casename,outdir,casename)
    elif SINGLE:
        controlsh="%s count -m 31 -s 100000 -t 10 -o %s/case_specific_kmers/%s.jf -F 2 -C <(gunzip -c %s)\n%s dump -c %s/case_specific_kmers/%s.jf | sort -k 1 -S 3000M --parallel 10|pigz -p 10 -c > %s/case_specific_kmers/%s.txt.gz\n"%(jfdir,outdir,controlname,wild1,jfdir,outdir,controlname,outdir,controlname)
        casesh="%s count -L 5 -m 31 -s 100000 -t 10 -o %s/case_specific_kmers/%s.jf -F 2 -C -L %s <(gunzip -c %s)\n%s dump -c %s/case_specific_kmers/%s.jf | sort -k 1 -S 3000M --parallel 10|pigz -p 10 -c > %s/case_specific_kmers/%s.txt.gz\n"%(jfdir,outdir,casename,support,mutant1,jfdir,outdir,casename,outdir,casename)
    o1=open(outdir+'/control_jf.sh','w')
    o2=open(outdir+'/case_jf.sh','w')
    o1.write(controlsh)
    o2.write(casesh)
    o1.close()
    o2.close()

def runjf(x):
    if x=='case':fpath='%s/case_specific_kmers/%s.txt.gz'%(outdir,casename)
    else:fpath='%s/case_specific_kmers/%s.txt.gz'%(outdir,controlname)
    if os.path.isfile(fpath) and os.path.getsize(fpath) > 1000:
        return True
    else:os.system('bash %s/%s_jf.sh'%(outdir,x))

def variant_call():
    fpath=outdir+'/variant_result/SNV_alignments.txt'
    if os.path.isfile(fpath) and os.path.getsize(fpath)>1000:return True
    os.system("mkdir %s/merged_contigs %s/contig_pair %s/contig_unpair %s/variant_result"%(outdir,outdir,outdir,outdir))
    print ("python3 variant_call.py %s %s/case_specific_kmers/%s.txt.gz %s/case_specific_kmers/%s.txt.gz %s %s %s %s %s %s %s"%(threads,outdir,casename,outdir,controlname,wild1,wild2,lowdepth,cutoff,support,nb_kmers,distance))
    os.system("python3 variant_call.py %s %s/case_specific_kmers/%s.txt.gz %s/case_specific_kmers/%s.txt.gz %s %s %s %s %s %s %s"%(threads,outdir,casename,outdir,controlname,rawwild1,rawwild2,lowdepth,cutoff,support,nb_kmers,distance))

def generate_vcf():
    fpath=outdir+'/variant_result/SNV_alignments.vcf'
    if os.path.isfile(fpath) and os.path.getsize(fpath)>1000:return True
    print ('python3 createVCF.py %s %s %s %s %s'%(os.path.realpath(outdir+'/variant_result/SNV_alignments.txt'),cutoff,lowdepth,threads,genome))
    os.system('python3 createVCF.py %s %s %s %s %s'%(os.path.realpath(outdir+'/variant_result/SNV_alignments.txt'),cutoff,lowdepth,threads,genome))
    #if os.path.exists(fpath):subprocess.call(r'''rm %s %s %s %s'''%(wild1,wild2,mutant1,mutant2),shell=True)
def compressfile(samid):
    if os.path.exists('%s/fa.gz'%outdir+'/'+samid):return True
    os.system('cat %s/x*%s.fa|gzip > %s.fa.gz;rm %s/x*%s*'%(outdir,samid,outdir+'/'+samid,outdir,samid))
if __name__ =='__main__':
    os.system('rm -r %s/variant_result'%outdir)
    if rmLQB=='cutadapt':
        [wild1,wild2,mutant1,mutant2]=[outdir+'/'+i.split('/')[-1] for i in [rawwild1,rawwild2,rawmutant1,rawmutant2]]
        pool=Pool(4)
        pool.starmap(trimLBQ,zip([rawmutant1,rawwild1],[rawmutant2,rawwild2]))
        pool.close()
        pool.join()
    else:
        [wild1,wild2,mutant1,mutant2]=[outdir+'/'+i.split('/')[-1].replace('fq.gz','fa.gz').replace('fastq.gz','fa.gz') for i in [rawwild1,rawwild2,rawmutant1,rawmutant2]]
        for fi in [rawwild1,rawwild2,rawmutant1,rawmutant2]:
            samid=fi.split('/')[-1].replace('.fq.gz','').replace('.fastq.gz','')
            if int(os.popen('ls %s/*fa.gz|wc -c'%outdir).readline().strip())>0:continue
            subprocess.call(r'''split -l %s -d --additional-suffix=%s <(zcat %s);mv x*%s %s'''%(20000000,samid,fi,samid,outdir),executable='/bin/bash',shell=True)
            fileidxs=[i.strip() for i in os.popen('ls %s/x*%s'%(outdir,samid)).readlines()]
            pool=Pool(int(threads))
            pool.map(OverrideN, fileidxs)
            pool.close()
            pool.join()
        pool=Pool(4)
        pool.map(compressfile,[i.split('/')[-1].replace('.fq.gz','').replace('.fastq.gz','') for i in [rawwild1,rawwild2,rawmutant1,rawmutant2]])
        pool.close()
        pool.join()
    print (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),'generate kmers')
    generate_kmers()
    pool=Pool(2)
    pool.map(runjf, ['case','control'])
    pool.close()
    pool.join()
    print (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),'calling variants')
    variant_call()
    print (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),'generate VCF')
    generate_vcf()
