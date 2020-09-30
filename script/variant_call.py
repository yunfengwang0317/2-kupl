import os,sys,time,subprocess,re,gzip
from math import ceil
from tqdm import *
import pandas as pd
import numpy as np
from tempfile import *
import scipy.stats as stats
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from difflib import SequenceMatcher
from pydna.assembly import Assembly
from pydna.dseqrecord import Dseqrecord
from Bio.Seq import Seq
from collections import defaultdict
from functools import partial
from multiprocessing import Process,Pool,Manager,Value
from fsplit.filesplit import FileSplit
def dist(s1,s2,start=0):
    if len(s1)!=len(s2):raise ValueError('undefined for sequences of unequal length')
    hd=0
    for e1,e2 in zip(s1[start:],s2[start:]):
        if e1!=e2:hd+=1
        if hd>1:return 9
    return hd
    #return sum(chr1!=chr2 for chr1,chr2 in zip(s1,s2))

def sizeof_fmt(num, suffix='B'):
    for unit in ['','Ki','Mi','Gi','Ti','Pi','Ei','Zi']:
        if abs(num) < 1024.0:
            return "%3.1f %s%s" % (num, unit, suffix)
        num /= 1024.0
    return "%.1f %s%s" % (num, 'Yi', suffix)

def createindex(lock,fi):
    global pair_T,pair_N
    idx_head,idx_tail=defaultdict(lambda: defaultdict(list)),defaultdict(lambda: defaultdict(list))
    dat=pd.read_csv(fi,header=None,index_col=None,sep=' ')#the 1st line in kmercount table is considered to be the head by kmerFilter
    subwildkmers = dict(zip(dat.iloc[:,0],dat.iloc[:,1]))
    for p in subwildkmers.keys():
        idx_head[p[:15]]['ref'].append(p)
        idx_tail[p[-15:]]['ref'].append(p)
    for q in spekmer:
        idx_head[q[:15]]['mut'].append(q)
        idx_tail[q[-15:]]['mut'].append(q)
    subpair_T,subpair_N=pairkmers(idx_head,idx_tail)
    lock.acquire()
    pair_N+=subpair_N
    pair_T+=subpair_T
    lock.release()
    idx_head.clear()
    idx_tail.clear()
def pairkmers(idx_head,idx_tail): 
    subpair_T,subpair_N=[],[]
    for key in tqdm(idx_head.keys()):
        if idx_head[key]['ref']==[] or idx_head[key]['mut']==[]:continue
        for q in idx_head[key]['mut']:
            for j in idx_head[key]['ref']:
                if q==j:break#second check if kmer is common in both T&N
                if dist(q,j,15)==1:
                    subpair_T.append(q)
                    subpair_N.append(j)
    for key in tqdm(idx_tail.keys()):
        if idx_tail[key]['ref']==[] or idx_tail[key]['mut']==[]:continue
        for q in idx_tail[key]['mut']:
            for j in idx_tail[key]['ref']:
                if q==j:break
                if dist(q,j,0)==1:
                    subpair_T.append(q)
                    subpair_N.append(j)
    return (subpair_T,subpair_N)

def contig2kmer(contig,rev=False):
    kmerlst=[]
    if rev==False:contig=contig
    else:contig=str(Seq(contig).reverse_complement())
    for i in range(len(contig)-30):
        kmerlst.append(contig[i:(i+31)])
    return kmerlst
def maxoverlap(s1,s2):
    d = SequenceMatcher(None,s1,s2)
    pos_a, pos_b, size = d.find_longest_match(0, len(s1), 0, len(s2))
    return (pos_a,pos_b,size)

def assemDNA(lst):
    res = lst[0]+'N'*100
    for tarkmer in lst[1:]:
        pos_s,pos_t,size = maxoverlap(res, tarkmer)
        if size<10:continue
        res=res.replace(res[pos_s:pos_s+31-pos_t],tarkmer)
    return (res.strip('N'))

def CAP(contig,reads):
    randseed=np.random.randint(1000,1000000)
    fa=NamedTemporaryFile(delete=False)
    out=open(fa.name,'w')
    for i in reads:
        if contig in i:return ('NA')
        else:
            out.write('>%s\n'%i)
            out.write(i+'\n')
    out.close()
    subprocess.call("cap3 %s -x %s > /dev/null"%(fa.name,randseed),shell=True)
    if os.path.exists('%s.%s.contigs'%(fa.name,randseed)):infered_contigs_nb=int(os.popen("grep Contig %s.%s.contigs|wc -l"%(fa.name,randseed)).readline().strip())
    else:
        print (contig,reads,' failed to infer ref')
        infered_contigs_nb=0
    if infered_contigs_nb==0:return('')
    putative_ref = {}
    for line in open(r'%s.%s.contigs'%(fa.name,randseed)):
        if line.startswith('>'):
            key = line.split(' ')[0]    
            putative_ref[key] = '' 
        else:
            putative_ref[key]+=line.strip()
    ref_qual = {}
    for line in open(r'%s.%s.contigs.qual'%(fa.name,randseed)):
        if line.startswith('>'):
            key = line.split(' ')[0] 
            ref_qual[key] = []
        else:
            ref_qual[key]+=line.split()
    ref=[]
    bestref,bestscore='',0
    for k,v in ref_qual.items():
        score=np.mean(np.array(ref_qual[k]).astype(int))
        if score>bestscore:
            bestscore=score
            bestref=putative_ref[k]
        if score>50:ref.append(putative_ref[k])
    if len(ref)==0:#use the putative ref with largest score
        ref=[bestref]
    os.system('rm %s.%s*'%(fa.name,randseed))
    if len(ref)==1:putative_ref=(ref[0])
    elif len(ref)>=2:
        putative_ref=bestref
    plus,minus=maxoverlap(contig,putative_ref)[2],maxoverlap(str(Seq(contig).reverse_complement()),putative_ref)[2]
    if plus>minus:return putative_ref
    else:return str(Seq(putative_ref).reverse_complement())

def nestDict(set1,set2):
    idx=defaultdict(list)
    keys_h,keys_t=[],[]
    for q in set1:
        keys_h.append(q[:min(int(len(q)/2),60)])
        keys_t.append(q[max(-60,-int(len(q)/2)):])
    for p in set2:
        for k in keys_h:
            if len(idx[k])>50:continue
            if k in p:
                idx[k].append(p[max(0,p.index(k)-15):p.index(k)+len(k)*2+15])
                break
        for k in keys_t:
            if len(idx[k])>50:continue
            for k in p:
                idx[k].append(p[max(p.index(k)-len(k),0):p.index(k)+len(k)+15])
    return idx

def contig_read_hash(fi,cor_reads):
    with open(fi)as f:
        contigs=list(map(lambda x:x.strip(),f))
    contig_reads=defaultdict(list)
    with open(cor_reads) as f:
        while True:
            line=f.readline()
            if not line:break
            if line.startswith('>')is False:continue
            seq=f.readline().strip()
            querys=re.findall('([ATCG]+)=',line)
            if len(querys)==0:continue
            for query in querys:
                if len(contig_reads[query])>10:
                    continue
                query=query.split('=')[0]#useless
                if query not in contigs:continue
                seq_s1,query_s1,len_overlap1=maxoverlap(seq,query)
                seqrv=str(Seq(seq).reverse_complement())
                seq_s2,query_s2,len_overlap2=maxoverlap(seqrv,query)
                if len_overlap1>len_overlap2:seq,seq_s,query_s,len_overlap=seq,seq_s1,query_s1,len_overlap1
                else:seq,seq_s,query_s,len_overlap=seqrv,seq_s2,query_s2,len_overlap2
                if len_overlap<15:continue
                contig_reads[query].append(seq[max(0,seq_s-query_s-30):min(len(seq),seq_s+len(query)-query_s+30)])
    return (contig_reads)

def infer_ref(line):
    contig=line[0]
    kmers=contig2kmer(contig)
    sp_case,cov_case,sp_control,cov_control,refpairs=contig_sp_cov.loc[contig].tolist()
    refpairs=refpairs.split(',')
    if len(refpairs)==1:return (refpairs[0],sp_case,cov_case,sp_control,cov_control)
    try:
        refseq=Assembly([Dseqrecord(i) for i in refpairs],limit=15).assemble_linear(max_nodes=3)[0].seq.watson
        if maxoverlap(contig,refseq)[2]<15:refseq=assemDNA(refpairs)#for sake of low complexity sequences
    except:refseq=assemDNA(refpairs)
    if maxoverlap(contig,refseq)[2]<15:refseq=str(Seq(refseq).reverse_complement())
    return (refseq,sp_case,cov_case,sp_control,cov_control)

def infer_ref_unpair(line,unpair_reads_dict):
    contig,refseq=line[0],''
    kmers=contig2kmer(contig)
    sp_case,cov_case,sp_control,cov_control,refpairs=contig_sp_cov.loc[contig].tolist()
    refpairs=refpairs.split(',')
    related_reads=unpair_reads_dict[contig]
    refseq='NA'
    if len(refpairs)>2:#indels should have no more than 2 paired refs.(head and tail)
        try:
            refseq=Assembly([Dseqrecord(i) for i in refpairs],limit=15).assemble_linear(max_nodes=3)[0].seq.watson
        except:refseq=assemDNA(refpairs)
    if len(related_reads)>lowdepth/2 and len(refseq)<len(contig):refseq=CAP(contig,related_reads)
    if maxoverlap(contig,refseq)[2]<15:refseq=str(Seq(refseq).reverse_complement())
    return (refseq,sp_case,cov_case,sp_control,cov_control)

def ana_contigs(fi,paired=True):
    if fi.split('/')[-1].startswith('x00'):contigs=pd.read_csv(fi,header=0,index_col=None,sep='\t')
    else:contigs=pd.read_csv(fi,header=None,index_col=None,sep='\t')
    contigs.columns=['contig']
    if paired==False:unpair_reads_dict=contig_read_hash(fi,cor_reads)
    a,b,c,d,e=[],[],[],[],[]
    for i in trange(contigs.shape[0]):
        if paired==False:aa,bb,cc,dd,ee=infer_ref_unpair(contigs.loc[i],unpair_reads_dict)
        else:aa,bb,cc,dd,ee=infer_ref(contigs.loc[i])
        a.append(aa)
        b.append(bb)
        c.append(cc)
        d.append(dd)
        e.append(ee)
    contigs['putative_ref'],contigs['sp_case'],contigs['cov_case'],contigs['sp_control'],contigs['cov_control']=a,b,c,d,e#zip(*contigs.apply(infer_ref_unpair,1,args=(unpair_reads_dict,)))
    if paired==False:contigs.to_csv('%s/contig_unpair/result%s.csv'%(outdir,fi.split('/')[-1]),header=False,index=False,sep='\t')
    else:contigs.to_csv('%s/contig_pair/result%s.csv'%(outdir,fi.split('/')[-1]),header=False,index=False,sep='\t')

def filter_unpaired_contigs(fi):
    #fileter low depth contig_unpaired
    seed=NamedTemporaryFile(delete=True).name.split('/')[-1]
    out=open('%s/contig_unpair/passedcontig_%s.fa'%(outdir,seed),'w')
    out2=open('%s/contig_unpair/contigs_unpaired_%s'%(outdir,seed),'w')
    weird_contig=open('%s/contig_unpair/FailedToInferRef_%s.txt'%(outdir,seed),'w')
    mutkmerpool=kmerpair.index.tolist()
    contig_unpair=pd.read_csv(fi,header=None,index_col=None,sep='\t')
    if contig_unpair.shape[1]==4:contig_unpair.columns=['nb_kmer','contig','tag','Pvalue']
    elif contig_unpair.shape[1]==1:contig_unpair.columns=['contig']
    for contig in tqdm(contig_unpair.contig.tolist()):
        headtailkmers=[contig[:31],contig[-31:]]
        refkmers=kmerpair.reindex(headtailkmers)[1].dropna().tolist()
        if len(refkmers)<2:
            weird_contig.write(contig+'\n')
            continue
        out.write('>%s\n%s\n'%(contig,contig))
        out2.write(contig+'\n')
    out2.close()
    out.close()
    weird_contig.close()
def usedkmers(fi):
    tag=fi.split('/')[-1].strip('.txt.gz')
    subprocess.call(r"""less %s/contig_pair/usedkmers|awk '{print ">contig_"NR"\n"$1}' > %s/contig_pair/usedkmers_%s.fa"""%(outdir,outdir,tag),shell=True)
    subprocess.call(r"""jellyfish query -s %s/contig_pair/usedkmers_%s.fa %s -o %s/contig_pair/usedkmers_%s"""%(outdir,tag,fi.replace('.txt.gz','.jf'),outdir,tag),shell=True,executable='/bin/bash')
    kmers=pd.read_csv('%s/contig_pair/usedkmers_%s'%(outdir,fi.split('/')[-1].replace('.txt.gz','')),header=None,index_col=None,sep=' ')
    kmers_rv=pd.DataFrame({0:[str(Seq(i).reverse_complement()) for i in kmers[0]],1:kmers[1]})
    kmers=pd.concat([kmers,kmers_rv]).drop_duplicates()
    kmers.index=kmers[0]
    del kmers[0]
    kmers.to_csv('%s/contig_pair/usedkmers_%s'%(outdir,fi.split('/')[-1].replace('.txt.gz','')),header=False,index=True,sep=' ')

def OnlyKeepMaxRef():
    kmerpair=pd.read_csv('%s/contig_pair/kmerpair.csv'%outdir,header=None,index_col=None,sep='\t')
    kmerpair.columns=['mut','wild']
    wildkmercount=pd.read_csv('%s/contig_pair/usedkmers_%s'%(outdir,kmerfile_N.split('/')[-1].replace('.txt.gz','')),header=None,index_col=None,sep=' ')
    wildkmercount.columns=['wild','count']
    kmerpair_wildcount=pd.merge(kmerpair,wildkmercount,left_on='wild',right_on='wild',how='left')
    kmerpair_wildcount=kmerpair_wildcount.sort_values('count',ascending=False).groupby('mut').first()
    kmerpair_wildcount.to_csv('%s/contig_pair/kmerpair.csv'%outdir,header=False,index=True,sep='\t')

def shrink():
    contigs=pd.read_csv('%s/merged_contigs/contigs_allspkmers'%outdir,header=0,index_col=None,sep='\t')
    length=contigs.contig.apply(lambda x:len(x),1)
    contigs=contigs[(length>=31+nb_kmers) & (length<100)]
    contigs.to_csv('%s/merged_contigs/contigs_allspkmers'%outdir,header=True,index=False,sep='\t')

def comm(param):
    if param=='12':
        #print (r'''comm -12 <(zcat %s|awk '{if($2>%s){print $1}}') <(zcat %s |awk '{if($2>%s){print $1}}') > %s/case_specific_kmers/shared_kmers'''%(kmerfile_T,support,kmerfile_N,support,outdir))
        subprocess.call(r'''comm -12 <(zcat %s|awk '{if($2>%s){print $1}}') <(zcat %s |awk '{if($2>%s){print $1}}') > %s/case_specific_kmers/shared_kmers'''%(kmerfile_T,support,kmerfile_N,support,outdir),shell=True,executable='/bin/bash')
    elif param=='23':
        #print (r'''comm -23 <(zcat %s|awk '{if($2>%s){print $1}}') <(zcat %s |awk '{if($2>0){print $1}}') > %s/case_specific_kmers/specific_kmer'''%(kmerfile_T,support,kmerfile_N,outdir))
        subprocess.call(r'''comm -23 <(zcat %s|awk '{if($2>%s){print $1}}') <(zcat %s |awk '{if($2>0){print $1}}') > %s/case_specific_kmers/specific_kmer'''%(kmerfile_T,support,kmerfile_N,outdir),shell=True,executable='/bin/bash')
    elif param=='homo':
        subprocess.call(r'''zcat %s|awk '{if($2>%s){print $1}}' > %s/case_specific_kmers/shared_kmers'''%(kmerfile_N,lowdepth,outdir),shell=True,executable='/bin/bash')#adaptable to homo variant

def Cal_sp_cov(contigs):
    col1,col2=[],[]
    for contig in contigs:
        kmers=contig2kmer(contig)
        col1+=[contig]*len(kmers)
        col2+=kmers
    df=pd.DataFrame({'contig':col1,'kmers':col2})
    df['refs']=kmerpair.reindex(df.kmers)[1].tolist()
    df['sp_T']=mutkmers.reindex(df.kmers)[1].tolist()
    df['allel2_T']=mutkmers.reindex(df.refs)[1].tolist()
    df['sp_N']=wildkmers.reindex(df.kmers)[1].tolist()
    df['allel2_N']=wildkmers.reindex(df.refs)[1].tolist()
    rawdf=df.dropna().copy()
    df=df.groupby('contig').median()
    df['cov_T']=df.sp_T+df.allel2_T
    df['cov_N']=df.sp_N+df.allel2_N
    df=df[['sp_T','cov_T','sp_N','cov_N']].dropna().astype(int)
    df=df[(df.sp_T>=support)&(df.cov_T>=lowdepth)]
    df['refpairs']=rawdf.groupby('contig')['refs'].apply(lambda x:','.join(x))
    df_rv=df.copy()
    df_rv.index=[str(Seq(i).reverse_complement()) for i in df.index]
    df=pd.concat([df,df_rv])
    df=df.loc[~df.index.duplicated(keep='first')]
    return df

def RemoveFP_via_control(contigs):
    ext_kmers=defaultdict(list)
    for contig in contigs:
        if len(contig)>60:continue
        for c in 'ATCG':
            ext_kmers[contig].append(c+contig[:30])
            ext_kmers[contig].append(contig[-30:]+c)
    fa=NamedTemporaryFile(delete=False)
    ext_kmers_count=NamedTemporaryFile(delete=False)
    out=open(fa.name+'.fa','w')
    for i_ in ext_kmers.values():
        for i in i_:
            out.write('>%s\n'%i)
            out.write(i+'\n')
    out.close()
    subprocess.call(r"""jellyfish query -s %s.fa %s/case_specific_kmers/%s -o %s/variant_result/ext_kmers_count"""%(fa.name,outdir,kmerfile_N.split('/')[-1].replace('.txt.gz','.jf'),outdir),shell=True,executable='/bin/bash')
    
if __name__ == '__main__':
    threads,kmerfile_T,kmerfile_N,wild1,wild2,lowdepth,cutoff,support,nb_kmers,distance=sys.argv[1:]
    threads,lowdepth,cutoff,support,nb_kmers,distance=int(threads),int(lowdepth),float(cutoff),int(support),int(nb_kmers),int(distance)
    samid=kmerfile_T.split('/')[-1].replace('.txt.gz','')
    outdir=os.path.dirname(os.path.dirname(kmerfile_T))
    ################# extract case specific kmers #######################
    nb_kmers_eachthread=10000000
    fpath='%s/case_specific_kmers/shared_kmers_count'%outdir
    if os.path.exists(fpath) is False and os.path.exists('%s/variant_result/SNV_alignments.vcf'%outdir) is False:
        print (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),'extract specific kmers')
        pool=Pool(2)
        pool.map(comm,['23','homo'])# homo mode for calling homozygote and 12 mode for somatic
        pool.close()
        pool.join()
        if int(os.popen('wc -l %s/case_specific_kmers/specific_kmer'%outdir).readline().strip().split()[0])>50000000:os._exit(0)
        os.system("echo 'tag\tpvalue' > %s/case_specific_kmers/specific_kmer_fix"%outdir)
        subprocess.call(r"""awk '{print $1"\t0"}' %s/case_specific_kmers/specific_kmer >> %s/case_specific_kmers/specific_kmer_fix"""%(outdir,outdir),shell=True)
        subprocess.call(r"./mergeTags -k 31 -m 25 -n %s/case_specific_kmers/specific_kmer_fix 2>/dev/null|awk '{if($1>%s){print $0}}'|gzip -c > %s/merged_contigs/contigs.gz;gunzip -c %s/merged_contigs/contigs.gz > %s/merged_contigs/contigs_allspkmers"%(outdir,nb_kmers,outdir,outdir,outdir),shell=True)
        subprocess.call(r"""less %s/case_specific_kmers/shared_kmers|awk '{print ">kmer"NR"\n"$1}' > %s/case_specific_kmers/shared_kmers.fa"""%(outdir,outdir),shell=True)
        subprocess.call(r"""jellyfish query -s %s/case_specific_kmers/shared_kmers.fa %s/case_specific_kmers/%s -o %s/case_specific_kmers/shared_kmers_count"""%(outdir,outdir,kmerfile_N.split('/')[-1].replace('.txt.gz','.jf'),outdir),shell=True)
    #shrink the contigs_allspkmers and remove useless kmers from the specific_kmer_fix
    if os.path.exists('%s/contig_pair/contigs_pairedspkmers'%outdir) is False:
        shrink()
        #print ("""rm %s/x*_N;split -l %s -d --additional-suffix=%s_N %s/case_specific_kmers/shared_kmers_count;mv x*%s_N %s"""%(outdir,min(10000000,nb_kmers_eachthread),samid,outdir,samid,outdir))
        os.system("""rm %s/x*_N;split -l %s -d --additional-suffix=%s_N %s/case_specific_kmers/shared_kmers_count;mv x*%s_N %s"""%(outdir,min(10000000,nb_kmers_eachthread),samid,outdir,samid,outdir))
        fileidxs=[i.strip() for i in os.popen('ls %s/x*%s_N'%(outdir,samid)).readlines()]
        print (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),'indexing and pairing kmers')
        with open('%s/case_specific_kmers/specific_kmer_fix'%outdir)as f:
            spekmer=list(map(lambda x:x.strip().split()[0],f))[1:]#the 1st line is tag\tpvalue
        spekmer=spekmer+[str(Seq(i).reverse_complement()) for i in spekmer]#add revcomp to the specific kmer list, in case kmerpair misses true pairs.remove and test
        with Manager() as manager:
        ################# pairing kmers from case and control #####################
            lock=manager.Lock()
            global pair_T,pair_N
            pair_T=manager.list()
            pair_N=manager.list()
            ncores=min(threads,len(fileidxs))
            pool=Pool(ncores)
            singlethread = partial(createindex, lock)
            pool.map(singlethread, fileidxs)
            pool.close()
            pool.join()
            '''
            print (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),'adding mut kmers into the hash')#what's the purpose of this part?
            idx_head,idx_tail=defaultdict(lambda: defaultdict(list)),defaultdict(lambda: defaultdict(list))
            for q in spekmer:
                idx_head[q[:15]]['mut'].append(q)
                idx_tail[q[-15:]]['mut'].append(q)
            subpair_T,subpair_N=pairkmers(idx_head,idx_tail)
            pair_T+=subpair_T
            pair_N+=subpair_N
            print (len(subpair_T))
            '''
            np.savetxt('%s/contig_pair/pair.txt'%outdir, np.stack([pair_T,pair_N]).T,delimiter="\t",fmt="%s")
            out=open('%s/contig_pair/kmer_T'%outdir,'w')
            kmerpair=pd.read_csv('%s/contig_pair/pair.txt'%outdir,header=None,index_col=None,sep='\t')
            pair_T_rv=[str(Seq(i).reverse_complement()) for i in pair_T]
            for p in set(pair_T).difference(set(pair_T_rv)):out.write(p+'\n')
            out.close()
            kmerpair_rv=pd.DataFrame({0:[str(Seq(i).reverse_complement()) for i in kmerpair[0]],1:[str(Seq(i).reverse_complement()) for i in kmerpair[1]]})
            kmerpair=pd.concat([kmerpair,kmerpair_rv]).drop_duplicates()
            kmerpair.index=kmerpair[0]
            del kmerpair[0]
            kmerpair.to_csv('%s/contig_pair/kmerpair.csv'%outdir,header=False,index=True,sep='\t')
            print (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),'shrink useful kmers from kmerpair')
            out=open('%s/contig_pair/usedkmers'%outdir,'w')
            for i in set(list(pair_T)+list(pair_N)):
                out.write(i+'\n')
            out.close()
            pool=Pool(2)
            pool.map(usedkmers, [kmerfile_T,kmerfile_N])
            pool.close()
            pool.join()
            OnlyKeepMaxRef()
            del pair_T_rv,spekmer,kmerpair_rv
    wildkmers=pd.read_csv('%s/contig_pair/usedkmers_%s'%(outdir,kmerfile_N.split('/')[-1].replace('.txt.gz','')),header=None,index_col=0,sep=' ')
    mutkmers=pd.read_csv('%s/contig_pair/usedkmers_%s'%(outdir,kmerfile_T.split('/')[-1].replace('.txt.gz','')),header=None,index_col=0,sep=' ')
    contig_all=pd.read_csv('%s/merged_contigs/contigs_allspkmers'%outdir,header=0,index_col=None,sep='\t')
    kmerpair=pd.read_csv('%s/contig_pair/kmerpair.csv'%outdir,header=None,index_col=0,sep='\t')
    contig_sp_cov=Cal_sp_cov(contig_all.contig)
    del wildkmers,mutkmers
    if os.path.exists('%s/contig_pair/contigs_pairedspkmers'%outdir) is False:
        print (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),'merging into paired contigs')
        os.system("echo 'tag\tpvalue' > %s/contig_pair/specific_kmer_fix"%outdir)
        subprocess.call(r"""awk '{print $1"\t0"}' %s/contig_pair/kmer_T >> %s/contig_pair/specific_kmer_fix"""%(outdir,outdir),shell=True)
        subprocess.call(r"./mergeTags -k 31 -m 25 -n %s/contig_pair/specific_kmer_fix 2>/dev/null|awk '{if($1>%s){print $0}}'|gzip -c > %s/contig_pair/contigs.gz;gunzip -c %s/contig_pair/contigs.gz |cut -f 2 > %s/contig_pair/contigs_pairedspkmers"%(outdir,nb_kmers,outdir,outdir,outdir),shell=True)
    contig_all=contig_all[contig_all.contig.isin(contig_sp_cov.index)]
    contig_pair=pd.read_csv('%s/contig_pair/contigs_pairedspkmers'%outdir,header=0,index_col=None,sep='\t')
    contig_pair=contig_pair[contig_pair.contig.isin(contig_sp_cov.index)]
    contig_pair_plus_rv=contig_pair.contig.tolist()+[str(Seq(i).reverse_complement()) for i in contig_pair.contig.tolist()]
    pairidx=contig_all.contig.isin(contig_pair_plus_rv)#leave too many unpaired contigs to bbduk step
    contig_pair=contig_all[pairidx]
    contig_pair.contig.to_csv('%s/contig_pair/contigs_pairedspkmers'%outdir,header=True,index=False,sep='\t')
    contig_unpair=contig_all[-pairidx]
    print ('orignial contigs: %s;paired contigs: %s; unpaired contigs: %s'%(contig_all.shape[0],contig_pair.shape[0],contig_unpair.shape[0]))

    for name, size in sorted(((name, sys.getsizeof(value)) for name, value in locals().items()),key= lambda x: -x[1])[:10]:
        if size > 100000:print("{:>30}: {:>8}".format(name, sizeof_fmt(size)))


    ################## call variants from paired contigs #######################
    nb_contigs_eachthread=ceil(contig_pair.shape[0]/threads)
    if nb_contigs_eachthread>0 and os.path.exists('%s/variant_result/SNV_results_pair.txt'%outdir) is False:
        print (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),'all threads for %s contig_paired are started'%contig_pair.shape[0])

        subprocess.call("""rm %s/x*;split -l %s -d --additional-suffix=pair_%s %s/contig_pair/contigs_pairedspkmers;mv x*pair_%s %s"""%(outdir,nb_contigs_eachthread,samid,outdir,samid,outdir),shell=True,executable='/bin/bash')
        fileidxs=[i.strip() for i in os.popen('ls %s/x*pair_%s|grep -v unpair'%(outdir,samid)).readlines()]
        print (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),'%s chunks are undergoing analysis'%len(fileidxs))
        pool=Pool(min(len(fileidxs),threads))
        pool.starmap(ana_contigs, zip(fileidxs,[True]*len(fileidxs)))
        pool.close()
        pool.join()
        subprocess.call(r'''cat %s/contig_pair/result* > %s/variant_result/SNV_results_pair.txt;rm %s/contig_pair/result*'''%(outdir,outdir,outdir),shell=True)
        print (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),'all threads for contig_paired are finished')
    ################## call variants from unpaired contigs ########################
    try: nb_unpair=contig_unpair.shape[0]
    except: nb_unpair=0
    if nb_unpair>0 and os.path.exists('%s/variant_result/SNV_results_unpair.txt'%outdir) is False:
        contig_unpair.to_csv('%s/contig_unpair/contigs_unpaired.pass'%outdir,header=False,index=False,sep='\t')
        print (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),'all threads for %s contig_unpaired are started'%(contig_unpair.shape[0]))
        nb_contigs_eachthread=ceil(int(os.popen('wc -l %s/contig_unpair/contigs_unpaired.pass'%outdir).readline().strip().split()[0])/threads)
        subprocess.call("""rm x*passed_%s;split -l %s -d --additional-suffix=passed_%s %s/contig_unpair/contigs_unpaired.pass;mv x*passed_%s %s"""%(samid,nb_contigs_eachthread,samid,outdir,samid,outdir),shell=True,executable='/bin/bash')
        fileidxs=[i.strip() for i in os.popen('ls %s/x*passed_%s'%(outdir,samid)).readlines()]
        print (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),'%s chunks of contig_unpaired undergo filtering'%len(fileidxs))
        pool=Pool(min(len(fileidxs),threads))
        pool.map(filter_unpaired_contigs, fileidxs)
        pool.close()
        pool.join()
        os.system("cat %s/contig_unpair/passedcontig_* > %s/contig_unpair/passedcontig.fa;echo 'contig' > %s/contig_unpair/contigs_unpaired;cat %s/contig_unpair/contigs_unpaired_* >> %s/contig_unpair/contigs_unpaired;cat %s/contig_unpair/FailedToInferRef_* > %s/contig_unpair/FailedToInferRef.txt"%(outdir,outdir,outdir,outdir,outdir,outdir,outdir)) 

    if os.path.exists('%s/variant_result/SNV_results_unpair.txt'%outdir) is False:
        nb_contigs_eachthread=ceil(int(os.popen('wc -l %s/contig_unpair/contigs_unpaired'%outdir).readline().strip().split()[0])/threads)
        print (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),'%s/%s unpaired contigs are dumped to bbduk'%(int(os.popen('wc -l %s/contig_unpair/contigs_unpaired'%outdir).readline().strip().split()[0])-1,contig_unpair.shape[0]))
        #retriev reads from fastq
        os.system("./bbduk.sh in=%s in2=%s ref=%s k=31 mm=f rcomp=t outm=%s fastawrap=500 rename=t hdist=%s speed=0 2>/dev/null"%(wild1,wild2,'%s/contig_unpair/passedcontig.fa'%outdir,'%s/contig_unpair/unpair_contigs_reads.fa'%outdir,distance))
        cor_reads='%s/contig_unpair/unpair_contigs_reads.fa'%outdir
        try:nb_weird_contig=int(os.popen('wc -l %s/contig_unpair/FailedToInferRef.txt'%outdir).readline().strip().split()[0])
        except:nb_weird_contig=0
        subprocess.call(r"""split -l %s -d --additional-suffix=unpair_%s %s/contig_unpair/contigs_unpaired;mv x*unpair_%s %s"""%(nb_contigs_eachthread,samid,outdir,samid,outdir),shell=True,executable='/bin/bash')
        fileidxs=[i.strip() for i in os.popen('ls %s/x*unpair_%s'%(outdir,samid)).readlines()]
        print (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),'%s chunks of contig_unpaired are started'%(len(fileidxs)))
        pool=Pool(min(len(fileidxs),threads))
        pool.starmap(ana_contigs, zip(fileidxs,[False]*len(fileidxs)))
        pool.close()
        pool.join()
        subprocess.call(r'''cat %s/contig_unpair/result* > %s/variant_result/SNV_results_unpair.txt;rm %s/contig_unpair/result* %s/contig_unpair/unpair_contigs_reads.fa'''%(outdir,outdir,outdir,outdir),shell=True)
        print (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),'all threads for contig_unpaired are finished')
    ##################### prepare final results ###################
    print (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),'prepare the final result')
    subprocess.call(r'''cat %s/variant_result/SNV_results_* > %s/variant_result/SNV_results.txt;rm %s/x* %s/contig_unpair/*_tmp*'''%(outdir,outdir,outdir,outdir),shell=True)
    result=pd.read_csv('%s/variant_result/SNV_results.txt'%outdir,header=None,index_col=None,sep='\t')
    result.columns=['contig','putative_ref','sp_case','cov_case','sp_control','cov_control']
    result=result[(result.sp_case>=support)&(result.cov_case>=lowdepth)&(result.cov_control>=lowdepth)&((result.sp_case>10)|(result.sp_case/result.cov_case>cutoff))&(-result.putative_ref.isna())].drop_duplicates()
    result=result.sort_values(['sp_case'],ascending=False)
    result.to_csv('%s/variant_result/SNV_report.txt'%outdir,header=True,index=False,sep='\t')
    with open('%s/variant_result/SNV_report.txt'%outdir)as f:
        snv=list(map(lambda x:x.strip().split(),f))
    RemoveFP_via_control([i[0] for i in snv[1:]])
    align=open('%s/variant_result/SNV_alignments.txt'%outdir,'w')
    for line in tqdm(snv[1:]):
        alignments=pairwise2.align.localms(line[0],line[1],5, -2.5, -8, -.4)
        try:
            #if len(alignments[0][0])<31:continue#the mutant contig is very different to the refseq and alignment failed
            align.write(format_alignment(*alignments[0]).rstrip()+'\t'+str(line[-4])+'\t'+str(line[-3])+'\t'+str(line[-2])+'\t'+str(line[-1])+'\n')
        except:
            print (line,alignments,'align failed')
    align.close()



