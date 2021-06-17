import pandas as pd
import numpy as np
import scipy.stats as stats
import os,re,sys
from collections import defaultdict
from collections import OrderedDict
from operator import itemgetter
from tempfile import *
from Bio.Seq import Seq
import datetime,time
from multiprocessing import Process,Pool,Manager,Value
SNVFILE=sys.argv[1]
cutoff=float(sys.argv[2])
lowdepth=int(sys.argv[3])
threads=int(sys.argv[4])
genome=sys.argv[5]

MUTREF,INFERED_REF,CONTIG,POS,REF,ALT,FILTER,Case,Control,INFO=[],[],[],[],[],[],[],[],[],[]
def extractrefandalt(refseq,contig):#this part of coding is not elegant.
    REF,ALT,POS_ref,POS_alt='','',[],[]
    try:head=len(re.findall('^-+',contig)[0])
    except:head=0
    try:tail=len(re.findall('-+$',contig)[0])
    except:tail=0
    for i in range(head,len(contig)-tail):
        if contig[i]!=refseq[i]:
            if POS_alt==[]:POS_alt=[i]
            if POS_ref==[]:POS_ref=[i]
            if contig[max(0,i-1)]!=refseq[max(0,i-1)] or REF=='':
                REF+=refseq[i]
                ALT+=contig[i]
            elif contig[max(0,i-1)]==refseq[max(0,i-1)]:
                REF+=','+refseq[i]
                ALT+=','+contig[i]
                POS_alt.append(i)
                POS_ref.append(i)
    if len(POS_alt)==0:return ('NA','NA','NA','NA','NA','NA')
    if REF.count('-')==len(REF):#insertion
        REF=refseq[POS_ref[0]-1]#to make our vcf compatible to snpEFF
        ALT=contig[POS_alt[0]-1]+ALT
    if ALT.count('-')==len(ALT):#deletion
        ALT=contig[POS_alt[0]-1]
        REF=refseq[POS_ref[0]-1]+REF
    if ',' in ALT and ('-' in ALT or '-' in REF):#multiple SNVs including indels
        ALT=ALT.split(',')
        REF=REF.split(',')
        for i in range(len(ALT)):
            if '-' in ALT[i]:#deletion
                ALT[i]=contig[POS_alt[i]-1]
                REF[i]=refseq[POS_ref[i]-1]+REF[i]
            if '-' in REF[i]:#insertion
                ALT[i]=contig[POS_alt[i]-1]+ALT[i]
                REF[i]=refseq[POS_ref[i]-1]
        ALT=','.join(ALT)
        REF=','.join(REF)
    if len(POS_alt)>1:
        POS_alt,POS_ref=np.array(POS_alt),np.array(POS_ref)
        POS_alt_tmp=abs(POS_alt-len(contig.replace('-',''))/2)
        POS_alt=POS_alt[np.where(POS_alt_tmp==np.amin(POS_alt_tmp))]
        POS_ref=POS_ref[np.where(POS_alt_tmp==np.amin(POS_alt_tmp))]
    return (REF,ALT,str(POS_ref[0]),str(POS_alt[0]-head),head,tail)#POS_alt and POS_ref are both index ranging from 0


def seqerror(x):
    if x.ALT*6 in x['#CHROM'][max(0,int(x.POS)-6):min(len(x['#CHROM']),int(x.POS)+7)]:return True
    else:return False

print (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),'prepare INFO')
tmp=NamedTemporaryFile(delete=False)
out=open(tmp.name,'w')
N=0
with open(SNVFILE)as f:
    pos,ref,contig=0,'',''
    for i,j in enumerate(f):
        if (i+1)%4==3:refseq=j.strip().split()[-1]
        if (i+1)%4==1:contig=j.strip().split()[-1]
        if (i+1)%4==0:
            stat=np.array(j.strip().split()[-4:])
            sp,DP,sp_N,DP_N=stat
            N+=1
            out.write('\t'.join([refseq,contig,sp,DP,sp_N,DP_N])+'\n')
out.close()
print ("split -l %s -d --additional-suffix=%s %s;mv x*%s %s"%(int(N/threads),tmp.name.split('/')[-1],tmp.name,tmp.name.split('/')[-1],os.path.dirname(SNVFILE)))
os.system("split -l %s -d --additional-suffix=%s %s;mv x*%s %s"%(int(N/threads),tmp.name.split('/')[-1],tmp.name,tmp.name.split('/')[-1],os.path.dirname(SNVFILE)))
def prepareINFO(fi):
    with open(fi)as f:
        dat=list(map(lambda x:x.strip().split(),f))
    MUTREF,INFERED_REF,CONTIG,POS,REF,ALT,FILTER,Case,Control,INFO=[],[],[],[],[],[],[],[],[],[]
    for line in dat:
        refseq,contig,sp,DP,sp_N,DP_N=line
        sp,DP,sp_N,DP_N=int(sp),int(DP),int(sp_N),int(DP_N)
        ref,alt,pos_r,pos_a,head,tail=extractrefandalt(refseq,contig)
        if pos_r=='NA' or pos_a=='NA' or len(contig)<31:continue
        mutref=refseq.replace(refseq[head:-tail],contig[head:-tail]).replace('-','')
        mutref=mutref[max(0,(int(pos_r)-25)):min(len(mutref),(int(pos_r)+26))]
        MUTREF.append(mutref)
        if sp>100:PhredP=50
        else:PhredP = round(-10*np.log10(stats.fisher_exact([[sp,DP], [sp_N, DP_N]])[1]),2)#abandon using fisher exact test and replace with regression prob
        FILTER.append('PASS')
        CONTIG.append(contig.replace('-',''))
        INFERED_REF.append(refseq.replace('-',''))
        ext=int(abs(int(pos_a)-len(contig.replace('-',''))/2))#min distance to the center
        if len(alt)>=2 and re.sub(',|-','',alt)[::-1]==re.sub(',|-','',ref):INFO.append('type=INV;length=%s;dist2cent=%s;PhredP=%s'%(len(contig.strip('-')),ext,PhredP))
        elif '-' not in alt and ',' in alt:INFO.append('type=MULTIPLE;length=%s;dist2cent=%s;PhredP=%s'%(len(contig.strip('-')),ext,PhredP))
        elif len(alt)<len(ref):INFO.append('type=DEL;length=%s;dist2cent=%s;PhredP=%s'%(len(contig.strip('-')),ext,PhredP))
        elif len(alt)>len(ref):INFO.append('type=INS;length=%s;dist2cent=%s;PhredP=%s'%(len(contig.strip('-')),ext,PhredP))
        elif len(alt)==len(ref):INFO.append('type=SNP;length=%s;dist2cent=%s;PhredP=%s'%(len(contig.strip('-')),ext,PhredP))
        else:INFO.append('type=COMPLEX;length=%s;dist2cent=%s;PhredP=%s'%(len(contig.strip('-')),ext,PhredP))
        Case.append('%s:%s,%s:%s'%(DP,sp,DP-sp,round(sp/DP,3)))
        Control.append('%s:%s,%s:%s'%(DP_N,sp_N,DP_N-sp_N,round(sp_N/DP_N,3)))
        POS.append(pos_a)
        REF.append(ref)
        ALT.append(alt)
    example_dict = {'#CHROM':CONTIG, 'POS':POS, 'REF':REF, 'ALT':ALT,'INFERED_REF':INFERED_REF,'MUTREF':MUTREF}
    df = pd.DataFrame(example_dict)
    df['ID'] = ""
    df['QUAL'] = 40
    df['FILTER']=FILTER
    df['INFO'] = INFO
    df['FORMAT']='DP:AD:VAF'
    df['CASE']=Case
    df['CONTROL']=Control
    df = df[['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'CASE','CONTROL','INFERED_REF']] 
    df.loc[df.INFO.apply(lambda x:'COMPLEX' in x),'FILTER']='FAIL'
    df['MAF_CASE']=df.CASE.apply(lambda x:float(x.split(':')[-1]))
    df['case_cov']=df.CASE.apply(lambda x:int(x.split(':')[0]))
    df.loc[df.apply(seqerror,1),'FILTER']='FAIL'#we take the base quality into account so this function is useless. Which turns out to be not!!! Still enabled
    df.loc[(df.ALT.apply(lambda x:len(x.split(',')))>2)|(df.case_cov<lowdepth)|(df.CONTROL.apply(lambda x:int(x.split(':')[0]))<lowdepth)|((df.INFO.apply(lambda x:float(x.split(';')[3].lstrip('PhredP='))<5)&(df.case_cov*df.MAF_CASE<3))),'FILTER']='FAIL'
    df.to_csv(output_VCF,sep='\t',mode='a',header=False,index=False)

def calpos(x):
    #remind x.POS is the index of alt position
    #remind cigger and x.POS serve for the seqplus instead of the input contig
    softclip_left=re.findall('^(\d+)[S,H]\d+M',x.cigger,re.S)
    if len(softclip_left)==0:sc_l=0
    else:sc_l=int(softclip_left[0])
    indel=re.findall('(\d+)M\d+[I,D](\d+)M$',x.cigger,re.S)
    if len(indel)!=0:
        indel=indel[0]
        if x.CONTIG==x.seqplus:return x.POS_s+int(indel[0])
        else:return x.POS_s+int(indel[1])
    if x.CONTIG==x.seqplus:return x.POS_s+x.POS-sc_l
    else:return x.POS_s+len(x.CONTIG)-x.POS-1-sc_l
        

def align2genome(df):
    tmpfile=NamedTemporaryFile(delete=True)
    out=open('%s.fa'%tmpfile.name,'w')
    for contig in df.CONTIG.tolist():
        out.write(">%s\n"%contig)
        out.write(contig+'\n')
    out.close()
    #Preliminary screening using gsnap
    if os.path.exists(os.path.dirname(genome)+'/gsnap')is False:#create gsnap index
        os.system("gmap_build -D %s %s -d gsnap"%(os.path.dirname(genome),genome))
    print ('gsnap -t 4 -A sam -N 1 -D %s -d gsnap -w 50000 %s.fa 2>/dev/null|samtools view -bh > %s/contigs.bam'%(os.path.dirname(genome),tmpfile.name,os.path.dirname(SNVFILE)))
    os.system('gsnap -t 4 -A sam -N 1 -D %s -d gsnap -w 50000 %s.fa 2>/dev/null|samtools view -bh > %s/contigs.bam'%(os.path.dirname(genome),tmpfile.name,os.path.dirname(SNVFILE)))
    tmpfile=NamedTemporaryFile(delete=True)
    os.system("samtools view %s/contigs.bam|awk '{print $1,$3,$4,$6,$10,$15}' > %s"%(os.path.dirname(SNVFILE),tmpfile.name))
    gsnapres=pd.read_csv(tmpfile.name,header=None,index_col=None,sep=' ')
    gsnapres=gsnapres[gsnapres.iloc[:,1].apply(lambda x:x.startswith('chr'))].groupby(0).first()
    gsnapres.to_csv(tmpfile.name,header=False,index=True,sep='\t')
    anno_gsnap=pd.read_csv(tmpfile.name,header=None,sep='\t')
    anno_gsnap.columns=['CONTIG','#CHROM','POS_s','cigger','seqplus','NM']#mapped contig is rev of the original contig
    mapped_contigs=anno_gsnap.CONTIG.tolist()+[str(Seq(i).reverse_complement()) for i in anno_gsnap.CONTIG.tolist()]
    fusion_gsnap=anno_gsnap[anno_gsnap.cigger.apply(lambda x:np.array(re.findall('(\d+)[S,D,I,N,H]',x,re.S)).astype(int).sum()>5)]
    indel_gsnap=anno_gsnap[anno_gsnap.cigger.apply(lambda x:np.array(re.findall('(\d+)[D,I]',x,re.S)).astype(int).sum()>0)]
    SV_gsnap=fusion_gsnap.append(indel_gsnap).drop_duplicates()
    #SNV_gsnap=anno_gsnap[-anno_gsnap['CONTIG'].isin(SV_gsnap['CONTIG'])]
    SNV_gsnap=anno_gsnap[anno_gsnap.NM!='NM:i:0']
    #SNV=df[df['#CHROM'].isin(SNV_gsnap['CONTIG'])]#only consider convincing SNVs
    SNV=df#consider all SNVs
    SNV=pd.merge(SNV,SNV_gsnap[['#CHROM','CONTIG','cigger','POS_s','seqplus']],left_on='CONTIG',right_on='CONTIG',how='inner')
    SNV.insert(0,'#CHROM',SNV['#CHROM_y'])
    SNV.POS=SNV.apply(calpos,1)
    del SNV['#CHROM_x'], SNV['#CHROM_y'], SNV['seqplus'], SNV['cigger'],SNV['POS_s']
    SV=df[df['#CHROM'].isin(SV_gsnap.CONTIG)]
    print ('potential SV events:',SV.shape)
    #secondary screening using blat
    tmpfile=NamedTemporaryFile(delete=True)
    out=open('%s.fa'%tmpfile.name,'w')
    for contig in SV.CONTIG.tolist():
        out.write(">%s\n"%contig)
        out.write(contig+'\n')
    out.close()
    os.system('pblat -threads=20 %s %s.fa -out=blast8 %s.txt'%(genome,tmpfile.name,tmpfile.name))
    blatres=pd.read_csv('%s.txt'%tmpfile.name,header=None,index_col=None,sep='\t')
    blatres.columns=['CONTIG', 'Subject', 'identity', 'length', 'mismatches', 'gap', 'q. start', 'q. end', 's. start', 's. end', 'e-value', 'bitscore']
    blatres=blatres[(blatres.length>15) & (blatres.Subject.apply(lambda x:x.startswith('chr')))]
    blatres.to_csv('%s/blat_result.txt'%os.path.dirname(SNVFILE),header=True,index=False,sep='\t')
    nomapdat=SV[-SV.CONTIG.isin(blatres.iloc[:,0])]
    #nomapdat.to_csv('%s/unmapped.vcf'%os.path.dirname(SNVFILE),header=True,index=False,sep='\t')
    df[-df['#CHROM'].isin(mapped_contigs)].to_csv('%s/unmapped.vcf'%os.path.dirname(SNVFILE),header=True,index=False,sep='\t')
    TOPTWO=blatres.groupby('CONTIG').head(2).sort_values(['CONTIG','q. start'])
    fusion_contig=[]
    for contig in set(TOPTWO.CONTIG):
        if TOPTWO[TOPTWO.CONTIG==contig].shape[0]!=2:continue
        twoparts=TOPTWO[TOPTWO.CONTIG==contig]
        dist=abs(twoparts['q. start'].tolist()[1]-twoparts['q. end'].tolist()[0])
        if dist<4:fusion_contig.append(contig)
    TOPTWO[TOPTWO.CONTIG.isin(fusion_contig)].to_csv('%s/rearrangement.txt'%os.path.dirname(SNVFILE),header=True,index=False,sep='\t')
    indels=TOPTWO[-TOPTWO.CONTIG.isin(fusion_contig)].groupby('CONTIG').first().reset_index()
    indels.loc[(indels.gap>0)].to_csv('%s/indels.txt'%os.path.dirname(SNVFILE),header=True,index=False,sep='\t')
    return SNV

def removeFP(vcf):
    ext_kmers_nonZero=pd.read_csv('%s/ext_kmers_count'%os.path.dirname(SNVFILE),header=None,sep=' ')
    ext_kmers_nonZero=ext_kmers_nonZero[ext_kmers_nonZero[1]>0]
    ext_kmers_nonZero=ext_kmers_nonZero[0].tolist()+ext_kmers_nonZero[0].apply(lambda x:str(Seq(x).reverse_complement())).tolist()
    FP=[]
    cs_count_N=vcf.CONTROL.apply(lambda x:int(x.split(':')[1].split(',')[0]))
    vcf.loc[cs_count_N>0,'FILTER']='FAIL'
    for i in range(vcf.shape[0]):
        contig,pos=vcf.iloc[i,0],vcf.iloc[i,1]
        if len(contig)>60:continue 
        if len(contig)-pos<30:
            for c in 'ATCG':
                if contig[-30:]+c in ext_kmers_nonZero:
                    FP.append(contig)
                    break
        elif pos<30:
            for c in 'ATCG':
                if c+contig[:30] in ext_kmers_nonZero:
                    FP.append(contig)
                    break
    vcf.loc[vcf['#CHROM'].isin(FP),'FILTER']='FAIL'
    #os.system('rm %s/ext_kmers_count'%os.path.dirname(SNVFILE))
    return vcf
header = """##fileformat=VCFv4.1
##filedate=%s
##source=mutcaller_V2
##SAMPLE=%s
##FILTER=<ID=MULTIPLE,Description="Mapping type : .">
##INFO=<ID=type,Number=1,Type=String,Description="SNP, INS, DEL, INV, MULTIPLE or COMPLEX.">
##INFO=<ID=length,Number=1,Type=Integer,Description="length of the contig">
##INFO=<ID=dist2cent,Number=1,Type=Integer,Description="distance between variant pos and the center of the contig">
##INFO=<ID=PhredP,Number=1,Type=Float,Description="Fisher exact test PhredP">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Estimated depth in each library">
##FORMAT=<ID=AD,Number=2,Type=Integer,Description="Support reads of each allele by library">
"""%(datetime.date.today(),SNVFILE.split('/')[-2])
if __name__ == '__main__':
    output_VCF = SNVFILE.replace('.txt','.vcf')
    with open(output_VCF, 'w') as vcf:
        vcf.write(header)
    fileidxs=[i.strip() for i in os.popen('ls %s/x*%s'%(os.path.dirname(SNVFILE),tmp.name.split('/')[-1])).readlines()]
    pool=Pool(min(len(fileidxs),threads))
    pool.map(prepareINFO, fileidxs)
    pool.close()
    pool.join()
    print (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),'INFO finished')
    df=pd.read_csv(output_VCF,header=None,index_col=None,sep='\t',skiprows=range(11))
    df.columns=['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','CASE','CONTROL','INFERED_REF','MAF_CASE','case_cov']
    sp=df.CASE.apply(lambda x:int(x.split(':')[1].split(',')[0]))
    df['sp']=sp
    df=df.sort_values('sp',ascending=False)
    del df['MAF_CASE'], df['case_cov'], df['sp']
    df=removeFP(df)
    with open(output_VCF,'w')as vcf:
        vcf.write(header)
    df[df.FILTER=='PASS'].to_csv(output_VCF,mode='a',header=True,index=False,sep='\t')
    print ('%d of %d records are predicted to be false positive'%(df[df.FILTER=='FAIL'].shape[0],df.shape[0]))
    df[df.FILTER=='FAIL'].to_csv('%s/Not_convincing_contigs.vcf'%os.path.dirname(SNVFILE),sep='\t',header=True,index=False)


if sys.argv[5]!='None':
    genome=sys.argv[5]
    fa=open(os.path.dirname(SNVFILE)+'/contigs.fa','w')
    for key in df['#CHROM'].tolist():
        fa.write('>%s\n'%key)
        fa.write('%s\n'%key)
    fa.close()
    output_VCF=SNVFILE.replace('.txt','_ref.vcf')
    with open(output_VCF, 'w') as vcf:
        vcf.write(header)
    df=pd.read_csv(SNVFILE.replace('.txt','.vcf'),header=0,sep='\t',skiprows=list(range(11)))
    df['CONTIG']=df['#CHROM']
    print (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),'mapping to the genome')
    map2ref=align2genome(df)
    print (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),'mapping finished')
    combineddf=map2ref
    sp=combineddf.CASE.apply(lambda x:int(x.split(':')[1].split(',')[0]))
    combineddf['sp']=sp
    combineddf=combineddf.sort_values('sp',ascending=False)
    del combineddf['sp']
    #combineddf=combineddf.sort_values(['#CHROM','POS'],ascending=True)
    combineddf[combineddf.FILTER=='PASS'].to_csv(output_VCF,sep='\t',mode='a',index=False)
    combineddf.loc[combineddf.FILTER=='FAIL',].to_csv('%s/Not_convincing_contigs.vcf'%os.path.dirname(SNVFILE),sep='\t',header=True,index=False)
os.system('rm %s/contigs.fa %s/blat_result.txt %s/SNV_results.txt %s/x*'%(os.path.dirname(SNVFILE),os.path.dirname(SNVFILE),os.path.dirname(SNVFILE),os.path.dirname(SNVFILE)))
