#!/usr/bin/env python2.7
#--coding:utf-8--


"""
Onfly_mapping_dna_sra.py
"""


__author__="CAO Yaqiang"
__date__="2014-09-19"
__modified__=""
__email__="caoyaqiang0410@gmail.com"



import glob,os,time,commands
from datetime import datetime
from joblib import Parallel,delayed

from Biolibs.rel.General.logger import getlogger
date=time.strftime(' %Y-%m-%d',time.localtime( time.time())) 
logger = getlogger( fn=os.getcwd()+"/"+date.strip()+"_"+os.path.basename(__file__)+".log" )



def prepare_sras(  ):
    sras = glob.glob( "../1.SRAs/*.sra" )
    nsras = {  }
    ns = [  ]
    for s in sras:
        if "_" not in s:
            ns.append( [ s ] )
            continue
        g = s.split( "/" )[ -1 ].split( "_" )[ 0 ]
        if g not in nsras:
            nsras[ g ] = [  ]
        nsras[ g ] += [ s ]
    ns.extend( nsras.values(  ) )
    return ns



def sra_dump( sras ):
    """
    sras: list of sra files
    """
    fn = sras[ 0 ].split( "/" )[ -1 ].split( "." )[ 0 ]
    if "_" in fn:
        fn = fn.split( "_" )[ 0  ]
    fq = fn + ".fastq"
    if os.path.exists( fq ):
        return fq
    for s in sras:
        cmd="fastq-dump --split-3 %s"%s
        print cmd
        os.system( cmd )
    if len( sras ) == 1:
        return fq
    else:
        fqs = glob.glob( "%s*.fastq"%fn )
        cmd1 = "cat %s > %s"%(" ".join( fqs ),fq  )
        cmd2 = "rm %s"%( " ".join( fqs ) )
        cmds = [ cmd1, cmd2 ]
        for c in cmds:
            print c
            os.system( c )
        return fq



def bowtie_mapping( fq ):
    hg38Ind = "/picb/molsysbio/usr/caoyaqiang/1.Projects/17.Alu/1.GenomeReference/1.hg38/3.BowtieTophatIndex/1.Bowtie1/hg38"
    sample = fq.split(  )[ 0 ]
    sample = os.path.splitext( os.path.basename( sample ) )[ 0 ]
    logger.info( "Start mapping %s.\n"%sample )
    os.mkdir( sample )
    sam = sample+"/"+sample+".sam"
    bam = sample+"/"+sample+".bam"
    doBowtie="/picb/molsysbio/usr/caoyaqiang/1.Projects/17.Alu/0.Tools/bowtie-1.1.0/bowtie -p 4 -q --phred33-quals -v 0 -y -k 1 --chunkmbs 200 -S --best --seed 123 %s %s %s"%( hg38Ind,fq,sam )
    samview = "samtools view -S %s -b -o %s"%( sam,bam )
    samsort = "samtools sort %s %s/%s"%( bam,sample,sample )
    samindex = "samtools index %s %s/%s.bai  "%( bam,sample,sample )
    rmsam = "rm %s"%( sam ) 
    cmds = [doBowtie,samview,samsort,samindex,rmsam]
    for c in cmds:
        print c
        status,output=commands.getstatusoutput(c )
        if "bowtie" in c:
            c = "FLAG_A "+c
            logger.info( c )
            #trim with "Warning"
            output = output.split( "\n" )
            output = [t for t in output if not t.startswith( "Warning" ) ]
            output = "\n".join( output )
            logger.info( sample+":\n"+output+"\n" )
    return bam
 

def rm_dup_unmap( bam ):
    fd = os.path.splitext( bam )[ 0 ]
    sam = fd + ".sam"
    bai = fd + ".bai"
    bed = fd + ".bed"
    tdf = fd + ".tdf"
    tmpbam = fd + "2.bam"
    tmpbai = fd + "2.bai"
    #get head
    cmd = "samtools view -H {}".format( bam )
    status,head = commands.getstatusoutput( cmd )
    with open( sam,"w" ) as f:
        f.write( head+"\n" )
    #remove pcr duplicates
    rmdup = "samtools rmdup -s {} {}".format( bam,tmpbam )
    samindex1 = "samtools index {} {}".format( tmpbam,tmpbai )
    #remove unmaped reads
    rmunmaped = "samtools view -F 4 {} >> {}".format( tmpbam,sam )
    #convert sam to bam and index
    samview = "samtools view -S {} -b -o {}".format( sam,bam )
    samindex2 = "samtools index {} {}".format( bam,bai )
    rm = "rm {} {} {}".format( sam,tmpbam,tmpbai )
    #convert bam to bed, keep no-redundance
    bam2bed = "bamToBed -i {} > {}".format( bam,bed )
    cmds = [ rmdup,samindex1,rmunmaped,samview,samindex2,rm,bam2bed ]
    for cmd in cmds:
        print cmd 
        status,output = commands.getstatusoutput( cmd )
    status,output = commands.getstatusoutput( "samtools flagstat %s"%bam )
    logger.info( "%s has been removed duplicates and unmapped reads, basic statistic: \n%s"%(bam,"FLAG_B "+"\n".join( output.split( "\n" )[ :3 ] )) )
    return bed


def trimBed( bed ):
    #trim the 4th column of raw bed file to thin it
    n = bed+".2"
    with open( n,"w" ) as f:
        for line in open( bed ):
            line = line.split( "\n" )[ 0 ].split( "\t" )
            try:
                line[ 3 ] = line[ 5 ]
            except:
                report = "FLAG_C an error line in %s: %s"( bed,"\t".join( line ) )
                continue
            line = "\t".join( line )+"\n"
            f.write( line )
    os.system( "mv %s %s"%( n,bed ) )



def mapping( sras ):
    fn = sras[ 0 ].split( "/" )[ -1 ].split( "." )[ 0 ]
    if "_" in fn:
        fn = fn.split( "_" )[ 0  ]
    bed = fn + ".bed"
    if os.path.exists( bed ) or os.path.exists( bed+".bz2" ):
        print bed,"already generated."
        return 
    fq = sra_dump( sras )
    bam = bowtie_mapping( fq )
    bed = rm_dup_unmap( bam )
    trimBed( bed )
    fd = bam.split( "/" )[ 0 ]
    cmd1 = "mv %s ./"%bed
    cmd2 = "rm -r %s"%fd
    cmd3 = "rm %s"%fq
    cmds = [ cmd1, cmd2,cmd3 ]
    for c in cmds:
        print c
        os.system( c )


def bzip( bed ):
    cmd = "bzip2 -z %s"%bed
    print cmd
    os.system( cmd )


def main(  ):
    #sras = prepare_sras(  )
    #Parallel( n_jobs=4 )( delayed( mapping )( f ) for f in sras )
    beds = glob.glob( "*.bed" )
    Parallel( n_jobs=10 )( delayed( bzip )( f ) for f in beds )


if __name__=='__main__':
    start_time=datetime.now(  )
    main()
    elapsed=datetime.now(  )-start_time
    print "The process is done"
    print "Time used:",elapsed
