#!/usr/bin/env python2.7
#!--coding:utf-8 --

"""
2014-11-17: Modified as to remove repeats like (TCGA)n
2015-01-29: Modified as using the GENCODE v21.
2015-02-02: Also caculating GC content. 
2015-02-03: Add phostCons score as conservation score.
"""

import copy,os,random
import numpy as np
import HTSeq
import pandas as pd
from Bio.SeqUtils import GC
from Bio import SeqIO
from joblib import Parallel,delayed
import bx.bbi.bigwig_file


def get_TSS_TES( gtf ):
    """
    """
    regions = { }
    for g in HTSeq.GFF_Reader( gtf ):
        gid = g.attr[ "gene_id" ]
        if gid not in regions:
            regions[ gid ] = g.iv
        else:
            if g.iv.start < regions[ gid ].start:
                regions[ gid ].start = g.iv.start
            if g.iv.end > regions[ gid ].end:
                regions[ gid ].end = g.iv.end
    return regions



def get_gene_model( gtf=None ):
    regions = get_TSS_TES( gtf )
    bgmodel = HTSeq.GenomicArrayOfSets( "auto", stranded = True )     
    for key,iv in regions.items(  ):
        bgmodel[ iv ] += key
    return bgmodel 


def is_contained( iv,bgmodel ):
    gs = set(  )
    for ivb,valueb in bgmodel[iv].steps():
        gs.update( valueb )
    if len( gs ) > 0:
        return 1
    else:
        return 0


def get_distal_assign( iv,bgmodel, proximal=10*1000.0, distal=100*1000.0 ):
    #default value
    ano = "Intergenic"
    #first proximal.upstream
    niv = copy.deepcopy( iv )
    niv.end += proximal
    if is_contained( niv,bgmodel ):
        ano = "Proximal.upstream"
        return ano
    #then proximal.downstream
    niv = copy.deepcopy( iv )
    if niv.start - proximal > 0:
        niv.start -= proximal
    else:
        niv.start = 0
    if is_contained( niv,bgmodel  ):
        ano = "Proximal.downstream"
        return ano
    #then distal.upstream 
    niv = copy.deepcopy( iv )
    niv.end += distal
    if is_contained( niv,bgmodel ):
        ano = "Distal.upstream"
        return ano
    #then distal.downstream 
    niv = copy.deepcopy( iv )
    if niv.start - distal > 0:
        niv.start -= distal
    else:
        niv.start = 0
    if is_contained( niv,bgmodel ):
        ano = "Distal.downstream"
        return ano
    return ano 


def parse_rep_locus( fnIn,fnOut,bgmodel,proximal=10*1000.0, distal=100*1000.0 ):
    header = 1
    with open( fnOut,"w" ) as f:
        for line in open( fnIn ):
            if header == 1:
                header = 0
                continue
            line = line.split( "\n" )[ 0 ].split( "\t" ) 
            if len( line ) < 7:
                continue
            ano = line[ 7  ]
            if ano == "NA":
                continue
            rep = line[ 0 ]
            if ")n" in rep:
                continue
            rep = rep.split( "-" )[ :-1 ]
            rep = "-".join( rep )
            chrom = line[ 1 ]
            start = int( line[ 2 ] )
            end = int( line[ 3 ] )
            strand = line[ 4 ] 
            ano = ano.split( "(" )[ 0 ].strip(  )
            dis = int( line[ 9  ] )
            if ano == "promoter-TSS":
                ano = "Proximal.upstream"
            if ano == "TTS":
                ano = "Proximal.downstream"
            if ano == "Intergenic" :
                iv = HTSeq.GenomicInterval( chrom,start,end,strand  )
                ano = get_distal_assign( iv,bgmodel ) 
            line = [ chrom,start,end,strand,rep,ano ]
            line = "\t".join( map(str,line) ) + "\n"
            f.write( line )


def anotate_rep( fnOut ):
    gtf = "genecode.vM4.gtf"
    bgmodel = get_gene_model( gtf=gtf )
    f = "mm10_repeatmasker_ano.txt"
    parse_rep_locus( f,fnOut,bgmodel )


def getGC( line ):
    line = line.split( "\n" )[ 0 ].split( "\t" ) 
    tmp = "test%s.fa"%random.random(  ) 
    cmd = "twoBitToFa /picb/molsysbio/usr/caoyaqiang/1.Projects/17.Alu/1.GenomeReference/1.hg38/1.hg38_Sequence/hg38.2bit:{chr}:{start}-{end} {tmp}".format( chr=line[ 0 ],start=line[ 1 ], end=line[ 2 ] ,tmp=tmp)
    os.system( cmd )
    for r in SeqIO.parse( open( tmp ),"fasta" ):
        gc = GC( r.seq )
        gc = gc/100
    os.system( "rm %s"%tmp )
    line += [ str( gc ) ]
    line = "\t".join( line ) + "\n"
    return line
        


def addGCinfo( fn ):
    with open( fn+"2","w" ) as f:
        for line in open( fn ):
            nline = getGC( line )
            print nline
            f.write( nline )
    #data = open( fn ).read(  ).split( "\n" )
    #data = Parallel( n_jobs=45 )( delayed( getGC )( line ) for line in data )
    #with open( fn+"2","w" ) as f:
    #    f.write( "".join( data ) )
    

def getConservation( bwh,line ):
    line = line.split( "\n" )[ 0 ].split( "\t" ) 
    d = bwh.get_as_array( line[ 0 ],int( line[ 1 ] ), int( line[ 2 ] ) )
    if d == None or len( d ) == 0:
        d = np.nan
    else:
        d = d[ ~np.isnan( d ) ]
        if len( d ) == 0:
            d = np.nan
        else:
            d = d.mean(  )
    line += [ str( d ) ]
    line = "\t".join( line ) + "\n"
    return line


def addConservationScore( fn ):
    bigWig = "hg38.phyloP7way.bw"
    #bigWig = "hg38.phastCons7way.bw"
    bwh = bx.bbi.bigwig_file.BigWigFile(open(bigWig,"rb"))
    data = open( fn ).read(  ).split( "\n" )
    #data = Parallel( n_jobs=20 )( delayed( getConservation )( bwh,line ) for line in data[ :100 ] )
    nd = [  ]
    for line in data[ :100 ]:
        nd += [ getConservation( bwh,line ) ]
    with open( fn+"3","w" ) as f:
        f.write( "".join( data ) )


def main(  ):
    fn = "mm10_rep_locus.txt"
    anotate_rep( fn )
    #addGCinfo( fn )
    #addConservationScore( fn )



if __name__ == "__main__":
    main(  )
