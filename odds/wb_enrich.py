#!/usr/bin/env python2.7
#--coding:utf-8--

"""
wb_enrich.py
"""


__author__="CAO Yaqiang"
__date__="2014-07-20"
__modified__=""
__email__="caoyaqiang0410@gmail.com"



#systematic library

#3rd library
import numpy as np
from Orange.bio import go
from Orange.bio.utils import stats
import pylab 
import prettyplotlib as ppl 


#own library

#global setting
ROOT_DIR = "/picb/molsysbio/usr/caoyaqiang/1.Projects/15.WD_miRNA/0.GenomeReference/5.WormBase_WS243_Ontology_Phenotype"


def wb_go_enrichment( gids,bg=None ):
    """
    WormBase GO annotation enrichment analysis.
    """
    obo = ROOT_DIR +"/gene_ontology.WS243.obo"
    gaf = ROOT_DIR + "/gene_association.WS243.wb.c_elegans"
    ontology = go.Ontology( obo )
    ontology.set_slims_subset( "goslim_generic" )
    annotations = go.Annotations( gaf,ontology=ontology )
    #res = annotations.get_enriched_terms( gids,reference=bg,slims_only=True,prob=stats.Binomial(), use_fdr=True,)
    res = annotations.get_enriched_terms( gids,slims_only=True,prob=stats.Hypergeometric(), use_fdr=False,)
    rids = [  ]
    for r in res.keys(  ):
        a = ontology[ r ]
        if a.namespace != "biological_process":
            del res[ r ]
    return res,ontology,annotations



def wb_phenotype_enrichment( gids,bg=None ):
    """
    WormBase phenotype annotation enrichment analysis.
    """
    obo = ROOT_DIR +"/phenotype_ontology.WS243.obo"
    gaf = ROOT_DIR + "/phenotype_association.WS243.wb"
    ontology = go.Ontology( obo )
    ontology.set_slims_subset( "phenotype_slim_wb" )
    annotations = go.Annotations( gaf,ontology=ontology )
    #res = annotations.get_enriched_terms( gids,reference=bg,slims_only=True,prob=stats.Binomial(), use_fdr=True,)
    res = annotations.get_enriched_terms( gids,slims_only=True,prob=stats.Hypergeometric(), use_fdr=False,)
    return res,ontology,annotations


 

def wb_other_enrichment(  gids,term="disease", bg=None):
    """
    WormBase disease,development annotation enrichment analysis.
    """
    if term == "disease":
        obo = ROOT_DIR +"/disease_ontology.WS243.obo"
        gaf = ROOT_DIR + "/disease_association.WS243.wb"
    if term == "development":
        obo = ROOT_DIR +"/development_ontology.WS243.obo"
        gaf = ROOT_DIR + "/development_association.WS243.wb"
    else:
        return None
    ontology = go.Ontology( obo )
    annotations = go.Annotations( gaf,ontology=ontology )
    res = annotations.get_enriched_terms( gids,reference=bg,slims_only=False,prob=stats.Binomial(), use_fdr=True,)
    return res,ontology,annotations


def enrich_plot( res,pre="test" ,x="-log2( p )", t="Enriched GO Terms" ):
    terms,ps = [  ],[  ]
    for r in res:
        terms.append("%s(%s)"%( r[ 1 ],r[ 0 ] ) )
        ps.append( r[ -1 ]  )
    fig = pylab.figure( figsize=( 16, int(8/30.0*len(terms))+2 ))
    ax = fig.add_subplot( 1,1,1 )
    pos = np.array(range( len( terms ) )) 
    ppl.barh( pos,ps,yticklabels=terms,ax=ax,grid="x",annotate=False)
    pylab.subplots_adjust( left=0.6,bottom=0.2,top=0.8 )
    ax.set_ylim( [0,len(terms)] )
    pylab.xlabel( x,fontsize=16 )
    pylab.title( t,fontsize=30)
    fig.tight_layout()
    pylab.savefig( pre+".svg",dpi=1000, bbox_inches='tight' )


def enrich( gids,term="go", bg=None,pcut=0.05,gcut=3,tcut=None,pre="test" ):
    if term == "go":
        res,ontology,annotations = wb_go_enrichment( gids,bg=bg)
    else:
        if term == "phenotype":
            res,ontology,annotations = wb_phenotype_enrichment( gids,bg=bg)
        else:
            res,ontology,annotations = wb_other_enrichment( gids,term=term,bg=bg)
    if res == None:
        return 
    results = [  ]
    for r in  res.keys(  ) :
        genes = res[ r ][ 0 ]
        for j in range( len( genes ) ):
            genes[ j ] = annotations.alias_mapper[ genes[ j ] ]
        #Id, Term, genes, all genes in this term, log2(fdr) 
        nr = [ r,ontology[ r ].name,",".join( list( set( genes ) ) ),len( annotations.get_all_genes(r) ), 0.0-np.log2( res[ r ][ 1 ] ) ]
        if res[ r ][ 1 ] > pcut or len( genes ) < gcut :
            continue
        results.append( nr )
    results = sorted( results,key=lambda x: x[ -1 ] )
    if len( results ) > 0:
        if term == "go":
            t = "Enriched GO Terms"
        else:
            t = "Enriched %s Terms"%term
        enrich_plot( results,pre=pre,x="-log2( FDR )", t=t )
    with open( pre+".txt","w" ) as f:
        for r in results:
            f.write( "\t".join( map( str,r ) ) + "\n" )


def test(  ):
    gids = [ "daf-2","daf-16","age-1","sir-2","let-7","mir-71","pax-3","ins-7","bus-4","C18E9.10","cec-3" ]
    enrich( gids,term="go",pre="go_test" )
    enrich( gids,term="phenotype",pre="phenotype_test" )


if __name__ == "__main__":
    test(  )
