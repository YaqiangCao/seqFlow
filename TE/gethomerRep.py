columns = "id,chr,start,end,strand,name,class/family"
with open( "UCSC_RMSK_HOMER.txt","w" ) as f:
    f.write( "\t".join( columns.split( "," ) )+"\n" )
    for i,line in enumerate( open( "UCSC_RMSK.txt" )):
        if i == 0:
            continue
        line = line.split( "\n" )[ 0 ].split( "\t" )
        nline = [ i ] + line[ :-2 ] + [ line[ -2 ] + "/" + line[ -1 ] ]
        nline = map( str,nline )
        f.write( "\t".join( nline ) + "\n" )
