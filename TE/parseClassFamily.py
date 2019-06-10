cs = {  }
fs = {  }
header = 1
for line in open( "UCSC_mm10_RMSK.txt" ):
    if header == 1:
        header = 0
        continue
    line = line.split( "\n" )[ 0 ].split( "\t" ) 
    repName,repClass,repFamily = line[ 4: ]
    if repClass not in cs:
        cs[ repClass ] = set(  )
    cs[ repClass ].add( repFamily )
    if repFamily not in fs:
        fs[ repFamily ] = set(  )
    fs[ repFamily ].add( repName )
with open( "repClass.gmt","w" ) as f:
    for key in cs.keys(  ):
        line = key + "\t" +",".join( cs[ key ] ) + "\n"
        f.write( line )
with open( "repFamily.gmt","w" ) as f:
    for key in fs.keys(  ):
        line = key + "\t" +",".join( fs[ key ] ) + "\n"
        f.write( line )

    
