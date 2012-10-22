library( biomaRt )

# Usage:
print("Usage: Rscript loadDomainAnnotationsFromUniprot.R path/2/accession_per_line.txt path/2/uniprot_domain_annotations_output.tbl")

# Read uniprot accessions
input.args <- commandArgs( trailingOnly = TRUE )
accs <- as.character( read.table( input.args[[1]] )$V1 )

# Download Uniprot domain annotations:
uni.mart <- useDataset( "uniprot", mart=useMart( "unimart" ) )

no.batches <- ceiling( length(accs) / 100 )

annos <- data.frame()

for ( i in 1 : no.batches ) {
  start.ind <- if ( i == 1 ) 1 else ( i - 1 ) * 100
  stop.ind <- i * 100 - 1
  accs.batch <- accs[ start.ind : stop.ind ]
  
  annos <- rbind( annos, 
    getBM( c( "accession", "interpro_id" ), filters=c( "accession" ), values=accs.batch, mart=uni.mart )
  )
}

# Write output:
# print( annos )
write.table( annos, file=input.args[[2]], row.names=F, col.names=F)
