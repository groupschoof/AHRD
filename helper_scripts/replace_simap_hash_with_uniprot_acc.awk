# USAGE: In the project's root directory execute
# awk -f helper_scripts/replace_simap_hash_with_uniprot_acc.awk test/resources/simap_iprscn_result.txt > path/2/outputfile.txt
## check the location of the proteina and feature files and specify before running the script 

BEGIN{
  FS="\t"
	while((getline < "./test/resources/simap_hash_to_uniprot_accessions.txt") > 0){ 
		if(parr[$1]){
			parr[$1] = (parr[$1] "," $2)
		} 
		else {
			parr[$1] = $2
      # print parr[$1]
		}
    # print $1
	}
} 
{
  if(parr[$1]) {
    simap_hash = $1 # copy $1, because later sub will
    sub($1,"",$0)   # make retrieval of $1 impossible
    split(parr[simap_hash],p,",")
    for(i in p) {
      print p[i] "\t" $0
    }
  } else {
    print "ERROR! Could not find any Uniprot-Accession for SIMAP-Hash " $1
  }
}
