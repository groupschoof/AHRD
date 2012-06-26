BEGIN{FS="\t"}
{
	while((getline < "proteins") > 0){
		if(parr[$2]){
			parr[$2] = (parr[$2] "," $1)
		} 
		else {
			parr[$2] = $1
		}
	}

	while((getline < "features_all.txt") > 0){
		if(iarr[$1]){
			iarr[$1] = (iarr[$1] "," $12)
		} 
		else {
			iarr[$1] = $12
		}
	}
} 
END {
	for(p in parr){
		# print p " -> " parr[p] " ~> " iarr[parr[p]]
		print p "\t" iarr[parr[p]] > "IDmapping.txt"
	}
}
