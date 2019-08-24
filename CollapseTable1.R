collapse_table1 <- function(){

	message('collapsing unit-level t1 into one level up t2.’)

	t1$t1Technique[t1$t1Technique=="PR-OH / sucrose"] <- "PROH / sucrose"
	temp_data= t1[, c("t1ID", "t1Outcome", "t1Grade", "t1Expansion", "t1ICM", "t1Troph", "t1CellStage", "t1Cell", "t1Technique")]
	temp_data <- temp_data[order(temp_data$t1ID),]

	#get unique values for each field of interest
	fields <- list()
	fields[[1]] <- c("Replaced", "Ongoing", "Cryoed", "Discard", "MFR", "MFC")
	fields[[2]] <- c("1", "2", "3", "4") 
	fields[[3]] <- c("0", "1", "2", "3", "4", "5", "6")
	fields[[4]] <- c("A", "B", "C", "D", "E", "F", "None")
	fields[[5]] <- c("A", "B", "C", "D", "E")
	fields[[6]] <- c("Day 1", "Day 2", "Day 3", "Day 4", "Day 5", "Day 6", "Day 7", "Day 8", "oocyte")
	fields[[7]] <- c("0_1", "2_3", "4_7", "8_16", "Morula", "Blast")
	fields[[8]] <- c("Vitrification", "PROH / sucrose", "glycerol / sucrose", "Slow Freeze")

	names(fields) <- c("t1Outcome", "t1EmbryoGrade", "t1Expansion", "t1ICM", "t1Troph", "t1CellStage", "t1Cell", "t1Technique")

	# Create new field names based on unique values
	field.names <- list()
	for (i in 1:length(colnames(temp_data[,!colnames(temp_data)%in%c("t1ID")]))) {
		field.names[[i]] <- paste0(colnames(temp_data[,!colnames(temp_data)%in%c("t1ID")])[i],na.omit(gsub(' ','\\', fields[[i]][fields[[i]]!=" "]))) 
		}
	field.names <- unlist(field.names)


	fun <- function(temp_data, fields){
		sapply(1:length(temp_data[,!colnames(temp_data)%in%c("t1ID")]), function(x){	   
			sapply(1:length(fields[[x]]), function(y) {	
				if(grepl("-", fields[[x]][[y]])){
					as.integer(temp_data[,names(fields)[x]] %in% eval(parse(text=(c(sub("-", ":", fields[[x]][[y]]))))) ) } 
				else{  as.integer(temp_data[,names(fields)[x]] %in% c(unlist(strsplit(fields[[x]][[y]], split="_"))) )  }		
				}	)	       
			}   )  
		}
	
	out <- sfLapply(split(temp_data,temp_data$t1ID), fun, fields)

	
	out2 <- t(sapply(1:length(out), function(x){
				unlist(	lapply(1:length(out[[x]]), function(y){
				if(is.matrix(out[[x]][[y]])){
				apply(out[[x]][[y]], 2, sum) }	
				else(out[[x]][[y]])	
				}  )   )  
			} )
		) 

	out2 <- as.data.frame(out2)
	out2$t1ID <- names(out)
	colnames(out2) <- c(field.names, "t1ID")

	t1 <<- out2

	rm(temp_data, out, out2)
         
}