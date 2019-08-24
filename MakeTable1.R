get_t1 <- function(){


    ## Merge t1 with t2, t3, t4, CMXMed:
    t1 <- merge(x=t1, y=tb_Stim[,c('t2StimDate','t2StimID','t2CycleID','t2MedProcedure','t2Event','t2Day')], by.x="t1StimID", by.y="t2StimID", all.x=T)
    t1 <- merge(x=t1, y=tb_t3[,c('t3__Event','t3__EventID')], by.x='t2Event', by.y="t3__EventID", all.x=T)
    t1 <- merge(x=t1, y=tb_t4[,c('t4MedProcedure','t4MedProcedureID')], by.x='t2MedProcedure', by.y="t4MedProcedureID", all.x=T)
    t1 <- merge(t1,tb_CMXMed[,c("t2CMXMedStart","t2CMXMedStartDate","MedInstructions__CycleID")],by.x='t2CycleID',by.y='MedInstructions__CycleID',all.x=T)
    t1 <- merge(t1, tb_CycStimMedProcedure[,c('t2CycleID','t2MedProcedure.Gnd','t2MedProcedure.GndDate')], by='t2CycleID',all.x=T)

    ## cleanup the data
    t1$t2Event <- NULL
    t1$t2MedProcedure <- NULL
    t1 <- t1[order(t1$t2CycleID,t1$t2StimDate,t1$t1StimID),]
    follicle.count.fields <-  c('t1RLessThan11','t1R11To13','t1R14To15','t1R16To17','t1R18To19','t1RGreaterThan19','t1LLessThan11','t1L11To13','t1L14To15','t1L16To17','t1L18To19','t1LGreaterThan19')
    for(field in follicle.count.fields) t1[,field] <- as.numeric(t1[,field])
    for(field in c('t3__Event','t4MedProcedure')) t1[,field] <- as.character(t1[,field])
   


    collapseDay3StimsToCycle <- function(dd){
        ## This function will be called for getting baseline AFC+Endo values. 
        ## It serves to:
        ## 1. take data for a given Cycle
        ## 2. subset the data depending according to a hierarchy of which Stim data is most reliable
        ## 3. send that data to collapseStimsToCycle() to return the Cycle's collapsed AFC+Endo values.
        ## On Srg and Srg+1 days, this function is skipped. Instead, data split by CycleID goes to collapseStimsToCycle() directly.

        # Check if Gnd Medprocedure is available first
        if (any(dd$t4MedProcedure %in% 'Gnd',na.rm=T)){
            tt = dd[dd$t4_MedProcedure %in% 'Gnd',]
            lb = as.POSIXct(tt$t2StimDate[1])-ddays(3)
            ub = as.POSIXct(tt$t2StimDate[1])+ddays(1)
            tt = tt[tt$t2StimDate>=lb & tt$t2StimDate<=ub & !is.na(tt$t2StimDate),]

            if (dim(tt)[1]>0){
                h = collapseStimsToCycle(tt, datefield="t2StimDate")

                # If no value within -3+1 days of Gnd, go to CMXMed  
            } else if (any(dd$t2CMXMedStart==1,na.rm=T)){
                tt = dd[dd$t2CMXMedStart==1,]
                lb = as.POSIXct(tt$t2CMXMedStartDate[1])-ddays(3)
                ub = as.POSIXct(tt$t2CMXMedStartDate[1])+ddays(1)
                tt = tt[tt$t2StimDate>=lb & tt$t2StimDate<=ub & !is.na(tt$t2StimDate),]

                if (dim(tt)[1]>0){
                    h = collapseStimsToCycle(tt, datefield="t2CMXMedStartDate")

                    # If no value within -3+1 days of CMXMed, go to StimDay 0-4
                } else if (dim(dd[(dd$t2Day>=0 & dd$t2Day<=4) & !is.na(dd$t2Day),])[1]>0){
                    h = collapseStimsToCycle(dd, datefield="t2Day")

                    # No value found using any logic above  
                } else {
                    h = NULL
                }

                # If Gnd and MedStart are false, check if values in Stim_Day0-4   
            } else if (dim(dd[(dd$t2Day>=0 & dd$t2Day<=4) & !is.na(dd$t2Day),])[1]>0){
                h = collapseStimsToCycle(dd, datefield="t2Day")

                # No value found using any logic above  
            } else {
                h = NULL
            }

            # Otherwise, check for MedStart  
        } else if (any(dd$t2CMXMedStart==1,na.rm=T)){
            tt = dd[dd$t2CMXMedStart==1,]
            lb = as.POSIXct(tt$t2CMXMedStartDate[1])-ddays(3)
            ub = as.POSIXct(tt$t2CMXMedStartDate[1])+ddays(1)
            tt = tt[tt$t2StimDate>=lb & tt$t2StimDate<=ub & !is.na(tt$t2StimDate),]

            if (dim(tt)[1]>0){
                h = collapseStimsToCycle(tt, datefield="t2CMXMedStartDate")

                # If no value found within -3 +1 days of CMXMed, go to StimDay 0-4
            } else if (dim(dd[(dd$t2Day>=0 & dd$t2Day<=4) & !is.na(dd$t2Day),])[1]>0){
                h = collapseStimsToCycle(dd, datefield="t2Day")

                # No value found using any logic above  
            } else {
                h = NULL
            }

            # If Gnd and MedStart are false, check if values in Stim_Day0-4  
        } else if (dim(dd[(dd$t2Day>=0 & dd$t2Day<=4) & !is.na(dd$t2Day),])[1]>0){
            h = collapseStimsToCycle(dd, datefield="t2Day")

            # No value found using any logic above  
        } else {
            h = NULL
        }

        return(h)
    }



    collapseStimsToCycle <- function(x,datefield) {
        ## This function takes data split by CycleID, calculates the Cycle's follicle counts, BAFC_GND, EndoType, EndoThick,
        ## and returns them as a single row of Cycle-level data.
        ## -- "datafield" is the field on which to split the Sum of L and R follicle counts, 
        ## -- before getting overall BAFC_GND for baseline data.


        ## fcn to get the sum of all the R and L counts (on a given day).
      if(any(x$t3__Event %in% c('Srg','NSrg'))){
        SumRLCountsPerDay <- function(afc) return( ifelse(any(!is.na(ifelse(any(!is.na(afc)), colSums(afc,na.rm=T), rep(NA, ncol(afc))))), sum(ifelse(any(!is.na(afc)), colSums(afc,na.rm=T), rep(NA, ncol(afc))),na.rm=T), NA) )
      }                                           
      else{  SumRLCountsPerDay <- function(afc) return( sum(colSums(afc,na.rm=T),na.rm=T))}

        if(nrow(x)>0) {

            ## 1. BAFC_GND -----
            ## For baseline measurements, split by datefield:
            if(length(unique(x[[datefield]])>1) & datefield!="") {
                ## BAFC is a vector with the sum of the R/L for each day
                BAFC <- sapply(split(x[,follicle.count.fields],x[[datefield]]), SumRLCountsPerDay) 
            } else {
                BAFC <- SumRLCountsPerDay(x[,follicle.count.fields])
            }

            ## if any days had non-NA BAFC, BAFC_GND will be the Max, else NA
            BAFC_GND <- ifelse(any(!is.na(BAFC) & length(BAFC>0)), max(BAFC,na.rm=T), NA)                
            BAFC = SumRLCountsPerDay(x[,follicle.count.fields]) 

            ## 2. For remaining fields: FolliclesGtr14mm, EndoThick, EndoType -----
            ## get column sums for the follicle.count.fields (convert to numerics first)  
            t2CycleID <- x$t2CycleID[1]
            x1 <- x[,which(names(x) %in% c('t1RLessThan11', 't1R11To13', 't1R14To15', 't1R16To17', 't1R18To19', 't1RGreaterThan19', 
                                           't1LLessThan11', 't1L11To13', 't1L14To15', 't1L16To17', 't1L18To19', 't1LGreaterThan19'))]
            for(i in 1:length(x1)) x1[,i] <- as.numeric(x1[,i]) 
            x1 <- colSums(x1)            
            n = names(x1)

            ## get overall total Follicle and Follicle Count > 14mm
            RTotal <- c('t1RLessThan11', 't1R11To13','t1R14To15','t1R16To17','t1R18To19','t1RGreaterThan19')
            LTotal <- c('t1LLessThan11', 't1L11To13','t1L14To15','t1L16To17','t1L18To19','t1LGreaterThan19')
            RGreater14 <- c('t1R14To15','t1R16To17','t1R18To19','t1RGreaterThan19')
            LGreater14 <- c('t1L14To15','t1L16To17','t1L18To19','t1LGreaterThan19')
           
            if(any(x$t3__Event %in% c('Srg','NSrg'))){
            RFolliclesTotal <- ifelse(any(!is.na(sum(x1[RTotal],na.rm=T))), sum(x1[RTotal],na.rm=T), NA)
            LFolliclesTotal <- ifelse(any(!is.na(sum(x1[LTotal],na.rm=T))), sum(x1[LTotal],na.rm=T), NA)
            FolliclesTotal <- ifelse(is.na(RFolliclesTotal)&is.na(LFolliclesTotal), NA, RFolliclesTotal+LFolliclesTotal)
              
            RFolliclesGtr14mm <- ifelse(any(!is.na(sum(x1[RGreater14],na.rm=T))), sum(x1[RGreater14],na.rm=T), NA)
            LFolliclesGtr14mm <- ifelse(any(!is.na(sum(x1[LGreater14],na.rm=T))), sum(x1[LGreater14],na.rm=T), NA)
            FolliclesGtr14mm <- ifelse(is.na(RFolliclesGtr14mm)&is.na(LFolliclesGtr14mm), NA, RFolliclesGtr14mm+LFolliclesGtr14mm)}
            
            else{
            RFolliclesTotal <- sum(x1[RTotal],na.rm=T)
            LFolliclesTotal <- sum(x1[LTotal],na.rm=T)
            FolliclesTotal <- RFolliclesTotal+LFolliclesTotal 
              
            RFolliclesGtr14mm <- sum(x1[RGreater14],na.rm=T)
            LFolliclesGtr14mm <- sum(x1[LGreater14],na.rm=T)
            FolliclesGtr14mm <- RFolliclesGtr14mm+LFolliclesGtr14mm }

            ## get EndoType, EndoThick from the latest value for this cycle
            type = x[,which(names(x)=='t1EndoType')]
            thick = x[,which(names(x)=='t1EndoThick')]
            EndoType = ifelse(length(type[type!=0])==0, NA, tail(type[type!=0],1))
            EndoThick = ifelse(length(thick[thick!=0])==0, NA, tail(thick[thick!=0],1))

            ## combine output values into a named-vector and return it.
            out <- c(t2CycleID, x1, RFolliclesTotal, LFolliclesTotal, FolliclesTotal, RFolliclesGtr14mm, LFolliclesGtr14mm, FolliclesGtr14mm, EndoType, EndoThick, BAFC_GND)
            names(out) = c('t2CycleID', n, 't1RFolliclesTotal', 't1LFolliclesTotal', 't1FolliclesTotal','t1RFolliclesGtr14mm', 't1LFolliclesGtr14mm', 't1FolliclesGtr14mm', 't1EndoType', 't1EndoThick', 't1BAFC_GND')
            out <- as.list(out)

        } else {
            out <- NULL
        }
        return(out)
    }


    cleanlist <- function(l) {
        drops <- sapply(l,is.null)
        out <- l[!drops]
        out <- as.data.frame(rbindlist(out))
        return(out)
    }

    sfExport(list=list('collapseStimsToCycle','follicle.count.fields'))
    x.baseline <- t1[order(t1$t2CycleID),]
    if(nrow(x.baseline)>0) {
        x.baseline <- cleanlist(sfLapply(split(x.baseline,x.baseline$t2CycleID), collapseDay3StimsToCycle))
    } 


    ## --------------------
    ## 2. Srg data (data on Surge day)
    x.srg <- t1[t1$t3__Event %in% c('Srg','NSrg'),]
    if(nrow(x.srg)>0) {
        x.srg <- cleanlist(sfLapply(split(x.srg,x.srg$t2CycleID), collapseStimsToCycle, datefield=''))
    } 


    ## --------------------    
    ## 3. Srg+1 data:
    ## 3a. start by removing Cycles that have NO Srg/NSrg, these won't be used here.
    srg.ids <- unique(t1$t2CycleID[t1$t3__Event %in% c('Srg','NSrg')])
    x.srgp1 <- t1[t1$t2CycleID %in% srg.ids,]

    if(nrow(x.srgp1)>0) {
        ## 3b. use the function 'get.srgp1' to extract the Stim-row for Srg+1:
        xx <- sfLapply(split(x.srgp1,x.srgp1$t2CycleID),
                       get.srgp1,cid.field='t2CycleID',srg.field='t3__Event')
        dxx <- as.data.frame(rbindlist(xx))
        x.srgp1 <- dxx[!is.na(dxx$t2CycleID),]

        ## 3c. proceed with the collapseStimsToCycle() function as before:
        x.srgp1 <- cleanlist(sfLapply(split(x.srgp1,x.srgp1$t2CycleID), collapseStimsToCycle, datefield=''))
    } 

    names(x.srg)[2:length(x.srg)] <- paste0(names(x.srg)[2:length(x.srg)],'_srg')
    names(x.srgp1)[2:length(x.srgp1)] <- paste0(names(x.srgp1)[2:length(x.srgp1)],'_srgp1')


    ## -----------------------
    ## 4. merge out 3 subsets:
    out <- merge(x.baseline,x.srg,by='t2CycleID',all=T)
    out <- merge(out,x.srgp1,by='t2CycleID',all=T)   
    to.num <- names(out)
    for(name in to.num) out[,name] <- as.numeric(out[,name]) ## everything here should be numeric


    ## ----------------------
    ## 5. Re-label EndoTypes according to legend provided by Product. 
    recodeEndoType <- function(x) {
        if(length(x)>0) {
            x[!(x %in% c(1:3))] <- NA
            x[x==1] <- "Late proliferative"
            x[x==2] <- "Early secretory"
            x[x==3] <- "Mid-late secretory"
        } else {
            x <- c()   
        }
        return(x)
    }

    out$t1EndoType <- recodeEndoType(out$t1EndoType)
    out$t1EndoType_srg <- recodeEndoType(out$t1EndoType_srg)
    out$t1EndoType_srgp1 <- recodeEndoType(out$t1EndoType_srgp1)

    out$t1FolliclesGtr14mm <- NULL
    out$t1BAFC_GND_srg <- NULL
    out$t1BAFC_GND_srgp1 <- NULL
    ## finished 
    tb_CycStimUSS <<- out
}

