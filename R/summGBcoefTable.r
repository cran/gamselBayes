########## R function: summGBcoefTable ##########

# For obtaining the linear coefficient summary table 
# for summary.gamselBayes().

# Last changed: 08 NOV 2021

summGBcoefTable <- function(namesPredsLin,credLev,meanVec,lowVec,uppVec,addMFVBwarn)
{
   if (is.null(namesPredsLin))
      coefTable <- NULL

   if (!is.null(namesPredsLin))
   {
      # Make a summary table:
   
      coefTable <- data.frame(meanVec,lowVec,uppVec)

      # Enforce alphabetical ordering on the summary table rows:

      alphabetOrder <- order(namesPredsLin)
      namesPredsLin <- namesPredsLin[alphabetOrder]
      coefTable <- coefTable[alphabetOrder,]

      # Padd the row dimension names with some spaces on the left:

      dimNamesOne <- paste("  ",namesPredsLin)

      if (addMFVBwarn)  # Add warning for MFVB accuracy limitations:
      {
         coefTable <- rbind(coefTable,c("","",""),c("WARNING:","method=\"MFVB\" credible intervals",""),
                         c("","have accuracy limitations when ",""),
                         c("","family=\"binomial\".\ \ \ \ \ \ \ \ \ \ \ \ \ ",""))
      
         dimNamesOne <- c(dimNamesOne,c("",">",">>",">>>"))
      } 
      dimnames(coefTable)[[1]] <- dimNamesOne
      credLevName <- paste("  ",as.character(round(100*credLev)),"% credible",sep="")
      dimnames(coefTable)[[2]] <- c("posterior mean",credLevName,"interval")  
   }
 
   # Return the coefficient table:

   return(coefTable)
}

############ summGBcoefTable ############