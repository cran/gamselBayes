########## R function: neatCharVecPrint ##########

# For neatly printing a vector of character strings.

# Last changed: 19 JAN 2022

neatCharVecPrint <- function(stringsVec,numLeftPadSpaces,
                             charsPerLine)
{
   iString <- 0
   allFinished <- FALSE
   while (!allFinished)
   {
      currLineFinished <- FALSE
      currLineCharCnt <- 0
      for (iSpace in 1:numLeftPadSpaces) cat(" ")
      while (!currLineFinished)
      {
         iString <- iString + 1
         if (iString>length(stringsVec))
         {
            allFinished <- TRUE
            currLineFinished <- TRUE
         }
         if (iString<=length(stringsVec))     
         {
            strCurr <- stringsVec[iString]
            currLineCharCnt <- currLineCharCnt + nchar(strCurr) + 1
            if (currLineCharCnt<=charsPerLine)
            {
               cat(strCurr) ; cat(" ")
            }
            if (currLineCharCnt>charsPerLine)
            {
               currLineFinished <- TRUE
               currLineCharCnt <- 0
               iString <- iString - 1
               cat("\n")
            }
         }
      }
   }
   cat("\n")

   return(invisible())
}

############ neatCharVecPrint ############