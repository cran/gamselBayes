.onAttach <- function(libname, pkgname)
   packageStartupMessage("gamselBayes 2.0 loaded.\nCopyright V.X. He and M.P. Wand 2023.\nFor details on the use of gamselBayes, issue the command:\ngamselBayesVignette()")

.onUnload <- function(libpath)
    library.dynam.unload("gamselBayes",libpath)
