.onAttach <- function(libname, pkgname) {
    RTver <- read.dcf(file=system.file("DESCRIPTION", package=pkgname),
                      fields="Version")
    cat(paste(pkgname, RTver, "\n"))
}

.onLoad<-function(lib,pkg) require(methods)

