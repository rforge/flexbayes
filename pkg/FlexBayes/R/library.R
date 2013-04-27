.onLoad <- function(libname, pkgname)
{
  dll <- library.dynam(pkgname, package = pkgname, lib.loc = libname)
  invisible()
}


