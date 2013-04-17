\name{getBUGSDirectory}
\alias{getBUGSDirectory}
\title{Get Location of WinBUGS Program}
\description{
Gets the location of the WinBUGS Program.
The location can be supplied as an argument, set as an option, or set as
an environment variable.
The function will check if the specified location exists and if it
contains the WinBUGS program.
}
\usage{
getBUGSDirectory(dir="c:/Program Files/WinBUGS14", program="WinBUGS14.exe")
}
\arguments{
  \item{dir}{
a character string specifying the location where WinBUGS is installed.
The default value will be used only if the option, \code{bugs.directory},
is not set and if the environment variable, \code{BUGS_HOME}, is not
set.
}
  \item{program}{
a character string specifying the name of the program that is
expected to be found in the WinBUGS installation directory.
}
}
\details{
This getBUGSDirectory allows the WinBUGS to be installed in locations
other than the default (\code{c:/Program Files/WinBUGS14}).
The alternative location can be specified as the option \code{bugs.directory}
with the \code{options()} function,
e.g. \code{options(bugs.directory = "c:/apps/bugs")}.
Or the alternative location can be specified as the environment variable,
\code{BUGS_HOME}.  This can be set before starting Spotfire S+, or set using the
\code{setenv} function.  The order of selection for the location is:
supplied dir argument, as an option, as an environment variable,
and finally, the default value of \code{dir}.

The function will report an error if the determined location does not
exist or if the determined location does not contain a file with the
name \code{program}.
}
\value{
a character string containing the full path to the
WinBUGS installation directory.
}
\author{
TIBCO Software Inc.
}
\examples{
\dontrun{
getBUGSDirectory()  # the default location

setenv(BUGS_HOME="d:/Program Files/WinBUGS14")
getBUGSDirectory()  # should be the above if WinBUGS is really there

options(bugs.directory="c:/apps/WinBUGS")
getBUGSDirectory()  # should be "c:/apps/WinBUGS" if WinBUGS is really there
                    # option overrides the environment variable set above
}
}
