"getBUGSDirectory" <- 
function(dir = "c:/Program Files/WinBUGS14", program = "WinBUGS14.exe")
{
	BUGSDir <- dir
	if(missing(dir)) {
		BUGSDir <- getOption("bugs.directory")
		if(is.null(BUGSDir)) {
			BUGSDir <- getenv("BUGS_HOME")
			if(BUGSDir == "") {
				BUGSDir <- dir
			}
		}
	}
	if(!is.dir(BUGSDir)) {
		stop(paste("The directory:", BUGSDir, "does not exist"))
	}
	if(!file.exists(paste(BUGSDir, program, sep = "/"))) {
		stop(paste("The directory:", BUGSDir, "does not contain", 
			program))
	}
	BUGSDir
}

"getOpenBUGSDirectory" <- 
function(dir = "c:/Program Files/OpenBUGS", program = "winbugs.exe")
{
	BUGSDir <- dir
	if(missing(dir)) {
		BUGSDir <- getOption("bugs.directory")
		if(is.null(BUGSDir)) {
			BUGSDir <- getenv("BUGS_HOME")
			if(BUGSDir == "") {
				BUGSDir <- dir
			}
		}
	}
	if(!is.dir(BUGSDir)) {
		stop(paste("The directory:", BUGSDir, "does not exist"))
	}
	if(!file.exists(paste(BUGSDir, program, sep = "/"))) {
		stop(paste("The directory:", BUGSDir, "does not contain", 
			program))
	}
	BUGSDir
}
