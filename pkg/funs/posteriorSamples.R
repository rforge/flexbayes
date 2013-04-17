# Run an MCMC and obtain the posterior samples of the parameters
posteriorSamples <- function (data, inits = NULL, parametersToSave,
		model = "model.txt", nChains = 1, nIter = 10000, 
	  nThin = max(1, floor(nIter / 10000)),
	  nBurnin = 1000 * nThin,
    DIC = TRUE, 
    workingDirectory = NULL, bugsDirectory = getBUGSDirectory(),
    debug = FALSE, engine = "OpenBUGS",
    overRelax = FALSE,
    disabledUpdaters = NULL){

	# Set the current working directory
	if( !is.null(workingDirectory) ) {
		oldwd <- getwd()
		setwd(workingDirectory)
		on.exit(setwd(oldwd))
	}
		
	# If the data has been specified by filename, give an error--it must be specified
	# as an S-PLUS object
	if( is.character( data ) )
		stop ("Data cannot be specified in a file")
		
	if (engine == "WinBUGS"){
 		result <- posteriorSamplesFromWinBUGS (data = data, inits = inits, parametersToSave = parametersToSave,
 			model = model, nChains = nChains, nIter = nIter, nBurnin = nBurnin, nThin = nThin, DIC = DIC,
 			workingDirectory = workingDirectory, bugsDirectory = bugsDirectory, debug = debug) 
 		sims <- result[[1]]
 		DIC.value <- result[[2]]
	} else if (engine == "OpenBUGS"){
    result <- posteriorSamplesFromOpenBUGS (data = data, inits = inits, parametersToSave = parametersToSave,
 			model = model, nChains = nChains, nIter = nIter, nBurnin = nBurnin, nThin = nThin, DIC = DIC,
 			workingDirectory = workingDirectory, overRelax = overRelax, disabledUpdaters = disabledUpdaters) 
 		sims <- result[[1]]
 		DIC.value <- result[[2]]
	} else if (engine == "BatchBUGS"){
    result <- posteriorSamplesFromBatchBUGS (data = data, inits = inits, parametersToSave = parametersToSave,
 			model = model, nChains = nChains, nIter = nIter, nBurnin = nBurnin, nThin = nThin, DIC = DIC,
 			workingDirectory = workingDirectory, bugsDirectory = bugsDirectory, overRelax = overRelax,
 			disabledUpdaters = disabledUpdaters) 
 		sims <- result[[1]]
 		DIC.value <- result[[2]]
	} else {
		stop("There is no Bayesian software engine of that name available in FlexBayes")
	}
	# Put the deviance parameter first
	is.deviance.param <- (varnames(sims) == "deviance")
	param.order <- c(which(is.deviance.param), which(!is.deviance.param))
	sims <- sims [,param.order]

	# If there is only one parameter, change the mcmc objects from vectors to matrices
	if (nvar(sims) == 1){
		mcpar <- attr(sims[[1]], "mcpar")
		for (i in (1:nchain(sims))){
			sims[[i]] <- matrix(as.vector(sims[[i]]), ncol=1)
			dimnames(sims[[i]]) <- list(NULL, c(parametersToSave[1]))
			attr(sims[[i]], "mcpar") <- mcpar
			attr(sims[[i]], "class") <- "mcmc"
		}
	}

	# Create the posterior object.  
	post <- posterior (sims = sims, DIC = DIC.value, call = match.call()) 
	
	return (post)
}

posteriorSamplesFromWinBUGS <- function(data, inits = NULL, parametersToSave,
		model = "model.txt", nChains = 1, nIter = 10000, 
		nBurnin = 1000,
	  nThin = max(1, floor(nIter / 10000)),
    DIC = TRUE, 
    workingDirectory = NULL, bugsDirectory = getBUGSDirectory(),
    debug = FALSE){

	loadR2WinBUGS()

  # Run the MCMC using R2WinBUGS
  nIterWinBUGS <- nIter * nThin + nBurnin  # In WinBUGS the number of iterations 
  																				 # includes the thinning and burnin
  
  # If the model has been specified within S-PLUS as a function, write the model to a 
	# file.
	if (is.character(model)){
		model.file <- model
	} else {
		write.model(model, con = "model.txt")
		model.file <- "model.txt"
	}

  if( is.null(bugsDirectory) ){
  	sims <- bugs (data = data, inits=inits, parameters.to.save = parametersToSave, 
			model.file = model.file, n.chains = nChains, 
			n.iter = nIterWinBUGS, n.burnin = nBurnin, n.thin = nThin, 
			DIC = DIC, codaPkg = T,
			working.directory = workingDirectory, debug = debug)
	} else {
  	sims <- bugs (data = data, inits=inits, parameters.to.save = parametersToSave, 
			model.file = model.file, n.chains = nChains, 
			n.iter = nIterWinBUGS, n.burnin = nBurnin, n.thin = nThin, 
			DIC = DIC, codaPkg = T, bugs.directory = bugsDirectory,
			working.directory = workingDirectory, debug = debug)
  }
	
	# Obtain the parameter posterior samples from the WinBUGS output
	loadCoda()
	sims <- read.bugs (sims, quiet = T)	

	for (i in 1:nChains){  # Fix the start, end, and thin attributes in the mcmc objects.  R2WinBUGS
												 # currently specifies this incorrectly.
	  mcmcAttr <- attr(sims[[i]], "mcpar") 
	  start <- mcmcAttr[1]
	  end <- mcmcAttr[2]
	  end <- (end-start) * nThin + (start-1) * nThin + 1
		start <- (start-1) * nThin + 1
		mcmcAttr[1] <- start
	  mcmcAttr[2] <- end
	  mcmcAttr[3] <- nThin
	  attr(sims[[i]], "mcpar") <- mcmcAttr
	}
	
	DIC.value <- NULL	
	# Obtain the DIC estimate from the WinBUGS output
  if (DIC){
	  DICfromLog <- bugs.log("log.txt")$DIC
	
	  if(any(is.na(DICfromLog))){
		  warning("WinBUGS has not calculated the DIC, either because it may not be appropriate for this model, or because the burnin was too short")
	  } else {
		  DIC.value <- DICfromLog
	  }
  }
  list(sims, DIC.value)
}

posteriorSamplesFromOpenBUGS <- function(data, inits = NULL, parametersToSave,
		model = "model.txt", nChains = 1, nIter = 10000, 
		nBurnin = 1000,
	  nThin = max(1, floor(nIter / 10000)),
    DIC = TRUE, 
    workingDirectory = NULL, overRelax = FALSE,
    disabledUpdaters = NULL){
	
	loadBRugs()
	options(BRugsVerbose = F)

  # If the model has been specified within S-PLUS as a function, write the model to a 
  # file.
  if (is.character(model)){
	  model.file <- model
  } else {
	  writeModel(model, con = "model.txt")
	  model.file <- "model.txt"
	}

	# check the model
	modelCheck("model.txt")
	
	# disable updaters that have been specified
	if (!is.null(disabledUpdaters)){
		if (!is.vector(disabledUpdaters))
			stop('The updaters to disable must be specified as a list of strings')
		nUpdaters <- length(disabledUpdaters)
		for (i in (1:nUpdaters)){
			updater <- disabledUpdaters[i]
			if (!is.character(updater))
				stop('The updaters to disable must be specified as a list of strings')
			modelDisable(updater)
		}
	}
	
	bugsData(data = data, fileName = "data.txt")
	modelData("data.txt")
	modelCompile(numChains = nChains)
	if ( !is.null(inits) ){
    initsFiles <- bugsInits(inits = inits, numChains = nChains)
    modelInits(initsFiles)
  }
	modelGenInits()
	modelUpdate(nBurnin, overRelax = overRelax)
	if (DIC){
    parametersToSave <- c(parametersToSave, "deviance")
		dicSet()
  }
  samplesSet(parametersToSave)
	modelUpdate(numUpdates = nIter, thin = nThin) 

	sims <- buildMCMC(parametersToSave)
	loadCoda()
	for (i in 1:nChains){  # Fix the start, end, and thin attributes in the mcmc objects.  BRugs 
												 # currently specifies this incorrectly.
	  mcmcAttr <- attr(sims[[i]], "mcpar") 
	  start <- mcmcAttr[1]
	  end <- mcmcAttr[2]
	  end <- (end - start) * nThin + start
	  mcmcAttr[2] <- end
	  mcmcAttr[3] <- nThin
	  attr(sims[[i]], "mcpar") <- mcmcAttr
	}
	# Obtain the DIC estimate
	DIC.value <- NULL
  if (DIC){
	  DICval <- dicStats()
	  if( any(is.missing(DICval)) || any(is.na(DICval)) ){
		  warning("OpenBUGS has not calculated the DIC, either because it may not be appropriate for this model, or because the burnin was too short")
	  } else {
		  DIC.value <- DICval
	  }
  }
  list(sims, DIC.value)
}

posteriorSamplesFromBatchBUGS <- function(data, inits = NULL, parametersToSave,
		model = "model.txt", nChains = 1, nIter = 10000, 
		nBurnin = 1000,
	  nThin = max(1, floor(nIter / 10000)),
    DIC = TRUE, 
    workingDirectory = NULL, bugsDirectory = getOpenBUGSDirectory(), overRelax = FALSE,
    disabledUpdaters = NULL){
	
	# Change the working directory to bugsDirectory
	currDir <- getwd()
	setwd( bugsDirectory )
	on.exit( setwd( currDir ) )

  # Put the startup commands in the text for the script file
  scriptTxt <- "modelDisplay(\'log\')\n"
  scriptTxt <- paste( scriptTxt, "modelPrecision(6)\n", sep = "" )
	
  # If the model has been specified within S-PLUS as a function, write the model to a 
  # file.
	loadBRugs()
	options(BRugsVerbose = F)
  if (is.character(model)){
	  model.file <- model
  } else {
	  writeModel( model, con = "model.txt" )
	  model.file <- "model.txt"
	}

	# check the model
	scriptTxt <- paste( scriptTxt, "modelCheck(\"model.txt\")\n", sep = "" )
	
	# disable updaters that have been specified
	if (!is.null(disabledUpdaters)){
		if (!is.vector(disabledUpdaters))
			stop('The updaters to disable must be specified as a list of strings')
		nUpdaters <- length(disabledUpdaters)
		for (i in (1:nUpdaters)){
			updater <- disabledUpdaters[i]
			if (!is.character(updater))
				stop('The updaters to disable must be specified as a list of strings')
			scriptTxt <- paste( scriptTxt, "modelDisable(", updater, ")\n", sep = "" )
		}
	}

	bugsData(data = data, fileName = "data.txt")
	scriptTxt <- paste( scriptTxt, "modelData(\"data.txt\")\n", sep = "" )
	scriptTxt <- paste( scriptTxt, "modelCompile(", nChains, ")\n", sep = "")
	if ( !is.null(inits) ){
    initsFiles <- bugsInits(inits = inits, numChains = nChains)
    for( i in (1:nChains) )
      scriptTxt <- paste( scriptTxt, "modelInits(", initsFiles[i], ")\n", sep = "" )
  }
	scriptTxt <- paste( scriptTxt, "modelGenInits()\n", sep = "" )
	scriptTxt <- paste( scriptTxt, "modelUpdate(", nBurnin, ", ", overRelax, ")\n", sep = "")
	if (DIC){
    parametersToSave <- c(parametersToSave, "deviance")
		scriptTxt <- paste( scriptTxt, "dicSet()\n", sep = "" )
  }
  for( i in (1:length(parametersToSave)) )
    scriptTxt <- paste( scriptTxt, "samplesSet(", parametersToSave[i], ")\n", sep = "" )
	scriptTxt <- paste( scriptTxt, "modelUpdate(", nIter, ", ", nThin, ", ", overRelax, ")\n", sep = "" )
	scriptTxt <- paste( scriptTxt, "samplesCoda(\"*\", \'samps\')\n", sep = "" )
	
	# add in the closing commands to the text for the script file
	scriptTxt <- paste( scriptTxt, "modelExternalize(\'restart\')\n", sep = "" )
  scriptTxt <- paste( scriptTxt, "modelSaveLog(\'log.txt\')\n", sep = "" )
  scriptTxt <- paste( scriptTxt, "modelQuit()\n", sep = "" )

  writeLines(scriptTxt, con = "script.txt")

  bugs.command <- paste( bugsDirectory, "/winbugs.exe", sep = "" )
  if( !file.exists( bugs.command ) )
    stop( paste( "The BUGS executable does not exist at: ", bugs.command ) )
  bugs.command <- paste( bugs.command, "/PAR \"script.txt\" /HEADLESS" )
  system(bugs.command)
    
	# Obtain the parameter posterior samples from the BUGS output
	loadCoda()
	loadR2WinBUGS()
	sims <- list()
	for( i in (1:nChains) ){
	  thisChainSims <- read.coda( paste( "sampsCODAchain", i, ".txt", sep = "" ), 
	                              "sampsCODAindex.txt", quiet = T )
	  sims <- c( sims, list( thisChainSims ) )
	}

  sims <- mcmc.list( sims )
	for (i in 1:nChains){  # Fix the start, end, and thin attributes in the mcmc objects.  BRugs 
												 # currently specifies this incorrectly.
	  mcmcAttr <- attr(sims[[i]], "mcpar") 
	  start <- mcmcAttr[1]
	  end <- mcmcAttr[2]
	  end <- (end - start) * nThin + start
	  mcmcAttr[2] <- end
	  mcmcAttr[3] <- nThin
	  attr(sims[[i]], "mcpar") <- mcmcAttr
	}
	# Obtain the DIC estimate
	DIC.value <- NULL
#  if (DIC){
#	  DICval <- dicStats()
#	  if( any(is.missing(DICval)) || any(is.na(DICval)) ){
#		  warning("OpenBUGS has not calculated the DIC, either because it may not be appropriate for this model, or because the burnin was too short")
#	  } else {
#		  DIC.value <- DICval
#	  }
#  }
  list( sims, DIC.value )
}

