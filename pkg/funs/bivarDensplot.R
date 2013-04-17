
bivarDensplot <- function(post, maxVars = 5, cuts = 8, xylimProb = 0.98, lwd = 1.5, col = 1, ...)
{

  xylim <- c( (1-xylimProb)/2, 1 - (1-xylimProb)/2 )

	##
	##	The univariate n (UVN) cannot be equal to the
	##	the multivariate n (MVN).
	##

	UVN <- 100
	MVN <- 50

  if( nvar( post ) > maxVars ){
    post <- post[,1:maxVars]
  }
  x <- getSamples(post)
  xNames <-  varnames(post)

	##
	##	kde2d.Mass is a modified version of the function
	##	kde2d found in the MASS library of Venables and
	##	Ripley.
	##

	kde2d.Mass <- function(x, y, h, n = 50, lims = c(range(x),
		range(y)))
	{
		nx <- length(x)
		if(length(y) != nx)
			stop("Data vectors must be the same length")
		gx <- seq(lims[1], lims[2], length = n)
		gy <- seq(lims[3], lims[4], length = n)
		if(missing(h))
			h <- 2*c(bandwidth.nrd(x), bandwidth.nrd(y))
		h <- h/4
		# for S's bandwidth scale
		ax <- outer(gx, x, "-")/h[1]
		ay <- outer(gy, y, "-")/h[2]
		z <- matrix(dnorm(ax), n, nx) %*% t(matrix(dnorm(ay), n, nx)) /
					(nx * h[1] * h[2])
		data.frame(expand.grid(x = gx, y = gy), z = as.vector(z))
	}

	p <- dim(x)[2]
	df <- NULL
	for(j in 1:p) {
		for(i in 1:j) {

			if(i == j) {
				temp <- density( x[,i], n = UVN, from = quantile( x[,i], xylim[1] ), 
				  to = quantile( x[,i], xylim[2] ), width = "nrd" )
				temp <- data.frame( x = temp$x, y = temp$y, z = temp$y, fac = xNames[j])
			}

			else{
				withinLimX <- ( x[,i] > quantile( x[,i], xylim[1] ) ) & 
				  ( x[,i] < quantile( x[,i], xylim[2] ) )
			  withinLimY <- ( x[,j] > quantile( x[,j], xylim[1] ) ) & 
				  ( x[,j] < quantile( x[,j], xylim[2] ) )
				withinLim <- withinLimX & withinLimY
			  restrX <- x[withinLim,i]
        restrY <- x[withinLim,j]
				temp <- data.frame( kde2d.Mass(restrX, restrY, n = MVN), fac = paste(xNames[j], xNames[i], sep = "~") )
		  }

			df <- rbind(df, temp)
		}
	}
	
	skipper <- matrix(F, p, p)
	skipper[row(skipper) > col(skipper)] <- T

	skipper <- as.vector(skipper)

	panel.special <- function(x, y, z, subscripts, at, col.regions,
		contour = F, region = T, n = 10, UVN = 100, ...)
	{

		if(length(x) == UVN) {
			z <- z[subscripts]
			lines(x, z, ... )
		}

		else {
			z <- z[subscripts]
			good <- !(is.na(x) | is.na(y) | is.na(z))
			x <- x[good]
			y <- y[good]
			z <- z[good]
			ux <- sort(unique(x))
			uy <- sort(unique(y))
			Z <- matrix(NA, length(ux), length(uy))
			Z[cbind(match(x, ux), match(y, uy))] <- z
			at <- seq(min(z), max(z), length = n)
			render.contour.trellis(ux, uy, Z, at, labels = F, ...)
		}
	}

	print(contourplot(z ~ x * y | fac,
		data = df,
		scales = list(x = "free", y = "free"),
		as.table = T,
		strip = function(...) strip.default(..., style = 1),
		panel = panel.special,
		col = col,
    lwd = lwd,
		xlab = '',
		ylab = '',
		layout = c(p, p, 1),
		skip = skipper,
		n = (cuts + 2),
		UVN = UVN))

  title( "Bivariate Posterior Densities", cex.main = 1, mar = c(5,4,2,2)+.1 )
	
	invisible(df)
}
