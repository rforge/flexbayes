
# Plots the estimated posterior density function for the parameters
densplot.posterior <- function(x, maxVars = 30, maxPerPage = 6, 
  	region = "credible", kernel.width = "nrd", 
  	level = 0.95, lwd = 2, col = 1, 
  	cex.axis = 0.5, cex = 0.7, xRelation = "free",
  	 yRelation = "free", ...) {

	region <- match.arg(region, c("none", "credible", "hpd"))

  if( ( ( maxPerPage %% 2 ) != 0 ) && ( maxPerPage > 1 ) )
    maxPerPage <- maxPerPage - ( maxPerPage %% 2 )

  mfrow.old <- par( "mfrow" )
  par( mfrow = c(1,1) )
  on.exit( par( mfrow = mfrow.old ) )

  if( nvar( x ) > maxVars ){
    x <- x[,1:maxVars]
  }
  
	panel.special <- switch(region,

		"none" = {
			function(x, y, print.p.value, level, ...) {
				panel.xyplot(x, y, ...)
				panel.abline(v = 0, lty = 2, lwd = 2)

				if(print.p.value) {
					pValue <- sum(y[x > 0]) / sum(y)

					if(x[y == max(y)][1] > 0) {
						corner <- c(1, 1)
						pValue <- 1 - pValue
						adj <- 1
					}

					else {
						corner <- c(0, 1)
						adj <- 0
					}

					key(text = paste("Bayes Factor:", round( pValue / ( 1 - pValue ), 4)),
						corner = corner,
						adj = adj,
						transparent = T)
			  }
			}
		},

		"credible" = {
      function( x, y, print.p.value, level, ... ) {
        level <- (1.0 - level) / 2.0
        level <- c(level, 1.0 - level)
        cdf <- cumsum(y)/sum(y)
        idx <- cdf > level[1] & cdf < level[2]

        cr.x <- x[idx]
        cr.x <- c(cr.x[1], cr.x, cr.x[length(cr.x)])
        polygon(cr.x, c(0.0, y[idx], 0.0), col = 3 )

        panel.xyplot(x, y, ...)
        panel.abline(v = 0, lty = 2, lwd = 2 )

				if ( print.p.value ) 
        {
	  			pValue <- sum(y[x > 0])/sum(y)

				  if ( x[ y == max( y ) ][ 1 ] > 0 ) 
          {
	    			corner <- c(1, 1)
	    			pValue <- 1-pValue
	    			adj <- 1
	  			}
	  			else 
          {
	    			corner <- c(0, 1)
	    			adj <- 0
	  			}

				  key( text = paste( "Posterior odds that beta > 0:", round( pValue / ( 1 - pValue ), 4 ) ),
	    	    corner = corner,
	       		adj = adj,
	       		transparent = T )
	  	  }
		  }
		},

    "hpd" = {
	    function( x, y, print.p.value, level, ... ) {
	      y.sorted <- sort(y)
	      idx <- min( which( ( cumsum( y.sorted ) / sum( y ) ) >= (1 - level) ) )

	      hpd.level <- y.sorted[ idx ]
	 	    idx <- seq( 1, length( y ) )[ (y >= hpd.level) ]
	  	  hpd.x <- x[idx]
	      if ( length( hpd.x ) > 0 ){
	        lidx <- length( idx )
	        idx.p1 <- idx[ seq( 2, lidx ) ]
	        idx.m1 <- idx[ seq( 1, lidx - 1 ) ]
        	index.breaks <- seq( 1, lidx - 1 )[ (idx.p1 - idx.m1) > 1 ]
	        if ( length( index.breaks ) == 0 ){
  	   	    hpd.x <- c( hpd.x[1], hpd.x, hpd.x[ length( hpd.x ) ] )
		        hpd.y <- c( 0.0, y[ idx ], 0.0 )
		        #hpd.y[2] <- hpd.y[ length( hpd.y ) - 1 ] <- hpd.level
		        polygon(hpd.x, hpd.y, col = 3 )
		      } else {
	          index.breaks <- c( index.breaks, lidx )
        	  start.idx <- idx[1]
	          for ( i in seq( 1, length( index.breaks ) ) ){
	            idx.loop <- seq( start.idx, idx[ index.breaks[i] ] )
        	    hpd.x <- x[ idx.loop ]
	            hpd.x <- c( hpd.x[1], hpd.x, hpd.x[ length( hpd.x ) ] )
		          hpd.y <- c( 0.0, y[ idx.loop ], 0.0 )

		          polygon( hpd.x, hpd.y, col = 3 )
	            if ( i < length( index.breaks ) ){
	              start.idx <- idx[ index.breaks[i] + 1 ]
        	    }
	          }
	        }#end else index lenght breaks > 0
	      }  
				
		    panel.xyplot(x, y, lwd = 2, ...)
		    panel.abline(v = 0, lty = 2, lwd = 2)

		    if ( print.p.value ) {
		      pValue <- sum(y[x > 0])/sum(y)

		      if(x[y == max(y)][1] > 0) {
		        corner <- c(1, 1)
		        pValue <- 1-pValue
		        adj <- 1
		      } else {
		        corner <- c(0, 1)
		        adj <- 0
		      }

		      key( text = paste( "Bayes Factor:", round( pValue / ( 1 - pValue ), 4 ) ),
			     	  corner = corner, adj = adj, transparent = T )
	      }   
      }
	  }
	)

  #arrange layout
  plots.per.page <- maxPerPage
  pages <- 1
  total.number.plots <- nvar(x)
  if ( total.number.plots > maxPerPage )
  {
    pages <- ceiling( total.number.plots / maxPerPage )
    plots.per.page <- min( total.number.plots, maxPerPage )
  }
  
  split.rows <- ceiling( plots.per.page / 2 )
  if( plots.per.page == 1 )
    split.cols <- 1
  else
    split.cols <- 2
  
  layout = c( split.cols, split.rows, 1 )
  
  myPlots <- vector( "list", 2 )
  last.plot <- 0
  for ( plot.page in seq( 1, pages ) ){
    first.plot <- last.plot + 1
    last.plot <- min( last.plot + plots.per.page, total.number.plots )
    
    df <- NULL
    for( i in first.plot:last.plot ){
      param.density <- density( getSamples(x=x, param=i), width = kernel.width )
    
      xVec <- param.density$x
      yVec <- param.density$y
      names <- rep( varnames(x)[i], length( param.density$x ) )
      temp <- data.frame( x = xVec, y = yVec, names = names )
      df <- rbind(df, temp)
    }
    print( xyplot( y ~ x | names,
            data = df,
            panel = panel.special,
            print.p.value = F,
            level = level,
            as.table = T,
            type = "l",
            col = col,
            lwd = lwd,
#            scales = list(x = "free", y = "free"),
            par.strip.text = list( cex = cex ),
            scales = list( cex.axis = cex.axis, x = xRelation,
              y = yRelation ),
            between = list(x = 1, y = 1),
            strip = function(...) strip.default(..., style = 1),
            xlab = "",
            ylab = "",
            layout = layout ) )
      title( "Parameter Posterior Densities", cex.main = 1, mar = c(5,4,2,2)+.1 )
  }

  return()
}

