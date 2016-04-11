area.BRB <- function(x = NULL, start.date = NULL, end.date = NULL, hab = NULL, iso = 95, plot.it = FALSE, t = "UD", vv = NULL) {

	# Subset
	start.date <- as.POSIXct(start.date)
	end.date <- as.POSIXct(end.date)

	x.sub <- gdltraj(x, min = start.date, max = end.date, type = "POSIXct")

	# Calculate home range
	ud.BRB <- BRB(x.sub, D = vv, Tmax = 90*60,
                Lmin = 5, hmin = 50, type = t, filtershort = FALSE,
                habitat = hab, grid = hab, b = FALSE, extent = 1,
                maxt = 90*60, radius = 150)

	if(plot.it == TRUE){
		image(ud.BRB, col = topo.colors(100))
	}

  po <- list()
  n <- list()

	for(i in 1:length(iso)){

    # polygons
    po[[paste("hr", iso[i], sep = "")]] <- getverticeshr(ud.BRB,
                                                         percent = iso[i])


    proj4string(po[[i]]) <- CRS("+proj=utm +zone=16 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")

    # ndvi

    n[[paste("ndvi", iso[i], sep = "")]] <- mean(ndvi[po[[i]], drop = TRUE]$ndvi,
                                                 na.rm=TRUE)
	}

	return(list(hr = po, ndvi = n, ud = ud.BRB))
}