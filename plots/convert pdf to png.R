#-----------------------------------------------------------------------------------------------------------------------------------------------
# can't easily insert pdf images into a google doc, so convert pdfs to png instead
# Requires Imagemagick and ghostscript to be installed locally
#-----------------------------------------------------------------------------------------------------------------------------------------------
pdfs <- list.files('plots',pattern='.pdf',full.names=T)
pngs <- sub('.pdf','.png',pdfs)
for(n in 1:length(pdfs)){
	command <- paste("magick convert -density 300 -flatten",pdfs[n],pngs[n])
	system(command)
	print(pdfs[n])
	}
#-----------------------------------------------------------------------------------------------------------------------------------------------

