plotImage <- function(img) {
	type=NULL
	if (length(dim(img)) == 2 && is.null(type)) type <- "grey"
	if (length(dim(img)) == 3 && is.null(type)) type <- "rgb"
	ncol=dim(img)[1]
	nrow=dim(img)[2]
	imgdim <- c(ncol, nrow, if (type == "rgb") 3 else NULL)
	img <- array(img, dim=imgdim)
	attr(img, "type") <- type
	colvec <- switch(attr(img, "type"),grey=grey(img),rgb=rgb(img[,,1], img[,,2], img[,,3]))
	colors <- unique(colvec)
	colmat <- array(match(colvec, colors), dim=dim(img)[1:2])
	image(x = 0:(dim(colmat)[2]), y=0:(dim(colmat)[1]),z = t(colmat[nrow(colmat):1, ]), col=colors,xlab="", ylab="", axes=FALSE, asp=1)
}