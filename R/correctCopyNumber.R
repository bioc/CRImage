correctCopyNumber <-
function(arr="Sample1", chr=NULL, p=NULL, z=NULL, min.value=-5) {
	require(DNAcopy)
require(aCGH)
	
# Modification of mergeLevels, in aCGH package to fix a bug in ansari.test
	MymergeLevels <- 
	function (vecObs, vecPred, pv.thres = 1e-04, ansari.sign = 0.05,
			  thresMin = 0.05, thresMax = 0.5, verbose = 0, scale = TRUE)
	{
		if (thresMin > thresMax) {
			cat("Error, thresMax should be equal to or larger than thresMin\n")
			return()
		}
		thresAbs = thresMin
		sq <- numeric()
		j = 0
		ansari = numeric()
		lv = numeric()
		flag = 0
		if (thresMin == thresMax) {
			flag = 2
		}
		else {
			l.step <- signif((thresMax - thresMin)/10, 1)
			s.step <- signif((thresMax - thresMin)/200, 1)
		}
		while (1) {
			if (verbose >= 1) {
				cat("\nCurrent thresAbs: ", thresAbs, "\n")
			}
			j = j + 1
			sq[j] <- thresAbs
			vecPredNow = vecPred
			mnNow = unique(vecPred)
			mnNow = mnNow[!is.na(mnNow)]
			cont = 0
			while (cont == 0 & length(mnNow) > 1) {
				mnNow = sort(mnNow)
				n <- length(mnNow)
				if (verbose >= 2) {
					cat("\r", n, ":", length(unique(vecPred)), "\t")
				}
				if (scale) {
					d <- (2 * 2^mnNow)[-n] - (2 * 2^mnNow)[-1]
				}
				else {
					d <- (mnNow)[-n] - (mnNow)[-1]
				}
				dst <- cbind(abs(d)[order(abs(d))], (2:n)[order(abs(d))],
							 (1:(n - 1))[order(abs(d))])
				for (i in 1:nrow(dst)) {
					cont = 1
					#Henrik#
					out = aCGH:::combine.func(diff = dst[i, 1], vecObs,
											  vecPredNow, mnNow, mn1 = mnNow[dst[i, 2]],
											  mn2 = mnNow[dst[i, 3]], pv.thres = pv.thres,
											  thresAbs = if (scale) {
											  2 * 2^thresAbs - 2
											  }
											  else {
											  thresAbs
											  })
					####
					#out = aCGH:::combine.func(diff = dst[i, 1], vecObs,
#						  vecPredNow, mnNow, mn1 = mnNow[dst[i, 2]],
#											  mn2 = mnNow[dst[i, 3]], pv.thres = pv.thres,
#											  thresAbs = if (scale) {
#											  2 * 2^thresAbs - 2
#											  }
#											  else {
#											  thresAbs
#											  })
					if (out$pv > pv.thres) {
						cont = 0
						vecPredNow = out$vecPredNow
						mnNow = out$mnNow
						break
					}
				}
			}
			ansari[j] = Myansari.test(sort(vecObs - vecPredNow), sort(vecObs -
																	  vecPred))$p.value
			if (is.na(ansari[j])) {
				ansari[j] = 0
			}
			lv[j] = length(mnNow)
			if (flag == 2) {
				break
			}
			if (ansari[j] < ansari.sign) {
				flag = 1
			}
			if (flag) {
				if (ansari[j] > ansari.sign | thresAbs == thresMin) {
					break
				}
				else {
					thresAbs = signif(thresAbs - s.step, 3)
					if (thresAbs <= thresMin) {
						thresAbs = thresMin
					}
				}
			}
			else {
				thresAbs = thresAbs + l.step
			}
			if (thresAbs >= thresMax) {
				thresAbs = thresMax
				flag = 2
			}
		}
		return(list(vecMerged = vecPredNow, mnNow = mnNow, sq = sq,
					ansari = ansari))
	}
	
	
	Myansari.test <- function(x, ...) UseMethod("Myansari.test")
	
# altered with respect to bug report 2252
# and some additional improvements
	
	Myansari.test.default <-
	function (x, y, alternative = c("two.sided", "less", "greater"), 
			  exact = NULL, conf.int = FALSE, conf.level = 0.95, ...) 
	{
		alternative <- match.arg(alternative)
		if (conf.int) {
			if (!((length(conf.level) == 1) && is.finite(conf.level) && 
				  (conf.level > 0) && (conf.level < 1))) 
            stop("'conf.level' must be a single number between 0 and 1")
		}
		DNAME <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
		x <- x[complete.cases(x)]
		y <- y[complete.cases(y)]
		m <- length(x)
		if (m < 1) 
        stop("not enough 'x' observations")
		n <- length(y)
		if (n < 1) 
        stop("not enough 'y' observations")
		N <- m + n
		r <- rank(c(x, y))
		STATISTIC <- sum(pmin(r, N - r + 1)[seq_along(x)])
		TIES <- (length(r) != length(unique(r)))
		if (is.null(exact)) 
        exact <- ((m < 50) && (n < 50))
		if (exact && !TIES) {
			#modified by Henrik
			R_pansari=NULL
			#	
			pansari <- function(q, m, n) {
				.C(R_pansari, as.integer(length(q)), p = as.double(q), 
				   as.integer(m), as.integer(n))$p
			}
			PVAL <- switch(alternative, two.sided = {
						   if (STATISTIC > ((m + 1)^2%/%4 + ((m * n)%/%2)/2)) 
						   p <- 1 - pansari(STATISTIC - 1, m, n)
						   else p <- pansari(STATISTIC, m, n)
						   min(2 * p, 1)
						   }, less = 1 - pansari(STATISTIC - 1, m, n), greater = pansari(STATISTIC, 
																						 m, n))
			if (conf.int) {
				#modified by Henrik
				R_qansari=NULL
				#
				qansari <- function(p, m, n) {
					.C(R_qansari, as.integer(length(p)), q = as.double(p), 
					   as.integer(m), as.integer(n))$q
				}
				alpha <- 1 - conf.level
				x <- sort(x)
				y <- sort(y)
				ab <- function(sig) {
					rab <- rank(c(x/sig, y))
					sum(pmin(rab, N - rab + 1)[seq_along(x)])
				}
				ratio <- outer(x, y, "/")
				aratio <- ratio[ratio >= 0]
				sigma <- sort(aratio)
				cci <- function(alpha) {
					u <- absigma - qansari(alpha/2, m, n)
					l <- absigma - qansari(1 - alpha/2, m, n)
					uci <- NULL
					lci <- NULL
					if (length(u[u >= 0]) == 0 || length(l[l > 0]) == 
						0) {
						warning("samples differ in location: cannot compute confidence set, returning NA")
						return(c(NA, NA))
					}
					if (is.null(uci)) {
						u[u < 0] <- NA
						uci <- min(sigma[which(u == min(u, na.rm = TRUE))])
					}
					if (is.null(lci)) {
						l[l <= 0] <- NA
						lci <- max(sigma[which(l == min(l, na.rm = TRUE))])
					}
					if (uci > lci) {
						l <- absigma - qansari(alpha/2, m, n)
						u <- absigma - qansari(1 - alpha/2, m, n)
						u[u < 0] <- NA
						uci <- min(sigma[which(u == min(u, na.rm = TRUE))])
						l[l <= 0] <- NA
						lci <- max(sigma[which(l == min(l, na.rm = TRUE))])
					}
					c(uci, lci)
				}
				cint <- if (length(sigma) < 1) {
					warning("cannot compute confidence set, returning NA")
					c(NA, NA)
				}
				else {
					absigma <- sapply(sigma + c(diff(sigma)/2, sigma[length(sigma)] * 
												1.01), ab)
					switch(alternative, two.sided = {
						   cci(alpha)
						   }, greater = {
						   c(cci(alpha * 2)[1], Inf)
						   }, less = {
						   c(0, cci(alpha * 2)[2])
						   })
				}
				attr(cint, "conf.level") <- conf.level
				u <- absigma - qansari(0.5, m, n)
				sgr <- sigma[u <= 0]
				if (length(sgr) == 0) 
                sgr <- NA
				else sgr <- max(sgr)
				sle <- sigma[u > 0]
				if (length(sle) == 0) 
                sle <- NA
				else sle <- min(sle)
				ESTIMATE <- mean(c(sle, sgr))
			}
		}
		else {
			EVEN <- ((N%%2) == 0)
			normalize <- function(s, r, TIES, m = length(x), n = length(y)) {
				m <- as.double(m)
				n <- as.double(n)
				z <- if (EVEN) 
                s - m * (N + 2)/4
				else s - m * (N + 1)^2/(4 * N)
				if (!TIES) {
					SIGMA <- if (EVEN) 
					sqrt((m * n * (N + 2) * (N - 2))/(48 * (N - 
															1)))
					else sqrt((m * n * (N + 1) * (3 + N^2))/(48 * 
															 N^2))
				}
				else {
					r <- rle(sort(pmin(r, N - r + 1)))
					SIGMA <- if (EVEN) 
					sqrt(m * n * (16 * sum(r$lengths * r$values^2) - 
								  N * (N + 2)^2)/(16 * N * (N - 1)))
					else sqrt(m * n * (16 * N * sum(r$lengths * r$values^2) - 
									   (N + 1)^4)/(16 * N^2 * (N - 1)))
				}
				z/SIGMA
			}
			p <- pnorm(normalize(STATISTIC, r, TIES))
			PVAL <- switch(alternative, two.sided = 2 * min(p, 1 - 
															p), less = 1 - p, greater = p)
			if (conf.int && !exact) {
				alpha <- 1 - conf.level
				ab2 <- function(sig, zq) {
					r <- rank(c(x/sig, y))
					s <- sum(pmin(r, N - r + 1)[seq_along(x)])
					TIES <- (length(r) != length(unique(r)))
					normalize(s, r, TIES, length(x), length(y)) - 
					zq
				}
				srangepos <- NULL
				srangeneg <- NULL
				if (length(x[x > 0]) && length(y[y > 0])) 
                srangepos <- c(min(x[x > 0], na.rm = TRUE)/max(y[y > 
															   0], na.rm = TRUE), max(x[x > 0], na.rm = TRUE)/min(y[y > 
																												  0], na.rm = TRUE))
				if (length(x[x <= 0]) && length(y[y < 0])) 
                srangeneg <- c(min(x[x <= 0], na.rm = TRUE)/max(y[y < 
																0], na.rm = TRUE), max(x[x <= 0], na.rm = TRUE)/min(y[y < 
																													0], na.rm = TRUE))
				if (any(is.infinite(c(srangepos, srangeneg)))) {
					warning("cannot compute asymptotic confidence set or estimator")
					conf.int <- FALSE
				}
				else {
					ccia <- function(alpha) {
						statu <- ab2(srange[1], zq = qnorm(alpha/2))
						statl <- ab2(srange[2], zq = qnorm(alpha/2, 
														   lower.tail = FALSE))
						if (statu > 0 || statl < 0) {
							warning("samples differ in location: cannot compute confidence set, returning NA")
							return(c(NA, NA))
						}
						u <- uniroot(ab2, srange, tol = 1e-04, zq = qnorm(alpha/2))$root
						l <- uniroot(ab2, srange, tol = 1e-04, zq = qnorm(alpha/2, 
																		  lower.tail = FALSE))$root
						sort(c(u, l))
					}
					srange <- range(c(srangepos, srangeneg), na.rm = FALSE)
					cint <- switch(alternative, two.sided = {
								   ccia(alpha)
								   }, greater = {
								   c(ccia(alpha * 2)[1], Inf)
								   }, less = {
								   c(0, ccia(alpha * 2)[2])
								   })
					attr(cint, "conf.level") <- conf.level
					statu <- ab2(srange[1], zq = 0)
					statl <- ab2(srange[2], zq = 0)
					if (statu > 0 || statl < 0) {
						ESTIMATE <- NA
						warning("cannot compute estimate, returning NA")
					}
					else ESTIMATE <- uniroot(ab2, srange, tol = 1e-04, 
											 zq = 0)$root
				}
			}
			if (exact && TIES) {
				warning("cannot compute exact p-value with ties")
				if (conf.int) 
                warning("cannot compute exact confidence intervals with ties")
			}
		}
		names(STATISTIC) <- "AB"
		RVAL <- list(statistic = STATISTIC, p.value = PVAL, null.value = c("ratio of scales" = 1), 
					 alternative = alternative, method = "Ansari-Bradley test", 
					 data.name = DNAME)
		if (conf.int) 
        RVAL <- c(RVAL, list(conf.int = cint, estimate = c("ratio of scales" = ESTIMATE)))
		class(RVAL) <- "htest"
		return(RVAL)
	}
	
	
##   function to estimate corrected intensities for each allele
	getCorrected <- function(LRR, BAF, p) {
		x <- 2 ^ (LRR+1)
		y <- tan((pi/2)*BAF)
		Sa <- (x + p - 1 - y*(1-p)) / (p*(y+1))
		Sb <- x/p + 2 - 2/p - Sa
		Tot <- Sa + Sb
		Tot[Tot<0] <- 0
		Sa <- pmax(0, Sa)
		Sb <- pmax(0, Sb)
		old.Sa <- Sa
		old.Sb <- Sb
		Sb <- pmin(old.Sa, old.Sb)
		Sa <- pmax (old.Sa, old.Sb)
		Sa <- Tot - Sb
		data.frame(Sa, Sb, LRR.tum=log2(Sa+Sb) - 1, BAF.tum=2/pi*atan(Sb/Sa))
	}
	if (is.null(z)) {
		stop("z must be a data.frame with columns Name, Chr, Pos, LRR and BAF and with no replicate probes\n")
	}
	
	if (!is.numeric(z[,2])) {
		stop("Chromosome must be numeric\n")
	}
	
	if (!is.null(chr)) {
		sub.z <- z[z$Chr == chr,]
	} else {
		sub.z <- z
	}
##   Segment the LRR
	sub.z <- sub.z[order(sub.z[,2], sub.z[,3]),]
	CN <- CNA(genomdat=as.matrix(sub.z[,4]), chrom=sub.z[,2], maploc=sub.z[,3], data.type="logratio", sampleid=arr)
	cat("\n...Segmenting Sample...\n")
	sm.CN <- smooth.CNA(CN)
	y <- segment(sm.CN)
	y.merged <- MymergeLevels(CN[,3], rep(y$output$seg.mean, y$output$num.mark))
	DNA <- NULL
	for (chr in unique(sub.z$Chr)) {
		rle.merged <- rle(y.merged$vecMerged[y$data$chrom==chr])
		loc.start <- c(1, cumsum(rle.merged$length) + 1)
		loc.start <- loc.start[-length(loc.start)]
		loc.start <- sub.z[sub.z[,2]==chr,][loc.start,3]
		loc.end <- sub.z[sub.z[,2]==chr,][cumsum(rle.merged$length),3]
		num.mark <- rle.merged$lengths
		seg.mean <- rle.merged$values
		ID <- rep(arr, length(seg.mean))
		chrom <- rep(chr, length(seg.mean))
		tmp  <- data.frame(ID, chrom, loc.start, loc.end, num.mark, seg.mean)
		DNA <- rbind(DNA,tmp)
	}
	cat("\n...Estimating BAF...\n")
##   Adjust BAF
	BAF <- ifelse(sub.z[,5] <=0.5, sub.z[,5], 1-sub.z[,5])
	BAF.ID <- ifelse(sub.z[,5] <=0.5, 1, -1)
	count <- 1
	median.BAF <- rep(1, nrow(DNA))
	for (j in 1:nrow(DNA)) {
		ids <- which(sub.z[,3] >= DNA[j,'loc.start'] & sub.z[,3] <= DNA[j, 'loc.end'])
		sub.BAF <- BAF[ids]
		if (quantile(sub.BAF, 0.95, na.rm=TRUE)>0.45) median.BAF[count] <- min(0.5, median(sub.z[ids,5], na.rm=TRUE))
		else if (mean(sub.BAF<0.15, na.rm=TRUE) > 0.85) median.BAF[count] <- median(sub.BAF, na.rm=TRUE)
		else  median.BAF[count] <- median(sub.BAF[sub.BAF > 0.15], na.rm=TRUE)
		count <- count + 1
	}
	DNA$BAF <- median.BAF
	DNA$num.BAF <- apply(DNA, 1, function(x) length(which(sub.z[,3]>=as.numeric(x[3]) & sub.z[,3] <= as.numeric(x[4]))))
	cat("\n...Correcting for cellularity...\n")
	DNA <- cbind(DNA, getCorrected(DNA$seg.mean, DNA$BAF, p))
	DNA$LRR.tum[is.infinite(DNA$LRR.tum)] <- min.value
##   transform back BAF values
	mean.BAF <- rep(DNA$BAF.tum, DNA$num.mark)
	mean.BAF[is.na(mean.BAF)] <- 0.5
	new.BAF <- sub.z[,5]
	new.BAF[new.BAF> 0.15 & new.BAF < 0.85 & mean.BAF < 0.4 & BAF.ID==1] <-
    new.BAF[new.BAF> 0.15 & new.BAF < 0.85 & mean.BAF < 0.4 & BAF.ID==1] +
	rep(DNA$BAF.tum - DNA$BAF, DNA$num.mark)[new.BAF> 0.15 & new.BAF < 0.85 & mean.BAF < 0.4 & BAF.ID==1]
	new.BAF[new.BAF> 0.15 & new.BAF < 0.85 & mean.BAF < 0.4 & BAF.ID==-1] <-
    new.BAF[new.BAF> 0.15 & new.BAF < 0.85 & mean.BAF < 0.4 & BAF.ID==-1] -
	rep(DNA$BAF.tum - DNA$BAF, DNA$num.mark)[new.BAF> 0.15 & new.BAF < 0.85 & mean.BAF < 0.4 & BAF.ID==-1]
	
	new.BAF[new.BAF> 0.15 & new.BAF < 0.85 & mean.BAF >= 0.4] <-
    new.BAF[new.BAF> 0.15 & new.BAF < 0.85 & mean.BAF >= 0.4] +
	sample(c(-1,1), sum(new.BAF> 0.15 & new.BAF < 0.85 & mean.BAF >= 0.4), replace=TRUE) * 
	rep(DNA$BAF.tum - DNA$BAF, DNA$num.mark)[new.BAF> 0.15 & new.BAF < 0.85 & mean.BAF >= 0.4]
    
	new.BAF <- pmax(0, pmin(new.BAF, 1))
	res <- list(x=data.frame(Chrom=sub.z[,2], Pos=sub.z[,3], Orig.LRR=sub.z[,4], Orig.BAF=sub.z[,5], Corr.LRR=sub.z[,4] + rep(DNA$LRR.tum - DNA$seg.mean, DNA$num.mark), Corr.BAF=new.BAF), seg=DNA)
	res
}

