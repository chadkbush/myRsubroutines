expand.dtfm <- function(x, na.strings = "NA", as.is = FALSE, dec = ".")
{
#	expand.dtfm(x, na.strings = "NA", as.is = FALSE, dec = ".")
#
# A function to expand a data.frame produced from a 2xk contigency table
# into a 'flattened' data structure with 'raw' bivariate data. This function
# was originally written by Mark Schwartz and posted in a discussion forum at
# https://stat.ethz.ch/pipermail/r-help/2006-October/115290.html 
  datfrm <- sapply(1:nrow(x), function(i) x[rep(i, each = x$Freq[i]), ], simplify = FALSE);
  datfrm <- subset(do.call("rbind", datfrm), select = -Freq);
  for (i in 1:ncol(datfrm))
  {
    datfrm[[i]] <- type.convert(as.character(datfrm[[i]]), na.strings = na.strings, as.is = as.is, dec = dec);
  }
  datfrm;
}
