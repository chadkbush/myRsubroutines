likelihoodRatioTest <- function(x, alpha=0.05, likelihood=TRUE, returnValue=FALSE) 
{
#
#	likelihoodRatioTest()
#	
#	Syntax: 	likelihoodRatioTest(x, alpha=0.05, likelihood=TRUE, returnValue=FALSE)
#	Purpose:	Performs either the Likelihood Ratio Chi-Squared (default) or 
#				the Pearson's Chi-Squared Test on a (j x k) Contingency Table
#	Arguments:
#					x		An object of class(table) of any dimensions.
#				alpha		Numeric. The alpha-level at which to perform the test. 
#							Default is 0.05.	
#		   likelihood		A logical. Default is TRUE. If FALSE, the function
#							performs Pearson's Chi-Squared rather than Likelihood 
#							Ratio Chi-Squared test.
#		  returnValue		A logical. At default, FALSE, the function prints results
#							to the screen but returns no value. If changed to TRUE,
#							the function returns an object of class list() containing
#							the test statistic, degrees of freedom and p-value.
#
#	Author: Chad K. Bush	Email: chadkbush at gmail dot com
#

	z <- as.vector(apply(x,1,sum)); 	# row margins
	y <- as.vector(apply(x,2,sum)); 	# column margins
	n <- sum(apply(x,1:2,sum));		# total observations
	d <- (nrow(x)-1)*(ncol(x)-1);	 	# degrees of freedom
	Obs <- as.vector(x);
	if(any(Obs == 0)) {
		message("\n Hm... Here's your table: \n")
		print(x);
		return(message("\n It looks like at least one cell has a value of 0. \n I can't do my job if you don't do yours... \n"));

	}
	Exp <- vector(mode='numeric', length=length(x));
	i <- 1; j <- 1; k <- 1; 			
	for(i in seq(along=x)) {			# Calculate expected frequencies
		Exp[i] <- ((z[j]*y[k])/n);
		if (j < nrow(x)) {
			j <- j + 1;
		}
		else if (j == nrow(x)) {
			j <- 1;
			k <- k + 1;
		}
		# cat("\n i = ",i," -- Expected = ", Exp[i], "\n");
	}
	if(any(Exp < 5)) {
		message("\n Here are the expected frequencies I've calculated: \n");
		print(matrix(Exp,nrow=2));
		message("\n I should warn you that since at least one cell has \n an expected frequency less than 5, my calculated \n statistics are unlikely to follow the Chi-Squared \n distribution... So you're probably going to have \n to do something else.");
	}
	if(likelihood) {
		L <- vector(mode='numeric', length=length(x));
		L2 <- vector(mode='numeric', length=length(x));
		i <- 1;
		for (i in seq(along=x)){
			L[i] <- Obs[i]*log(Obs[i]/Exp[i]);
			L2 <- 2*sum(L);
			# cat("\n i = ",i,"\t L = ", L[i], "\t L2 = ", L2);
		}
		p <- 1-pchisq(L2,d);
		# # Uncomment this if you like more "standard" 
		# # stars-for-significance notation:
		# if (magnitude > 4) {
		# 	magnitude <- 4;
		# }
		sig <- c('ns');
		if ((p < alpha)) {
			sig <- paste0(rep('*', times=ceiling(abs(log10(p)))), collapse='');
		}
		if(!returnValue) {
			cat("\n\tLikelihood Ratio Chi-Squared Test\n");
			cat("Likelihood Ratio^2 Statistic =", L2, " df = ", d, " p = ", p, " ", sig, "\n\n");
		}
		else if(returnValue) {
			cat("\n\tLikelihood Ratio Chi-Squared Test\n");
			cat("Likelihood Ratio^2 Statistic =", L2, " df = ", d, " p = ", p, " ", sig, "\n\n");
			return(list(L2=L2,df=d,p=p));
		}
	}
	else if (!likelihood) {
 		Chi <- vector(mode='numeric', length=length(x));
 		Chi2 <- vector(mode='numeric', length=1);
		i <- 1;
		for (i in seq(along=x)){
			Chi[i] <- (((Obs[i]-Exp[i])^2)/Exp[i]);
			Chi2 <- sum(Chi);
			# cat("\n i = ",i,"\t Chi = ", Chi[i], "\t Chi2 = ", Chi2);
		}
		p <- 1-pchisq(Chi2,d);
		magnitude <- ceiling(abs(log10(p)));
		# # Uncomment this if you like more "standard" 
		# # stars-for-significance notation:
		# if (magnitude > 4) {
		# 	magnitude <- 4;
		# }
		sig <- c('ns');
		if ((p < alpha)) {
			sig <- paste0(rep('*', times=magnitude), collapse='');
		}
		if (!returnValue) {
			cat("\n\tPearson's Chi-Squared Test\n");
			cat("Chi-Squared Statistic =", Chi2, " df = ", d, " p = ", p, " ", sig,"\n\n");
		}
		else if (returnValue) {
			cat("\n\tPearson's Chi-Squared Test\n");
			cat("Chi-Squared Statistic =", Chi2, " df = ", d, " p = ", p, " ", sig,"\n\n");
			return(list(Chi2=Chi2,df=d,p=p));
		}
	}
}
