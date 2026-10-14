#' Variance ratio estimation for individual-variant analysis with a sparse GRM
#'
#' The \code{Variance_Ratio_Estimation} function estimates the variance ratio used by \code{\link{Individual_Analysis}}
#' with \code{variance_method = "adaptive_ratio"}, for a null model fitted with a sparse genetic relatedness matrix (GRM) in \code{kins},
#' and stores it in the null model object: \code{obj_nullmodel <- Variance_Ratio_Estimation(genofile, obj_nullmodel)},
#' then \code{Individual_Analysis(..., obj_nullmodel = obj_nullmodel, variance_method = "adaptive_ratio")}.
#' For each calibration variant, the ratio is the exact score test variance under the null model,
#' \eqn{g'P g} with \eqn{P = \Sigma^{-1} - \Sigma^{-1}X(X'\Sigma^{-1}X)^{-1}X'\Sigma^{-1}},
#' divided by the score test variance without the GRM, \eqn{V_0 = g'Dg - g'DX(X'DX)^{-1}X'Dg},
#' where \eqn{D} is \eqn{\mu(1-\mu)} for a dichotomous phenotype and \eqn{1/\tau_e} for a quantitative phenotype.
#' The variance ratio is the mean of these ratios; calibration variants with non-positive or non-finite variances are skipped.
#' Calibration variants are drawn from a candidate queue that is fixed before any phenotype-specific check:
#' up to max(1000, \code{variants_num}) variants passing variant-level QC (\code{QC_label} and \code{variant_type}) are
#' selected at random using \code{seed} and put in random order. Candidates are examined in that order, only as needed.
#' Their minor allele count and minor allele frequency (MAF) in the samples of the null model are evaluated lazily,
#' and only candidates with minor allele count of at least \code{mac_cutoff} and MAF of at least 0.01 are used for the
#' calibration: the variance ratio is calibrated on common variants and applied to all variants tested by
#' \code{\link{Individual_Analysis}}. Additional candidates are examined only when needed to obtain the requested
#' number of valid calibration variants.
#' As in SAIGE, the estimation starts with \code{variants_num} variants and adds 10 variants at a time until the
#' coefficient of variation of the ratios, calculated as sd/mean/(number of variants), is at most \code{ratio_cv_cutoff},
#' or no more calibration variants are available.
#' If the variance component estimate of the GRM is 0, \eqn{V_0} is the exact variance and the variance ratio is 1
#' (no calibration variants are used; \code{Individual_Analysis} does not need this calibration in that case).
#' Null models fitted without \code{kins} are not supported.
#' @param genofile an object of opened annotated GDS (aGDS) file of an autosome.
#' @param obj_nullmodel an object from fitting the null model with a sparse GRM in \code{kins}, which is the output from \code{\link{fit_nullmodel}} function.
#' @param mac_cutoff the cutoff of minimum minor allele count of calibration variants (default = 20).
#' Calibration variants also need a minor allele frequency of at least 0.01.
#' @param QC_label channel name of the QC label in the GDS/aGDS file (default = "annotation/filter").
#' @param variant_type type of variant included in the analysis. Choices include "variant", "SNV", or "Indel" (default = "variant").
#' @param geno_missing_imputation method of handling missing genotypes. Either "mean" or "minor" (default = "mean").
#' @param variants_num the initial number of calibration variants (default = 30).
#' @param ratio_cv_cutoff the threshold of the coefficient of variation of the variance ratio estimate (default = 0.001).
#' @param seed a numerical value to set the seed for selecting calibration variants (default = 7590).
#' @return The null model object \code{obj_nullmodel} with the variance ratio calibration added as \code{obj_nullmodel$variance_ratio},
#' a list containing the variance ratio (\code{ratio_v0}), the number of calibration variants (\code{n_markers})
#' and the coefficient of variation of the variance ratio estimate (\code{cv_v0}).
#' It is used internally by \code{\link{Individual_Analysis}} with \code{variance_method = "adaptive_ratio"}.
#' @references Zhou, W., et al. (2018). Efficiently controlling for case-control imbalance and sample relatedness
#' in large-scale genetic association studies. \emph{Nature Genetics}, \emph{50}(9), 1335-1341.
#' (\href{https://doi.org/10.1038/s41588-018-0184-y}{pub})
#' @export

Variance_Ratio_Estimation <- function(genofile,obj_nullmodel,mac_cutoff=20,QC_label="annotation/filter",
                                      variant_type=c("variant","SNV","Indel"),geno_missing_imputation=c("mean","minor"),
                                      variants_num=30,ratio_cv_cutoff=0.001,seed=7590){

	## set to default if none provided
	variant_type <- match.arg(variant_type)
	geno_missing_imputation <- match.arg(geno_missing_imputation)

	## supported null models: single phenotype, one residual variance group, sparse GRM supplied in kins
	## (fit_nullmodel also sets sparse_kins to TRUE when kins is NULL)
	is.supported <- obj_nullmodel$sparse_kins && any(grepl("^kins", names(obj_nullmodel$theta))) &&
	                (obj_nullmodel$n.pheno == 1) && (obj_nullmodel$n.groups == 1)
	if(!is.supported)
	{
		stop("adaptive_ratio requires a single-phenotype null model fitted by fit_nullmodel with a sparse GRM in kins!")
	}

	## the working weights D are proportional to 1/theta[1]: a residual variance estimate of 0 (boundary fit) is not supported
	if(!(is.finite(obj_nullmodel$theta[1]) && (obj_nullmodel$theta[1] > 0)))
	{
		stop("adaptive_ratio requires a positive finite residual variance estimate (theta[1]) in obj_nullmodel, but it is ",obj_nullmodel$theta[1],"!")
	}

	## variance component estimate of 0: Sigma_i is D and the variance without the GRM is exact
	if(all(obj_nullmodel$theta[-1] == 0))
	{
		obj_nullmodel$variance_ratio <- list(ratio_v0=1,n_markers=0,cv_v0=0)
		return(obj_nullmodel)
	}

	## Null Model
	phenotype.id <- as.character(obj_nullmodel$id_include)
	samplesize <- length(phenotype.id)
	residuals.phenotype <- as.vector(obj_nullmodel$scaled.residuals)
	Sigma_i <- obj_nullmodel$Sigma_i
	Sigma_iX <- as.matrix(obj_nullmodel$Sigma_iX)
	cov <- obj_nullmodel$cov

	## working weights without the GRM: mu(1-mu) for binary, 1/tau_e for quantitative phenotypes
	X <- obj_nullmodel$X
	D <- obj_nullmodel$family$mu.eta(obj_nullmodel$linear.predictors)^2/obj_nullmodel$family$variance(obj_nullmodel$fitted.values)/obj_nullmodel$theta[1]
	## projection onto the covariates with weights D using the QR decomposition of D^{1/2}X: Q'D^{1/2} and the identity
	## replace X'D and (X'DX)^{-1}, as X'DX squares the condition number and can be numerically singular for few cases
	X_qr <- qr(sqrt(D)*X)
	XD <- t(qr.Q(X_qr)[,seq_len(X_qr$rank),drop=FALSE]*sqrt(D))
	XDX_inv <- diag(X_qr$rank)

	## get SNV id
	filter <- seqGetData(genofile, QC_label)
	if(variant_type=="variant")
	{
		SNVlist <- filter == "PASS"
	}

	if(variant_type=="SNV")
	{
		SNVlist <- (filter == "PASS") & isSNV(genofile)
	}

	if(variant_type=="Indel")
	{
		SNVlist <- (filter == "PASS") & (!isSNV(genofile))
	}

	variant.id <- seqGetData(genofile, "variant.id")
	SNV.id <- variant.id[SNVlist]

	## randomly selected candidate variants (in GDS order), then put in random order
	set.seed(seed)
	SNV.id <- SNV.id[sort(sample.int(length(SNV.id),min(length(SNV.id),max(1000,variants_num))))]
	SNV.id <- SNV.id[sample.int(length(SNV.id))]

	ratio <- c()
	num_tested <- 0
	num_target <- variants_num
	repeat
	{
		## draw calibration variants until num_target valid ratios are available or no candidates are left
		while((length(ratio) < num_target) && (num_tested < length(SNV.id)))
		{
			is.in <- (num_tested+1):min(num_tested+num_target-length(ratio),length(SNV.id))
			num_tested <- max(is.in)

			## phenotype-specific QC of these candidates in the samples of the null model:
			## minor allele count >= mac_cutoff and minor allele frequency >= 0.01 (calibration on common variants)
			seqSetFilter(genofile, variant.id = SNV.id[is.in], sample.id = phenotype.id)
			variant.id.in <- seqGetData(genofile, "variant.id")
			AF_AC_Missing <- seqGetAF_AC_Missing(genofile,minor=FALSE,parallel=FALSE)
			seqResetFilter(genofile)
			REF_AC <- AF_AC_Missing$ac
			ALT_AC <- 2*round(samplesize*(1-AF_AC_Missing$miss))-REF_AC
			MAC <- ifelse(REF_AC>=ALT_AC,ALT_AC,REF_AC)
			MAF <- pmin(AF_AC_Missing$af,1-AF_AC_Missing$af)
			is.mac <- which(!((MAC<mac_cutoff) | is.na(MAC) | (MAF<0.01) | is.na(MAF)))
			if(length(is.mac) == 0)
			{
				next
			}

			## Genotype imputation and flip
			Genotype_sp <- Genotype_sp_extraction(genofile,variant.id=variant.id.in[is.mac],
			                                      sample.id=phenotype.id,
			                                      REF_AF=AF_AC_Missing$af[is.mac],Missing_rate=AF_AC_Missing$miss[is.mac])
			Geno <- Genotype_sp$Geno
			MAF.in <- Genotype_sp$results_information$MAF

			if (geno_missing_imputation == "mean")
			{
				Geno <- na.replace.sp(Geno,m=2*MAF.in)
			}
			if (geno_missing_imputation == "minor")
			{
				Geno <- na.replace.sp(Geno,is_NA_to_Zero=TRUE)
			}

			## exact variance under the null model and variance without the GRM; invalid variances and ratios are skipped
			V_exact <- as.vector(Individual_Score_Test_sp(Geno, Sigma_i, Sigma_iX, cov, residuals.phenotype)$Score_se)^2
			V_0 <- as.vector(Individual_Score_Test_sp_ratio(Geno, D, XD, XDX_inv, residuals.phenotype, 1)$Score_se)^2
			is.valid <- is.finite(V_exact) & is.finite(V_0) & (V_exact > 0) & (V_0 > 0)
			ratio_in <- V_exact[is.valid]/V_0[is.valid]
			ratio <- c(ratio,ratio_in[is.finite(ratio_in) & (ratio_in > 0)])
		}

		if(length(ratio) < max(variants_num,2))
		{
			stop(paste0("Fewer than ",max(variants_num,2)," valid calibration variants with minor allele count >= ",mac_cutoff," and minor allele frequency >= 0.01!"))
		}

		cv_v0 <- sd(ratio)/mean(ratio)/length(ratio)
		if((cv_v0 <= ratio_cv_cutoff) || (num_tested == length(SNV.id)))
		{
			break
		}
		num_target <- num_target + 10
	}

	obj_nullmodel$variance_ratio <- list(ratio_v0=mean(ratio),n_markers=length(ratio),cv_v0=cv_v0)
	return(obj_nullmodel)
}

