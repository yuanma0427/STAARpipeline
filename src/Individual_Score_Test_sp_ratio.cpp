// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;

// [[Rcpp::export]]
List Individual_Score_Test_sp_ratio(arma::sp_mat G, arma::vec D, arma::mat XD, arma::mat XDX_inv, arma::vec residuals, double ratio)
{
	int i;

	// number of markers
	int p = G.n_cols;

	// Uscore
	arma::rowvec Uscore = trans(residuals)*G;
	// log(p-value)
	arma::vec pvalue_log;
	pvalue_log.zeros(p);

	// SE of Uscore
	arma::vec Uscore_se;
	Uscore_se.zeros(p);
	// Effect size estimation
	arma::vec Est;
	Est.zeros(p);
	// SE of Effect size estimation
	arma::vec Est_se;
	Est_se.zeros(p);

	double test_stat = 0;

	int q = XD.n_rows;

	arma::mat XDG;
	XDG.zeros(q,p);
	XDG = XD*G;

	for(i = 0; i < p; i++)
	{
		// variance without the GRM: g'Dg - (X'Dg)'(X'DX)^{-1}(X'Dg), scaled by the variance ratio
		double Cov_ii = 0;
		for(arma::sp_mat::const_col_iterator it = G.begin_col(i); it != G.end_col(i); ++it)
		{
			Cov_ii = Cov_ii + D(it.row())*(*it)*(*it);
		}

		if (q > 0)
		{
			arma::vec XDG_i = XDG.col(i);
			Cov_ii = Cov_ii - arma::as_scalar(trans(XDG_i)*XDX_inv*XDG_i);
		}

		Cov_ii = ratio*Cov_ii;

		if (Cov_ii == 0)
		{
			pvalue_log(i) = 0;
			Uscore_se(i) = 0;
			Est(i) = 0;
			Est_se(i) = 0;
		}
		else
		{
			test_stat = pow(Uscore(i),2)/Cov_ii;
			pvalue_log(i) = -R::pchisq(test_stat,1,false,true);

			Uscore_se(i) = sqrt(Cov_ii);
			Est(i) = Uscore(i)/Cov_ii;
			Est_se(i) = 1/Uscore_se(i);
		}

	}

	return List::create(Named("Score") = trans(Uscore), Named("Score_se") = Uscore_se, Named("pvalue_log") = pvalue_log, Named("Est") = Est, Named("Est_se") = Est_se);
}

