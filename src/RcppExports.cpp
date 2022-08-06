// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// computeECC
NumericVector computeECC(NumericMatrix D, int maxhomDim, NumericVector scaleSeq);
RcppExport SEXP _TDAvec_computeECC(SEXP DSEXP, SEXP maxhomDimSEXP, SEXP scaleSeqSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type D(DSEXP);
    Rcpp::traits::input_parameter< int >::type maxhomDim(maxhomDimSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type scaleSeq(scaleSeqSEXP);
    rcpp_result_gen = Rcpp::wrap(computeECC(D, maxhomDim, scaleSeq));
    return rcpp_result_gen;
END_RCPP
}
// computePES
NumericVector computePES(NumericMatrix D, int homDim, NumericVector scaleSeq);
RcppExport SEXP _TDAvec_computePES(SEXP DSEXP, SEXP homDimSEXP, SEXP scaleSeqSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type D(DSEXP);
    Rcpp::traits::input_parameter< int >::type homDim(homDimSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type scaleSeq(scaleSeqSEXP);
    rcpp_result_gen = Rcpp::wrap(computePES(D, homDim, scaleSeq));
    return rcpp_result_gen;
END_RCPP
}
// computePI
NumericVector computePI(NumericMatrix D, int homDim, int res, double sigma, double minB, double maxB, double minP, double maxP);
RcppExport SEXP _TDAvec_computePI(SEXP DSEXP, SEXP homDimSEXP, SEXP resSEXP, SEXP sigmaSEXP, SEXP minBSEXP, SEXP maxBSEXP, SEXP minPSEXP, SEXP maxPSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type D(DSEXP);
    Rcpp::traits::input_parameter< int >::type homDim(homDimSEXP);
    Rcpp::traits::input_parameter< int >::type res(resSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type minB(minBSEXP);
    Rcpp::traits::input_parameter< double >::type maxB(maxBSEXP);
    Rcpp::traits::input_parameter< double >::type minP(minPSEXP);
    Rcpp::traits::input_parameter< double >::type maxP(maxPSEXP);
    rcpp_result_gen = Rcpp::wrap(computePI(D, homDim, res, sigma, minB, maxB, minP, maxP));
    return rcpp_result_gen;
END_RCPP
}
// computePL
NumericVector computePL(NumericMatrix D, int homDim, NumericVector scaleSeq,int k);
RcppExport SEXP _TDAvec_computePL(SEXP DSEXP, SEXP homDimSEXP, SEXP scaleSeqSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type D(DSEXP);
    Rcpp::traits::input_parameter< int >::type homDim(homDimSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type scaleSeq(scaleSeqSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(computePL(D, homDim, scaleSeq,k));
    return rcpp_result_gen;
END_RCPP
}
// computePS
NumericVector computePS(NumericMatrix D, int homDim, NumericVector scaleSeq,int p);
RcppExport SEXP _TDAvec_computePS(SEXP DSEXP, SEXP homDimSEXP, SEXP scaleSeqSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type D(DSEXP);
    Rcpp::traits::input_parameter< int >::type homDim(homDimSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type scaleSeq(scaleSeqSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(computePS(D, homDim, scaleSeq, p));
    return rcpp_result_gen;
END_RCPP
}
// computeVAB
NumericVector computeVAB(NumericMatrix D, int homDim, NumericVector scaleSeq);
RcppExport SEXP _TDAvec_computeVAB(SEXP DSEXP, SEXP homDimSEXP, SEXP scaleSeqSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type D(DSEXP);
    Rcpp::traits::input_parameter< int >::type homDim(homDimSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type scaleSeq(scaleSeqSEXP);
    rcpp_result_gen = Rcpp::wrap(computeVAB(D, homDim, scaleSeq));
    return rcpp_result_gen;
END_RCPP
}
// computeVPB
NumericVector computeVPB(NumericMatrix D, int homDim, NumericVector xSeq, NumericVector ySeq, double tau);
RcppExport SEXP _TDAvec_computeVPB(SEXP DSEXP, SEXP homDimSEXP, SEXP xSeqSEXP, SEXP ySeqSEXP, SEXP tauSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type D(DSEXP);
    Rcpp::traits::input_parameter< int >::type homDim(homDimSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type xSeq(xSeqSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ySeq(ySeqSEXP);
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    rcpp_result_gen = Rcpp::wrap(computeVPB(D, homDim, xSeq, ySeq, tau));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_TDAvec_computeECC", (DL_FUNC) &_TDAvec_computeECC, 3},
    {"_TDAvec_computePES", (DL_FUNC) &_TDAvec_computePES, 3},
    {"_TDAvec_computePI", (DL_FUNC) &_TDAvec_computePI, 8},
    {"_TDAvec_computePL", (DL_FUNC) &_TDAvec_computePL, 4},
    {"_TDAvec_computePS", (DL_FUNC) &_TDAvec_computePS, 4},
    {"_TDAvec_computeVAB", (DL_FUNC) &_TDAvec_computeVAB, 3},
    {"_TDAvec_computeVPB", (DL_FUNC) &_TDAvec_computeVPB, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_TDAvec(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
