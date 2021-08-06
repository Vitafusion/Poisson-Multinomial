// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/PMD.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <string>
#include <set>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// pmn_mdfft_arma
arma::vec pmn_mdfft_arma(int nnt, arma::mat pp, arma::vec nn_vec, arma::vec l_vec, arma::vec cn_vec);
static SEXP _PMD_pmn_mdfft_arma_try(SEXP nntSEXP, SEXP ppSEXP, SEXP nn_vecSEXP, SEXP l_vecSEXP, SEXP cn_vecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< int >::type nnt(nntSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type pp(ppSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type nn_vec(nn_vecSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type l_vec(l_vecSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type cn_vec(cn_vecSEXP);
    rcpp_result_gen = Rcpp::wrap(pmn_mdfft_arma(nnt, pp, nn_vec, l_vec, cn_vec));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _PMD_pmn_mdfft_arma(SEXP nntSEXP, SEXP ppSEXP, SEXP nn_vecSEXP, SEXP l_vecSEXP, SEXP cn_vecSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_PMD_pmn_mdfft_arma_try(nntSEXP, ppSEXP, nn_vecSEXP, l_vecSEXP, cn_vecSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// rpmd_arma
arma::vec rpmd_arma(arma::mat pp);
static SEXP _PMD_rpmd_arma_try(SEXP ppSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< arma::mat >::type pp(ppSEXP);
    rcpp_result_gen = Rcpp::wrap(rpmd_arma(pp));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _PMD_rpmd_arma(SEXP ppSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_PMD_rpmd_arma_try(ppSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// pmd_simulation_singlepoint
double pmd_simulation_singlepoint(arma::mat pp, arma::vec x_vec, int t);
static SEXP _PMD_pmd_simulation_singlepoint_try(SEXP ppSEXP, SEXP x_vecSEXP, SEXP tSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< arma::mat >::type pp(ppSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type x_vec(x_vecSEXP);
    Rcpp::traits::input_parameter< int >::type t(tSEXP);
    rcpp_result_gen = Rcpp::wrap(pmd_simulation_singlepoint(pp, x_vec, t));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _PMD_pmd_simulation_singlepoint(SEXP ppSEXP, SEXP x_vecSEXP, SEXP tSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_PMD_pmd_simulation_singlepoint_try(ppSEXP, x_vecSEXP, tSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// pmd_simulation_allpoints
arma::vec pmd_simulation_allpoints(arma::mat pp, int nnt, arma::vec l_vec, arma::vec cn_vec, int t);
static SEXP _PMD_pmd_simulation_allpoints_try(SEXP ppSEXP, SEXP nntSEXP, SEXP l_vecSEXP, SEXP cn_vecSEXP, SEXP tSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< arma::mat >::type pp(ppSEXP);
    Rcpp::traits::input_parameter< int >::type nnt(nntSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type l_vec(l_vecSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type cn_vec(cn_vecSEXP);
    Rcpp::traits::input_parameter< int >::type t(tSEXP);
    rcpp_result_gen = Rcpp::wrap(pmd_simulation_allpoints(pp, nnt, l_vec, cn_vec, t));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _PMD_pmd_simulation_allpoints(SEXP ppSEXP, SEXP nntSEXP, SEXP l_vecSEXP, SEXP cn_vecSEXP, SEXP tSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_PMD_pmd_simulation_allpoints_try(ppSEXP, nntSEXP, l_vecSEXP, cn_vecSEXP, tSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}

// validate (ensure exported C++ functions exist before calling them)
static int _PMD_RcppExport_validate(const char* sig) { 
    static std::set<std::string> signatures;
    if (signatures.empty()) {
        signatures.insert("arma::vec(*pmn_mdfft_arma)(int,arma::mat,arma::vec,arma::vec,arma::vec)");
        signatures.insert("arma::vec(*rpmd_arma)(arma::mat)");
        signatures.insert("double(*pmd_simulation_singlepoint)(arma::mat,arma::vec,int)");
        signatures.insert("arma::vec(*pmd_simulation_allpoints)(arma::mat,int,arma::vec,arma::vec,int)");
    }
    return signatures.find(sig) != signatures.end();
}

// registerCCallable (register entry points for exported C++ functions)
RcppExport SEXP _PMD_RcppExport_registerCCallable() { 
    R_RegisterCCallable("PMD", "_PMD_pmn_mdfft_arma", (DL_FUNC)_PMD_pmn_mdfft_arma_try);
    R_RegisterCCallable("PMD", "_PMD_rpmd_arma", (DL_FUNC)_PMD_rpmd_arma_try);
    R_RegisterCCallable("PMD", "_PMD_pmd_simulation_singlepoint", (DL_FUNC)_PMD_pmd_simulation_singlepoint_try);
    R_RegisterCCallable("PMD", "_PMD_pmd_simulation_allpoints", (DL_FUNC)_PMD_pmd_simulation_allpoints_try);
    R_RegisterCCallable("PMD", "_PMD_RcppExport_validate", (DL_FUNC)_PMD_RcppExport_validate);
    return R_NilValue;
}

static const R_CallMethodDef CallEntries[] = {
    {"_PMD_pmn_mdfft_arma", (DL_FUNC) &_PMD_pmn_mdfft_arma, 5},
    {"_PMD_rpmd_arma", (DL_FUNC) &_PMD_rpmd_arma, 1},
    {"_PMD_pmd_simulation_singlepoint", (DL_FUNC) &_PMD_pmd_simulation_singlepoint, 3},
    {"_PMD_pmd_simulation_allpoints", (DL_FUNC) &_PMD_pmd_simulation_allpoints, 5},
    {"_PMD_RcppExport_registerCCallable", (DL_FUNC) &_PMD_RcppExport_registerCCallable, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_PMD(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
