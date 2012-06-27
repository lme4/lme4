//
// Nelder_Mead.cpp: implementation of Nelder-Mead optimization algorithm
//
// Based on the files nldrmd.h, nldrmd.c, stop.c and nlopt-util.h from NLopt 2.2.4 
// Steven G. Johnson, The NLopt nonlinear-optimization package,
// http://ab-initio.mit.edu/nlopt
//
// Original implementation Copyright (C) 2007-2011 Massachusetts Institute of Technology
// Modifications Copyright (C) 2011 Douglas Bates, Martin Maechler and Ben Bolker
//
// This file is part of lme4.

#include "optimizer.h"

namespace optimizer {
    using     std::invalid_argument;
    using     std::runtime_error;

    typedef VectorXd::Scalar  Scalar;
    typedef VectorXd::Index   Index;

/** 
 * Determine if two values are approximately equal relative to floating-point precision
 * 
 * @param a first value to compare
 * @param b second value to compare
 * 
 * @return true if a and b are approximately equal else false
 */
    static bool close(const Scalar& a, const Scalar& b) {
	return (std::abs(a - b) <= 1e-13 * (std::abs(a) + std::abs(b)));
    }

/** 
 *
 * 
 * @param lb lower bounds
 * @param ub upper bounds
 * @param xstep initial step sizes
 * @param x initial parameter vector
 * @param f value of function at initial parameter vector
 */
    Nelder_Mead::Nelder_Mead(const VectorXd& lb, const VectorXd& ub,
			     const VectorXd& xstep, const VectorXd& x, 
			     const nl_stop& stp) 
	: d_lb(    lb),
	  d_ub(    ub),
	  d_xstep( xstep),
	  d_x(     x),
	  d_n(     x.size()),
	  d_pts(   d_n, d_n + 1),
	  d_vals(  d_n + 1),
	  d_c(     d_n),
	  d_xcur(  d_n),
	  d_xeval( x),
	  d_minf(  std::numeric_limits<Scalar>::infinity()),
	  d_stage( nm_restart),
	  d_stop(  stp),
	  d_verb(  10)
    {
	d_stop.setForce_stop( 0); // BMB: this was undefined, bad news on Win64
	if (!d_n || d_lb.size() != d_n ||
	    d_ub.size() != d_n || d_xstep.size() != d_n)
	    throw invalid_argument("dimension mismatch");
	if (((d_x - d_lb).array() < 0).any() || ((d_ub - d_x).array() < 0).any())
	    throw std::invalid_argument("initial x is not a feasible point");
	d_stop.resetEvals();
	init_pos = 0;
	for (int i = 0; i <= d_n; ++i) d_vals[i] = std::numeric_limits<double>::min();
	d_pts = d_x.replicate(1, d_n + 1);
	for (Index i = 0; i < d_n; ++i) { // generate and check the initial positions
	    Index j(i + 1);
	    d_pts(i, j) += d_xstep[i];
	    if (d_pts(i, j) > d_ub[i]) {
		d_pts(i, j) =  (d_ub[i] - d_x[i] > std::abs(d_xstep[i]) * 0.1) ? d_ub[i]
				// ub is too close to pt, go in other direction
		    : d_x[i] - std::abs(d_xstep[i]);
	    }
	    if (d_pts(i,j) < d_lb[i]) {
		if (d_x[i] - d_lb[i] > std::abs(d_xstep[i]) * 0.1) d_pts(i,j) = d_lb[i];
		else { // lb is too close to pt, go in other direction
		    d_pts(i,j) = d_x[i] + std::abs(d_xstep[i]);
		    if (d_pts(i, j) > d_ub[i]) // go towards farther of lb, ub */
			d_pts(i, j) = 0.5 * ((d_ub[i] - d_x[i] > d_x[i] - d_lb[i] ?
					      d_ub[i] : d_lb[i]) + d_x[i]);
		}
	    }
	    if (close(d_pts(i,j), d_x[i]))
		throw std::invalid_argument("cannot generate feasible simplex");
	}
    }
    
/** 
 * Install the function value at d_xeval
 * 
 * @param f value of function at d_xeval
 * 
 * @return status
 */
    nm_status Nelder_Mead::newf(const Scalar& f) {
	d_stop.incrEvals();
	if (d_verb > 0 && (d_stop.ev() % d_verb) == 0)
	    Rcpp::Rcout << "(NM) " << d_stop.ev() << ": " << "f = " << value() << " at " << d_x.adjoint() << std::endl;
	if (d_stop.forced()) {
	    if (d_verb==1) {
		Rcpp::Rcout << "(NM) stop_forced" << std::endl;
	    }
	    return nm_forced;
	}
	if (f < d_minf) {
	    d_minf = f;
	    d_x = d_xeval;	// save the value generating current minimum
	    if (d_minf < d_stop.minfMax()) {
		if (d_verb==1) {
		    Rcpp::Rcout << "(NM) nm_minf_max: " << d_minf << ", " 
				<< d_stop.minfMax() << ", " << d_x << std::endl;
		}
		return nm_minf_max;
	    }
	}
	if (d_stop.evals()) {
	    if (d_verb==1) {
		Rcpp::Rcout << "(NM) nm_evals" << std::endl;
	    }
	    return nm_evals;
	}
	if (init_pos <= d_n) {
	    if (d_verb==1) {
		Rcpp::Rcout << "(NM) init_pos <= d_n" << std::endl;
	    }
	    return init(f);
	}
	switch (d_stage) {
	case nm_restart:      return restart(f);
	case nm_postreflect:  return postreflect(f);
	case nm_postexpand:   return postexpand(f);
	case nm_postcontract: return postcontract(f);
	}
	return nm_active;	// -Wall
    }

/** 
 * Initialization of d_vals from the positions in d_pts;
 * 
 * @param f function value
 * 
 * @return status
 */
    nm_status Nelder_Mead::init(const Scalar& f) {
	if (init_pos > d_n) throw std::runtime_error("init called after n evaluations");
	d_vals[init_pos++] = f;
	if (init_pos > d_n) return restart(f);
	d_xeval = d_pts.col(init_pos);
	return nm_active;
    }
    
/** 
 * Recompute the high/low function values (d_fh and d_fl) and indices
 * (d_ih and d_il) plus the centroid of the n-1 simplex opposite the high
 * vertex.  Check if the simplex has collapsed.  If not, attempt a reflection.
 * 
 * @param f function value
 * 
 * @return status
 */
    nm_status Nelder_Mead::restart(const Scalar& f) {
	d_fl = d_vals.minCoeff(&d_il);
	d_fh = d_vals.maxCoeff(&d_ih);
	d_c = (d_pts.rowwise().sum() - d_pts.col(d_ih)) / d_n; // compute centroid

	// check for x convergence by calculating the maximum absolute
	// deviation from the centroid for each coordinate in the simplex
	if (d_stop.x(VectorXd::Constant(d_n, 0.),
		     (d_pts.colwise() - d_c).array().abs().rowwise().maxCoeff()))
	    return nm_xcvg;

	if (!reflectpt(d_xcur, d_c, alpha, d_pts.col(d_ih))) return nm_xcvg;
	d_xeval = d_xcur;
	d_stage = nm_postreflect;
	return nm_active;
    }

    nm_status Nelder_Mead::postreflect(const Scalar& f) {
	// Rcpp::Rcout << "postreflect: ";
	if (f < d_fl) {	// new best point, try to expand
	    if (!reflectpt(d_xeval, d_c, gamm, d_pts.col(d_ih))) return nm_xcvg;
	    // Rcpp::Rcout << "New best point" << std::endl;
	    d_stage = nm_postexpand;
	    f_old = f;
	    return nm_active;
	}
	if (f < d_fh) {		// accept new point
	    // Rcpp::Rcout << "Accept new point" << std::endl;
	    d_vals[d_ih] = f;
	    d_pts.col(d_ih) = d_xeval;
	    return restart(f);
	}
				// new worst point, contract
	// Rcpp::Rcout << "New worst point" << std::endl;
	if (!reflectpt(d_xcur, d_c, d_fh <= f ? -beta : beta, d_pts.col(d_ih))) return nm_xcvg;
	f_old = f;
	d_xeval = d_xcur;
	d_stage = nm_postcontract;
	return nm_active;
    }

    nm_status Nelder_Mead::postexpand(const Scalar& f) {
	if (f < d_vals[d_ih]) { // expanding improved
	    // Rcpp::Rcout << "successful expand" << std::endl;
	    d_pts.col(d_ih) = d_xeval;
	    d_vals[d_ih]    = f;
	} else {
	    // Rcpp::Rcout << "unsuccessful expand" << std::endl;
	    d_pts.col(d_ih) = d_xcur;
	    d_vals[d_ih]    = f_old;
	}
	return restart(f);
    }

    nm_status Nelder_Mead::postcontract(const Scalar& f) {
	if (f < f_old && f < d_fh) {
	    // Rcpp::Rcout << "successful contraction:" << std::endl;
	    d_pts.col(d_ih) = d_xeval;
	    d_vals[d_ih] = f;
	    return restart(f);
	}
	// Rcpp::Rcout << "unsuccessful contraction, shrink simplex" << std::endl;
	for (Index i = 0; i <= d_n; ++i) {
	    if (i != d_il) {
		if (!reflectpt(d_xeval, d_pts.col(d_il), -delta, d_pts.col(i))) return nm_xcvg;
		d_pts.col(i) = d_xeval;
	    }
	}
	init_pos = 0;
	d_xeval = d_pts.col(0);
	return nm_active;
    }

/* Perform the reflection xnew = c + scale * (c - xold),
   returning 0 if xnew == c or xnew == xold (coincident points), 1 otherwise.

   The reflected point xnew is "pinned" to the lower and upper bounds
   (lb and ub), as suggested by J. A. Richardson and J. L. Kuester,
   "The complex method for constrained optimization," Commun. ACM
   16(8), 487-489 (1973).  This is probably a suboptimal way to handle
   bound constraints, but I don't know a better way.  The main danger
   with this is that the simplex might collapse into a
   lower-dimensional hyperplane; this danger can be ameliorated by
   restarting (as in subplex), however. */
    bool Nelder_Mead::reflectpt(VectorXd& xnew, const VectorXd& c,
				const Scalar& scale, const VectorXd& xold) {
	xnew = c + scale * (c - xold);
	bool equalc = true, equalold = true;
	for (Index i = 0; i < d_n; ++i) {
	    Scalar newx = std::min(std::max(xnew[i], d_lb[i]), d_ub[i]);
	    equalc = equalc && close(newx, d_c[i]);
	    equalold = equalold && close(newx, xold[i]);
	    xnew[i] = newx;
	}
	return !(equalc || equalold);
    }

    nl_stop::nl_stop(const VectorXd& xtol)
	: xtol_abs( xtol),
	  maxeval(  300),
	  minf_max( std::numeric_limits<Scalar>::min()),
	  ftol_rel( 1e-15),
	  xtol_rel( 1e-7) {
    }
    
    bool nl_stop::x(const VectorXd& x, const VectorXd& oldx) const {
	for (Index i = 0; i < x.size(); ++i)
	    if (!relstop(oldx[i], x[i], xtol_rel, xtol_abs[i])) return false;
	return true;
    }	     

    bool nl_stop::dx(const VectorXd& x, const VectorXd& dx) const {
	for (Index i = 0; i < x.size(); ++i)
	    if (!relstop(x[i] - dx[i], x[i], xtol_rel, xtol_abs[i])) return false;
	return true;
    }

    bool nl_stop::xs(const VectorXd& xs, const VectorXd& oldxs,
		     const VectorXd& scale_min, const VectorXd& scale_max) const {
	for (Index i = 0; i < xs.size(); ++i)
	    if (relstop(sc(oldxs[i], scale_min[i], scale_max[i]), 
			sc(xs[i], scale_min[i], scale_max[i]),
			xtol_rel, xtol_abs[i]))
		return true;
	return false;
    }

    Golden::Golden(const Scalar& lower, const Scalar& upper)
	: d_lower(lower), d_upper(upper) {
	if (lower >= upper)
	    throw invalid_argument("lower >= upper");
	d_invratio = 2./(1. + std::sqrt(5.));
	double range = upper - lower;
	d_x[0] = lower + range * (1. - d_invratio);
	d_x[1] = lower + range * d_invratio;
	d_init = true;
	d_ll   = true;
    }

    void Golden::newf(const Scalar& fv) {
	Rcpp::Rcout << "f = " << fv << " at x = " << xeval() << std::endl;
	d_f[d_ll ? 0 : 1] = fv;
	if (d_init) {
	    d_init = false;
	    d_ll   = false;
	    return;
	}
	if (d_f[0] > d_f[1]) {	// discard left portion of interval
	    d_lower = d_x[0];
	    d_x[0]  = d_x[1];
	    d_f[0]  = d_f[1];
	    d_x[1]  = d_lower + (d_upper - d_lower) * d_invratio;
	    d_ll    = false;
	} else {
	    d_upper = d_x[1];
	    d_x[1]  = d_x[0];
	    d_f[1]  = d_f[0];
	    d_x[0]  = d_lower + (d_upper - d_lower) * (1 - d_invratio);
	    d_ll    = true;
	}
    }
    
}

