// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
//
// Nelder_Mead.h: NLopt's Nelder-Mead optimizer, modified to use Eigen
//
// Copyright (C)       2011 Douglas Bates, Martin Maechler and Ben Bolker
//
// This file is part of lme4.

#ifndef LME4_NELDER_MEAD_H
#define LME4_NELDER_MEAD_H

#include <RcppEigen.h>

namespace optimizer {
    using Eigen::MatrixXd;
    using Eigen::VectorXd;
    using Eigen::VectorXi;

    typedef VectorXd::Scalar  Scalar;
    typedef VectorXd::Index   Index;

    class nl_stop {
    private:			// utilities
	bool relstop(const Scalar& vold, const Scalar& vnew,
		     const Scalar& reltol, const Scalar& abstol) const;
	Scalar sc(const Scalar& x, const Scalar& smin, const Scalar& smax) const {
	    return smin + x * (smax - smin);
	}
    protected:
	const VectorXd xtol_abs;
	unsigned n, nevals, maxeval;
	Scalar minf_max, ftol_rel, ftol_abs, xtol_rel;
	bool force_stop;
    public:
	nl_stop(const VectorXd&); // constructor

	void       incrEvals() {nevals++;}
	void      resetEvals() {nevals = 0;}	// setters
	void     setFtol_rel(const Scalar& ftr) {ftol_rel = ftr;}
	void     setFtol_abs(const Scalar& fta) {ftol_abs = fta;}
	void   setForce_stop(const bool& stp) {force_stop = stp;}
	void     setMinf_max(const Scalar& mm) {minf_max = mm;}
	void     set_Maxeval(const unsigned int& mm) {maxeval = mm;}
	int      get_Maxeval() const {return maxeval;}

	bool     f(const Scalar& f, const Scalar& oldf) const { // convergence checking
	    return (f <= minf_max || ftol(f, oldf));
	}
	bool  ftol(const Scalar& f, const Scalar& oldf) const {
	    return relstop(oldf, f, ftol_rel, ftol_abs);
	}
	bool         x(const VectorXd& x, const VectorXd& oldx) const;
	bool        dx(const VectorXd& x, const VectorXd& dx) const;
	bool        xs(const VectorXd& xs, const VectorXd& oldxs,
		       const VectorXd& scale_min, const VectorXd& scale_max) const;
	bool     evals() const {return maxeval > 0 && nevals > maxeval;}
	bool    forced() const {return force_stop;}

	int         ev() const {return nevals;}

	Scalar minfMax() const {return minf_max;}
    };

    inline bool nl_stop::relstop(const Scalar& vold, const Scalar& vnew,
			  const Scalar& reltol, const Scalar& abstol) const {
	if (std::abs(vold) == std::numeric_limits<Scalar>::infinity()) return false;
	return std::abs(vnew - vold) < abstol
	    || std::abs(vnew - vold) < reltol * (std::abs(vnew) + std::abs(vold)) * 0.5
	    || (reltol > 0 && vnew == vold);
    }
    
    enum nm_status {nm_active, nm_x0notfeasible, nm_nofeasible, nm_forced, nm_minf_max,
		    nm_evals, nm_fcvg, nm_xcvg};

    enum nm_stage {nm_restart, nm_postreflect, nm_postexpand, nm_postcontract};

    /* heuristic "strategy" constants: */
    static const double alpha = 1, beta = 0.5, gamm = 2, delta = 0.5;

    class Nelder_Mead {
    private:
	Scalar           f_old;
	Index         init_pos;
	nm_status         init(const Scalar&);
	nm_status      restart(const Scalar&);
	bool         reflectpt(VectorXd&, const VectorXd&, const Scalar&, const VectorXd&);
	nm_status  postreflect(const Scalar&);
	nm_status   postexpand(const Scalar&);
	nm_status postcontract(const Scalar&);
    protected:
	const VectorXd     d_lb; /*<< lower bounds  */
	const VectorXd     d_ub; /*<< upper bounds  */
	const VectorXd  d_xstep; /*<< initial step sizes  */
	VectorXd            d_x; /*<< initial value and optimum */
	Index              d_ih; /**< index in d_vals of largest value */
	Index              d_il; /**< index in d_vals of smallest value */
	Index               d_n; /**< size of parameter vector */
	MatrixXd          d_pts; /*<< points */
	VectorXd         d_vals; /*<< function values */
	VectorXd            d_c; /*<< centroid  */
	VectorXd         d_xcur; /*<< current x */
	VectorXd        d_xeval; /*<< x at which next evaluation is requested */
	Scalar             d_fl, d_fh, d_minf;
	nm_status        d_stat;
	nm_stage        d_stage;
	nl_stop          d_stop;
	Index            d_verb; /**< verbosity, if > 0 results are displayed every d_verb evaluations */
    public:
	Nelder_Mead(const VectorXd&, const VectorXd&, const VectorXd&, const VectorXd&,
		    const nl_stop&);

	const MatrixXd&   pts() const {return d_pts;}

	const VectorXd&    lb() const {return d_lb;}
	const VectorXd&    ub() const {return d_ub;}
	const VectorXd&  vals() const {return d_vals;}
	const VectorXd& xstep() const {return d_xstep;}
	const VectorXd& xeval() const {return d_xeval;}
	const VectorXd&  xpos() const {return d_x;}
	Index              ih() const {return d_ih;}
	Index              il() const {return d_il;}
	Index           evals() const {return d_stop.ev();}

	Scalar          value() const {return d_minf;}
	nm_status        newf(const Scalar&);

	void    setForce_stop(const bool& stp)        {d_stop.setForce_stop(stp);}
	void      setFtol_abs(const Scalar& fta)      {d_stop.setFtol_abs(fta);}
	void      setFtol_rel(const Scalar& ftr)      {d_stop.setFtol_rel(ftr);}
	void      set_Maxeval(const unsigned int& mm) {d_stop.set_Maxeval(mm);}
	void      set_Iprint(const int& ip)           {d_verb = ip;}
	void      setMinf_max(const Scalar& mm)       {d_stop.setMinf_max(mm);}
    };

    class Golden {
    protected:
	Scalar           d_invratio, d_lower, d_upper;
	Eigen::Vector2d  d_x, d_f;
	bool             d_init, d_ll;
    public:
	Golden(const Scalar&, const Scalar&);
	void    newf(const Scalar&);
	Scalar xeval() const {return d_x[d_ll ? 0 : 1];}
	Scalar value() const {return d_f[0];}
	Scalar  xpos() const {return d_x[0];}
    };
}

#endif // LME4_NELDER_MEAD_H
