/*! \file epot_bicgstabsolver.cpp
 *  \brief BiCGSTAB matrix solver for electric potential problem
 */

/* Copyright (c) 2005-2013 Taneli Kalvas. All rights reserved.
 *
 * You can redistribute this software and/or modify it under the terms
 * of the GNU General Public License as published by the Free Software
 * Foundation; either version 2 of the License, or (at your option)
 * any later version.
 * 
 * This library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this library (file "COPYING" included in the package);
 * if not, write to the Free Software Foundation, Inc., 51 Franklin
 * Street, Fifth Floor, Boston, MA 02110-1301 USA
 * 
 * If you have questions about your rights to use or distribute this
 * software, please contact Berkeley Lab's Technology Transfer
 * Department at TTD@lbl.gov. Other questions, comments and bug
 * reports should be sent directly to the author via email at
 * taneli.kalvas@jyu.fi.
 * 
 * NOTICE. This software was developed under partial funding from the
 * U.S.  Department of Energy.  As such, the U.S. Government has been
 * granted for itself and others acting on its behalf a paid-up,
 * nonexclusive, irrevocable, worldwide license in the Software to
 * reproduce, prepare derivative works, and perform publicly and
 * display publicly.  Beginning five (5) years after the date
 * permission to assert copyright is obtained from the U.S. Department
 * of Energy, and subject to any subsequent five (5) year renewals,
 * the U.S. Government is granted for itself and others acting on its
 * behalf a paid-up, nonexclusive, irrevocable, worldwide license in
 * the Software to reproduce, prepare derivative works, distribute
 * copies to the public, perform publicly and display publicly, and to
 * permit others to do so.
 */

#include <chrono>
#include <limits>
#include "compmath.hpp"
#include "hbio.hpp"
#include "epot_bicgstabsolver.hpp"
#include "statusprint.hpp"
#include "ibsimu.hpp"

using Clock = std::chrono::steady_clock;

// All of these are reset to 0 at the start of every subsolve() call (see
// its top) so the "Timing summary" printed at the end of subsolve()
// reflects that one Poisson solve only, not an accumulation across every
// major cycle of the whole program run.
static double time_spmv = 0.0;
static double time_prec = 0.0;
static double time_dot = 0.0;
static double time_norm = 0.0;
static double time_axpy = 0.0;
static double time_gmg_prepare = 0.0;
static double time_gmg_construct = 0.0;

EpotBiCGSTABSolver::EpotBiCGSTABSolver( Geometry &geom, 
					double eps, 
					uint32_t imax,
					double newton_eps,
					uint32_t newton_imax,
					bool gnewton )
    : EpotMatrixSolver(geom), _eps(eps), _imax(imax), _iter(0), 
      _res(0.0), _err(0.0), _newton_res(0.0), _newton_res_norm(0.0), 
      _newton_step(0.0), _newton_step_norm(0.0), _gnewton(gnewton), 
      _newton_eps(newton_eps), _newton_imax(newton_imax),
      _epot(0), _callback(0), _callback_nonlinear(0)
{
    if( eps <= 0.0 || newton_eps <= 0.0 )
        throw( ErrorDim( ERROR_LOCATION, "invalid accuracy request" ) );

    // Default preconditioner: ILU0.
    //
    // ILU0 gives a better-conditioned system per BiCGSTAB iteration
    // than the alternatives below, at the cost of being a serial
    // (sequential forward/backward substitution) operation that scales
    // poorly with problem size -- it is kept as the default for
    // unchanged convergence behavior, not because it's the fastest
    // option at high resolution.
    //
    // A parallel alternative is available and generally faster overall
    // for high-resolution 3D problems (fewer serial bottlenecks per
    // BiCGSTAB iteration, at the cost of needing more iterations since
    // it's a weaker preconditioner per application):
    //   solver.set_preconditioner( RBSOR_Precond() );
    // RBSOR_Precond(nsweeps=4, w=1.4) performs red-black Gauss-Seidel
    // sweep pairs, parallelised across ibsimu.get_thread_count() threads.
    _pc = new ILU0_Precond();
}


EpotBiCGSTABSolver::EpotBiCGSTABSolver( Geometry &geom, std::istream &s )
    : EpotMatrixSolver(geom,s)
{
    _pc = new ILU0_Precond();
    throw( ErrorUnimplemented( ERROR_LOCATION ) );

}


void EpotBiCGSTABSolver::save( std::ostream &s ) const
{
    throw( ErrorUnimplemented( ERROR_LOCATION ) );
}


EpotBiCGSTABSolver::~EpotBiCGSTABSolver()
{
    if( _pc )
	delete _pc;
}


void EpotBiCGSTABSolver::set_preconditioner( Precond &pc )
{
    if( _pc )
	delete _pc;
    _pc = pc.copy();
}


void EpotBiCGSTABSolver::set_gnewton( bool enable ) 
{
    _gnewton = enable;
}


void EpotBiCGSTABSolver::set_eps( double eps ) 
{
    if( eps <= 0.0 )
        throw( ErrorDim( ERROR_LOCATION, "invalid accuracy request" ) );
    _eps = eps;
}


void EpotBiCGSTABSolver::set_imax( uint32_t imax ) 
{
    _imax = imax;
}


void EpotBiCGSTABSolver::set_newton_imax( uint32_t newton_imax ) 
{
    _newton_imax = newton_imax;
}


void EpotBiCGSTABSolver::set_newton_eps( double newton_eps ) 
{
    if( newton_eps <= 0.0 )
        throw( ErrorDim( ERROR_LOCATION, "invalid accuracy request" ) );
    _newton_eps = newton_eps;
}


double EpotBiCGSTABSolver::get_scaled_residual( void ) const
{
    return( _res );
}


double EpotBiCGSTABSolver::get_error_estimate( void ) const
{
    if( linear() )
	return( _err );
    else {
	double maxsize = _geom.size().max();
	return( maxsize * _newton_res / 1e7 );
    }
}


uint32_t EpotBiCGSTABSolver::get_iter( void ) const
{
    return( _iter );
}


double EpotBiCGSTABSolver::get_newton_residual( void ) const
{
    return( _newton_res );
}


double EpotBiCGSTABSolver::get_newton_residual_norm( void ) const
{
    return( _newton_res_norm );
}


double EpotBiCGSTABSolver::get_newton_step( void ) const
{
    return( _newton_step );
}


double EpotBiCGSTABSolver::get_newton_step_norm( void ) const
{
    return( _newton_step_norm );
}


void EpotBiCGSTABSolver::reset_problem( void )
{
    reset_matrix();
    _pc->clear();
}


void EpotBiCGSTABSolver::bicgstab( const Matrix &mat, const Vector &rhs, Vector &sol,
				   const Precond &pc )
{
   
    // Checks
    if( mat.columns() != mat.rows() )
	throw( ErrorDim( ERROR_LOCATION, "matrix not square" ) );
    if( mat.rows() != rhs.size() )
	throw( ErrorDim( ERROR_LOCATION, "matrix dimension does not match vector" ) );

    double omega = 0, alpha = 0, beta, rho_1, rho_2 = 0;
    Vector p, phat, s, shat, t, v;
    double maxsize = _geom.size().max();
    double errscale = maxsize*maxsize;
	auto t0 = Clock::now();
    double norm_rhs = norm2(rhs);
	auto t1 = Clock::now();
    time_norm += std::chrono::duration<double>(t1-t0).count();
    if( sol.size() != mat.columns() ) {
	sol.resize( mat.columns() );
	sol.clear();
    }

	 t0 = Clock::now();
    Vector r = mat * sol;
	 t1 = Clock::now();
	time_spmv += std::chrono::duration<double>(t1-t0).count();

    r = rhs - r;
    Vector rtilde = r;

    if( norm_rhs == 0.0 )
	norm_rhs = 1;
	 t0 = Clock::now();
    _res = norm2(r) / norm_rhs;
	 t1 = Clock::now();
    time_norm += std::chrono::duration<double>(t1-t0).count();
    _err = errscale*_res;
    if( _err <= _eps )
	return;

    StatusPrint sp;
    if( ibsimu.get_message_threshold(MSG_VERBOSE) && ibsimu.output_is_cout() ) {
	std::stringstream ss;
	ss << "  " << std::setw(5) << 0 << " " << std::scientific << std::setw(20) << _err;
	sp.print( ss.str() );
    }

    uint32_t i = 0; // Local iteration counter
    while( _iter < _imax ) {
	 t0 = Clock::now();
	rho_1 = dot_prod( rtilde, r );
	 t1 = Clock::now();
    time_dot += std::chrono::duration<double>(t1-t0).count();
	if( rho_1 == 0 )
	    break;
	if( i == 0 )
	    p = r;
	else {
		 t0 = Clock::now();
	    beta = (rho_1/rho_2) * (alpha/omega);
	    p = r + beta * (p - omega * v);
		 t1 = Clock::now();
    	time_axpy += std::chrono::duration<double>(t1-t0).count();
	}
	 t0 = Clock::now();
	pc.solve( phat, p );
	 t1 = Clock::now();
    time_prec += std::chrono::duration<double>(t1-t0).count();
	 t0 = Clock::now();
	v = mat * phat;
	 t1 = Clock::now();
    time_spmv += std::chrono::duration<double>(t1-t0).count();
	 t0 = Clock::now();
	alpha = rho_1 / dot_prod( rtilde, v );
	 t1 = Clock::now();
    time_dot += std::chrono::duration<double>(t1-t0).count();
	 t0 = Clock::now();
	s = r - alpha * v;
	 t1 = Clock::now();
    time_axpy += std::chrono::duration<double>(t1-t0).count();
	i++;
	_iter++;
	 t0 = Clock::now();
	_res = norm2(s)/norm_rhs;
	 t1 = Clock::now();
    time_norm += std::chrono::duration<double>(t1-t0).count();
	_err = errscale*_res;
	if( _err < _eps ) {
		 t0 = Clock::now();
	    sol += alpha * phat;
		 t1 = Clock::now();
    	time_axpy += std::chrono::duration<double>(t1-t0).count();
	    break;
	}
	 t0 = Clock::now();
	pc.solve( shat, s );
	 t1 = Clock::now();
    time_prec += std::chrono::duration<double>(t1-t0).count();
	 t0 = Clock::now();
	t = mat * shat;
	 t1 = Clock::now();
    time_spmv += std::chrono::duration<double>(t1-t0).count();
	 t0 = Clock::now();
	omega = dot_prod( t, s ) / dot_prod( t, t );
	 t1 = Clock::now();
    time_dot += std::chrono::duration<double>(t1-t0).count();
	 t0 = Clock::now();
	sol += alpha * phat + omega * shat;
	r = s - omega * t;
	 t1 = Clock::now();
    time_axpy += std::chrono::duration<double>(t1-t0).count();

	rho_2 = rho_1;
	 t0 = Clock::now();
	_res = norm2(r) / norm_rhs;
	 t1 = Clock::now();
    time_norm += std::chrono::duration<double>(t1-t0).count();
	_err = errscale*_res;
	if( _err < _eps )
	    break;
	if( comp_isnan(_res) || omega == 0 )
	    throw( Error( ERROR_LOCATION, "convergence failure" ) );

	if( ibsimu.get_message_threshold(MSG_VERBOSE) && ibsimu.output_is_cout() ) {
	    std::stringstream ss;
	    ss << "  " << std::setw(5) << i << " " << std::scientific << std::setw(20) << _err;
	    sp.print( ss.str() );
	}

	if( _callback ) {
	    set_solution( *_epot, sol );
	    _callback();
	}
    }

    // Force print last
    if( ibsimu.get_message_threshold(MSG_VERBOSE) && ibsimu.output_is_cout() ) {
	std::stringstream ss;
	ss << "  " << std::setw(5) << i << " " << std::scientific << std::setw(20) << _err;
	sp.print( ss.str(), true );
    }
    
    return;
}


void EpotBiCGSTABSolver::subsolve( MeshScalarField &epot, const MeshScalarField &scharge )
{
    _epot = &epot;

    // Zero every per-solve timing accumulator -- see their doc
    // comments/declarations. Must happen before preprocess()/build_mat_vec()
    // get a chance to run, since the linear-build/dielectric-cache timers
    // are accumulated from inside those.
    time_spmv = time_prec = time_dot = time_norm = time_axpy = 0.0;
    time_gmg_prepare = time_gmg_construct = 0.0;
    reset_build_timing();

    // Preprocess and set starting guess
    preprocess( epot, scharge );

	GMG_Precond *gmg = dynamic_cast<GMG_Precond*>( _pc );
	if( gmg ) {
		std::vector<int32_t> node_map;
		build_node_map( node_map );   // EpotMatrixSolver::build_node_map(), see earlier
		gmg->set_node_map( node_map );
	}
	
    Vector X;
    set_initial_guess( epot, X );
    double maxsize = _geom.size().max();

    if( linear() ) {

	ibsimu.message( 1 ) << "Using BiCGSTAB-" << _pc->typestring() << " solver ("
			    << "imax = " << _imax
			    << ", eps = " << _eps
			    << ", threads = " << ibsimu.get_thread_count()
			    << ")\n";
	ibsimu.flush();
	    
	// Fetch matrix form of problem
	const CRowMatrix *A;
	const Vector *B;
	get_vecmat( &A, &B );

	auto t0 = Clock::now();
	if( !_pc->is_prepared() )
	    _pc->prepare( *A );
	auto t1 = Clock::now();
	time_gmg_prepare += std::chrono::duration<double>(t1-t0).count();
	_pc->construct( *A );
	auto t2 = Clock::now();
	time_gmg_construct += std::chrono::duration<double>(t2-t1).count();

	_iter = 0;
        bicgstab( *A, *B, X, *_pc );
	ibsimu.message( 1 ) << "\n";

	if( _iter == _imax )
	    ibsimu.message( 1 ) << "Maximum number of iteration rounds done.\n";
	ibsimu.message( 1 ) << "error estimate = " << _err << "\n";
	ibsimu.message( 1 ) << "iterations = " << _iter << "\n";
	ibsimu.flush();

    } else {

	int32_t a;
	uint32_t iter_old = 0;
	double error_estimate = 0.0;
        const CRowMatrix *J;
        const Vector *R;
        Vector dX;

	_newton_res = 0.0;
	_newton_step = 0.0;
	_iter = 0;
	if( _gnewton ) {

	    ibsimu.message( 1 ) << "Using Globally convergent Newton-Raphson BiCGSTAB-" 
				<< _pc->typestring() << " solver ("
				<< "imax = " << _imax
				<< ", eps = " << _eps
				<< ", newton_imax = " << _newton_imax
				<< ", newton_eps = " << _newton_eps
				<< ", threads = " << ibsimu.get_thread_count()
				<< ")\n";
	    ibsimu.message( 1 ) << std::setw(5)  << "Round" << " " 
				<< std::setw(8)  << "Iter" << " " 
				<< std::setw(14) << "Step size" << " " 
				<< std::setw(14) << "Step fac" << " " 
				<< std::setw(14) << "Residual" << " "
				<< std::setw(14) << "Err estimate" << "\n";
	    ibsimu.flush();
	    
	    // Globally convergent Newton-Raphson
            Vector Xold( X.size() );

            // First jacobian and residual
            get_resjac( &J, &R, X );
            double f = ssqr( *R );

	    {
		auto t0 = Clock::now();
		if( !_pc->is_prepared() )
		    _pc->prepare( *J );
		auto t1 = Clock::now();
		time_gmg_prepare += std::chrono::duration<double>(t1-t0).count();
	    }

	    for( a = 0; a < (int)_newton_imax; a++ ) {

                // Calculate dX = J^{-1}*R
		{
		    auto t0 = Clock::now();
		    _pc->construct( *J );
		    auto t1 = Clock::now();
		    time_gmg_construct += std::chrono::duration<double>(t1-t0).count();
		}
                dX.clear();
		iter_old = _iter;
                bicgstab( *J, *R, dX, *_pc );

                // Search for acceptable step for which residual decreases
                double t = 2.0;
                double fold = f;
                Xold = X;
                while( f >= fold ) {

                    t *= 0.5;
                    X = Xold - t*dX;
                    get_resjac( &J, &R, X );
                    f = ssqr( *R );
                    if( t <= std::numeric_limits<double>::epsilon() )
                        break;
                }

                // Check for convergence
                _newton_res = max_abs( *R );
                _newton_res_norm = norm2( *R );
                _newton_step = t*max_abs( dX );
                _newton_step_norm = t*norm2( dX );
		error_estimate = maxsize * _newton_res / 1e7;

		ibsimu.message( 1 ) << "\r  " 
				    << std::setw(5)  << a << " " 
				    << std::setw(8)  << _iter-iter_old << " " 
				    << std::setw(14) << _newton_step << " " 
				    << std::setw(14) << t << " " 
				    << std::setw(14) << _newton_res << " "
				    << std::setw(14) << error_estimate << "\n";
		ibsimu.flush();
                
		if( _callback_nonlinear ) {
		    set_solution( *_epot, X );
		    _callback_nonlinear();
		}

		if( error_estimate < _newton_eps || _iter >= _imax )
                    break;
            }

	} else {

	    ibsimu.message( 1 ) << "Using Newton-Raphson BiCGSTAB-" << _pc->typestring() << " solver ("
				<< "imax = " << _imax
				<< ", eps = " << _eps
				<< ", newton_imax = " << _newton_imax
				<< ", threads = " << ibsimu.get_thread_count()
				<< ")\n";
	    ibsimu.message( 1 ) << std::setw(5)  << "Round" << " " 
				<< std::setw(8)  << "Iter" << " " 
				<< std::setw(14) << "Step size" << " " 
				<< std::setw(14) << "Residual" << " "
				<< std::setw(14) << "Err estimate" << "\n";
	    ibsimu.flush();

	    for( a = 0; a < (int)_newton_imax; a++ ) {

		// Calculate dX = J^{-1}*R
		get_resjac( &J, &R, X );
		auto t0 = Clock::now();
		if( !_pc->is_prepared() )
		    _pc->prepare( *J );
		auto t1 = Clock::now();
		time_gmg_prepare += std::chrono::duration<double>(t1-t0).count();
		_pc->construct( *J );
		auto t2 = Clock::now();
		time_gmg_construct += std::chrono::duration<double>(t2-t1).count();

		dX.clear();
		iter_old = _iter;
		bicgstab( *J, *R, dX, *_pc );

		// Take step
		X -= dX;

		// Check for convergence
		_newton_res = max_abs( *R );
                _newton_res_norm = norm2( *R );
		_newton_step = max_abs( dX );
                _newton_step_norm = norm2( dX );
		error_estimate = maxsize * _newton_res / 1e7;

		ibsimu.message( 1 ) << "\r  " 
				    << std::setw(5) << a << " " 
				    << std::setw(8) << _iter-iter_old << " " 
				    << std::setw(14) << _newton_step << " " 
				    << std::setw(14) << _newton_res << " "
				    << std::setw(14) << error_estimate << "\n";
		ibsimu.flush();
		
		if( _callback_nonlinear ) {
		    set_solution( *_epot, X );
		    _callback_nonlinear();
		}

		if( error_estimate < _newton_eps || _iter >= _imax )
		    break;
	    }
        }

	if( error_estimate < _newton_eps )
	    ibsimu.message( 1 ) << "Newton-Raphson converged\n";
	else if( _iter >= _imax )
	    ibsimu.message( 1 ) << "Maximum number of BiCGSTAB iterations\n";
	else
	    ibsimu.message( 1 ) << "Maximum number of Newton-Raphson iterations\n";
	
	ibsimu.message( 1 ) << "total iterations = " << _iter << "\n";

    }

    // Per-solve timing breakdown -- every value here was reset to 0 at
    // the top of this subsolve() call (see there), so this reflects
    // only this Poisson solve, not an accumulation across the whole
    // program run. Matrix-build figures come from EpotMatrixSolver
    // (this class' base) via build_mat_vec()'s own timing.
    ibsimu.message(1) << "\nTiming summary:\n";
    ibsimu.message(1)
    << "  Linear matrix build   : " << time_linbuild() << " s\n"
    << "  Nonlinear (plasma) update: " << time_nonlin() << " s\n"
    << "  GMG prepare           : " << time_gmg_prepare << " s\n"
    << "  GMG construct         : " << time_gmg_construct << " s\n"
    << "  SpMV                  : " << time_spmv << " s\n"
    << "  Preconditioner apply  : " << time_prec << " s\n"
    << "  Dot products          : " << time_dot << " s\n"
    << "  Norms                 : " << time_norm << " s\n"
    << "  Vector ops            : " << time_axpy << " s\n";

    // Postprocess and set solution
    set_solution( epot, X );
    postprocess();
}


void EpotBiCGSTABSolver::set_analysis_callback( void (*func)(void) )
{
    _callback = func;
}

    
void EpotBiCGSTABSolver::set_analysis_callback_nonlinear( void (*func)(void) )
{
    _callback_nonlinear = func;
}

    
void EpotBiCGSTABSolver::debug_print( std::ostream &os ) const
{
    EpotMatrixSolver::debug_print( os );
    os << "**EpotBiCGSTABSolver\n";


}
