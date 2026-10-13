/*! \file sheath.cpp
 *  \brief Planar plasma sheath: Bohm entry, Boltzmann electrons, cold ions
 */

/* This file is part of the IBSimu library.
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
 */

#include <cmath>
#include <sstream>
#include "sheath.hpp"
#include "constants.hpp"
#include "error.hpp"


double sheath_bohm_speed( double Te, double mass_amu, double charge_state )
{
    if( Te <= 0.0 || mass_amu <= 0.0 || charge_state <= 0.0 )
	throw( Error( ERROR_LOCATION, "Te, mass and charge state must be positive" ) );
    return( sqrt( charge_state*CHARGE_E*Te / (mass_amu*MASS_U) ) );
}


double sheath_edge_density( double J, double Te, double mass_amu, double charge_state,
			    double bohm_factor )
{
    if( J < 0.0 )
	throw( Error( ERROR_LOCATION, "current density must not be negative" ) );
    if( bohm_factor <= 0.0 )
	throw( Error( ERROR_LOCATION, "bohm_factor must be positive" ) );
    // Entry speed at kinetic energy bohm_factor*Te is the Bohm speed times
    // sqrt(2 f), so f = 0.5 reproduces the plain Bohm speed exactly.
    double v_entry = sheath_bohm_speed( Te, mass_amu, charge_state )
	* sqrt( 2.0*bohm_factor );
    return( J / (charge_state*CHARGE_E*v_entry) );
}


double sheath_debye_length( double Te, double n )
{
    if( Te <= 0.0 || n <= 0.0 )
	throw( Error( ERROR_LOCATION, "Te and density must be positive" ) );
    return( sqrt( EPSILON0*Te / (CHARGE_E*n) ) );
}


double sheath_G( double xi, double bohm_factor )
{
    // The Bohm criterion is a MINIMUM: f = 0.5 is marginal, f > 0.5 enters
    // supersonically. Below it G_f goes negative near the edge and no
    // monotonic sheath exists, so refuse rather than clamp.
    if( bohm_factor < 0.5 ) {
	std::ostringstream ss;
	ss << "bohm_factor = " << bohm_factor << " is below 0.5, violating the "
	   << "Bohm criterion: G would be negative near the sheath edge";
	throw( Error( ERROR_LOCATION, ss.str() ) );
    }
    // ErrorRange is for integer indexing; a continuous out-of-domain
    // value gets a plain Error with the offending number in the message.
    if( xi < -bohm_factor ) {
	std::ostringstream ss;
	ss << "xi = " << xi << " is below -bohm_factor = " << -bohm_factor
	   << ", where 1 + xi/f turns negative";
	throw( Error( ERROR_LOCATION, ss.str() ) );
    }
    return( 2.0*bohm_factor*sqrt( 1.0 + xi/bohm_factor )
	    - 2.0*bohm_factor - 1.0 + exp( -xi ) );
}


double sheath_field_factor( double xi, double bohm_factor )
{
    double g = sheath_G( xi, bohm_factor );
    if( g < 0.0 )
	g = 0.0;                      // rounding at xi ~ 0, where G -> xi^3/3
    return( sqrt( 2.0*g ) );
}


/* Integrand of the position quadrature, 1/sqrt(2 G(xi)).
 *
 * Diverges as sqrt(3/2) xi^(-3/2) at the sheath edge. sheath_distance()
 * refuses xi_from <= 0 for that reason, so the singularity is never inside
 * the interval being integrated.
 */
static double sheath_integrand( double xi, double bohm_factor )
{
    double g = sheath_G( xi, bohm_factor );
    if( g <= 0.0 )
	return( HUGE_VAL );
    return( 1.0 / sqrt( 2.0*g ) );
}


/* Adaptive Simpson. The integrand is smooth and slowly varying over most
 * of the range (it decays as (2 xi)^(-1/4) for large xi) but steepens
 * sharply as xi_from approaches the sheath edge, so a fixed rule would
 * either waste effort at the flat end or under-resolve the steep one.
 */
static double sheath_simpson( double a, double b, double fa, double fm, double fb,
			      double whole, double eps, int depth,
			      double bohm_factor )
{
    double m = 0.5*(a+b);
    double lm = 0.5*(a+m);
    double rm = 0.5*(m+b);
    double flm = sheath_integrand( lm, bohm_factor );
    double frm = sheath_integrand( rm, bohm_factor );

    double left  = (m-a)/6.0 * (fa + 4.0*flm + fm);
    double right = (b-m)/6.0 * (fm + 4.0*frm + fb);
    double delta = left + right - whole;

    if( depth <= 0 || fabs( delta ) <= 15.0*eps )
	return( left + right + delta/15.0 );

    return( sheath_simpson( a, m, fa, flm, fm, left, 0.5*eps, depth-1, bohm_factor ) +
	    sheath_simpson( m, b, fm, frm, fb, right, 0.5*eps, depth-1, bohm_factor ) );
}


double sheath_distance( double xi_from, double xi_to, double bohm_factor )
{
    if( xi_from <= 0.0 )
	throw( Error( ERROR_LOCATION, "xi_from must be positive: the distance "
		      "measured from the sheath edge (xi = 0) is infinite" ) );
    if( xi_to < xi_from )
	throw( Error( ERROR_LOCATION, "xi_to must not be less than xi_from" ) );
    if( xi_to == xi_from )
	return( 0.0 );

    double fa = sheath_integrand( xi_from, bohm_factor );
    double fb = sheath_integrand( xi_to, bohm_factor );
    double m  = 0.5*(xi_from + xi_to);
    double fm = sheath_integrand( m, bohm_factor );
    double whole = (xi_to - xi_from)/6.0 * (fa + 4.0*fm + fb);

    // Relative target: the result spans a wide range (a fraction of a
    // Debye length up to thousands of them for a high-voltage gap).
    double eps = 1.0e-10 * (whole > 0.0 ? whole : 1.0);

    return( sheath_simpson( xi_from, xi_to, fa, fm, fb, whole, eps, 50, bohm_factor ) );
}


double sheath_wall_xi_metal( double V_plasma, double Te )
{
    if( Te <= 0.0 )
	throw( Error( ERROR_LOCATION, "Te must be positive" ) );
    if( V_plasma <= 0.0 )
	throw( Error( ERROR_LOCATION, "plasma potential must be above the wall" ) );
    return( V_plasma / Te );
}


double sheath_wall_xi_dielectric( double V_plasma, double Te, double lambda_D,
				  double eps_r, double d, double bohm_factor )
{
    if( Te <= 0.0 || lambda_D <= 0.0 || eps_r <= 0.0 || d <= 0.0 )
	throw( Error( ERROR_LOCATION, "Te, lambda_D, eps_r and thickness must "
		      "be positive" ) );
    if( V_plasma <= 0.0 )
	throw( Error( ERROR_LOCATION, "plasma potential must be above the "
		      "grounded backing" ) );

    // f(xi) = sheath field available - field the slab demands.
    // f(0) < 0 (no sheath field, full drop across the slab) and f grows
    // with xi, so the root is bracketed by (0, V_plasma/Te).
    double lo = 1.0e-12;
    double hi = (V_plasma/Te) * (1.0 - 1.0e-12);

    double f_lo = Te*sheath_field_factor( lo, bohm_factor )/lambda_D
	- eps_r*(V_plasma - lo*Te)/d;
    double f_hi = Te*sheath_field_factor( hi, bohm_factor )/lambda_D
	- eps_r*(V_plasma - hi*Te)/d;
    if( f_lo*f_hi > 0.0 )
	throw( Error( ERROR_LOCATION, "no dielectric surface potential brackets "
		      "the sheath/slab field balance" ) );

    for( int i = 0; i < 200; i++ ) {
	double mid = 0.5*(lo + hi);
	double f_mid = Te*sheath_field_factor( mid, bohm_factor )/lambda_D
	    - eps_r*(V_plasma - mid*Te)/d;
	if( f_lo*f_mid <= 0.0 ) {
	    hi = mid;
	    f_hi = f_mid;
	} else {
	    lo = mid;
	    f_lo = f_mid;
	}
	if( hi - lo < 1.0e-14*(1.0 + fabs( hi )) )
	    break;
    }

    return( 0.5*(lo + hi) );
}
