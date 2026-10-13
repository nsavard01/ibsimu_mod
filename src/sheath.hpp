/*! \file sheath.hpp
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

#ifndef SHEATH_HPP
#define SHEATH_HPP 1


/* Planar, collisionless sheath with cold ions entering at the Bohm speed
 * and Boltzmann electrons. With
 *
 *     xi   = (phi_s - phi) / Te      drop below the sheath edge, in Te
 *     zeta = x / lambda_Ds           distance in sheath-edge Debye lengths
 *     n_i  = n_s / sqrt(1 + 2 xi)    flux + energy conservation
 *     n_e  = n_s exp(-xi)            Boltzmann
 *
 * Poisson reduces to
 *
 *     xi'' = 1/sqrt(1 + 2 xi) - exp(-xi)
 *
 * whose first integral is (xi')^2 = 2 G(xi), so position follows from a
 * single quadrature -- see sheath_distance().
 *
 * These routines were checked against an independent shooting solution of
 * the same ODE (integrate backwards from the wall, bisecting on the wall
 * field until the solution decays into the flat quasineutral region). That
 * reference uses neither the first integral nor the quadrature, and its
 * wall field agrees with sqrt(2 G(xi_wall)) to 12 significant figures.
 */


/*! \brief Bohm speed sqrt(q e Te / M) in m/s.
 *
 *  \a Te in volts, \a mass_amu in atomic mass units, \a charge_state in
 *  units of the elementary charge. Each species has its own Bohm speed,
 *  so a multi-species plasma must call this per species.
 */
double sheath_bohm_speed( double Te, double mass_amu, double charge_state = 1.0 );


/*! \brief Sheath-edge number density J / (q e v_B) in 1/m3.
 *
 *  The density a current density \a J implies at the sheath edge if the
 *  species enters there at its Bohm speed. This is what makes the density
 *  a DERIVED quantity: given J and Te it is not free.
 */
double sheath_edge_density( double J, double Te, double mass_amu,
			    double charge_state = 1.0,
			    double bohm_factor = 0.5 );


/*! \brief Debye length sqrt(eps0 Te / (e n)) in m, with \a Te in volts. */
double sheath_debye_length( double Te, double n );


/*! \brief G(xi) = sqrt(1 + 2 xi) + exp(-xi) - 2, the first integral.
 *
 *  (xi')^2 = 2 G(xi). G vanishes as xi^3/3 at the sheath edge, which is
 *  why the edge is approached asymptotically and is a convention rather
 *  than a point.
 */
/*  \a bohm_factor is the ion kinetic energy AT THE SHEATH EDGE in units
 *  of Te. The Bohm criterion makes 0.5 the marginal (and default) value;
 *  entering supersonically means f > 0.5, which the presheath permits.
 *
 *  The entry speed enters G: with n_i = n_s/sqrt(1 + xi/f),
 *
 *      G_f(xi) = 2 f sqrt(1 + xi/f) - 2 f - 1 + exp(-xi)
 *
 *  which reduces to sqrt(1+2xi) + exp(-xi) - 2 exactly at f = 0.5. Near
 *  the edge G_f ~ xi^2 (1/2 - 1/(4f)): the quadratic term cancels only at
 *  f = 0.5 (leaving xi^3/3), and for f < 0.5 it is NEGATIVE -- the Bohm
 *  criterion violated, with no monotonic sheath solution -- so f < 0.5
 *  throws rather than being silently clamped.
 */
double sheath_G( double xi, double bohm_factor = 0.5 );


/*! \brief Normalised field xi' = sqrt(2 G(xi)) at drop \a xi.
 *
 *  The physical field is Te * sheath_field_factor(xi) / lambda_Ds in V/m.
 */
double sheath_field_factor( double xi, double bohm_factor = 0.5 );


/*! \brief Distance from drop \a xi_from to drop \a xi_to, in Debye lengths.
 *
 *  The quadrature
 *
 *      zeta = integral from xi_from to xi_to of dt / sqrt(2 G(t))
 *
 *  evaluated by adaptive Simpson. Only DIFFERENCES of position are
 *  physical: any common lower limit cancels, so this takes two drops and
 *  needs no sheath-edge convention.
 *
 *  \a xi_from must be strictly positive -- the integrand diverges as
 *  xi^(-3/2) at the sheath edge, so the distance measured from xi = 0 is
 *  infinite. Throws if xi_from <= 0 or xi_to < xi_from.
 */
double sheath_distance( double xi_from, double xi_to,
			double bohm_factor = 0.5 );


/*! \brief Drop to a grounded conductor, in units of Te: simply V_plasma/Te. */
double sheath_wall_xi_metal( double V_plasma, double Te );


/*! \brief Drop to the surface of a dielectric slab backing onto ground.
 *
 *  A charge-free slab of thickness \a d and relative permittivity
 *  \a eps_r in front of a grounded conductor carries a linear potential
 *  fall, and D is continuous at its surface, so the sheath field must
 *  supply what the slab demands:
 *
 *      Te sqrt(2 G(xi)) / lambda_D  =  eps_r (V_plasma - xi Te) / d
 *
 *  Solved for xi by bisection. Note the surface generally sits DEEPER
 *  than the floating potential would put it: current balance alone gives
 *  about 2.84 Te for hydrogen, but that field cannot always carry the
 *  drop the slab needs.
 *
 *  \a Te and \a V_plasma in volts, \a lambda_D and \a d in metres.
 */
double sheath_wall_xi_dielectric( double V_plasma, double Te, double lambda_D,
				  double eps_r, double d,
				  double bohm_factor = 0.5 );


#endif
