/*! \file parallel_util.hpp
 *  \brief Minimal pthread-based parallel-for utility for matrix solvers
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

#ifndef PARALLEL_UTIL_HPP
#define PARALLEL_UTIL_HPP 1


/* Design note
 * -----------
 * IBSimu already ships a thread abstraction, Scheduler<Solv,Prob,Err>
 * (scheduler.hpp), used by e.g. ParticleIterator. That class is a
 * producer-consumer work queue intended for dynamic streams of
 * independent, heterogeneous, long-running jobs (one job per particle,
 * jobs of very different individual cost, problems possibly added
 * while consumers are running). Spinning it up involves creating a new
 * Scheduler, a vector of per-thread solver objects and a fresh set of
 * pthreads for every batch of problems.
 *
 * That is the wrong cost model for what is needed inside BiCGSTAB and
 * a red-black SOR preconditioner: a fixed-size array of N independent,
 * equal-cost updates (one sparse row, or one mesh node) that has to be
 * repeated every single linear-solver iteration (often hundreds of
 * times per nonlinear Newton step, itself repeated many times per
 * Vlasov round). Re-creating a Scheduler and its pthreads on every
 * BiCGSTAB iteration or every red/black colour sweep would dwarf the
 * actual floating point work for anything but very large problems.
 *
 * What is needed instead is the classic "parallel for": split a known
 * range [0,n) into ibsimu.get_thread_count() contiguous chunks, run a
 * worker per chunk, join, return -- i.e. a single barrier per call,
 * implemented directly with pthread_create()/pthread_join() and no
 * persistent thread pool. This keeps the dependency surface identical
 * to what the rest of the library already uses (pthreads only, no
 * OpenMP, no other threading library) while matching the actual access
 * pattern of sparse matrix-vector products, dot products and red-black
 * sweeps.
 */


#include <pthread.h>
#include <vector>
#include <cstdint>


/*! \brief Internal data passed to each worker thread by parallel_for().
 *
 *  \tparam Func is any callable (function pointer, functor or lambda)
 *  with signature void(uint32_t a, uint32_t b), called once per chunk
 *  with the half-open index range [a,b) that the calling thread should
 *  process.
 */
template <class Func>
struct ParallelForChunk {
    const Func *func;
    uint32_t a;
    uint32_t b;
};


template <class Func>
void *parallel_for_entry( void *data )
{
    ParallelForChunk<Func> *chunk = (ParallelForChunk<Func> *)data;
    (*chunk->func)( chunk->a, chunk->b );
    return( NULL );
}


/*! \brief Run \a func( a, b ) once per chunk, splitting [0,n) into
 *  \a nthreads contiguous, roughly equal chunks and running each chunk
 *  in its own pthread. Blocks until all chunks have completed.
 *
 *  If \a n is zero or \a nthreads is 1, \a func is called once,
 *  in-line, for the whole range, without spawning any thread, which
 *  keeps the small-problem case (low resolution geometries, where
 *  thread creation overhead would dominate) as fast as the original
 *  serial code.
 */
template <class Func>
void parallel_for( uint32_t n, const Func &func, uint32_t nthreads )
{
    if( nthreads < 1 )
        nthreads = 1;
    if( n == 0 )
        return;
    if( nthreads > n )
        nthreads = n;

    if( nthreads == 1 || n < 10000) {
        func( 0, n );
        return;
    }

    std::vector<pthread_t> threads( nthreads );
    std::vector<ParallelForChunk<Func> > chunks( nthreads );

    uint32_t base  = n / nthreads;
    uint32_t extra = n % nthreads;
    uint32_t pos = 0;
    for( uint32_t t = 0; t < nthreads; t++ ) {
        uint32_t csize = base + (t < extra ? 1 : 0);
        chunks[t].func = &func;
        chunks[t].a    = pos;
        chunks[t].b    = pos+csize;
        pos += csize;
    }

    // Launch threads 1..nthreads-1; run chunk 0 on the calling thread
    // to avoid paying for one extra thread creation/join.
    for( uint32_t t = 1; t < nthreads; t++ )
        pthread_create( &threads[t], NULL, parallel_for_entry<Func>, (void *)&chunks[t] );

    func( chunks[0].a, chunks[0].b );

    for( uint32_t t = 1; t < nthreads; t++ )
        pthread_join( threads[t], NULL );
}


#endif
