// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2010, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================

#ifndef SEQAN_HEADER_BASIC_H
#define SEQAN_HEADER_BASIC_H

//____________________________________________________________________________
// prerequisites

#include <seqan/platform.h>

//#include <cstring>
#ifdef PLATFORM_WINDOWS
#include <limits>   // limits include file exists only for g++ >= 3.0
#endif

#include <cstddef>  // size_t
#include <cstdio>   // FILE, basic_debug
#include <cstdlib>  // posix_memalign
#include <ctime>
#include <cmath>
#include <iterator>
#include <algorithm>
#include <cstring>  // memset, memcpy
#include <string>   // basic_profile
#ifdef PLATFORM_WINDOWS
#include <malloc.h> // _aligned_malloc
#endif  // PLATFORM_WINDOWS

#define SEQAN_NAMESPACE_MAIN seqan

//____________________________________________________________________________

#include <seqan/basic/basic_forwards.h>
#ifdef SEQAN_SWITCH_USE_FORWARDS
#include <seqan/basic/basic_generated_forwards.h>
#endif

// --------------------------------------------------------------------------
// Support Code
// --------------------------------------------------------------------------
//
// This code is accessed through macros.  No classes, tags, or functions are
// defined that are to be used directly outside these files.

#include <seqan/basic/basic_testing.h>
#include <seqan/basic/basic_debug.h>
#include <seqan/basic/basic_testing.h>  // new, better debug
#include <seqan/basic/basic_profile.h>
#include <seqan/basic/basic_parallelism.h>  // include after basic_testing.h!

// --------------------------------------------------------------------------
// Metaprogramming Support Code
// --------------------------------------------------------------------------
//
// This code lies the foundation for metaprogramming in SeqAn.

#include <seqan/basic/basic_metaprogramming.h>
#include <seqan/basic/basic_metaprogramming_enable_if.h>
#include <seqan/basic/basic_tag.h>
#include <seqan/basic/basic_type.h>

// --------------------------------------------------------------------------
// Concept Checking
// --------------------------------------------------------------------------

#include <seqan/basic/boost_preprocessor_subset.h>
#include <seqan/basic/concept_checking.h>
#include <seqan/basic/basic_concepts.h>

// ... and some utility code for computing logarithms and such.
#include <seqan/basic/basic_math.h>
#include <seqan/basic/basic_logvalue.h>

// --------------------------------------------------------------------------
// SeqAn Foundation Code
// --------------------------------------------------------------------------
//
// This code contains the "language extension" code in SeqAn: Construction,
// conversion, comparison and the transport (assign/set/move) functions.

#include <seqan/basic/basic_transport.h>
#include <seqan/basic/basic_converter.h>
#include <seqan/basic/basic_compare.h>

// --------------------------------------------------------------------------
// Iterators
// --------------------------------------------------------------------------

#include <seqan/basic/iterator_interface.h>
#include <seqan/basic/iterator_base.h>
#include <seqan/basic/iterator_adaptor.h>
#include <seqan/basic/iterator_position.h>
#include <seqan/basic/iterator_adapt_std.h>

// --------------------------------------------------------------------------
// Construction / Destruction Code
// --------------------------------------------------------------------------
//
// Contains functions wrapping construction and destuction of objects.

#include <seqan/basic/basic_construct_destruct.h>

// --------------------------------------------------------------------------
// Holders
// --------------------------------------------------------------------------

#include <seqan/basic/holder_base.h>
#include <seqan/basic/holder_tristate.h>
#include <seqan/basic/holder_simple.h>

#include <seqan/basic/basic_host.h>

// --------------------------------------------------------------------------
// Allocators
// --------------------------------------------------------------------------

#include <seqan/basic/allocator_interface.h>
#include <seqan/basic/allocator_simple.h>
#include <seqan/basic/allocator_singlepool.h>
#include <seqan/basic/allocator_multipool.h>
#include <seqan/basic/allocator_chunkpool.h>
                      
#include <seqan/basic/allocator_to_std.h>

// --------------------------------------------------------------------------
// Proxies
// --------------------------------------------------------------------------

#include <seqan/basic/proxy_base.h>
#include <seqan/basic/proxy_iterator.h>

#include <seqan/basic/basic_pointer.h>  // TODO(holtgrew): Belongs in sequence module.

// --------------------------------------------------------------------------
// Alphabets
// --------------------------------------------------------------------------

#include <seqan/basic/alphabet_storage.h>
#include <seqan/basic/alphabet_bio.h>
#include <seqan/basic/alphabet_math.h>
#include <seqan/basic/alphabet_simple.h>
#include <seqan/basic/alphabet_qualities.h>
#include <seqan/basic/alphabet_adapt_builtins.h>

#include <seqan/basic/basic_profchar.h>

#include <seqan/basic/basic_sse2.h>

#include <seqan/basic/basic_alphabet_simple_tabs.h>
#include <seqan/basic/basic_alphabet_simple.h>

#include <seqan/basic/alphabet_concept.h>
#include <seqan/basic/container_concept.h>

// --------------------------------------------------------------------------
// Aggregate Types: Pairs, Triples, Tuples.
// --------------------------------------------------------------------------

#include <seqan/basic/aggregate_concept.h>
#include <seqan/basic/pair_base.h>
#include <seqan/basic/pair_packed.h>
#include <seqan/basic/pair_bit_compressed.h>
#include <seqan/basic/triple_base.h>
#include <seqan/basic/triple_packed.h>
#include <seqan/basic/tuple_base.h>
#include <seqan/basic/tuple_bit_compressed.h>

// --------------------------------------------------------------------------
// Misc Code
// --------------------------------------------------------------------------

// TODO(holtgrew): Does this actually belong here?
#include <seqan/basic/basic_volatile_ptr.h>

// --------------------------------------------------------------------------
// Concept tests
// --------------------------------------------------------------------------

#include <seqan/basic/basic_test_concepts.h>

#endif //#ifndef SEQAN_HEADER_...
