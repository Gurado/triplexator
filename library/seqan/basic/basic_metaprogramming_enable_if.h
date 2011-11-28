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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// "enable if" functionality.
// ==========================================================================

// SEQAN_NO_GENERATED_FORWARDS

#ifndef SEQAN_BASIC_BASIC_METAPROGRAMMING_ENABLE_IF_H_
#define SEQAN_BASIC_BASIC_METAPROGRAMMING_ENABLE_IF_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction EnableIf
// ----------------------------------------------------------------------------

//
// Example for enable-if function: 
//
//    template <typename TContainer>
//    typename EnableIf<
//        IsContainer<TContainer>,                  // 1st arg: enable-if condition
//        typename Size<TContainer>::Type >::Type   // 2nd arg: return type
//    length(TContainer & cont) 
//    {
//        SEQAN_CONCEPT_ASSERT((ContainerConcept<TContainer>));
//        return end(cont) - begin(cont);
//    }
    
template <bool b, typename T>
struct EnableIf;

template <bool b, typename T = void>
struct EnableIf
{
    typedef T Type;
};

template <typename T>
struct EnableIf<false, T> {};

// ----------------------------------------------------------------------------
// Metafunction DisableIf
// ----------------------------------------------------------------------------

template <bool b, typename T>
struct DisableIf;

template <bool b, typename T = void>
struct DisableIf
{
    typedef T Type;
};

template <typename T>
struct DisableIf<true, T> {};

// ----------------------------------------------------------------------------
// Metafunction EnableIf2
// ----------------------------------------------------------------------------

template <typename TCondition, typename T>
struct EnableIf2;

template <typename TCondition, typename T = void>
struct EnableIf2
{
    typedef T Type;
};

template <typename T>
struct EnableIf2<False, T> {};

// ----------------------------------------------------------------------------
// Metafunction DisableIf2
// ----------------------------------------------------------------------------

template <typename TCondition, typename T>
struct DisableIf2;

template <typename TCondition, typename T = void>
struct DisableIf2
{
    typedef T Type;
};

template <typename T>
struct DisableIf2<True, T> {};

}  // namespace seqan

// ============================================================================
// Macros
// ============================================================================

//
// Example for enable-if constructor: 
//
//    Rational(T const & n, SEQAN_CTOR_ENABLE_IF( IsInteger<T> )) :  // macro must be extra c'tor argument
//        num(n), den(1)
//    { 
//      (void)dummy;    // necessary to avoid unused warning
//    }
//

#define SEQAN_CTOR_ENABLE_IF(cond) typename EnableIf<cond::VALUE>::Type * dummy = 0

//
// Example for enable-if function (with macro): 
//
//    template <typename TContainer>
//    SEQAN_FUNC_ENABLE_IF(
//        IsContainer<TContainer>, 
//        typename Size<TContainer>::Type)
//    length(TContainer & cont) 
//    {
//        SEQAN_CONCEPT_ASSERT((ContainerConcept<TContainer>));
//        return end(cont) - begin(cont);
//    }
//

#define SEQAN_FUNC_ENABLE_IF(cond, retVal) typename EnableIf<cond::VALUE, retVal>::Type

// ============================================================================
// Functions
// ============================================================================

#endif  // #ifndef SEQAN_BASIC_BASIC_METAPROGRAMMING_ENABLE_IF_H_
