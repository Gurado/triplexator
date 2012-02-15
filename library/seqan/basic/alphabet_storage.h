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
// Author: Andreas Gogol-Doering <andreas.doering@mdc-berlin.de>
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Storage-related alphabet interface part.  This means both
// construction type (simple, non-simple) and storage size.
// ==========================================================================

#ifndef SEQAN_BASIC_ALPHABET_STORAGE_H_
#define SEQAN_BASIC_ALPHABET_STORAGE_H_

#include <float.h>

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// TODO(holtgrew): Remove Ascii and Unicode alias. Also see #849.
typedef char Ascii;
//typedef unsigned char Byte;  // TODO(holtgrew): Disabling, remove together with Ascii and Unicode with #849
typedef wchar_t Unicode;

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction BitsPerValue
// ----------------------------------------------------------------------------

/**
.Metafunction.BitsPerValue:
..cat:Basic
..summary:Number of bits needed to store a value.
..signature:BitsPerValue<T>::VALUE
..param.T:A class.
..returns.param.VALUE:Number of bits needed to store $T$.
...default:$sizeof<T> * 8$
..see:Metafunction.ValueSize
..include:seqan/basic.h
*/

template <typename TValue>
struct BitsPerValue
{
    enum { VALUE = sizeof(TValue) * 8 };
};

template <typename TValue>
struct BitsPerValue<TValue const> : public BitsPerValue<TValue>
{};

// ----------------------------------------------------------------------------
// Metafunction ValueSize
// ----------------------------------------------------------------------------

/**
.Metafunction.ValueSize:
..cat:Basic
..summary:Number of different values a value type object can have.
..signature:ValueSize<T>::VALUE
..param.T:A class.
..returns.param.VALUE:Value size of $T$.
..remarks
...text:This function is only defined for integral types like $unsigned int$, $double$ or @Spec.Dna@.
..see:Metafunction.Value
..include:seqan/basic.h
*/

template <typename T>
struct ValueSize
{
    enum { VALUE = 1 << BitsPerValue<T>::VALUE };
};
template <typename TValue>
struct ValueSize<TValue const>
        : public ValueSize<TValue>
{};

// TODO(holtgrew): What is this used for?
template <typename TValue> 
struct InternalValueSize_
        : public ValueSize<TValue>
{};

// ----------------------------------------------------------------------------
// Metafunction BytesPerValue
// ----------------------------------------------------------------------------

/**
.Metafunction.BytesPerValue:
..cat:Basic
..summary:Number of bytes needed to store a value.
..signature:BytesPerValue<T>::VALUE
..param.T:A class.
..returns.param.VALUE:Number of bytes needed to store $T$.
...default:$BitsPerValue / 8$, rounded up. For built-in types, this is the same as $sizeof(T)$.
..see:Metafunction.ValueSize
..see:Metafunction.BitsPerValue
..include:seqan/basic.h
*/

template <typename TValue>
struct BytesPerValue
{
    enum { VALUE = (BitsPerValue<TValue>::VALUE + 7) / 8 };
};

// ----------------------------------------------------------------------------
// Metafunction IntegralForValue
// ----------------------------------------------------------------------------

/**
.Metafunction.IntegralForValue:
..cat:Basic
..summary:Returns an itegral type that provides sufficient space to store a value.
..signature:IntegralForValue<T>::Type
..param.T:A class.
..returns.param.Type:An integral type that can store $T$ values.
..remarks:The type is the smallest unsigned integral type that has a size of at least @Metafunction.BytesPerValue@ bytes.
...tableheader:bytes|integral type
...table:1|$unsigned char$
...table:2|$unsigned short$
...table:3|$unsigned int$
...table:4|$unsigned int$
...table:5 and above|$__int64$
..remarks:Note that the returned integral type cannot store $T$ values, if $T$ takes more than 8 bytes, 
    since there exists no integral type that provides sufficient space to store types of this size.
..see:Metafunction.ValueSize
..see:Metafunction.BitsPerValue
..see:Metafunction.BytesPerValue
..include:seqan/basic.h
*/

template <int SIZE>
struct IntegralForValueImpl_
{
    typedef __int64 Type;
};

template <>
struct IntegralForValueImpl_<1>
{
    typedef unsigned char Type;
};

template <>
struct IntegralForValueImpl_<2>
{
    typedef unsigned short Type;
};

template <>
struct IntegralForValueImpl_<3>
{
    typedef unsigned int Type;
};

template <>
struct IntegralForValueImpl_<4>
{
    typedef unsigned int Type;
};

template <typename TValue>
struct IntegralForValue
        : IntegralForValueImpl_<BytesPerValue<TValue>::VALUE>
{};

// ============================================================================
// Functions
// ============================================================================

}  // namespace seqan

#endif  // #ifndef SEQAN_BASIC_ALPHABET_STORAGE_H_
