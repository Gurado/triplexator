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
// Author: Andres Gogol-Döring <andreas.doering@mdc-berlin.de>
// ==========================================================================
// Adaptions of builting types such as bool, int, but also "builtin-level"
// user defined types such as wchar_t, __int64, __uint64 to the alphabet
// concepts they are in.
// ==========================================================================

#ifndef SEQAN_BASIC_ALPHABET_ADAPT_BUILTINS_H_
#define SEQAN_BASIC_ALPHABET_ADAPT_BUILTINS_H_

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
// Metafunction IsSimple
// ----------------------------------------------------------------------------

template <> struct IsSimple_<bool> { typedef True Type; };
template <> struct IsSimple_<char> { typedef True Type; };

template <> struct IsSimple_<unsigned char> { typedef True Type; };
template <> struct IsSimple_<unsigned short> { typedef True Type; };
template <> struct IsSimple_<unsigned int> { typedef True Type; };
template <> struct IsSimple_<unsigned long> { typedef True Type; };

template <> struct IsSimple_<signed char> { typedef True Type; };
template <> struct IsSimple_<signed short> { typedef True Type; };
template <> struct IsSimple_<signed int> { typedef True Type; };
template <> struct IsSimple_<signed long> { typedef True Type; };

template <> struct IsSimple_<float> { typedef True Type; };
template <> struct IsSimple_<double> { typedef True Type; };
template <> struct IsSimple_<long double> { typedef True Type; };

// user defined types (re-specializations are allowed here)
template <> struct IsSimple<wchar_t> { typedef True Type; };
template <> struct IsSimple<__int64> { typedef True Type; };
template <> struct IsSimple<__uint64> { typedef True Type; };

// ----------------------------------------------------------------------------
// Metafunction BitsPerValue
// ----------------------------------------------------------------------------

template <> struct BitsPerValue<bool> { enum { VALUE = 1 }; };

/**
.Metafunction.IsCharType
..cat:Alphabets
..summary:Return whether the argument is $char$, $wchar_t$, $char const$, or $wchar_t const$.
..signature:IsCharType<T>::Type
..signature:IsCharType<T>::VALUE
..param.T:Type to check type of.
..remarks:This metafunction is used to enable and disable templated adaptions of arrays to sequences for builtin character types only.
..remarks:The return value is $True$/$true$ for $char$, $wchar_t$, $char const$, and $wchar_t const$.
..include:seqan/sequence.h
*/

// TODO(holtgrew): Write tests for this.

template <typename T>
struct IsCharType;

template <typename T>
struct IsCharType
{
    typedef False Type;
    enum { VALUE = 0 };
};

template <typename T>
struct IsCharType<T const>
    : IsCharType<T> {};

template <>
struct IsCharType<char>
{
    typedef True Type;
    enum { VALUE = 1 };
};

template <>
struct IsCharType<wchar_t>
{
    typedef True Type;
    enum { VALUE = 1 };
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function gapValueImpl()
// ----------------------------------------------------------------------------

// TODO(holtgrew): Why references? Profiles are faster then...

inline char const &
gapValueImpl(char *)
{
    SEQAN_CHECKPOINT;
    static char const _gap = '-';
    return _gap;
}

inline char const &
gapValueImpl(char const *)
{
    SEQAN_CHECKPOINT;
    static char const _gap = '-';
    return _gap;
}

// ----------------------------------------------------------------------------
// Function unknownValueImpl()
// ----------------------------------------------------------------------------

inline char const &
unknownValueImpl(char *)
{
    SEQAN_CHECKPOINT;
    static char const _unknown = 'N';
    return _unknown;
}

inline char const &
unknownValueImpl(char const *)
{
    SEQAN_CHECKPOINT;
    static char const _unknown = 'N';
    return _unknown;
}

// ----------------------------------------------------------------------------
// Function supremumValueImpl()
// ----------------------------------------------------------------------------

template <typename T>
inline T const &
supremumValueImpl(T *)
{
    // TODO(holtgrew): We probably do not want a default specialization.
    SEQAN_CHECKPOINT;
    return MaxValue<T>::VALUE;
}

inline long double const &
supremumValueImpl(long double *)
{
    SEQAN_CHECKPOINT;
#ifdef PLATFORM_WINDOWS
    static long double const _value = ::std::numeric_limits<long double>::infinity( );
#else
    static long double const _value = 1.7976931348623157e+308;
#endif
    return _value;
}

inline double const &
supremumValueImpl(double *)
{
    SEQAN_CHECKPOINT;
#ifdef PLATFORM_WINDOWS
    static double const _value = ::std::numeric_limits<double>::infinity( );
#else
    static double const _value = 1.7976931348623157e+308;
#endif
    return _value;
}
inline float const &
supremumValueImpl(float *)
{
    SEQAN_CHECKPOINT;
#ifdef PLATFORM_WINDOWS
    static float const _value = ::std::numeric_limits<float>::infinity( );
#else
    static float const _value = 3.40282347e+38F;
#endif
    return _value;
}

// ----------------------------------------------------------------------------
// Function infimumValueImpl()
// ----------------------------------------------------------------------------

template <typename T>
inline T const &
infimumValueImpl(T *)
{
    // TODO(holtgrew): We probably do not want a default specialization.
    SEQAN_CHECKPOINT;
    return MinValue<T>::VALUE;
}

inline float const &
infimumValueImpl(float *)
{
    SEQAN_CHECKPOINT;
#ifdef PLATFORM_WINDOWS
    static float const _value = -::std::numeric_limits<float>::infinity( );
#else
    static float const _value = -3.40282347e+38F;
#endif
    return _value;
}

inline double const &
infimumValueImpl(double *)
{
    SEQAN_CHECKPOINT;
#ifdef PLATFORM_WINDOWS
    static double const _value = -::std::numeric_limits<double>::infinity( );
#else
    static double const _value = -1.7976931348623157e+308;
#endif
    return _value;
}

inline long double const &
infimumValueImpl(long double *)
{
    SEQAN_CHECKPOINT;
#ifdef PLATFORM_WINDOWS
    static long double const _value = -::std::numeric_limits<long double>::infinity( );
#else
    static long double const _value = -1.7976931348623157e+308;
#endif
    return _value;
}

// TODO(holtgrew): The following functions were in basic_alphabet_traits.h but commented out because of infimumValue->minValue, supremumValue->maxValue change. Eventually, they might die out since infimumValue/supremumValue are only worth keeping for floating point numbers where there are infinity values.

/*
//////////////////////////////////////////////////////////////////////////////
// char 
//////////////////////////////////////////////////////////////////////////////

inline char const &
supremumValueImpl(char *)
{
SEQAN_CHECKPOINT
    static char const _value = (char) 127;
    return _value;
}
inline char const &
infimumValueImpl(char *)
{
SEQAN_CHECKPOINT
    static char const _value = (char) -128;
    return _value;
}

//////////////////////////////////////////////////////////////////////////////
// signed char 
//////////////////////////////////////////////////////////////////////////////

inline signed char const &
supremumValueImpl(signed char *)
{
SEQAN_CHECKPOINT
    static signed char const _value = 127;
    return _value;
}
inline signed char const &
infimumValueImpl(signed char *)
{
SEQAN_CHECKPOINT
    static signed char const _value = -128;
    return _value;
}

//////////////////////////////////////////////////////////////////////////////
// unsigned char 
//////////////////////////////////////////////////////////////////////////////

inline unsigned char const &
supremumValueImpl(unsigned char *)
{
SEQAN_CHECKPOINT
    static unsigned char const _value = 255;
    return _value;
}
inline unsigned char const &
infimumValueImpl(unsigned char *)
{
SEQAN_CHECKPOINT
    static unsigned char const _value = 0;
    return _value;
}

//////////////////////////////////////////////////////////////////////////////
// wchar_t
//////////////////////////////////////////////////////////////////////////////

inline wchar_t const &
supremumValueImpl(wchar_t *)
{
SEQAN_CHECKPOINT
    static wchar_t const _value = 1UL << (BitsPerValue<wchar_t>::VALUE) - 1;
    return _value;
}
inline wchar_t const &
infimumValueImpl(wchar_t *)
{
SEQAN_CHECKPOINT
    static wchar_t const _value = 0;
    return _value;
}

//////////////////////////////////////////////////////////////////////////////
// signed short 
//////////////////////////////////////////////////////////////////////////////

inline signed short const &
supremumValueImpl(signed short *)
{
SEQAN_CHECKPOINT
    static signed short const _value = (((1 << (BitsPerValue<signed short>::VALUE - 2)) - 1) << 1) + 1;
    return _value;
}
inline signed short const &
infimumValueImpl(signed short *dummy)
{
SEQAN_CHECKPOINT
    static signed short const _value = -supremumValueImpl(dummy) - 1;
    return _value;
}

//////////////////////////////////////////////////////////////////////////////
// unsigned short 
//////////////////////////////////////////////////////////////////////////////

inline unsigned short const &
supremumValueImpl(unsigned short *)
{
SEQAN_CHECKPOINT
    static unsigned short const _value = (((1 << (BitsPerValue<unsigned short>::VALUE - 1)) - 1) << 1) + 1;
    return _value;
}
inline unsigned short const &
infimumValueImpl(unsigned short *)
{
SEQAN_CHECKPOINT
    static unsigned short const _value = 0;
    return _value;
}

//////////////////////////////////////////////////////////////////////////////
// signed int 
//////////////////////////////////////////////////////////////////////////////

inline signed int const &
supremumValueImpl(signed int *)
{
SEQAN_CHECKPOINT
    static signed int const _value = (((1 << (BitsPerValue<signed int>::VALUE - 2)) - 1) << 1) + 1;
    return _value;
}

inline signed int const &
infimumValueImpl(signed int *dummy)
{
SEQAN_CHECKPOINT
    static signed int const _value = -supremumValueImpl(dummy) - 1;
    return _value;
}

//////////////////////////////////////////////////////////////////////////////
// unsigned int 
//////////////////////////////////////////////////////////////////////////////

inline unsigned int const &
supremumValueImpl(unsigned int *)
{
SEQAN_CHECKPOINT
    static unsigned int const _value = ~0ul;
    return _value;
}

inline unsigned int const &
infimumValueImpl(unsigned int *)
{
SEQAN_CHECKPOINT
    static unsigned int const _value = 0;
    return _value;
}

//////////////////////////////////////////////////////////////////////////////
// signed long
//////////////////////////////////////////////////////////////////////////////

inline signed long const &
supremumValueImpl(signed long *)
{
SEQAN_CHECKPOINT
    static signed long const _value = (((1 << (BitsPerValue<signed long>::VALUE - 2)) - 1) << 1) + 1;
    return _value;
}

inline signed long const &
infimumValueImpl(signed long *dummy)
{
SEQAN_CHECKPOINT
    static signed long const _value = -supremumValueImpl(dummy) - 1;
    return _value;
}

//////////////////////////////////////////////////////////////////////////////
// unsigned long
//////////////////////////////////////////////////////////////////////////////

inline unsigned long const &
supremumValueImpl(unsigned long *)
{
SEQAN_CHECKPOINT
    static unsigned long const _value = ~0ul;
    return _value;
}

inline unsigned long const &
infimumValueImpl(unsigned long *)
{
SEQAN_CHECKPOINT
    static unsigned long const _value = 0;
    return _value;
}

//////////////////////////////////////////////////////////////////////////////
// signed 64bit int (cannot use long long <- no ISO C++)
//////////////////////////////////////////////////////////////////////////////

inline __int64 const &
supremumValueImpl(__int64 *)
{
SEQAN_CHECKPOINT
    static __int64 const _value = ((((__int64)1 << (BitsPerValue<__int64>::VALUE - 2)) - 1) << 1) + 1;
    return _value;
}

inline __int64 const &
infimumValueImpl(__int64 *dummy)
{
SEQAN_CHECKPOINT
    static __int64 const _value = -supremumValueImpl(dummy) - 1;
    return _value;
}
*/

}  // namespace seqan

#endif  // #ifndef SEQAN_BASIC_ALPHABET_ADAPT_BUILTINS_H_
