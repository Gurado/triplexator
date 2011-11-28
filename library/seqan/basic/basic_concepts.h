// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2011, Knut Reinert, FU Berlin
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
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================
// Basic concept definitions, e.g. DefaultConstructible, Comparable, ...
// ==========================================================================

// SEQAN_NO_GENERATED_FORWARDS

#ifndef CORE_INCLUDE_SEQAN_BASIC_BASIC_CONCEPTS_H_
#define CORE_INCLUDE_SEQAN_BASIC_BASIC_CONCEPTS_H_

namespace seqan {

// ---------------------------------------------------------------------------
// ==> boost/concept_check.hpp <==
// ---------------------------------------------------------------------------

// (C) Copyright Jeremy Siek 2000.
// Copyright 2002 The Trustees of Indiana University.
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)


// ============================================================================
// Assignment Concepts
// ============================================================================

SEQAN_CONCEPT(DefaultConstructible,(TT))
{
	SEQAN_CONCEPT_USAGE(DefaultConstructible) {
		TT a;               // require default constructor
		ignoreUnusedVariableWarning(a);
	}
};

SEQAN_CONCEPT(Assignable,(TT))
{
	SEQAN_CONCEPT_USAGE(Assignable) {
#if !defined(_ITERATOR_) // back_insert_iterator broken for VC++ STL
		a = b;             // require assignment operator
#endif
		constConstraints(b);
	}
private:
	void constConstraints(const TT& x) {
#if !defined(_ITERATOR_) // back_insert_iterator broken for VC++ STL
		a = x;              // const required for argument to assignment
#else
		ignoreUnusedVariableWarning(x);
#endif
	}
private:
	TT a;
	TT b;
};

SEQAN_CONCEPT(CopyConstructible,(TT))
{
	SEQAN_CONCEPT_USAGE(CopyConstructible) {
		TT a(b);            // require copy constructor
		TT* ptr = &a;       // require address of operator
		constConstraints(a);
		ignoreUnusedVariableWarning(ptr);
	}
private:
	void constConstraints(const TT& a) {
		TT c(a);            // require const copy constructor
		const TT* ptr = &a; // require const address of operator
		ignoreUnusedVariableWarning(c);
		ignoreUnusedVariableWarning(ptr);
	}
	TT b;
};


// ============================================================================
// Relation Concepts
// ============================================================================

// The C++ standard requirements for many concepts talk about return
// types that must be "convertible to bool".  The problem with this
// requirement is that it leaves the door open for evil proxies that
// define things like operator|| with strange return types.  Two
// possible solutions are:
// 1) require the return type to be exactly bool
// 2) stay with convertible to bool, and also
//    specify stuff about all the logical operators.
// For now we just test for convertible to bool.
template <class TT>
void requireBooleanExpr(const TT& t) {
	bool x = t;
	ignoreUnusedVariableWarning(x);
}

SEQAN_CONCEPT(EqualityComparable,(TT))
{
	SEQAN_CONCEPT_USAGE(EqualityComparable) {
		requireBooleanExpr(a == b);
		requireBooleanExpr(a != b);
	}
private:
	TT a, b;
};

SEQAN_CONCEPT(LessThanComparable,(TT))
{
	SEQAN_CONCEPT_USAGE(LessThanComparable) {
		requireBooleanExpr(a < b);
	}
private:
	TT a, b;
};

// This is equivalent to SGI STL's LessThanComparable.
SEQAN_CONCEPT(Comparable,(TT))
{
	SEQAN_CONCEPT_USAGE(Comparable) {
		requireBooleanExpr(a < b);
		requireBooleanExpr(a > b);
		requireBooleanExpr(a <= b);
		requireBooleanExpr(a >= b);
	}
private:
	TT a, b;
};


// ============================================================================
// Forwards
// ============================================================================

SEQAN_CONCEPT(IntegerConcept, (T));
SEQAN_CONCEPT(SignedIntegerConcept, (T));
SEQAN_CONCEPT(UnsignedIntegerConcept, (T));

// ============================================================================
// Test fulfilled concepts
// ============================================================================

template <typename T>
struct Is< SignedIntegerConcept<T> >
{
	typedef
		// Explicitely unsigned.
		typename If< IsSameType<T, signed char>::VALUE,     True,
		typename If< IsSameType<T, short>::VALUE,           True,
		typename If< IsSameType<T, int>::VALUE,             True,
		typename If< IsSameType<T, long>::VALUE,            True,
		typename If< IsSameType<T, __int64>::VALUE,         True,
		False
		>::Type>::Type>::Type>::Type>::Type Type;
		enum { VALUE = Type::VALUE };
};

template <typename T>
struct Is< UnsignedIntegerConcept<T> >
{
	typedef
		// Explicitely unsigned.
		typename If< IsSameType<T, unsigned char>::VALUE,   True,
		typename If< IsSameType<T, unsigned short>::VALUE,  True,
		typename If< IsSameType<T, unsigned int>::VALUE,    True,
		typename If< IsSameType<T, unsigned long>::VALUE,   True,
		typename If< IsSameType<T, __uint64>::VALUE,        True,
		False
		>::Type>::Type>::Type>::Type>::Type Type;
		enum { VALUE = Type::VALUE };
};

template <typename T>
struct Is< IntegerConcept<T> >
{
	typedef
		// char is not necessarily equal to signed or unsigned char.
		typename If< IsSameType<T, char>::VALUE,                True,
		// Is T a signed integer?
		typename If< Is< SignedIntegerConcept<T> >::VALUE,      True,
		// Is T an unsigned integer?
		typename If< Is< UnsignedIntegerConcept<T> >::VALUE,    True,
		False
		>::Type>::Type>::Type Type;
		enum { VALUE = Type::VALUE };
};

template <typename T>
struct Is< SignedIntegerConcept<T const> > : Is< SignedIntegerConcept<T> > {};

template <typename T>
struct Is< UnsignedIntegerConcept<T const> > : Is< UnsignedIntegerConcept<T> > {};

template <typename T>
struct Is< IntegerConcept<T const> > : Is< IntegerConcept<T> > {};


/**
.Metafunction.IsSignedInteger:
..cat:Basic
..summary:Tests for a type to be of signed integral value.
..signature:IsSignedInteger<T>::Type
..param.T:Type that is tested.
..returns.param.Type:@Tag.Logical Values.True@, if $T$ is a signed integral type, @Tag.Logical Values.False@ otherwise.
...default:@Tag.Logical Values.False@
..status:deprecated
..include:seqan/basic.h
 */

/**
.Metafunction.IsUnsignedInteger:
..cat:Basic
..summary:Tests for a type to be of unsigned integral value.
..signature:IsUnsignedInteger<T>::Type
..param.T:Type that is tested.
..returns.param.Type:@Tag.Logical Values.True@, if $T$ is an unsigned integral type, @Tag.Logical Values.False@ otherwise.
...default:@Tag.Logical Values.False@
..status:deprecated
..include:seqan/basic.h
 */


/**
.Metafunction.IsInteger:
..cat:Basic
..summary:Tests for a type to be of integral value.
..signature:IsInteger<T>::Type
..param.T:Type that is tested.
..returns.param.Type:@Tag.Logical Values.True@, if $T$ is an ingegral type, @Tag.Logical Values.False@ otherwise.
...default:@Tag.Logical Values.False@
..status:deprecated
..include:seqan/basic.h
 */

// deprecation wrappers
template <typename T>
struct IsSignedInteger : Is< SignedIntegerConcept<T> > {};
template <typename T>
struct IsUnsignedInteger : Is< UnsignedIntegerConcept<T> > {};
template <typename T>
struct IsInteger : Is< IntegerConcept<T> > {};

template <typename T>
struct IsIntegral : IsInteger<T> {};

// ============================================================================
// Concepts for integers
// ============================================================================

// Generally there are two ways to check concepts:
//
// 1. Define expressions and rules a passing type must fulfill.
// 2. Let the concept check fail by default and only let it pass for the types
//    that definitely fulfill the concept

// We try variant 1, as it lets the user to define his/her own integer types
// without the need to specialize all kinds of Integer/SignedInteger/UnsignedInteger concepts.


SEQAN_CONCEPT(IntegerConcept, (TValue)) :
	Comparable<TValue>
{
    TValue val, val2;
    
    SEQAN_CONCEPT_USAGE(IntegerConcept)
    {
        val  = 0u;
        val2 = 1u;
		
        val2 = val + 1u;
        val2 = val + val;
        val2 += val;
        val2 += 1u;
        val2 = val++;
        val2 = ++val;

        val2 = val - val;
        val2 = val - 1u;
        val2 -= val;
        val2 -= 1u;
        val2 = val--;
        val2 = --val;
        
        val2 = val * val;
        val2 = val * 1u;
        val2 *= val;
        val2 *= 1u;

        val2 = val / val;
        val2 = val / 1u;
        val2 /= val;
        val2 /= 1u;
        
        val2 = val << val;
        val2 = val << 1;
        val2 <<= val;
        val2 <<= 1;

        val2 = val >> val;
        val2 = val >> 1;
        val2 >>= val;
        val2 >>= 1;
        
		SEQAN_STATIC_ASSERT_MSG(static_cast<TValue>(0u) < static_cast<TValue>(1u), "Integer has wrong order.");
    }
};

// an integer that must have a sign
SEQAN_CONCEPT(SignedIntegerConcept, (TValue)) :
    IntegerConcept<TValue>
{
    TValue val;
    
    SEQAN_CONCEPT_USAGE(SignedIntegerConcept)
    {
        val = -1;
        
		val = val - val;
		val = val + 1;
		val = val - 1;

        val = val / 2;
        
		SEQAN_STATIC_ASSERT_MSG(static_cast<TValue>(-1) < static_cast<TValue>(0), "Signed integer is either not signed or has wrong order.");
    }
};

// an integer that mustn't have a sign
SEQAN_CONCEPT(UnsignedIntegerConcept, (TValue)) :
    IntegerConcept<TValue>
{
    TValue val;
    
    SEQAN_CONCEPT_USAGE(UnsignedIntegerConcept)
    {
		SEQAN_STATIC_ASSERT_MSG(static_cast<TValue>(0) < static_cast<TValue>(-1), "Unsigned integer is either signed or has wrong order.");
    }
};


// This would be variant 2 (disabled for now):
/*
SEQAN_CONCEPT(IntegerConcept, (TValue))
{
    SEQAN_CONCEPT_USAGE(IntegerConcept)
    {
        x.error_type_must_be_an_integer_type();     // for the sake of readability we break the coding style rule here.
    }

private:
    T x;
};

template <> struct IntegerConcept<char> {};
template <> struct IntegerConcept<signed char> {};
template <> struct IntegerConcept<unsigned char> {};
template <> struct IntegerConcept<short> {};
template <> struct IntegerConcept<unsigned short> {};
template <> struct IntegerConcept<int> {};
template <> struct IntegerConcept<unsigned int> {};
template <> struct IntegerConcept<long> {};
template <> struct IntegerConcept<unsigned long> {};
//template <> struct IntegerConcept<__int64> {};
//template <> struct IntegerConcept<__uint64> {};

SEQAN_CONCEPT(SignedIntegerConcept, (TValue))
{
    SEQAN_CONCEPT_USAGE(SignedIntegerConcept)
    {
        x.error_type_must_be_a_signed_integer_type();
    }

private:
    T x;
};

template <> struct SignedIntegerConcept<char> {};
template <> struct SignedIntegerConcept<short> {};
template <> struct SignedIntegerConcept<int> {};
template <> struct SignedIntegerConcept<long> {};
//template <> struct SignedIntegerConcept<__int64> {};

SEQAN_CONCEPT(UnignedIntegerConcept, (TValue))
{
    SEQAN_CONCEPT_USAGE(UnignedIntegerConcept)
    {
        x.error_type_must_be_an_unsigned_integer_type();
    }

private:
    T x;
};

template <> struct UnignedIntegerConcept<unsigned char> {};
template <> struct UnignedIntegerConcept<unsigned short> {};
template <> struct UnignedIntegerConcept<unsigned int> {};
template <> struct UnignedIntegerConcept<unsigned long> {};
//template <> struct UnignedIntegerConcept<__uint64> {};
*/


}  // namespace seqan

#endif  // #ifndef CORE_INCLUDE_SEQAN_BASIC_BASIC_CONCEPTS_H_
