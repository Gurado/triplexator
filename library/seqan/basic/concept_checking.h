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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// A minimal subset of the Boost Concept Checking Library.  A lot of the code
// in the BCCL deals with support of non-conforming compilers and we cut this
// away.  The code here has been adjusted to work with the compilers supported
// by SeqAn and be as simple as possible while still creating useful compiler
// errors.
// ==========================================================================

// SEQAN_NO_GENERATED_FORWARDS

#ifndef CORE_INCLUDE_SEQAN_BASIC_CONCEPT_CHECKING_H_
#define CORE_INCLUDE_SEQAN_BASIC_CONCEPT_CHECKING_H_

namespace seqan {

// ---------------------------------------------------------------------------
// ==> boost/static_assert.hpp <==
// ---------------------------------------------------------------------------

//  (C) Copyright John Maddock 2000.
//  Use, modification and distribution are subject to the 
//  Boost Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  See http://www.boost.org/libs/static_assert for documentation.

#ifdef SEQAN_CXX11_STANDARD
#  define SEQAN_STATIC_ASSERT_MSG( B, Msg ) static_assert(B, Msg)
#else
#  define SEQAN_STATIC_ASSERT_MSG( B, Msg ) SEQAN_STATIC_ASSERT( B )
#endif

//
// If the compiler issues warnings about old C style casts,
// then enable this:
//
//#if defined(__GNUC__) && ((__GNUC__ > 3) || ((__GNUC__ == 3) && (__GNUC_MINOR__ >= 4)))
//#  define BOOST_STATIC_ASSERT_BOOL_CAST( x ) ((x) == 0 ? false : true)
//#else
#  define SEQAN_STATIC_ASSERT_BOOL_CAST(x) (bool)(x)
//#endif

#ifdef SEQAN_CXX11_STANDARD
#  define SEQAN_STATIC_ASSERT( B ) static_assert(B, #B)
#else

// HP aCC cannot deal with missing names for template value parameters
template <bool x> struct STATIC_ASSERTION_FAILURE;

template <> struct STATIC_ASSERTION_FAILURE<true> { enum { value = 1 }; };

// HP aCC cannot deal with missing names for template value parameters
template<int x> struct static_assert_test{};

//
// Implicit instantiation requires that all member declarations be
// instantiated, but that the definitions are *not* instantiated.
//
// It's not particularly clear how this applies to enum's or typedefs;
// both are described as declarations [7.1.3] and [7.2] in the standard,
// however some compilers use "delayed evaluation" of one or more of
// these when implicitly instantiating templates.  We use typedef declarations
// by default, but try defining SEQAN_USE_ENUM_STATIC_ASSERT if the enum
// version gets better results from your compiler...
//
// Implementation:
// Both of these versions rely on sizeof(incomplete_type) generating an error
// message containing the name of the incomplete type.  We use
// "STATIC_ASSERTION_FAILURE" as the type name here to generate
// an eye catching error message.  The result of the sizeof expression is either
// used as an enum initialiser, or as a template argument depending which version
// is in use...
// Note that the argument to the assert is explicitly cast to bool using old-
// style casts: too many compilers currently have problems with static_cast
// when used inside integral constant expressions.
//
//#if !defined(SEQAN_BUGGY_INTEGRAL_CONSTANT_EXPRESSIONS)
/*
#if defined(SEQAN_MSVC) && (SEQAN_MSVC < 1300)
// __LINE__ macro broken when -ZI is used see Q199057
// fortunately MSVC ignores duplicate typedef's.
#define SEQAN_STATIC_ASSERT( B ) \
   typedef static_assert_test<\
      sizeof(STATIC_ASSERTION_FAILURE< (bool)( B ) >)\
      > seqan_static_assert_typedef_
#elif defined(SEQAN_MSVC)
*/
#if defined(PLATFORM_WINDOWS_VS)
#define SEQAN_STATIC_ASSERT( B ) \
   typedef static_assert_test<\
      sizeof(STATIC_ASSERTION_FAILURE< SEQAN_STATIC_ASSERT_BOOL_CAST ( B ) >)>\
         SEQAN_JOIN(seqan_static_assert_typedef_, __COUNTER__)
/*
#elif defined(SEQAN_INTEL_CXX_VERSION) || defined(SEQAN_SA_GCC_WORKAROUND)
// agurt 15/sep/02: a special care is needed to force Intel C++ issue an error 
// instead of warning in case of failure
# define SEQAN_STATIC_ASSERT( B ) \
    typedef char SEQAN_JOIN(seqan_static_assert_typedef_, __LINE__) \
        [ STATIC_ASSERTION_FAILURE< SEQAN_STATIC_ASSERT_BOOL_CAST( B ) >::value ]
#elif defined(__sgi)
// special version for SGI MIPSpro compiler
#define SEQAN_STATIC_ASSERT( B ) \
   SEQAN_STATIC_CONSTANT(bool, \
     SEQAN_JOIN(boost_static_assert_test_, __LINE__) = ( B )); \
   typedef static_assert_test<\
     sizeof(STATIC_ASSERTION_FAILURE< \
       SEQAN_JOIN(boost_static_assert_test_, __LINE__) >)>\
         SEQAN_JOIN(seqan_static_assert_typedef_, __LINE__)
#elif SEQAN_WORKAROUND(__MWERKS__, <= 0x3003)
// special version for CodeWarrior <= 8.x
#define SEQAN_STATIC_ASSERT( B ) \
   SEQAN_STATIC_CONSTANT(int, \
     SEQAN_JOIN(boost_static_assert_test_, __LINE__) = \
       sizeof(STATIC_ASSERTION_FAILURE< SEQAN_STATIC_ASSERT_BOOL_CAST( B ) >) )
*/
#else
// generic version
#define SEQAN_STATIC_ASSERT( B ) \
   typedef static_assert_test<\
      sizeof(STATIC_ASSERTION_FAILURE< SEQAN_STATIC_ASSERT_BOOL_CAST( B ) >)>\
         SEQAN_JOIN(seqan_static_assert_typedef_, __LINE__)
#endif
/*
#else
// alternative enum based implementation:
#define SEQAN_STATIC_ASSERT( B ) \
   enum { SEQAN_JOIN(boost_static_assert_enum_, __LINE__) \
      = sizeof(STATIC_ASSERTION_FAILURE< (bool)( B ) >) }
#endif
*/
#endif

// ---------------------------------------------------------------------------
// ==> boost/parameter/aux_/paranthesized_type.hpp <==
// ---------------------------------------------------------------------------

// Copyright David Abrahams 2006. Distributed under the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt).

template <class UnaryFunctionPointer>
struct unaryfunptr_arg_type;

template <class Arg>
struct unaryfunptr_arg_type<void(*)(Arg)>
{
    typedef Arg type; 
};

template <>
struct unaryfunptr_arg_type<void(*)(void)>
{
    typedef void type;
};

// ---------------------------------------------------------------------------
// ==> boost/concept_check/general.hpp <==
// ---------------------------------------------------------------------------

// Copyright David Abrahams 2006. Distributed under the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt).

namespace concept_checking
{
template <void(*)()> struct instantiate {};
}

template <class ModelFn> struct concept_check_;

template <class Model>
void concept_check_failed()
{
    ((Model*)0)->~Model();
}

template <class Model>
struct concept_check
{
    concept_checking::instantiate<concept_check_failed<Model> > x;
    enum { instantiate = 1 };
};

template <class Model>
struct concept_check_<void(*)(Model)>
        : concept_check<Model>
{};

#  define SEQAN_CONCEPT_ASSERT_FN( ModelFnPtr )             \
    typedef ::seqan::detail::instantiate<          \
    &::seqan::requirement_<ModelFnPtr>::failed>    \
      SEQAN_PP_CAT(seqan_concept_check,__LINE__)

// ---------------------------------------------------------------------------
// ==> boost/concept/assert.hpp <==
// ---------------------------------------------------------------------------

// Copyright David Abrahams 2006. Distributed under the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt).

// Usage, in class or function context:
//     SEQAN_CONCEPT_ASSERT((UnaryFunctionConcept<F,bool,int>));
# define SEQAN_CONCEPT_ASSERT(ModelInParens) \
    SEQAN_CONCEPT_ASSERT_FN(void(*)ModelInParens)

// usage.hpp

template <class Model>
struct usage_requirements
{
    ~usage_requirements() { ((Model*)0)->~Model(); }
};

#define SEQAN_CONCEPT_USAGE(model)                                      \
    SEQAN_CONCEPT_ASSERT((seqan::usage_requirements<model>));           \
    ~model()

// ---------------------------------------------------------------------------
// ==> boost/concept/detail/has_constraints.hpp <==
// ---------------------------------------------------------------------------

// Copyright David Abrahams 2006. Distributed under the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt).

namespace detail {
  typedef char yes;
  typedef char (&no)[2];

  template <class Model, void (Model::*)()>
  struct wrap_constraints {};
    
  template <class Model>
  inline yes has_constraints_(Model*, wrap_constraints<Model,&Model::constraints>* = 0);
  inline no has_constraints_(...);
}

// This would be called "detail::has_constraints," but it has a strong
// tendency to show up in error messages.
template <class Model>
struct not_satisfied
{
    enum {value = sizeof( detail::has_constraints_((Model*)0) ) == sizeof(detail::yes) };
    typedef typename Eval<value>::Type Type;
};

// ---------------------------------------------------------------------------
// ==> boost/concept/detail/concept_def.hpp <==
// ---------------------------------------------------------------------------

// Copyright David Abrahams 2006. Distributed under the Boost
// Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
    
# define SEQAN_CONCEPT(name, params)                                            \
    template < SEQAN_PP_SEQ_FOR_EACH_I(SEQAN_CONCEPT_typename,~,params) >       \
    struct name                                                                
    
// Helper for SEQAN_CONCEPT, above.
# define SEQAN_CONCEPT_typename(r, ignored, index, t) \
    SEQAN_PP_COMMA_IF(index) typename t
    
// ---------------------------------------------------------------------------
// ==> boost/concept_check/detail/general.hpp <==
// ---------------------------------------------------------------------------

// Copyright David Abrahams 2006. Distributed under the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt).

template <class ModelFn>
struct requirement_;

namespace detail
{
  template <void(*)()> struct instantiate {};
}

template <class Model>
struct requirement
{
    static void failed() { ((Model*)0)->~Model(); }
};

struct failed {};

template <class Model>
struct requirement<failed ************ Model::************>
{
    static void failed() { ((Model*)0)->~Model(); }
};

template <class Model>
struct constraint
{
    static void failed() { ((Model*)0)->constraints(); }
};
  
template <class Model>
struct requirement_<void(*)(Model)>
        : If<not_satisfied<Model>::Type::VALUE, /* should be called "has_constraints", see above */
             constraint<Model>,
             requirement<failed ************ Model::************>
             >::Type
{};

#  define SEQAN_CONCEPT_ASSERT_FN( ModelFnPtr )             \
    typedef ::seqan::detail::instantiate<          \
    &::seqan::requirement_<ModelFnPtr>::failed>    \
      SEQAN_PP_CAT(seqan_concept_check,__LINE__)

// ---------------------------------------------------------------------------
// ==> boost/concept_check/detail/requires.hpp <==
// ---------------------------------------------------------------------------

// Copyright David Abrahams 2006. Distributed under the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt).

// Template for use in handwritten assertions
template <class Model, class More>
struct requires_ : More
{
    SEQAN_CONCEPT_ASSERT((Model));
};

// Template for use by macros, where models must be wrapped in parens.
// This isn't in namespace detail to keep extra cruft out of resulting
// error messages.

template <class ModelFn>
struct _requires_
{
    enum { value = 0 };
    SEQAN_CONCEPT_ASSERT_FN(ModelFn);
};

template <int check, class Result>
struct Requires_ : unaryfunptr_arg_type<Result>
{};

#  define SEQAN_CONCEPT_REQUIRES_(r,data,t) + (::seqan::_requires_<void(*)t>::value)

#if defined(NDEBUG)

# define SEQAN_CONCEPT_REQUIRES(models, result)                      \
    typename unaryfunptr_arg_type<void(*)result>::type

#else  // #if defined(NDEBUG)

# define SEQAN_CONCEPT_REQUIRES(models, result)                                        \
    typename ::seqan::Requires_<                                                       \
      (0 SEQAN_PP_SEQ_FOR_EACH(SEQAN_CONCEPT_REQUIRES_, ~, models)),                   \
      void(*)result                                                                 \
    >::type

#endif  // #if defined(NDEBUG)

// ---------------------------------------------------------------------------
// ==> boost/concept_check.hpp <==
// ---------------------------------------------------------------------------

//
// (C) Copyright Jeremy Siek 2000.
// Copyright 2002 The Trustees of Indiana University.
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)    //
// Backward compatibility
//

template <class Model>
inline void functionRequires(Model* = 0)
{
    SEQAN_CONCEPT_ASSERT((Model));
}
template <class T> inline void ignoreUnusedVariableWarning(T const&) {}


// ============================================================================
// Functions
// ============================================================================

template <typename T>
void sameType(T, T) { }


// ============================================================================
// Metafunctions
// ============================================================================

// test whether a concept is fulfilled (without concept checking)
template <typename T>
struct Is : False {};

    
}  // namespace seqan

#endif  // #ifndef CORE_INCLUDE_SEQAN_BASIC_CONCEPT_CHECKING_H_
