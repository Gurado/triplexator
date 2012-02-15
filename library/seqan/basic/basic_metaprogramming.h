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
// Code supporting metaprogramming.
// ==========================================================================

#ifndef SEQAN_BASIC_BASIC_METAPROGRAMMING_H_
#define SEQAN_BASIC_BASIC_METAPROGRAMMING_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/**
.Tag.Logical Values:
..cat:Metaprogramming
..summary:Tag that represents true and false.
..tag.True:The logical value "true".
..tag.False:The logical value "false".
..remarks:These tags also function as Metafunctions that return a boolean value $VALUE$ and themselves ($True$/$False$) as $Type$.
..example.text:Print the values of these tags/metafunctions.
..example.code:
std::cout << False::VALUE << std::endl;                         // 0
std::cout << True::VALUE << std::endl;                          // 1
std::cout << IsSameType<False,False::Type>::VALUE << std::endl; // 1
..include:seqan/basic.h
..see:Metafunction.Eval
..see:Metafunction.IsSameType
*/

struct True
{
    typedef True Type;
    enum { VALUE = true };
};

struct False
{
    typedef False Type;
    enum { VALUE = false };
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Eval
// ----------------------------------------------------------------------------

/**
.Metafunction.Eval
..cat:Metaprogramming
..summary:Convert from $bool$ values to types ($True$ and $False$).
..signature:Eval<b>::Type
..param.b:The boolean to evaluate.
...type:nolink:$bool$
..returns:Either @Tag.Logical Values.tag.True@ or @Tag.Logical Values.tag.False@, depending on $b$.
..include:seqan/basic.h
..example.text:Demonstrate the usage of $Eval$ to bridge between $bool$ values and the logical tags.
..example.code:
void printBoolType(True const &)
{
    std::cout << "true" << std::endl;
}

void printBoolType(False const &)
{
    std::cout << "false" << std::endl;
}

int main(int argc, char ** argv)
{
    using namespace seqan;

    printBoolType(Eval<true>::Type());   // => "true\n"
    printBoolType(Eval<false>::Type());  // => "false\n"
    return 0;
}
*/

template <bool b>
struct Eval: False {};

template <>
struct Eval<true>: True {};

// ----------------------------------------------------------------------------
// Metafunction Or
// ----------------------------------------------------------------------------

/**
.Metafunction.Or
..cat:Metaprogramming
..summary:Metaprogramming boolean "or" operator.
..signature:Or<B1, B2>::Type
..param.B1:Left-hand argument.
...type:Tag.Logical Values.tag.True
...type:Tag.Logical Values.tag.False
..param.B2:Right-hand argument.
...type:Tag.Logical Values.tag.True
...type:Tag.Logical Values.tag.False
..returns:One of @Tag.Logical Values.tag.True@ and @Tag.Logical Values.tag.False@, the result of logical or.
The arguments $B1$/$B2$ can either be @Tag.Logical Values.tag.True@/@Tag.Logical Values.tag.False@
or boolean metafunctions that return @Tag.Logical Values.tag.True@/@Tag.Logical Values.tag.False@.
..example.code:
Or<False,False>
Or<False,True>
Or<typename And<T1,T2>::Type,T3>
Or<And<T1,T2>,T3>
..include:seqan/basic.h
*/

template <typename TBool1, typename TBool2>
struct Or: Or<typename TBool1::Type, typename TBool2::Type> {}; 

template <>
struct Or<False, False>: False {};
template <>
struct Or<False, True>: True {};
template <>
struct Or<True, False>: True {};
template <>
struct Or<True, True>: True {};

// ----------------------------------------------------------------------------
// Metafunction And
// ----------------------------------------------------------------------------

/**
.Metafunction.And
..cat:Metaprogramming
..summary:Metaprogramming boolean "And" operator.
..signature:Or<B1, B2>::Type
..param.B1:Left-hand argument.
...type:Tag.Logical Values.tag.True
...type:Tag.Logical Values.tag.False
..param.B2:Right-hand argument.
...type:Tag.Logical Values.tag.True
...type:Tag.Logical Values.tag.False
..returns:One of @Tag.Logical Values.tag.True@ and @Tag.Logical Values.tag.False@, the result of logical and.
..include:seqan/basic.h
*/

template <typename TBool1, typename TBool2>
struct And: And<typename TBool1::Type, typename TBool2::Type> {}; 

template <>
struct And<False, False>: False {};
template <>
struct And<False, True>: False {};
template <>
struct And<True, False>: False {};
template <>
struct And<True, True>: True {};

// ----------------------------------------------------------------------------
// Metafunction If
// ----------------------------------------------------------------------------

/**
.Metafunction.If
..cat:Metaprogramming
..summary:Metaprogramming "if".
..signature:If<b, T1, T2>::Type
..param.b:The condition to evaluate.
...type:nolink:$bool$
..param.T1:Result if $b$.
..param.T2:Result if not $b$.
..returns:If $b$ is $true$ then $T1$, otherwise $T2$.
..include:seqan/basic.h
*/

//TODO: change bool Flag which is a boolean into typename TFlag which is True/False
template <bool Flag, class T1, class T2>
struct If
{
    typedef T1 Type;
};

template <class T1, class T2>
struct If<false, T1, T2>
{
    typedef T2 Type;
};

// ----------------------------------------------------------------------------
// Metafunction IsSameType
// ----------------------------------------------------------------------------

/**
.Metafunction.IsSameType
..cat:Metaprogramming
..summary:Metaprogramming type comparison.
..signature:Or<T1, T2>::Type
..signature:Or<T1, T2>::VALUE
..param.T1:Left-hand argument.
..param.T2:Right-hand argument.
..returns:@Tag.Logical Values.tag.True@/$true$ if $T1$ is the same as $T2$, otherwise @Tag.Logical Values.tag.False@/$false$.
..include:seqan/basic.h
*/

template <class Type1, class Type2>
struct IsSameType: False {};

template <class Type1>
struct IsSameType<Type1, Type1>: True {};

// ----------------------------------------------------------------------------
// Metafunction Switch;  Supporting Tags Case, NilCase.
// ----------------------------------------------------------------------------

/**
.Tag.NilCase
..cat:Metaprogramming
..summary:Metaprogramming $default:$ case expression.
..signature:NilCase
..remarks:$NilCase$ is returned by metafunction @Metafunction.Switch@ if no case matched.
..include:seqan/basic.h
..remarks:The documentation of @Metafunction.Switch@ gives an example.
..see:Tag.Case
..see:Metafunction.Switch

.Tag.Case
..cat:Metaprogramming
..summary:Metaprogramming $case$ expression.
..signature:Case<TAG[, NextCase]>
..param.TAG:
..param.NextCase:Optional next case of type @Tag.Case@, defaults to @Tag.NilCase@.
...type:Tag.NilCase
...type:Tag.Case
..include:seqan/basic.h
..remarks:The documentation of @Metafunction.Switch@ gives an example.
..see:Tag.NilCase
..see:Metafunction.Switch

.Metafunction.Switch
..cat:Metaprogramming
..summary:Metaprogramming $switch$ expression.
..signature:Switch<TAG, Case>::Type
..param.TAG:The case label.
...type:nolink:$int$
..param.Case:Cascade of $Case$ tags.
...type:Tag.Case
...type:Tag.NilCase
..returns:The selected type from the @Tag.Case@ cascade or @Tag.NilCase@.
..include:seqan/basic.h
..see:Tag.NilCase
..see:Tag.Case
..example.code:
int switchTest(Nothing const &) { return -1; }
int switchTest(False const &) { return 0; }
int switchTest(True const &) { return 1; }
int switchTest(NilCase const &) { return 2; }

template <int X>
struct SwitchTest
{
    typedef typename Switch<
        X,
        Case<-1, Nothing,
        Case<0, False,
        Case<1, True
        > > > >::Type Type;
};

SEQAN_DEFINE_TEST(test_metaprogramming_switch)
{
    typedef typename SwitchTest<-1>::Type T1;
    typedef typename SwitchTest< 0>::Type T2;
    typedef typename SwitchTest< 1>::Type T3;
    typedef typename SwitchTest< 2>::Type T4;

    std::cout << switchTest(T1()) << std::endl;  // => -1
    std::cout << switchTest(T2()) << std::endl;  // =>  0
    std::cout << switchTest(T3()) << std::endl;  // =>  1
    std::cout << switchTest(T4()) << std::endl;  // =>  2
 }
 */

const int DEFAULT_ = ~(~0u >> 1); // initialize with the smallest int

struct NilCase {};

template <int TAG, typename Type_, typename Next_ = NilCase>
struct Case
{
    enum { TAG_ = TAG };
    typedef Type_ Type;
    typedef Next_ Next;
};

template <int TAG, typename Case_>
struct Switch
{
    typedef typename Case_::Next NextCase_;
    enum
    {
        CASE_TAG_ = Case_::TAG_,
        FOUND_    = (CASE_TAG_ == TAG || CASE_TAG_ == DEFAULT_)
    };

    typedef typename If<FOUND_,
                        typename Case_::Type,
                        typename Switch<TAG, NextCase_>::Type
                        >::Type Type;
};

template <int TAG>
struct Switch<TAG, NilCase>
{
    typedef NilCase Type;
};

// ----------------------------------------------------------------------------
// Metafunction Loop
// ----------------------------------------------------------------------------

/**
.Metafunction.Loop:
..cat:Metaprogramming
..summary:Metafunction returning a function that iterates over a static integer range.
..signature:Loop<Worker, I>::run(Arg & arg)
..param.Worker:A worker $struct$. It has to implement the static (and preferably finline/inline) function $body$ that accepts two parameters. The first one will be a reference to $arg$, as given to $run()$.  The second will be the current value of the iterating variable.
..param.I:The upper limit for the iteration.
..param.arg:The argument to be passed into the workers' $body()$ function.
..remarks:The loop will go from 1 up to and including $I$.
..see:Metafunction.LoopReverse
..include:seqan/basic.h
..example.text:Print the values 1, 2, ..., $I-1$, $I$.
..example.code:
struct PrintWorker
{
    static inline void body(Nothing & arg, int I)
    {
        (void)arg;  // ignored
        printf("%d\n", I);
    }
};

Loop<PrintWorker, 10>::run(Nothing());
// This will print the numbers 1, 2, ..., 9, 10.
 */

// Example of a loop Worker class.  Could be removed, serves no
// purpose.
struct WorkerNothing
{
    template <typename Arg>
    static inline void body(Arg & /*arg*/, int /*I*/)
    {}
};

template <typename Worker, int I>
class Loop {
public:
    template <typename Arg>
    static inline void run(Arg & arg)
    {
        Loop<Worker, I - 1>::run(arg);
        Worker::body(arg, I);
    }
};

template <typename Worker>
class Loop<Worker, 0>
{
public:
    // end of loop
    template <typename Arg>
    static inline void run(Arg &)
    {}
};

// ----------------------------------------------------------------------------
// Metafunction LoopReverse
// ----------------------------------------------------------------------------

/**
.Metafunction.LoopReverse:
..cat:Metaprogramming
..summary:Metafunction returning a function that iterates over a static integer range in reverse order.
..signature:LoopReverse<Worker, I>::run(Arg & arg)
..param.Worker:A worker $struct$. It has to implement the static (and preferably finline/inline) function $body$ that accepts two parameters. The first one will be a reference to $arg$, as given to $run()$.  The second will be the current value of the iterating variable.
..param.I:The upper limit for the iteration.
..param.arg:The argument to be passed into the workers' $body()$ function.
..remarks:The loop will go from $I$ down to and including 1.
..include:seqan/basic.h
..see:Metafunction.Loop
..example.text:Print the values $I$, $I - 1$, ..., 2, 1.
..example.code:

struct PrintWorker
{
    static inline body(Nothing & arg, int I)
    {
        (void)arg;  // ignored
        printf("%d\n", I);
    }
};

Loop<PrintWorker, 10>::run(Nothing());
// This will print the numbers 1, 2, ..., 9, 10.
 */

template <typename Worker, int I>
class LoopReverse {
public:
    template <typename Arg>
    static inline void run(Arg & arg)
    {
        Worker::body(arg, I);
        LoopReverse<Worker, I - 1>::run(arg);
    }
};

template <typename Worker>
class LoopReverse<Worker, 0> {
  public:
    // end of loop
    template <typename Arg>
    static inline void run(Arg &) {}
};

// ----------------------------------------------------------------------------
// Metafunction Log2
// ----------------------------------------------------------------------------

/**
.Metafunction.Log2
..cat:Metaprogramming
..summary:Compute ceiled logarithm to base 2 using metaprogramming.
..signature:Log2<x>::VALUE
..param.x:The value to take the logarithm of.
...type:nolink:$__int64$
..returns:$ceil(log2(x))$.
..include:seqan/basic.h
 */

template <__int64 numerus>
struct Log2
{
    static const __uint64 VALUE = Log2<(numerus + 1) / 2>::VALUE + 1;		// ceil(log_2(n))
};

// Base cases.
template <> struct Log2<1> { static const __uint64 VALUE = 0; };
template <> struct Log2<0> { static const __uint64 VALUE = 0; };

// ----------------------------------------------------------------------------
// Metafunction Log2Floor
// ----------------------------------------------------------------------------

/**
.Metafunction.Log2Floor
..cat:Metaprogramming
..summary:Compute floored logarithm to base 2 using metaprogramming.
..signature:Log2<x>::VALUE
..param.x:The value to take the logarithm of.
...type:nolink:$__int64$
..returns:$floor(log2(x))$.
..include:seqan/basic.h
 */

template <__int64 numerus>
struct Log2Floor
{
    static const __uint64 VALUE = Log2Floor<numerus / 2>::VALUE + 1;		// floor(log_2(n))
};

// Base cases.
template <> struct Log2Floor<1> { static const __uint64 VALUE = 0; };
template <> struct Log2Floor<0> { static const __uint64 VALUE = 0; };

// ----------------------------------------------------------------------------
// Metafunction Power
// ----------------------------------------------------------------------------

/**
.Metafunction.Power
..cat:Metaprogramming
..summary:Compute power of a number.
..signature:Power<b, e>::VALUE
..param.b:The base.
...type:nolink:$__int64$
..param.e:The exponent.
...type:nolink:$__int64$
..returns:$b^e$
..include:seqan/basic.h
 */

template <__int64 base, __int64 exponent>
struct Power {
    static const __uint64 VALUE =
            Power<base, exponent / 2>::VALUE *
            Power<base, exponent - (exponent / 2)>::VALUE;
};

// Base cases.
template <__int64 base> struct Power<base, 1> { static const __uint64 VALUE = base; };
template <__int64 base> struct Power<base, 0> { static const __uint64 VALUE = 1; };

// ----------------------------------------------------------------------------
// Metafunction MakeUnsigned
// ----------------------------------------------------------------------------

/**
.Metafunction.MakeUnsigned:
..cat:Basic
..summary:Converts an integral value into an unsigned integral value.
..signature:MakeUnsigned<T>::Type
..param.T:Input integral type.
..returns.param.Type:A type without a sign of the same domain, e.g. $unsigned int$ for $T$ = $int$.
...default:$T$
..include:seqan/basic.h
 */

template <typename T>
struct MakeUnsigned
{
	typedef
		typename If< IsSameType<T, char>::VALUE,         unsigned char,
		typename If< IsSameType<T, signed char>::VALUE,  unsigned char,
		typename If< IsSameType<T, signed short>::VALUE, unsigned short,
		typename If< IsSameType<T, signed int>::VALUE,   unsigned int,
		typename If< IsSameType<T, signed long>::VALUE,  unsigned long,
		typename If< IsSameType<T, __int64>::VALUE,      __uint64, T
		>::Type>::Type>::Type>::Type>::Type>::Type Type;
};

template <typename T>
struct MakeUnsigned<T const>
{
	typedef typename MakeUnsigned<T>::Type const Type;
};

/**
.Internal.MakeUnsigned_:
..signature:MakeUnsigned_<T>
..status:deprecated, please use @Metafunction.MakeUnsigned@
..returns:$unsigned t$ if $T$ is not $unsigned t$, otherwise $T$.
*/
template <typename T>
struct MakeUnsigned_ : MakeUnsigned<T> {};

// ----------------------------------------------------------------------------
// Metafunction MakeSigned
// ----------------------------------------------------------------------------

// TODO(holtgrew): Make public, complete documentation.

/**
.Metafunction.MakeSigned:
..cat:Basic
..summary:Converts an integral value into a signed integral value.
..signature:MakeSigned<T>::Type
..param.T:Input integral type.
..returns.param.Type:A type with a sign of the same domain, e.g. $int$ for $T$ = $unsigned int$.
...default:$T$
..include:seqan/basic.h
..see:Metafunction.MakeUnsigned
 */

template <typename T>
struct MakeSigned
{
	typedef
		typename If< IsSameType<T, char>::VALUE,           signed char,
		typename If< IsSameType<T, unsigned char>::VALUE,  signed char,
		typename If< IsSameType<T, unsigned short>::VALUE, signed short,
		typename If< IsSameType<T, unsigned int>::VALUE,   signed int,
		typename If< IsSameType<T, unsigned long>::VALUE,  signed long,
		typename If< IsSameType<T, __uint64>::VALUE,       __int64, T
		>::Type>::Type>::Type>::Type>::Type>::Type Type;
};

template <typename T>
struct MakeSigned<T const>
{
	typedef typename MakeSigned<T>::Type const Type;
};

/**
.Internal.MakeSigned_:
..signature:MakeSigned_<T>
..status:deprecated, please use @Metafunction.MakeSigned@
..returns:$signed t$ if $T$ is not $signed t$, otherwise $T$.
*/
template <typename T>
struct MakeSigned_ : MakeSigned<T> {};

// ----------------------------------------------------------------------------
// Metafunction RemoveReference
// ----------------------------------------------------------------------------

/**
.Metafunction.RemoveReference:
..cat:Basic
..summary:Converts a (reference) type into the same type without reference.
..signature:RemoveReference<T>::Type
..param.T:Input type.
..returns.param.Type:A corresponding non-reference type, e.g. $int$ for $T$ = $&int$.
...default:$T$
..include:seqan/basic.h
..see:Metafunction.RemoveConst
*/

/**
.Internal.RemoveReference_:
..signature:RemoveReference_<T>
..status:deprecated, please use @Metafunction.RemoveReference@
..returns:$t$ if $T$ is $t &$, otherwise $T$.
*/

template <typename T>
struct RemoveReference
{
	typedef T Type;
};

template <typename T>
struct RemoveReference<T &> : RemoveReference<T> {};

template <typename T>
struct RemoveReference_ : RemoveReference<T> {};

// ----------------------------------------------------------------------------
// Metafunction RemoveConst
// ----------------------------------------------------------------------------

/**
.Metafunction.RemoveConst:
..cat:Basic
..summary:Converts a (const) type into the corresponding non-const type.
..signature:RemoveConst<T>::Type
..param.T:Input type.
..returns.param.Type:A corresponding non-const type, e.g. $int$ for $T$ = $const int$.
...default:$T$
..include:seqan/basic.h
*/

/**
.Internal.RemoveConst_:
..signature:RemoveConst_<T>
..status:deprecated, please use @Metafunction.RemoveConst@
..returns:$t$ if $T$ is $t const$, otherwise $T$.
*/

template <typename T>
struct RemoveConst
{
	typedef T Type;
};

template <typename T>
struct RemoveConst<T const> : public RemoveConst<T> {};

template <typename T>
struct RemoveConst<T &>
{
	typedef typename RemoveConst<T>::Type & Type;
};

template <typename T>
struct RemoveConst<T *>
{
	typedef typename RemoveConst<T>::Type * Type;
};

template <typename T, size_t I>
struct RemoveConst<T const [I]>
{
	typedef T * Type;
};

template <typename T>
struct RemoveConst_ : RemoveConst<T> {};

// ----------------------------------------------------------------------------
// Metafunction CopyConst_
// ----------------------------------------------------------------------------

// TODO(holtgrew): Make public, document.

// returns TTo const, if TFrom is const, TTo otherwise

template <typename TFrom, typename TTo>
struct CopyConst_
{
	typedef TTo Type;
};

template <typename TFrom, typename TTo>
struct CopyConst_<TFrom const, TTo>
{
	typedef TTo const Type;
};

// ----------------------------------------------------------------------------
// Metafunction IsConst_
// ----------------------------------------------------------------------------

// TODO(holtgrew): Make public, complete documentation.

/**
.Internal.IsConst_:
..signature:IsConst_<T>
..returns:@Tag.Logical Values.tag.True@ if $T$ is $t const$, otherwise @Tag.Logical Values.tag.False@.
*/

template <typename T>
struct IsConst_
{
	typedef False Type;
};

template <typename T>
struct IsConst_<T const>
{
    typedef True Type;
};

// ----------------------------------------------------------------------------
// Metafunction ClassIdentifier_
// ----------------------------------------------------------------------------

// TODO(holtgrew): Make public, complete documentation.

/**
.Internal.ClassIdentifier_:
..signature:void * ClassIdentifier_<T>::getID()
..returns:A void * that identifies $T$.
...text:The returned values of two calls of $getID$ are equal if and only if
the used type $T$ was the same.
 */

template <typename T>
struct ClassIdentifier_
{
	static inline void *
	getID()
	{
        SEQAN_CHECKPOINT;
		static bool _id_dummy;
		return &_id_dummy;
	}
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function memset()
// ----------------------------------------------------------------------------

/**
.Function.memset
..cat:Memory
..summary:An implementation of $memset$ with fixed number of bytes using Metaprogramming.
..signature:memset<SIZE>(ptr, c)
..signature:memset<SIZE, c>(ptr)
..param.SIZE:The number of bytes to set.
...type:nolink:$unsigned$
..param.ptr:Pointer to the data to set.
...type:nolink:$unsigned char *$
..param.c:The character to fill the memory with.
....type:nolink:$unsigned char$
..remarks:These functions can be completely unrolled and inlined by the compiler.
..include:seqan/basic.h
 */

// TODO(holtgrew): Does memset() really belong in this header? Used in find_myers_ukknonen.h, pump_lcp_core.h, pipe_sample.h, file_async

using ::std::memset;

// Implementation of memset() with fill size.

template <unsigned SIZE, bool direct>
struct MemsetWorker
{
    finline static
    void run(unsigned char * ptr, unsigned char c)
    {
        ::std::memset(ptr, c, SIZE);
    }
};

template <unsigned  SIZE>
struct MemsetWorker<SIZE, true>
{
    finline static
    void run(unsigned char* ptr, unsigned char c)
    {
        *((unsigned*)ptr) = ((unsigned)c << 24) + ((unsigned)c << 16) + ((unsigned)c << 8) + (unsigned)c;
        MemsetWorker<SIZE - 4, true>::run(ptr + 4, c);
    }
};

template <>
struct MemsetWorker<0, true>
{
    finline static void
    run(unsigned char*, unsigned char)
    {}
};

template <>
struct MemsetWorker<1, true>
{
    finline static
    void run(unsigned char* ptr, unsigned char c)
    {
        *ptr = c;
    }
};

template <>
struct MemsetWorker<2, true>
{
    finline static
    void run(unsigned char* ptr, unsigned char c)
    {
        *(unsigned short *)ptr = ((unsigned short)c << 8) + (unsigned short)c;
    }
};

template <>
struct MemsetWorker<3, true> {
    finline static
    void run(unsigned char* ptr, unsigned char c)
    {
        MemsetWorker<2, true>::run(ptr, c);
        MemsetWorker<1, true>::run(ptr + 2, c);
    }
};

template <unsigned SIZE>
finline void memset(void* ptr, unsigned char c)
{
    MemsetWorker<SIZE, SIZE <= 32>::run((unsigned char*)ptr, c);
}

// Implementation of memset() with fill value.

template <unsigned SIZE, bool direct, unsigned char c>
struct MemsetConstValueWorker
{
    finline static void run(unsigned char* ptr) { ::std::memset(ptr, c, SIZE); }
};

template <unsigned  SIZE, unsigned char c>
struct MemsetConstValueWorker<SIZE, true, c>
{
    finline static
    void run(unsigned char* ptr)
    {
        *((unsigned*)ptr) = ((unsigned)c << 24) + ((unsigned)c << 16) + ((unsigned)c << 8) + (unsigned)c;
        MemsetConstValueWorker<SIZE - 4, true, c>::run(ptr + 4);
    }
};

template <unsigned char c>
struct MemsetConstValueWorker<0, true, c>
{
    finline static
    void run(unsigned char*) {}
};

template <unsigned char c>
struct MemsetConstValueWorker<1, true, c>
{
    finline static
    void run(unsigned char* ptr) {
        *ptr = c;
    }
};

template <unsigned char c>
struct MemsetConstValueWorker<2, true, c>
{
    finline static
    void run(unsigned char* ptr)
    {
        *(unsigned short *)ptr = ((unsigned short)c << 8) + (unsigned short)c;
    }
};

template <unsigned char c>
struct MemsetConstValueWorker<3, true, c>
{
    finline static
    void run(unsigned char* ptr)
    {
        MemsetConstValueWorker<2, true, c>::run(ptr);
        MemsetConstValueWorker<1, true, c>::run(ptr + 2);
    }
};

template <unsigned SIZE, unsigned char c>
finline void
memset(void* ptr)
{
    MemsetConstValueWorker<SIZE, SIZE <= 32, c>::run((unsigned char*)ptr);
}

}  // namespace seqan

#endif  // #ifndef SEQAN_BASIC_BASIC_METAPROGRAMMING_H_
