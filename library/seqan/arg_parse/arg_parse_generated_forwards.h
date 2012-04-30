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


#ifndef SEQAN_HEADER_ARG_PARSE_GENERATED_FORWARDS_H 
#define SEQAN_HEADER_ARG_PARSE_GENERATED_FORWARDS_H 

//////////////////////////////////////////////////////////////////////////////
// NOTE: This file is automatically generated by build_forwards.py
//       Do not edit this file manually!
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
// CLASSES
//////////////////////////////////////////////////////////////////////////////

namespace seqan {

//____________________________________________________________________________
// ArgParseArgument

class ArgParseArgument;       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/arg_parse_argument.h"(95)
class ArgParseArgument;       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/arg_parse_type_support.h"(132)

//____________________________________________________________________________
// ArgParseOption

class ArgParseOption;       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/arg_parse_option.h"(90)
class ArgParseOption;       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/arg_parse_type_support.h"(116)
class ArgParseOption;       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/argument_parser.h"(61)

//____________________________________________________________________________
// ArgumentParser

class ArgumentParser;       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/argument_parser.h"(60)
class ArgumentParser;       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/argument_parser.h"(116)

//____________________________________________________________________________
// InvalidOptionException

class InvalidOptionException;       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/arg_parse_exceptions.h"(92)

//____________________________________________________________________________
// MissingArgumentException

class MissingArgumentException;       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/arg_parse_exceptions.h"(121)

//____________________________________________________________________________
// NotEnoughArgumentsException

class NotEnoughArgumentsException;       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/arg_parse_exceptions.h"(146)

//____________________________________________________________________________
// ParseException

class ParseException;       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/arg_parse_exceptions.h"(60)

} //namespace seqan


//////////////////////////////////////////////////////////////////////////////
// TYPEDEFS


//////////////////////////////////////////////////////////////////////////////
// FUNCTIONS
//////////////////////////////////////////////////////////////////////////////

namespace seqan {

//____________________________________________________________________________
// _addMinMaxRestrictions

inline void _addMinMaxRestrictions(std::vector<std::string> & restrictions, ArgParseOption const & opt);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/arg_parse_ctd_support.h"(135)

//____________________________________________________________________________
// _addUsage

void _addUsage(ToolDoc & toolDoc, ArgumentParser const & me);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/arg_parse_doc.h"(168)

//____________________________________________________________________________
// _addValidValuesRestrictions

inline void _addValidValuesRestrictions(std::vector<std::string> & restrictions, ArgParseOption const & opt);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/arg_parse_ctd_support.h"(161)

//____________________________________________________________________________
// _allArgumentsSet

inline bool _allArgumentsSet(ArgumentParser const & me);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/argument_parser.h"(461)

//____________________________________________________________________________
// _allRequiredSet

inline bool _allRequiredSet(ArgumentParser const & me);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/argument_parser.h"(448)

//____________________________________________________________________________
// _cast

template <typename TTarget, typename TString> inline TTarget _cast(TString const s);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/arg_parse_type_support.h"(73)

//____________________________________________________________________________
// _checkNumericArgument

template <typename TNumerical> inline void _checkNumericArgument(ArgParseArgument const & me, std::string const & value);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/arg_parse_argument.h"(566)

//____________________________________________________________________________
// _checkStringRestrictions

inline void _checkStringRestrictions(ArgParseArgument const & me, std::string const & value);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/arg_parse_argument.h"(586)

//____________________________________________________________________________
// _convertArgumentValue

inline bool _convertArgumentValue(bool & dst, ArgParseOption const & opt, std::string const & src);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/arg_parse_type_support.h"(120)
inline bool _convertArgumentValue(int & dst, ArgParseArgument const & opt, std::string const & src);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/arg_parse_type_support.h"(138)
inline bool _convertArgumentValue(unsigned int & dst, ArgParseArgument const & opt, std::string const & src);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/arg_parse_type_support.h"(146)
inline bool _convertArgumentValue(__int64 & dst, ArgParseArgument const & opt, std::string const & src);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/arg_parse_type_support.h"(154)
inline bool _convertArgumentValue(__uint64 & dst, ArgParseArgument const & opt, std::string const & src);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/arg_parse_type_support.h"(162)
inline bool _convertArgumentValue(float & dst, ArgParseArgument const & opt, std::string const & src);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/arg_parse_type_support.h"(170)
inline bool _convertArgumentValue(double & dst, ArgParseArgument const & opt, std::string const & src);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/arg_parse_type_support.h"(178)
template <typename TObject> inline bool _convertArgumentValue(TObject & dst, ArgParseArgument const & opt, std::string const & src);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/arg_parse_type_support.h"(187)

//____________________________________________________________________________
// _getOptionIndex

inline ArgumentParser::TOptionMapSize _getOptionIndex(ArgumentParser const & me, std::string const & name);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/argument_parser.h"(302)

//____________________________________________________________________________
// _includeInCTD

inline bool _includeInCTD(ArgParseOption const & opt);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/arg_parse_ctd_support.h"(192)

//____________________________________________________________________________
// _intervalAssert

template <typename TIntervalBorder> inline void _intervalAssert(const std::string minValueAsString, const std::string maxValueAsString);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/arg_parse_argument.h"(415)

//____________________________________________________________________________
// _isCastable

template <typename TTarget, typename TString> inline bool _isCastable(TString const s);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/arg_parse_type_support.h"(86)

//____________________________________________________________________________
// _isDouble

template <typename TString> inline bool _isDouble(TString const s);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/arg_parse_type_support.h"(98)

//____________________________________________________________________________
// _isInInterval

template <typename TTarget, typename TString> inline bool _isInInterval(TString value, TString lowerIntervalBound, TString upperIntervalBound);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/arg_parse_argument.h"(548)

//____________________________________________________________________________
// _isInt

template <typename TString> inline bool _isInt(TString const s);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/arg_parse_type_support.h"(108)

//____________________________________________________________________________
// _join

template <typename TValue> inline CharString _join(StringSet<TValue> const & v, CharString const & delimiter);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/arg_parse_ctd_support.h"(60)

//____________________________________________________________________________
// _parseAppName

inline void _parseAppName(ArgumentParser & parser, std::string const & candidate);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/arg_parse_doc.h"(74)

//____________________________________________________________________________
// _tryCast

template <typename TTarget, typename TString> inline bool _tryCast(TTarget & dest, TString const source);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/arg_parse_type_support.h"(62)

//____________________________________________________________________________
// _typeToString

inline std::string _typeToString(ArgParseArgument const & me);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/arg_parse_argument.h"(199)

//____________________________________________________________________________
// _writeOptName

template <typename TStream> inline void _writeOptName(TStream & target, ArgParseOption const & me);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/arg_parse_option.h"(301)

//____________________________________________________________________________
// _xmlEscape

template <typename TSequence> inline TSequence _xmlEscape(TSequence const & original);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/arg_parse_ctd_support.h"(88)
inline std::string _xmlEscape(std::string const & original);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/arg_parse_ctd_support.h"(109)

//____________________________________________________________________________
// addArgument

inline void addArgument(ArgumentParser & me, ArgParseArgument const & arg);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/argument_parser.h"(278)

//____________________________________________________________________________
// addDescription

inline void addDescription(ArgumentParser & me, std::string const & description);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/arg_parse_doc.h"(196)

//____________________________________________________________________________
// addLine

template <typename TString> inline void addLine(ArgumentParser & me, TString const & line);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/arg_parse_doc.h"(104)

//____________________________________________________________________________
// addOption

inline void addOption(ArgumentParser & me, ArgParseOption const & opt);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/argument_parser.h"(247)

//____________________________________________________________________________
// addSection

template <typename TString> inline void addSection(ArgumentParser & me, TString const & line);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/arg_parse_doc.h"(138)

//____________________________________________________________________________
// addText

inline void addText(ArgumentParser & me, std::string const & text);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/arg_parse_doc.h"(367)

//____________________________________________________________________________
// addTextSection

inline void addTextSection(ArgumentParser & me, std::string const & title);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/arg_parse_doc.h"(324)

//____________________________________________________________________________
// addTextSubSection

inline void addTextSubSection(ArgumentParser & me, std::string const & title);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/arg_parse_doc.h"(346)

//____________________________________________________________________________
// addUsageLine

inline void addUsageLine(ArgumentParser & me, std::string const & line);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/arg_parse_doc.h"(159)

//____________________________________________________________________________
// getAppName

inline CharString const & getAppName(ArgumentParser const & parser);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/arg_parse_doc.h"(65)

//____________________________________________________________________________
// getArgument

inline ArgParseArgument & getArgument(ArgumentParser & me, unsigned position);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/argument_parser.h"(408)
inline ArgParseArgument const & getArgument(ArgumentParser const & me, unsigned position);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/argument_parser.h"(415)

//____________________________________________________________________________
// getArgumentLabel

inline std::string const getArgumentLabel(ArgParseArgument const & me);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/arg_parse_argument.h"(371)
inline std::string const getArgumentLabel(ArgParseOption const & me);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/arg_parse_option.h"(288)

//____________________________________________________________________________
// getArgumentValue

inline std::string const & getArgumentValue(ArgParseArgument const & me, unsigned position);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/arg_parse_argument.h"(695)
inline std::string const & getArgumentValue(ArgParseArgument const & me);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/arg_parse_argument.h"(706)
template <typename TValue> inline bool getArgumentValue(TValue & value, ArgumentParser & me, unsigned argumentPosition, unsigned argNo);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/argument_parser.h"(581)
template <typename TValue> inline bool getArgumentValue(TValue & value, ArgumentParser & me, unsigned argumentPosition);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/argument_parser.h"(590)

//____________________________________________________________________________
// getArgumentValueCount

inline unsigned getArgumentValueCount(ArgumentParser const & me, unsigned argumentPosition);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/argument_parser.h"(552)

//____________________________________________________________________________
// getArgumentValues

inline std::vector<std::string> const & getArgumentValues(ArgParseArgument const & me);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/arg_parse_argument.h"(726)
inline std::vector<std::string> const & getArgumentValues(ArgumentParser & me, unsigned argumentPosition);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/argument_parser.h"(635)

//____________________________________________________________________________
// getOption

inline ArgParseOption & getOption(ArgumentParser & me, std::string const & name);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/argument_parser.h"(332)
inline ArgParseOption const & getOption(ArgumentParser const & me, std::string const & name);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/argument_parser.h"(338)

//____________________________________________________________________________
// getOptionValue

template <typename TValue> inline bool getOptionValue(TValue & val, ArgumentParser const & me, std::string const & name, unsigned argNo);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/argument_parser.h"(493)
template <typename TValue> inline bool getOptionValue(TValue & val, ArgumentParser const & me, std::string const & name);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/argument_parser.h"(502)

//____________________________________________________________________________
// getOptionValueCount

inline unsigned getOptionValueCount(ArgumentParser const & me, std::string const & name);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/argument_parser.h"(523)

//____________________________________________________________________________
// getOptionValues

inline std::vector<std::string> const & getOptionValues(ArgumentParser & me, std::string const & name);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/argument_parser.h"(612)

//____________________________________________________________________________
// getVersion

inline CharString const & getVersion(ArgumentParser const & me);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/arg_parse_doc.h"(281)

//____________________________________________________________________________
// hasOption

inline bool hasOption(ArgumentParser const & me, std::string const & name);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/argument_parser.h"(226)

//____________________________________________________________________________
// hasValue

inline bool hasValue(ArgParseArgument const & arg, unsigned position);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/arg_parse_argument.h"(751)
inline bool hasValue(ArgParseArgument const & arg);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/arg_parse_argument.h"(756)

//____________________________________________________________________________
// hideOption

inline void hideOption(ArgParseOption & me, bool hide );       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/arg_parse_option.h"(236)
inline void hideOption(ArgumentParser & me, std::string const & name, bool hide );       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/argument_parser.h"(386)

//____________________________________________________________________________
// isBooleanOption

inline bool isBooleanOption(ArgParseOption const & me);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/arg_parse_option.h"(195)

//____________________________________________________________________________
// isDoubleArgument

inline bool isDoubleArgument(ArgParseArgument const & me);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/arg_parse_argument.h"(310)

//____________________________________________________________________________
// isInputFileArgument

inline bool isInputFileArgument(ArgParseArgument const & me);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/arg_parse_argument.h"(330)

//____________________________________________________________________________
// isIntegerArgument

inline bool isIntegerArgument(ArgParseArgument const & me);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/arg_parse_argument.h"(290)

//____________________________________________________________________________
// isListArgument

inline bool isListArgument(ArgParseArgument const & me);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/arg_parse_argument.h"(248)

//____________________________________________________________________________
// isOutputFileArgument

inline bool isOutputFileArgument(ArgParseArgument const & me);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/arg_parse_argument.h"(351)

//____________________________________________________________________________
// isRequired

inline bool isRequired(ArgParseOption const & me);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/arg_parse_option.h"(256)

//____________________________________________________________________________
// isSet

inline bool isSet(ArgParseArgument const & me);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/arg_parse_argument.h"(776)
inline bool isSet(ArgumentParser const & me, std::string const & name);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/argument_parser.h"(438)

//____________________________________________________________________________
// isStringArgument

inline bool isStringArgument(ArgParseArgument const & me);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/arg_parse_argument.h"(268)
inline bool isStringArgument(ArgParseOption const & me);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/arg_parse_option.h"(175)

//____________________________________________________________________________
// isVisible

inline bool isVisible(ArgParseOption const & me);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/arg_parse_option.h"(215)

//____________________________________________________________________________
// numberOfAllowedValues

inline unsigned numberOfAllowedValues(ArgParseArgument const & me);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/arg_parse_argument.h"(796)

//____________________________________________________________________________
// operator<<

template <typename TStream> inline TStream & operator<<(TStream & target, ArgParseOption const & source);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/arg_parse_option.h"(345)

//____________________________________________________________________________
// parse

inline ArgumentParser::ParseResult parse(ArgumentParser & me, int argc, const char * argv[], std::ostream & outputStream, std::ostream & errorStream);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/arg_parse_parse.h"(68)
inline ArgumentParser::ParseResult parse(ArgumentParser & me, int argc, const char * argv[]);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/arg_parse_parse.h"(230)

//____________________________________________________________________________
// printHelp

inline void printHelp(ArgumentParser const & me, std::ostream & stream, CharString const & format);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/arg_parse_doc.h"(445)
inline void printHelp(ArgumentParser const & me, std::ostream & stream);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/arg_parse_doc.h"(518)

//____________________________________________________________________________
// printShortHelp

inline void printShortHelp(ArgumentParser const & me, std::ostream & stream);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/arg_parse_doc.h"(387)

//____________________________________________________________________________
// printVersion

inline void printVersion(ArgumentParser const & me, std::ostream & stream);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/arg_parse_doc.h"(418)
inline void printVersion(ArgumentParser const & me);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/arg_parse_doc.h"(423)

//____________________________________________________________________________
// setAppName

inline void setAppName(ArgumentParser & me, std::string const & name);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/arg_parse_doc.h"(217)

//____________________________________________________________________________
// setDate

inline void setDate(ArgumentParser & me, std::string const & date);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/arg_parse_doc.h"(302)

//____________________________________________________________________________
// setMaxValue

inline void setMaxValue(ArgParseArgument & me, const std::string maxValue);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/arg_parse_argument.h"(474)
inline void setMaxValue(ArgumentParser & me, std::string const & name, std::string const & _maxValue);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/argument_parser.h"(689)
inline void setMaxValue(ArgumentParser & me, unsigned argumentPosition, std::string const & _minValue);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/argument_parser.h"(697)

//____________________________________________________________________________
// setMinValue

inline void setMinValue(ArgParseArgument & me, const std::string minValue);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/arg_parse_argument.h"(440)
inline void setMinValue(ArgumentParser & me, std::string const & name, std::string const & _minValue);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/argument_parser.h"(658)
inline void setMinValue(ArgumentParser & me, unsigned argumentPosition, std::string const & _minValue);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/argument_parser.h"(666)

//____________________________________________________________________________
// setRequired

inline void setRequired(ArgParseOption & me, bool required);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/arg_parse_option.h"(277)
inline void setRequired(ArgumentParser & me, std::string const & name, bool required );       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/argument_parser.h"(362)

//____________________________________________________________________________
// setShortDescription

inline void setShortDescription(ArgumentParser & me, std::string const & description);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/arg_parse_doc.h"(238)

//____________________________________________________________________________
// setValidValues

inline void setValidValues(ArgParseArgument & me, std::vector<std::string> const & values);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/arg_parse_argument.h"(510)
inline void setValidValues(ArgParseArgument & me, std::string const & valuesString);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/arg_parse_argument.h"(518)
inline void setValidValues(ArgumentParser & me, std::string const & name, std::vector<std::string> const & values);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/argument_parser.h"(721)
inline void setValidValues(ArgumentParser & me, std::string const & name, std::string const & values);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/argument_parser.h"(729)
inline void setValidValues(ArgumentParser & me, unsigned argumentPosition, std::vector<std::string> const & values);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/argument_parser.h"(737)
inline void setValidValues(ArgumentParser & me, unsigned argumentPosition, std::string const & values);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/argument_parser.h"(745)

//____________________________________________________________________________
// setVersion

inline void setVersion(ArgumentParser & me, std::string const & versionString);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/arg_parse_doc.h"(259)

//____________________________________________________________________________
// throw

inline void _assignArgumentValue(ArgParseArgument & me, std::string const & value) throw (ParseException);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/arg_parse_argument.h"(651)

//____________________________________________________________________________
// toCString

inline char const * toCString(std::string const & me);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/arg_parse_type_support.h"(53)

//____________________________________________________________________________
// write

template <typename TStream> inline void write(TStream & target, ArgParseOption const & me);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/arg_parse_option.h"(330)

//____________________________________________________________________________
// writeCTD

inline void writeCTD(ArgumentParser const & me);       	// "/Users/f.buske/Documents/biosim/triplexes/progs/seqan/core/include/seqan/arg_parse/arg_parse_ctd_support.h"(212)

} //namespace seqan

#endif

