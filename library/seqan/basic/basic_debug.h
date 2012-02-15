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

#ifndef SEQAN_HEADER_BASIC_DEBUG_H
#define SEQAN_HEADER_BASIC_DEBUG_H

// TODO(holtgrew): Move to basic_testing.h

/* IOREV
 *
 * contains a TODO by holtgrew indicating that this should be merged into
 * basic_testing.h
 *
 * these functions use custom IO, i.e. FILE* fopen etc
 * might be ok for binary comparison
 *
 * not sure what the text-comparison is all about. If you want to compare
 * text, ignoring differences in EOL-style, you should not use binary
 * file reading at all, but text-mode reading which converts line-endings
 * on the fly anyhow, rendering this code useless (IMHO)
 *
 * code contains goto statements which could easily be replaced by clean code
 * (I didn't even know this is possible in c++
 */


// compare two files, do not translate linebreaks
inline bool 
_compareBinaryFiles(const char * file1, const char * file2)
{
//IOREV see above
    bool ret = false;

    FILE * fl1 = fopen(file1, "rb");
    if (!fl1) return ret;

    FILE * fl2 = fopen(file2, "rb");
    if (!fl2)
    {
        fclose(fl1);
        return ret;
    }

    while (!feof(fl1) && !feof(fl2))
    {
        if (fgetc(fl1) != fgetc(fl2)) goto End;
    }

    ret = feof(fl1) && feof(fl2);

End:
    fclose(fl2);
    fclose(fl1);

    return ret;

}
//____________________________________________________________________________

//one line break is either \r, \n, or \r\n.
//a single line break is skipped.
//the second line break is transformed into \n
inline void 
_compareTextFilesReadChar(FILE * fl, char & c, int & num_lb, bool & is_eof)
{
//IOREV see above
    num_lb = 0;
    is_eof = false;

    c = fgetc(fl);
    while ((c == '\r') || (c == '\n'))
    {
        ++num_lb;
        if (c == '\r')
        {
            c = fgetc(fl);
            if (feof(fl)) is_eof = true;
            else
            {
                if (c == '\n')
                {
                    c = fgetc(fl);
                    if (feof(fl)) is_eof = true;
                }
            }
        }
        else if (c == '\n')
        {
            c = fgetc(fl);
            if (feof(fl)) is_eof = true;
        }
    }
}

// compare two files, translate linebreaks
inline bool 
_compareTextFiles(const char * file1, const char * file2)
{
//IOREV see above
    FILE * fl1 = fopen(file1, "rb");
    if (!fl1) return false;

    FILE * fl2 = fopen(file2, "rb");
    if (!fl2)
    {
        fclose(fl1);
        return false;
    }

    bool ret = false;

    int num_lb1, num_lb2;
    bool is_eof1, is_eof2;
    char c1, c2;

    while (!feof(fl1) && !feof(fl2))
    {
        _compareTextFilesReadChar(fl1, c1, num_lb1, is_eof1);
        _compareTextFilesReadChar(fl2, c2, num_lb2, is_eof2);

        if (num_lb1 != num_lb2)
        {
            goto End;
        }
        if (is_eof1 ^ is_eof2)
        {
            goto End;
        }
        if (c1 != c2)
        {
            goto End;
        }
    }

    ret = feof(fl1) && feof(fl2);

End:
    fclose(fl2);
    fclose(fl1);

    return ret;

}

#endif //#ifndef SEQAN_HEADER_...
