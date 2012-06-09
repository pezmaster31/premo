// ***************************************************************************
// fastqreader.h (c) 2012 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 9 June 2012 (DB)
// ---------------------------------------------------------------------------
// FASTQ file reader
// ***************************************************************************

#ifndef FASTQREADER_H
#define FASTQREADER_H

#include <cstdio>
#include <string>
class Fastq;

class FastqReader {

    // ctor & dtor
    public:
        FastqReader(void);
        ~FastqReader(void);

    // FastqReader interface
    public:
        void close(void);
        std::string errorString(void) const;
        std::string filename(void) const;
        bool isOpen(void) const;
        bool open(const std::string& filename);
        bool readNext(Fastq* entry);

    // data members
    private:
        FILE*       m_stream;
        char*       m_buffer;
        size_t      m_bufferLength;
        std::string m_filename;
        std::string m_errorString;
};

#endif // FASTQREADER_H
