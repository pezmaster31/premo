// ***************************************************************************
// fastqreader.cpp (c) 2012 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 19 June 2012 (DB)
// ---------------------------------------------------------------------------
// FASTQ file reader
// ***************************************************************************

#include "fastqreader.h"
#include "fastq.h"
#include <bamtools/api/bamtools_global.h>
#include <cassert>
#include <cstring>
#include <sstream>
using namespace std;

// ------------------------
// static utility methods
// ------------------------

static
void chomp(char* s) {

    size_t length = strlen(s);
    if ( length == 0 )
        return;
    --length;

    while ( (s[length] == 10) || (s[length] == 13) ) {
        s[length] = 0;
        --length;
        if ( length < 0 )
            break;
    }
}

// ----------------------------
// FastqReader implementation
// ----------------------------

// these macros allow for compression-independent file I/O below
#define STREAM_CLEAR          ( m_isCompressed ? m_zStream = 0                               : m_stream = 0 )
#define STREAM_CLOSE          ( m_isCompressed ? gzclose(m_zStream)                          : fclose(m_stream) )
#define STREAM_EOF            ( m_isCompressed ? gzeof(m_zStream)                            : feof(m_stream) )
#define STREAM_GETC           ( m_isCompressed ? gzgetc(m_zStream)                           : fgetc(m_stream) )
#define STREAM_GETS           ( m_isCompressed ? gzgets(m_zStream, m_buffer, m_bufferLength) : fgets(m_buffer, m_bufferLength, m_stream) )
#define STREAM_OPEN(filename) ( m_isCompressed ? m_zStream = gzopen(filename, "rb")          : m_stream = fopen(filename, "rb") )
#define STREAM_UNGETC(ch)     ( m_isCompressed ? gzungetc(ch, m_zStream)                     : ungetc(ch, m_stream) )
#define STREAM                ( m_isCompressed ? m_zStream                                   : m_stream )

FastqReader::FastqReader(void)
    : m_stream(0)
    , m_zStream(0)
    , m_isCompressed(false)
    , m_buffer(0)
    , m_bufferLength(0)
{ }

FastqReader::~FastqReader(void) {
    close();
}

void FastqReader::close(void) {

    // close file stream
    if ( isOpen() ) {
        STREAM_CLOSE;
        STREAM_CLEAR;
    }

    // clean up allocated memory
    if ( m_buffer ) {
        delete[] m_buffer;
        m_buffer = 0;
        m_bufferLength = 0;
    }

    // clear any other file-dependent data
    m_filename.clear();
    m_isCompressed = false;
}

string FastqReader::errorString(void) const {
    return m_errorString;
}

string FastqReader::filename(void) const {
    return m_filename;
}

bool FastqReader::isOpen(void) const {
    return ( STREAM != 0 );
}

bool FastqReader::open(const string& filename) {

    // ensure clean slate
    close();

    // -----------------------------
    // check the compression state
    // -----------------------------

    FILE* checkStream = fopen(filename.c_str(), "rb");
    if ( checkStream == 0 ) {

        // if failed, set error & return failure
        m_errorString = "could not open input FASTQ file: ";
        m_errorString.append(filename);
        return false;
    }

    const uint16_t GZIP_MAGIC_NUMBER = 0x8b1f;
    uint16_t magicNumber = 0;
    const size_t numBytes = fread((char*)&magicNumber, sizeof(magicNumber), 1, checkStream);
    (void)numBytes;
    fclose(checkStream);

    m_isCompressed = ( magicNumber == GZIP_MAGIC_NUMBER );

    // ----------------------
    // attempt to open file
    // ----------------------

    STREAM_OPEN(filename.c_str());
    if ( !isOpen() ) {

        // if failed, set error & return failure
        m_errorString = "could not open input FASTQ file: ";
        m_errorString.append(filename);
        return false;
    }

    // create an input buffer
    m_bufferLength = 4096;
    m_buffer = new char[m_bufferLength]();

    // store filename & return success
    m_filename = filename;
    return true;
}

bool FastqReader::readNext(Fastq *entry) {

    // fail if unopened file
    if ( !isOpen() ) {
        m_errorString = "cannot read from unopened reader";
        return false;
    }

    // sanity checks
    assert(entry);
    assert(m_buffer);

    // read header
    char* result;

    result = STREAM_GETS;
    if ( STREAM_EOF ) {
        m_errorString = "could not read full FASTQ entry from file: ";
        m_errorString.append(m_filename);
        return false;
    }

    if ( m_buffer[0] != '@' ) {
        m_errorString = "malformed FASTQ entry - expected '@' in header, instead found: ";
        m_errorString.append(1, m_buffer[0]);
        return false;
    }
    chomp(m_buffer);
    entry->Header.assign(m_buffer);

    // read bases
    ostringstream sb("");
    while ( true ) {
        const char c = STREAM_GETC;
        STREAM_UNGETC(c);
        if ( c == '+' || STREAM_EOF )
            break;
        result = STREAM_GETS;
        chomp(m_buffer);
        sb << m_buffer;
    }
    entry->Bases.assign(m_buffer);
    const size_t numBases = entry->Bases.length();

    // read qualities
    result = STREAM_GETS;
    sb.str("");
    size_t numQualities = 0;
    while ( true ) {
        const char c = STREAM_GETC;
        STREAM_UNGETC(c);
        if ( STREAM_EOF )
            break;
        result = STREAM_GETS;
        chomp(m_buffer);
        numQualities += strlen(m_buffer);
        sb << m_buffer;
        if ( numQualities >= numBases )
            break;
    }
    entry->Qualities.assign(m_buffer);

    // sanity check
    if ( entry->Qualities.length() != entry->Bases.length() ) {
        m_errorString = "malformed FASTQ entry - the number of qualities does not match the number of bases";
        return false;
    }

    // return success
    return true;
}
