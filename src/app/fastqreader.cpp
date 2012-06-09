// ***************************************************************************
// fastqreader.cpp (c) 2012 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 8 June 2012 (DB)
// ---------------------------------------------------------------------------
// FASTQ file reader
// ***************************************************************************

#include "fastqreader.h"
#include "fastq.h"
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

FastqReader::FastqReader(void)
    : m_stream(0)
    , m_buffer(0)
    , m_bufferLength(0)
{ }

FastqReader::~FastqReader(void) {
    close();
}

void FastqReader::close(void) {

    // close file stream
    if ( m_stream ) {
        fclose(m_stream);
        m_stream = 0;
    }

    // clean up allocated memory
    if ( m_buffer ) {
        delete[] m_buffer;
        m_buffer = 0;
        m_bufferLength = 0;
    }
}

string FastqReader::errorString(void) const {
    return m_errorString;
}

bool FastqReader::isOpen(void) const {
    return ( m_stream != 0 );
}

bool FastqReader::open(const string& filename) {

    // ensure clean slate
    close();

    // attempt to open
    m_stream = fopen(filename.c_str(), "r");
    if ( m_stream == 0 ) {

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
    if ( m_stream == 0 ) {
        m_errorString = "cannot read from unopened reader";
        return false;
    }

    // sanity checks
    assert(entry);
    assert(m_buffer);

    // read header
    char* result;

    result = fgets(m_buffer, m_bufferLength, m_stream);
    if ( feof(m_stream) ) {
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
        const char c = fgetc(m_stream);
        ungetc(c, m_stream);
        if ( c == '+' || feof(m_stream) )
            break;
        result = fgets(m_buffer, m_bufferLength, m_stream);
        chomp(m_buffer);
        sb << m_buffer;
    }
    entry->Bases.assign(m_buffer);
    const size_t numBases = entry->Bases.length();

    // read qualities
    result = fgets(m_buffer, m_bufferLength, m_stream);
    sb.str("");
    size_t numQualities = 0;
    while ( true ) {
        const char c = fgetc(m_stream);
        ungetc(c, m_stream);
        if ( feof(m_stream) )
            break;
        result = fgets(m_buffer, m_bufferLength, m_stream);
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
