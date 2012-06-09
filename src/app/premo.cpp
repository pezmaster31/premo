// ***************************************************************************
// premo.cpp (c) 2012 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 7 June 2012 (DB)
// ---------------------------------------------------------------------------
// Main Premo workhorse
// ***************************************************************************

#include "premo.h"
#include "batch.h"
#include "fastqreader.h"
#include "options.h"
#include <iostream>
#include <sstream>
using namespace std;

// ------------------------
// static utility methods
// ------------------------

static inline
bool endsWith(const string& str, const string& query) {
    return ( str.find(query) == (str.length() - query.length()) );
}

// ----------------------
// Premo implementation
// ----------------------

Premo::Premo(const PremoSettings& settings)
    : m_settings(settings)
    , m_isFinished(false)
{
    // if no output file provided, use stdout for JSON output
    if ( !m_settings.HasOutputFilename )
        m_settings.OutputFilename = Options::StandardOut();

    // do any upfront stats prep ??
    initializeStats();
}

Premo::~Premo(void) { }

string Premo::errorString(void) const {
    return m_errorString;
}

void Premo::initializeStats(void) {

}

bool Premo::openInputFiles(void) {

    // open FASTQ input files for reading
    bool openedOk = true;
    openedOk &= m_reader1.open(m_settings.FastqFilename1);
    openedOk &= m_reader2.open(m_settings.FastqFilename2);

    // check for failures
    if ( !openedOk ) {

        // build error string
        m_errorString = "premo ERROR: could not open the following input FASTQ file(s):\n";
        if ( !m_reader1.isOpen() )
            m_errorString.append(m_settings.FastqFilename1);
        if ( !m_reader2.isOpen() )
            m_errorString.append(m_settings.FastqFilename2);

        // return failure
        return false;
    }

    // otherwise, opened OK
    else return true;
}

bool Premo::outputResults(void) {

    return true;
}

bool Premo::run(void) {

    // check that settings are valid
    if ( !validateSettings() )
        return false;

    // open our input files for reading (FastqReader dtor closes FASTQ file)
    if ( !openInputFiles() )
        return false;

    // main loop - batch processing
    int batchNumber = 0;
    while ( !m_isFinished ) {

        cerr << "Running batch: " << batchNumber << endl;

        Batch batch(batchNumber, m_settings, &m_reader1, &m_reader2);
        if ( batch.run() ) {

            // fetch results & update

            // check for end condition


        } else {
            m_errorString = "premo ERROR: batch failed - \n";
            m_errorString.append(batch.errorString());
            return false;
        }

        ++batchNumber;

        // FOR NOW
        if ( batchNumber == 10 )
            m_isFinished = true;
    }

    // output results
    if ( !outputResults() )
        return false;

    // return success
    return true;
}

bool Premo::validateSettings(void) {

    // -------------------------------
    // check for required parameters
    // -------------------------------

    stringstream missing("");
    bool hasMissing = false;

    if ( !m_settings.HasFastqFilename1 || m_settings.FastqFilename1.empty() ) {
        missing << endl << "\t-fq1 (FASTQ filename)";
        hasMissing = true;
    }

    if ( !m_settings.HasFastqFilename2 || m_settings.FastqFilename2.empty() ) {
        missing << endl << "\t-fq2 (FASTQ filename)";
        hasMissing = true;
    }

    if ( !m_settings.HasMosaikPath || m_settings.MosaikPath.empty() ) {
        missing << endl << "\t-mosaik (path/to/Mosaik/bin)";
        hasMissing = true;
    } else {

        // append dir separator if missing
        if ( !endsWith(m_settings.MosaikPath, "/") )
            m_settings.MosaikPath.append("/");
    }

    if ( !m_settings.HasReferenceArchiveFilename || m_settings.ReferenceArchiveFilename.empty() ) {
        missing << endl << "\t-ref (Mosaik reference archive)";
        hasMissing = true;
    }

    if ( !m_settings.HasScratchPath || m_settings.ScratchPath.empty() ) {
        missing << endl << "\t-tmp (scratch directory for generated files)";
        hasMissing = true;
    } else {

        // append dir separator if missing
        if ( !endsWith(m_settings.ScratchPath, "/") )
            m_settings.ScratchPath.append("/");
    }

    // -----------------------------------------
    // check other parameters for valid ranges
    // -----------------------------------------

    stringstream invalid("");
    bool hasInvalid = false;

    if ( m_settings.HasBatchSize && m_settings.BatchSize == 0 ) {
        invalid << endl << "\t-n cannot be zero";
        hasInvalid = true;
    }

    if ( m_settings.HasDeltaMean && m_settings.DeltaMean <= 0.0 ) {
        invalid << endl << "\t-delta-mean must be a positive, non-zero value";
        hasInvalid = true;
    }

    if ( m_settings.HasDeltaStdDev && m_settings.DeltaStdDev <= 0.0 ) {
        invalid << endl << "\t-delta-sd must be a positive, non-zero value";
        hasInvalid = true;
    }

    if ( m_settings.HasActSlope && m_settings.ActSlope <= 0.0 ) {
        invalid << endl << "\t-act-slope must be a positive, non-zero value";
        hasInvalid = true;
    }

    if ( m_settings.HasBandwidth && m_settings.Bandwidth <= 0.0 ) {
        invalid << endl << "\t-bw must be a positive, non-zero value";
        hasInvalid = true;
    }

    if ( m_settings.HasMhp && m_settings.Mhp == 0 ) {
        invalid << endl << "\t-mhp cannot be zero";
        hasInvalid = true;
    }

    if ( m_settings.HasMmp && (m_settings.Mmp < 0.0 || m_settings.Mmp > 1.0) ) {
        invalid << endl << "\t-mmp must be in the range [0.0 - 1.0]";
        hasInvalid = true;
    }

    // ---------------------------------------------------------------
    // set error string if anything missing/invalid
    // ---------------------------------------------------------------

    m_errorString.clear();

    if ( hasMissing ) {
        m_errorString.append("\nthe following parameters are missing:");
        m_errorString.append(missing.str());
    }

    if ( hasInvalid ) {
        m_errorString.append("\nthe following parameters are invalid:");
        m_errorString.append(invalid.str());
    }

    // --------------------------
    // return validation status
    // --------------------------
    return !( hasMissing || hasInvalid );
}
