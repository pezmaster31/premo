// ***************************************************************************
// premo.cpp (c) 2012 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 7 June 2012 (DB)
// ---------------------------------------------------------------------------
// Main Premo workhorse
// ***************************************************************************

#include "premo.h"
#include "options.h"
#include <sstream>
using namespace std;

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

}

bool Premo::outputResults(void) {

}

bool Premo::run(void) {

    // check that settings are valid
    if ( !validateSettings() )
        return false;

    // open our input files for reading (FastqReader dtor closes FASTQ file)
    if ( !openInputFiles() )
        return false;

    // initialize anything else ??

    // main loop - batch processing
    while ( !m_isFinished ) {

        // process batch
        // fire off batch                    -> reads N entries, write N entries to temp FASTQ file, MosaikBuild, MosaikAlign, parse BAM, calc stats
        // update stats with batch results
        // batch manages temp file destruction
        // check for end condition

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

    stringstream missing("the following parameters are missing: ");
    bool hasMissing = false;

    if ( !m_settings.HasFastqFilename1 || m_settings.FastqFilename1.empty() ) {
        missing << endl << "-fq1 (FASTQ filename)";
        hasMissing = true;
    }

    if ( !m_settings.HasFastqFilename2 || m_settings.FastqFilename2.empty() ) {
        missing << endl << "-fq2 (FASTQ filename)";
        hasMissing = true;
    }

    if ( !m_settings.HasMosaikPath || m_settings.MosaikPath.empty() ) {
        missing << endl << "-mosaik (path/to/Mosaik/bin)";
        hasMissing = true;
    }

    if ( !m_settings.HasReferenceArchiveFilename || m_settings.ReferenceArchiveFilename.empty() ) {
        missing << endl << "-ref (Mosaik reference archive)";
        hasMissing = true;
    }

    if ( !m_settings.HasScratchPath || m_settings.ScratchPath.empty() ) {
        missing << endl << "-tmp (scratch directory for generated files)";
        hasMissing = true;
    }

    // -----------------------------------------
    // check other parameters for valid ranges
    // -----------------------------------------

    stringstream invalid("the following parameters are invalid: ");
    bool hasInvalid = false;

    if ( m_settings.HasBatchSize && m_settings.BatchSize == 0 ) {
        invalid << endl << "-n cannot be zero";
        hasInvalid = true;
    }

    if ( m_settings.HasDeltaMean && m_settings.DeltaMean <= 0.0 ) {
        invalid << endl << "-delta-mean must be a positive, non-zero value";
        hasInvalid = true;
    }

    if ( m_settings.HasDeltaStdDev && m_settings.DeltaStdDev <= 0.0 ) {
        invalid << endl << "-delta-sd must be a positive, non-zero value";
        hasInvalid = true;
    }

    if ( m_settings.HasActSlope && m_settings.ActSlope <= 0.0 ) {
        invalid << endl << "-act-slope must be a positive, non-zero value";
        hasInvalid = true;
    }

    if ( m_settings.HasBandwidth && m_settings.Bandwidth <= 0.0 ) {
        invalid << endl << "-bw must be a positive, non-zero value";
        hasInvalid = true;
    }

    if ( m_settings.HasMhp && m_settings.Mhp == 0 ) {
        invalid << endl << "-mhp cannot be zero";
        hasInvalid = true;
    }

    if ( m_settings.HasMmp && (m_settings.Mmp < 0.0 || m_settings.Mmp > 1.0) ) {
        invalid << endl << "-mmp must be in the range [0.0 - 1.0]";
        hasInvalid = true;
    }

    // ---------------------------------------------------------------
    // set error string & return failure if anything missing/invalid
    // ---------------------------------------------------------------

    if ( hasMissing ) {
        m_errorString = missing.str();
        return false;
    }

    if ( hasInvalid ) {
        m_errorString = invalid.str();
        return false;
    }

    // -----------------------------
    // otherwise, everything is OK
    // -----------------------------
    return true;
}

