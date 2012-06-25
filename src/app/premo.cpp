// ***************************************************************************
// premo.cpp (c) 2012 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 24 June 2012 (DB)
// ---------------------------------------------------------------------------
// Main Premo workhorse
// ***************************************************************************

#include "batch.h"
#include "options.h"
#include "premo.h"
#include "stats.h"
#include "jsoncpp/json_value.h"
#include "jsoncpp/json_writer.h"
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
using namespace std;

// ------------------------
// static utility methods
// ------------------------

static
Json::Value containerStats(const vector<int>& container) {

    Json::Value result(Json::objectValue);
    result["count"] = static_cast<Json::UInt>(container.size());

    if ( !container.empty() ) {

        vector<int> c = container;
        sort(c.begin(), c.end());

        const Quartiles quartiles = calculateQuartiles(c);
        result["median"] = quartiles.Q2;
        result["Q1"] = quartiles.Q1;
        result["Q3"] = quartiles.Q3;
    }

    return result;
}

static
Json::Value resultToJson(const Result& result) {
    Json::Value json(Json::objectValue);
    json["fragment length"] = containerStats(result.FragmentLengths);
    json["read length"]     = containerStats(result.ReadLengths);
    return json;
}

static
bool isConverged(const vector<int>& previous,
                 const vector<int>& current,
                 const double cutoffDelta)
{
    // sort (a copy of) input containers (req'd for median calculation)
    vector<int> currentValues  = current;
    vector<int> previousValues = previous;
    sort(currentValues.begin(),  currentValues.end());
    sort(previousValues.begin(), previousValues.end());

    // calculate medians
    const double currentMedian  = calculateMedian(currentValues);
    const double previousMedian = calculateMedian(previousValues);

    // calculate difference between previous & current values
    const double diff = fabs( currentMedian - previousMedian );

    // calculate delta (ratio) from old values
    const double observedDelta = ( diff / previousMedian );

    // return whether observed delta is below requested cutoff
    return ( observedDelta <= cutoffDelta );
}

static
bool checkFinished(const Result& previousResult,
                   const Result& currentResult,
                   const PremoSettings& settings)
{
    return isConverged(previousResult.FragmentLengths,
                       currentResult.FragmentLengths,
                       settings.DeltaFragmentLength)  &&
           isConverged(previousResult.ReadLengths,
                       currentResult.ReadLengths,
                       settings.DeltaReadLength);
}

template<typename T>
void append(std::vector<T>& dest, const std::vector<T>& source) {
    dest.insert(dest.end(), source.begin(), source.end());
}

static inline
bool endsWith(const string& str, const string& query) {
    return ( str.find_last_of(query) == (str.length() - query.length()) );
}

// ----------------------
// Premo implementation
// ----------------------

Premo::Premo(const PremoSettings& settings)
    : m_settings(settings)
    , m_isFinished(false)
{ }

Premo::~Premo(void) { }

string Premo::errorString(void) const {
    return m_errorString;
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
        if ( !m_reader1.isOpen() ) m_errorString.append(m_settings.FastqFilename1);
        if ( !m_reader2.isOpen() ) m_errorString.append(m_settings.FastqFilename2);
        // return failure
        return false;
    }

    if ( m_settings.IsVerbose )
        cerr << "input FASTQ files opened OK" << endl;

    // otherwise, opened OK
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

        if ( m_settings.IsVerbose )
            cerr << "running batch: " << batchNumber << endl;

        // run batch
        Batch batch(batchNumber, &m_settings, &m_reader1, &m_reader2);
        if ( !batch.run() ) {

            // if batch failed, set error & return failure
            stringstream s("premo ERROR: batch ");
            s << batchNumber;
            s << " failed - " << endl;
            s << batch.errorString();
            m_errorString = s.str();
            return false;
        }

        // store batch results
        const Result result = batch.result();
        m_batchResults.push_back( result );

        // store previous result before adding batch data to "current" result
        const Result previousResult = m_currentResult;

        // add batch's data to current, overall result
        append(m_currentResult.FragmentLengths, result.FragmentLengths);
        append(m_currentResult.ReadLengths,     result.ReadLengths);

        // if not first batch, check to see if we're done
        if ( batchNumber > 0 )
            m_isFinished = checkFinished(previousResult, m_currentResult, m_settings);

        // increment our batch counter
        ++batchNumber;
    }

    // output results
    if ( !writeOutput() )
        return false;

    // if we get here, return success
    return true;
}

bool Premo::validateSettings(void) {

    // -------------------------------
    // check for required parameters
    // -------------------------------

    stringstream missing("");
    bool hasMissing = false;

    if ( !m_settings.HasAnnPeFilename || m_settings.AnnPeFilename.empty() ) {
        missing << endl << "\t-annpe (paired-end neural network filename)";
        hasMissing = true;
    }

    if ( !m_settings.HasAnnSeFilename || m_settings.AnnSeFilename.empty() ) {
        missing << endl << "\t-annse (single-end neural network filename)";
        hasMissing = true;
    }

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

        // append dir separator if missing from path
        if ( !endsWith(m_settings.MosaikPath, "/") )
            m_settings.MosaikPath.append("/");
    }

    if ( !m_settings.HasOutputFilename || m_settings.OutputFilename.empty() ) {
        missing << endl << "\t-out (output filename)";
        hasMissing = true;
    }

    if ( !m_settings.HasReferenceFilename || m_settings.ReferenceFilename.empty() ) {
        missing << endl << "\t-ref (Mosaik reference archive)";
        hasMissing = true;
    }

    if ( !m_settings.HasScratchPath || m_settings.ScratchPath.empty() ) {
        missing << endl << "\t-tmp (scratch directory for generated files)";
        hasMissing = true;
    } else {

        // append dir separator if missing from path
        if ( !endsWith(m_settings.ScratchPath, "/") )
            m_settings.ScratchPath.append("/");
    }

    if ( !m_settings.HasSeqTech || m_settings.SeqTech.empty() ) {
        missing << endl << "\t-st (sequencing technology)";
        hasMissing = true;
    }

    // -----------------------------------------
    // check other parameters for valid ranges
    // -----------------------------------------

    stringstream invalid("");
    bool hasInvalid = false;


    if ( m_settings.HasActSlope && m_settings.ActSlope <= 0.0 ) {
        invalid << endl << "\t-act-slope must be a positive, non-zero value";
        hasInvalid = true;
    }

    if ( m_settings.HasBatchSize && m_settings.BatchSize == 0 ) {
        invalid << endl << "\t-n cannot be zero";
        hasInvalid = true;
    }

    if ( m_settings.HasBwMultiplier && m_settings.BwMultiplier <= 0.0 ) {
        invalid << endl << "\t-bwm must be a positive, non-zero value";
        hasInvalid = true;
    }

    if ( m_settings.HasDeltaFragmentLength && m_settings.DeltaFragmentLength <= 0.0 ) {
        invalid << endl << "\t-delta-fl must be a positive, non-zero value";
        hasInvalid = true;
    }

    if ( m_settings.HasDeltaReadLength && m_settings.DeltaReadLength <= 0.0 ) {
        invalid << endl << "\t-delta-rl must be a positive, non-zero value";
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

    const bool settingsOk = !( hasMissing || hasInvalid );

    if ( settingsOk && m_settings.IsVerbose )
        cerr << "command-line settings OK" << endl;

    return settingsOk;
}

bool Premo::writeOutput(void) {

    Json::Value root(Json::objectValue);

    // ------------------------------
    // store top-level results
    // ------------------------------

    root["overall result"] = resultToJson(m_currentResult);

    // -------------------------
    // store per-batch results
    // -------------------------

    Json::Value batches(Json::arrayValue);
    vector<Result>::const_iterator batchIter = m_batchResults.begin();
    vector<Result>::const_iterator batchEnd  = m_batchResults.end();
    for ( ; batchIter != batchEnd; ++batchIter )
        batches.append( resultToJson(*batchIter) );

    root["batch results"] = batches;

    // ------------------------------
    // store settings used
    // ------------------------------

    Json::Value settings(Json::objectValue);
    settings["act intercept"]         = m_settings.ActIntercept;
    settings["act slope"]             = m_settings.ActSlope;
    settings["bandwidth multiplier"]  = m_settings.BwMultiplier;
    settings["batch size"]            = m_settings.BatchSize;
    settings["delta fragment length"] = m_settings.DeltaFragmentLength;
    settings["delta read length"]     = m_settings.DeltaReadLength;
    settings["mhp"]                   = m_settings.Mhp;
    settings["mmp"]                   = m_settings.Mmp;
    settings["seq tech"]              = m_settings.SeqTech;

    root["settings"] = settings;

    // -------------------------------
    // generate Mosaik parameter set
    // -------------------------------

    vector<int> fragmentLengths = m_currentResult.FragmentLengths;
    vector<int> readLengths     = m_currentResult.ReadLengths;
    sort(fragmentLengths.begin(), fragmentLengths.end());
    sort(readLengths.begin(),     readLengths.end());

    const double fragLengthMedian = calculateMedian(fragmentLengths);
    const double readLengthMedian = calculateMedian(readLengths);

    Json::Value mosaikAlignerParameters(Json::objectValue);
    mosaikAlignerParameters["-act"] = (m_settings.ActSlope * readLengthMedian) + m_settings.ActIntercept;
    mosaikAlignerParameters["-bw"]  = ceil( m_settings.BwMultiplier * readLengthMedian );
    mosaikAlignerParameters["-ls"]  = fragLengthMedian;
    mosaikAlignerParameters["-mhp"] = m_settings.Mhp;
    mosaikAlignerParameters["-mmp"] = m_settings.Mmp;
    mosaikAlignerParameters["-st"]  = m_settings.SeqTech;

    Json::Value parameters(Json::objectValue);
    parameters["MosaikAligner"] = mosaikAlignerParameters;

    root["parameters"] = parameters;

    // ---------------------------
    // write JSON to output file
    // ---------------------------

    // open stream on file
    ofstream outFile(m_settings.OutputFilename.c_str());
    if ( !outFile ) {
        m_errorString = "premo ERROR: could not open final output file: ";
        m_errorString.append(m_settings.OutputFilename);
        return false;
    }

    // write "pretty-printed" JSON contents to file
    Json::StyledStreamWriter writer("  ");
    writer.write(outFile, root);

    if ( m_settings.IsVerbose )
        cerr << "results written OK" << endl;

    // clean up & return success
    outFile.close();
    return true;
}
