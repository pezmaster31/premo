// ***************************************************************************
// batch.cpp (c) 2012 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 9 June 2012 (DB)
// ---------------------------------------------------------------------------
// Premo batch
// ***************************************************************************

#include "batch.h"
#include "fastq.h"
#include "fastqreader.h"
#include "fastqwriter.h"
#include "premo_settings.h"
#include "bamtools/api/BamReader.h"
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <sstream>
using namespace std;

// -----------------------------
// BatchResults implementation
// -----------------------------

BatchResults::BatchResults(void)
    : NumReads(0)
    , MeanReadLength(0.0)
    , MeanFragmentLength(0.0)
    , StdDevFragmentLength(0.0)
{ }

BatchResults::BatchResults(const BatchResults& other)
    : NumReads(other.NumReads)
    , MeanReadLength(other.MeanReadLength)
    , MeanFragmentLength(other.MeanFragmentLength)
    , StdDevFragmentLength(other.StdDevFragmentLength)
{ }

BatchResults::~BatchResults(void) { }

// ----------------------
// Batch implementation
// ----------------------

Batch::Batch(const int batchNumber,
             PremoSettings* settings,
             FastqReader* reader1,
             FastqReader* reader2)
    : m_settings(settings)
    , m_reader1(reader1)
    , m_reader2(reader2)
{
    // ----------------------------
    // set up generated filenames (scratchpath should already end in '/')
    // ----------------------------

    stringstream s;

    // mate1 FASTQ
    s.str(m_settings->ScratchPath);
    s << "premo_batch" << batchNumber << "_mate1.fq";
    m_generatedFastq1 = s.str();

    // mate2 FASTQ
    s.str(m_settings->ScratchPath);
    s << "premo_batch" << batchNumber << "_mate2.fq";
    m_generatedFastq2 = s.str();

    // Mosaik read archive
    s.str(m_settings->ScratchPath);
    s << "premo_batch" << batchNumber << "_reads.mkb";
    m_generatedReadArchive = s.str();

    // Mosaik alignment
    s.str(m_settings->ScratchPath);
    s << "premo_batch" << batchNumber << "_aligned.bam";
    m_generatedBam = s.str();
}

Batch::~Batch(void) {

    // auto-delete any generated files (unless requested otherwise)
    if ( !m_settings->IsKeepGeneratedFiles ) {
        remove(m_generatedFastq1.c_str());
        remove(m_generatedFastq2.c_str());
        remove(m_generatedReadArchive.c_str());
        remove(m_generatedBam.c_str());
    }
}

string Batch::errorString(void) const {
    return m_errorString;
}

bool Batch::generateTempFastqFiles(void) {

    // ------------------------------
    // open temp FASTQ output files
    // ------------------------------

    FastqWriter writer1;
    FastqWriter writer2;

    bool openedOk = true;
    openedOk &= writer1.open(m_generatedFastq1);
    openedOk &= writer2.open(m_generatedFastq2);

    // check for failures
    if ( !openedOk ) {

        // build error string
        m_errorString = "premo ERROR: could not create the following temp FASTQ file(s):\n";
        if ( !writer1.isOpen() )
            m_errorString.append(m_generatedFastq1);
        if ( !writer2.isOpen() )
            m_errorString.append(m_generatedFastq2);

        // return failure
        return false;
    }

    // -----------------------------------------------------------
    // copy next batch of FASTQ entries from input to temp files
    // -----------------------------------------------------------

    Fastq f1;
    Fastq f2;

    // iterate over requested number of entries
    for ( size_t i = 0; i < m_settings->BatchSize; ++i ) {

        // read from input FASTQ files
        if ( !m_reader1->readNext(&f1) ) {
            m_errorString = "premo ERROR: could not read from input FASTQ file: ";
            m_errorString.append(m_reader1->filename());
            return false;
        }
        if ( !m_reader2->readNext(&f2) ) {
            m_errorString = "premo ERROR: could not read from input FASTQ file: ";
            m_errorString.append(m_reader2->filename());
            return false;
        }

        // write to temp FASTQ files
        if ( !writer1.write(&f1) ) {
            m_errorString = "premo ERROR: could not write to temp FASTQ file: ";
            m_errorString.append(writer1.filename());
            return false;
        }
        if ( !writer2.write(&f2) ) {
            m_errorString = "premo ERROR: could not write to temp FASTQ file: ";
            m_errorString.append(writer2.filename());
            return false;
        }
    }

    // if we get here, all should be OK
    // cleanup & return success
    writer1.close();
    writer2.close();
    return true;
}

bool Batch::parseAlignmentForStats(void) {

    // open reader on new BAM alignment file
    BamTools::BamReader reader;
    if ( !reader.Open(m_generatedBam) ) {
        m_errorString = "premo ERROR: could not open generated BAM file: ";
        m_errorString.append(m_generatedBam);
        m_errorString.append(" to calculate stats");
        return false;
    }

    uint32_t numReads = 0;
    uint64_t numBases = 0;

    // plow through alignments
    BamTools::BamAlignment a;
    while ( reader.GetNextAlignmentCore(a) ) {

        ++numReads;
        numBases += a.Length;




    }
    reader.Close();

    // store stats for reporting later
    const double meanReadLength = ( numBases / numReads );

    // if we get here, all should be OK
    return true;
}

bool Batch::run(void) {

    // generate temp files
    if ( !generateTempFastqFiles() )
        return false;

    // run mosaik
    if ( !runMosaikPipeline() )
        return false;

    // parse BAM for stats
    if ( !parseAlignmentForStats() )
        return false;

    // if we get here, all should be OK
    return true;
}

bool Batch::runMosaikAlign(void) {

    // setup MosaikAlign command line
    stringstream commandStream("");
    commandStream << m_settings->MosaikPath << "MosaikAlign"   // (mosaikpath should already end in '/')
                  << " -ia "  << m_settings->ReferenceArchiveFilename
                  << " -in "  << m_generatedReadArchive
                  << " -out " << m_generatedBam
                  << " -mhp " << m_settings->Mhp
                  << " -mmp " << m_settings->Mmp;
    if ( !m_settings->IsVerbose )
        commandStream << " -quiet";

    // run MosaikAlign
    const string command = commandStream.str();
    const int result = system(command.c_str());
    return ( result == 0 );
}

bool Batch::runMosaikBuild(void) {

    // setup MosaikBuild command line
    stringstream commandStream("");
    commandStream << m_settings->MosaikPath << "MosaikBuild"   // (mosaikpath should already end in '/')
                  << " -q "   << m_generatedFastq1
                  << " -q2 "  << m_generatedFastq2
                  << " -out " << m_generatedReadArchive;
    if ( !m_settings->IsVerbose )
        commandStream << " -quiet";

    // run MosaikBuild
    const string command = commandStream.str();
    const int result = system(command.c_str());
    return ( result == 0 );
}

bool Batch::runMosaikPipeline(void) {

    // build mosaik archives for new batch FASTQ files
    if ( !runMosaikBuild() )
        return false;

    // align batch
    if ( !runMosaikAlign() )
        return false;

    // if we get here, all should be OK
    return true;
}
