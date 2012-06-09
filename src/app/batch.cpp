// ***************************************************************************
// batch.cpp (c) 2012 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 8 June 2012 (DB)
// ---------------------------------------------------------------------------
// Premo batch
// ***************************************************************************

#include "batch.h"
#include "fastq.h"
#include "fastqreader.h"
#include "premo_settings.h"
#include <cstdio>
#include <iostream>
#include <sstream>
using namespace std;

Batch::Batch(const int batchNumber,
             const PremoSettings& settings,
             FastqReader* reader1,
             FastqReader* reader2)
    : m_reader1(reader1)
    , m_reader2(reader2)
    , m_autoDeleteGeneratedFiles()
{
    // ----------------------------
    // set up generated filenames (scratchpath should already end in '/')
    // ----------------------------

    stringstream s;

    // mate1 FASTQ
    s.str(settings.ScratchPath);
    s << "premo_batch" << batchNumber << "_mate1.fq";
    m_generatedFastq1 = s.str();

    // mate2 FASTQ
    s.str(settings.ScratchPath);
    s << "premo_batch" << batchNumber << "_mate2.fq";
    m_generatedFastq2 = s.str();

    // Mosaik read archive
    s.str(settings.ScratchPath);
    s << "premo_batch" << batchNumber << "_reads.mkb";
    m_generatedReadArchive = s.str();

    // Mosaik alignment
    s.str(settings.ScratchPath);
    s << "premo_batch" << batchNumber << "_aligned.bam";
    m_generatedBam = s.str();
}

Batch::~Batch(void) {

    // delete any generated files
    remove(m_generatedFastq1.c_str());
    remove(m_generatedFastq2.c_str());
    remove(m_generatedReadArchive.c_str());
    remove(m_generatedBam.c_str());
}

string Batch::errorString(void) const {
    return m_errorString;
}

bool Batch::run(void) {

    Fastq f1;

    if ( m_reader1->readNext(&f1) ) {

        cerr << "Read FASTQ entry: "       << endl
             << "Header: " << f1.Header    << endl
             << "Bases : " << f1.Bases     << endl
             << "Quals : " << f1.Qualities << endl;
        return true;
    } else
        return false;


    // create writer1, writer2 for temp output
    // open tmp files in scratch dir
    //
    // Fastq f1, f2;
    // for ( int i = 0; i < BATCH_SIZE; ++i ) {
    //
    //     // read FASTQ entries
    //     if ( !m_reader1->readNext(&f1) ) {
    //         // report error
    //         // return failure
    //     }
    //     if ( !m_reader2->readNext(&f2) ) {
    //         // report error
    //         // return failure
    //     }
    //
    //     // write to tmp file
    //     writer1.write(f1);
    //     writer2.write(f2);
    // }

    // runMosaik(tmp1, tmp2);

    // MosaikBuild, MosaikAlign, parse BAM, calc stats

    return true;
}

// bool runMosaik(string& fn1, string& fn2) {
//     string readArchiveFilename = runMosaikBuild(fn1, fn2);
//     runMosaikAlign(readArchiveFilename);
//     return status
// }

// runMosaikBuild(string& fn1, string& fn2) {
//    string readArchiveFilename = "something";
//    system("@settings.MosaikPath/MosaikBuild... -q @fn1 -q2 @fn2 .... -out @readArchiveFilename");
//    return readArchiveFilename;
// }
