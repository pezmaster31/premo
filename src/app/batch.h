// ***************************************************************************
// batch.h (c) 2012 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 9 June 2012 (DB)
// ---------------------------------------------------------------------------
// Premo batch
// ***************************************************************************

#ifndef BATCH_H
#define BATCH_H

#include <string>
class FastqReader;
class PremoSettings;

struct BatchResults {

    // data members
    uint32_t NumReads;
    double   MeanReadLength;
    double   MeanFragmentLength;
    double   StdDevFragmentLength;

    // ctors & dtor
    BatchResults(void);
    BatchResults(const BatchResults& other);
    ~BatchResults(void);
};

class Batch {

    // ctor & dtor
    public:
        Batch(const int batchNumber,
              PremoSettings* settings,
              FastqReader* reader1,
              FastqReader* reader2);
        ~Batch(void);

    // Batch interface
    public:
        std::string errorString(void) const;
        bool run(void);

        // results() const;

    // internal methods
    private:
        bool generateTempFastqFiles(void);
        bool parseAlignmentForStats(void);
        bool runMosaikAlign(void);
        bool runMosaikBuild(void);
        bool runMosaikPipeline(void);

    // data members
    private:

        PremoSettings* m_settings; // copy, not owned
        FastqReader*   m_reader1;  // copy, not owned
        FastqReader*   m_reader2;  // copy, not owned

        std::string m_generatedFastq1;
        std::string m_generatedFastq2;
        std::string m_generatedReadArchive;
        std::string m_generatedBam;

        std::string m_errorString;
};

#endif // BATCH_H
