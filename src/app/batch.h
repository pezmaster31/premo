// ***************************************************************************
// batch.h (c) 2012 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 27 June 2012 (DB)
// ---------------------------------------------------------------------------
// Premo batch
// ***************************************************************************

#ifndef BATCH_H
#define BATCH_H

#include "result.h"
#include <string>
#include <vector>
class FastqReader;
class PremoSettings;

class Batch {

    public:
        enum RunStatus { Normal = 0  // batch processed normally, result available
                       , HitEOF      // batch hit EOF before processing settings.BatchSize reads,
                                     // result available but any further batch runs will return NoData
                       , NoData      // empty file (or starting from EOF) - NO result available
                       , Error       // any other error case - NO result available
                       };

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
        Result result(void) const;
        RunStatus run(void);

    // internal methods
    private:
        RunStatus generateTempFastqFiles(void);
        RunStatus parseAlignmentFile(void);
        RunStatus runMosaikAligner(void);
        RunStatus runMosaikBuild(void);
        RunStatus runMosaikPipeline(void);

    // data members
    private:

        // copies from main Premo app, not owned
        PremoSettings* m_settings;
        FastqReader*   m_reader1;
        FastqReader*   m_reader2;

        // store all possible generated filenames, for proper cleanup
        std::string m_generatedFastq1;
        std::string m_generatedFastq2;
        std::string m_generatedReadArchive;
        std::string m_generatedBamStub;
        std::string m_generatedBam;
        std::string m_generatedMosaikLog;
        std::string m_generatedMultipleBam;
        std::string m_generatedSpecialBam;
        std::string m_generatedStatFile;

        // our main result
        Result m_result;

        // error reporting
        std::string m_errorString;
};

#endif // BATCH_H
