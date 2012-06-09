// ***************************************************************************
// batch.h (c) 2012 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 8 June 2012 (DB)
// ---------------------------------------------------------------------------
// Premo batch
// ***************************************************************************

#ifndef BATCH_H
#define BATCH_H

#include <string>
class FastqReader;
class PremoSettings;

class Batch {

    // ctor & dtor
    public:
        Batch(const int batchNumber,
              const PremoSettings& settings,
              FastqReader* reader1,
              FastqReader* reader2);
        ~Batch(void);

    // Batch interface
    public:
        std::string errorString(void) const;
        bool run(void);

    // internal methods
    private:


    // data members
    private:

        FastqReader* m_reader1; // copy, not owned
        FastqReader* m_reader2; // copy, not owned

        bool m_autoDeleteGeneratedFiles;
        std::string m_generatedFastq1;
        std::string m_generatedFastq2;
        std::string m_generatedReadArchive;
        std::string m_generatedBam;

        std::string m_errorString;
};

#endif // BATCH_H
