// ***************************************************************************
// premo.h (c) 2012 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 6 June 2012 (DB)
// ---------------------------------------------------------------------------
// Main Premo workhorse
// ***************************************************************************

#ifndef PREMO_H
#define PREMO_H

#include "premo_settings.h"
#include <string>

class Premo {

    // ctor & dtor
    public:
        Premo(const PremoSettings& settings);
        ~Premo(void);

    // Premo interface
    public:
        std::string errorString(void) const;
        bool run(void);

    // internal methods
    private:
        void initializeStats(void);
        bool openInputFiles(void);
        bool outputResults(void);
        bool validateSettings(void);

    // data members
    private:
        PremoSettings m_settings;
        bool m_isFinished;




        std::string m_errorString;
};

#endif // PREMO_H
