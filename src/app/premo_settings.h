// ***************************************************************************
// premo_settings.h (c) 2012 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 7 June 2012 (DB)
// ---------------------------------------------------------------------------
// Premo app settings
// ***************************************************************************

#ifndef PREMO_SETTINGS_H
#define PREMO_SETTINGS_H

#include "bamtools/api/bamtools_global.h"  // for uint32_t type
#include <string>

// default numerical values (can override any/all from cmdline)
const uint32_t DEFAULT_ACTINTERCEPT = 13;
const double   DEFAULT_ACTSLOPE     = 0.2;
const double   DEFAULT_BANDWIDTH    = 2.5;
const uint32_t DEFAULT_BATCHSIZE    = 1000;
const double   DEFAULT_DELTAMEAN    = 0.01;
const double   DEFAULT_DELTASTDDEV  = 0.01;
const uint32_t DEFAULT_MHP          = 200;
const double   DEFAULT_MMP          = 0.15;

struct PremoSettings {

    // I/O flags
    bool HasFastqFilename1;
    bool HasFastqFilename2;
    bool HasMosaikPath;
    bool HasOutputFilename;
    bool HasReferenceArchiveFilename;
    bool HasScratchPath;
    bool IsVerbose;
    bool IsVersionRequested;

    // premo flags
    bool HasBatchSize;
    bool HasDeltaMean;
    bool HasDeltaStdDev;

    // mosaik "rule" flags
    bool HasActIntercept;
    bool HasActSlope;
    bool HasBandwidth;
    bool HasMhp;
    bool HasMmp;

    // I/O parameters
    std::string FastqFilename1;
    std::string FastqFilename2;
    std::string MosaikPath;
    std::string OutputFilename;
    std::string ReferenceArchiveFilename;
    std::string ScratchPath;

    // premo parameters
    uint32_t BatchSize;
    double   DeltaMean;
    double   DeltaStdDev;

    // mosaik "rule" parameters
    uint32_t ActIntercept;
    double   ActSlope;
    double   Bandwidth;
    uint32_t Mhp;
    double   Mmp;

    // ctors
    PremoSettings(void)
        : HasFastqFilename1(false)
        , HasFastqFilename2(false)
        , HasMosaikPath(false)
        , HasOutputFilename(false)
        , HasReferenceArchiveFilename(false)
        , HasScratchPath(false)
        , IsVerbose(false)
        , IsVersionRequested(false)
        , HasBatchSize(false)
        , HasDeltaMean(false)
        , HasDeltaStdDev(false)
        , HasActIntercept(false)
        , HasActSlope(false)
        , HasBandwidth(false)
        , HasMhp(false)
        , HasMmp(false)
        , FastqFilename1("")
        , FastqFilename2("")
        , MosaikPath("")
        , OutputFilename("")
        , ReferenceArchiveFilename("")
        , ScratchPath("")
        , BatchSize(DEFAULT_BATCHSIZE)
        , DeltaMean(DEFAULT_DELTAMEAN)
        , DeltaStdDev(DEFAULT_DELTASTDDEV)
        , ActIntercept(DEFAULT_ACTINTERCEPT)
        , ActSlope(DEFAULT_ACTSLOPE)
        , Bandwidth(DEFAULT_BANDWIDTH)
        , Mhp(DEFAULT_MHP)
        , Mmp(DEFAULT_MMP)
    { }

    PremoSettings(const PremoSettings& other)
        : HasFastqFilename1(other.HasFastqFilename1)
        , HasFastqFilename2(other.HasFastqFilename2)
        , HasMosaikPath(other.HasMosaikPath)
        , HasOutputFilename(other.HasOutputFilename)
        , HasReferenceArchiveFilename(other.HasReferenceArchiveFilename)
        , HasScratchPath(other.HasScratchPath)
        , IsVerbose(other.IsVerbose)
        , IsVersionRequested(other.IsVersionRequested)
        , HasBatchSize(other.HasBatchSize)
        , HasDeltaMean(other.HasDeltaMean)
        , HasDeltaStdDev(other.HasDeltaStdDev)
        , HasActIntercept(other.HasActIntercept)
        , HasActSlope(other.HasActSlope)
        , HasBandwidth(other.HasBandwidth)
        , HasMhp(other.HasMhp)
        , HasMmp(other.HasMmp)
        , FastqFilename1(other.FastqFilename1)
        , FastqFilename2(other.FastqFilename2)
        , MosaikPath(other.MosaikPath)
        , OutputFilename(other.OutputFilename)
        , ReferenceArchiveFilename(other.ReferenceArchiveFilename)
        , ScratchPath(other.ScratchPath)
        , BatchSize(other.BatchSize)
        , DeltaMean(other.DeltaMean)
        , DeltaStdDev(other.DeltaStdDev)
        , ActIntercept(other.ActIntercept)
        , ActSlope(other.ActSlope)
        , Bandwidth(other.Bandwidth)
        , Mhp(other.Mhp)
        , Mmp(other.Mmp)
    { }
};

#endif // PREMO_SETTINGS_H
