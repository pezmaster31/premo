// ***************************************************************************
// main.cpp (c) 2012 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 7 June 2012 (DB)
// ---------------------------------------------------------------------------
// Main entry point for the Premo app.
// ***************************************************************************

#include "options.h"
#include "premo.h"
#include "premo_settings.h"
#include "premo_version.h"
#include <iostream>
#include <string>
using namespace std;

static
void printVersion(void) {

    cerr << endl
         << "------------------------------" << endl
         << "premo v" << PREMO_VERSION_MAJOR << "." << PREMO_VERSION_MINOR << "." << PREMO_VERSION_BUILD << endl
         << "(c) 2012 Derek Barnett, Gabor Marth" << endl
         << "Boston College, Biology Dept." << endl
         << "------------------------------" << endl
         << endl;
}

int main(int argc, char* argv[]) {

    PremoSettings settings;

    // -------------------------------------------------------
    // command line parameters & help info
    // -------------------------------------------------------

    // set program details
    const string programName("premo");
    const string description("generates Mosaik parameters from sequence data using a bootstrapping heuristic");
    const string usage("-fq1 <file> -fq2 <file> -ref <file> -mosaik <path> -tmp <path> [...options...] > mosaik_params.json");
    Options::SetProgramInfo(programName, description, usage);

    // set up options
    const string FN("filename");
    const string DIR("directory");

    OptionGroup* IO_Opts = Options::CreateOptionGroup("Input & Output");

    const string fq1("input FASTQ file (mate 1)");
    const string fq2("input FASTQ file (mate 2)");
    const string mosaik("path/to/Mosaik/bin");
    const string out("destination for JSON output (writes to stdout if no option provided)");
    const string ref("MosaikBuild-generated reference archive");
    const string tmp("scratch directory for any generated files");

    Options::AddValueOption("-fq1",    FN,  fq1,    "", settings.HasFastqFilename1,           settings.FastqFilename1,           IO_Opts);
    Options::AddValueOption("-fq2",    FN,  fq2,    "", settings.HasFastqFilename2,           settings.FastqFilename2,           IO_Opts);
    Options::AddValueOption("-mosaik", DIR, mosaik, "", settings.HasMosaikPath,               settings.MosaikPath,               IO_Opts);
    Options::AddValueOption("-out",    FN,  out,    "", settings.HasOutputFilename,           settings.OutputFilename,           IO_Opts);
    Options::AddValueOption("-ref",    FN,  ref,    "", settings.HasReferenceArchiveFilename, settings.ReferenceArchiveFilename, IO_Opts);
    Options::AddValueOption("-tmp",    DIR, tmp,    "", settings.HasScratchPath,              settings.ScratchPath,              IO_Opts);

    Options::AddOption("-v", "verbose output (to stderr)", settings.IsVerbose, IO_Opts);
    Options::AddOption("-version", "show premo version information", settings.IsVersionRequested, IO_Opts);

    OptionGroup* PremoOpts = Options::CreateOptionGroup("Premo Bootstrapping Options");

    const string dm("delta mean - stop iterating when the mean fragment length differs by less than this value");
    const string ds("delta standard deviation - stop iterating when std. dev. of fragment length differs by less than this value");
    const string n("# sequences to use per Mosaik bootstrap iteration");

    Options::AddValueOption("-delta-mean", "double", dm, "", settings.HasDeltaMean,   settings.DeltaMean,   PremoOpts, DEFAULT_DELTAMEAN);
    Options::AddValueOption("-delta-sd",   "double", ds, "", settings.HasDeltaStdDev, settings.DeltaStdDev, PremoOpts, DEFAULT_DELTASTDDEV);
    Options::AddValueOption("-n",          "int",    n,  "", settings.HasBatchSize,   settings.BatchSize,   PremoOpts, DEFAULT_BATCHSIZE);

    OptionGroup* MosaikOpts = Options::CreateOptionGroup("Mosaik Parameter-Generation Options");

    const string actI("used to calculate alignment candidate threshold. Mosaik -act paramater will be ((ActSlope * ReadLength) + ActIntercept)");
    const string actS("used to calculate alignment candidate threshold. Mosaik -act paramater will be ((ActSlope * ReadLength) + ActIntercept)");
    const string bw("used to calculate Smith-Waterman bandwidth. Mosaik -bw parameter will be (BW-multiplier * MMP * ReadLength)");
    const string mhp("maximum hash positions. Mosaik -mhp parameter will be this value");
    const string mmp("mismatch percent. Mosaik -mmp parameter will be this value");

    Options::AddValueOption("-act-intercept", "int",    actI, "", settings.HasActIntercept, settings.ActIntercept, MosaikOpts, DEFAULT_ACTINTERCEPT);
    Options::AddValueOption("-act-slope",     "double", actS, "", settings.HasActSlope,     settings.ActSlope,     MosaikOpts, DEFAULT_ACTSLOPE);
    Options::AddValueOption("-bw",            "int",    bw,   "", settings.HasBandwidth,    settings.Bandwidth,    MosaikOpts, DEFAULT_BANDWIDTH);
    Options::AddValueOption("-mhp",           "int",    mhp,  "", settings.HasMhp,          settings.Mhp,          MosaikOpts, DEFAULT_MHP);
    Options::AddValueOption("-mmp",           "double", mmp,  "", settings.HasMmp,          settings.Mmp,          MosaikOpts, DEFAULT_MMP);

    // -------------------------------------------------------
    // parse command line
    // -------------------------------------------------------

    // options class will show help, if requested
    Options::Parse(argc, argv);

    // show version info, if requested
    if ( settings.IsVersionRequested ) {
        printVersion();
        return 0;
    }

    // -------------------------------------------------------
    // run Premo using settings
    // -------------------------------------------------------

    // initialize our PremoApp with cmdline settings
    Premo p(settings);

    // run PremoApp... if failed:
    if ( !p.run() ) {

        // print error & return failure
        cerr << "premo ERROR: " << p.errorString() << endl;
        return 1;
    }

    // otherwise return success
    return 0;
}
