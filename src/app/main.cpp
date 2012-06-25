// ***************************************************************************
// main.cpp (c) 2012 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 24 June 2012 (DB)
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
         << "(c) 2012 Derek Barnett" << endl
         << "Boston College, Biology Dept." << endl
         << "------------------------------" << endl
         << endl;
}

int main(int argc, char* argv[]) {

    // -------------------------------------------------------
    // command line parameters & help info
    // -------------------------------------------------------

    // set program details
    const string name("premo");
    const string description("\"pre-Mosaik\" application that generates MosaikAligner parameters "
                             "for paired-end sequencing data. Premo uses a bootstrapping heuristic "
                             "to estimate the overall read length & fragment length, running "
                             "Mosaik on samples from the input until it sees convergence on both of "
                             "these values. The resulting parameters, reported in JSON format, "
                             "should allow Mosaik to perform well on the full dataset.");
    const string usage("-annpe <file> "
                       "-annse <file> "
                       "-fq1 <file> "
                       "-fq2 <file> "
                       "-ref <file> "
                       "-jmp <file prefix> "
                       "-mosaik <dir> "
                       "-out <file> "
                       "-ref <file> "
                       "-st <technology> "
                       "[-tmp <dir>] "
                       "[options]" );
    Options::SetProgramInfo(name, description, usage);

    // hook up command-line options to our settings structure
    PremoSettings settings;

    OptionGroup* IO_Opts = Options::CreateOptionGroup("Input & Output");

    const string annpe("neural network filename (paired-end)");
    const string annse("neural network filename (single-end)");
    const string fq1("input FASTQ file (mate 1)");
    const string fq2("input FASTQ file (mate 2)");
    const string jump("stub for jump database files");
    const string keep("keep generated files (auto-deleted by default)");
    const string mosaik("/path/to/Mosaik/bin");
    const string out("output file (JSON). Contains generated Mosaik parameters & raw batch results");
    const string ref("MosaikBuild-generated reference archive");
    const string tmp("scratch directory for any generated files");
    const string verbose("verbose output (to stderr)");
    const string version("show version information");

    const string FN("filename");
    const string DIR("directory");

    Options::AddValueOption("-annpe",  FN,  annpe,  "", settings.HasAnnPeFilename,     settings.AnnPeFilename,     IO_Opts);
    Options::AddValueOption("-annse",  FN,  annse,  "", settings.HasAnnSeFilename,     settings.AnnSeFilename,     IO_Opts);
    Options::AddValueOption("-fq1",    FN,  fq1,    "", settings.HasFastqFilename1,    settings.FastqFilename1,    IO_Opts);
    Options::AddValueOption("-fq2",    FN,  fq2,    "", settings.HasFastqFilename2,    settings.FastqFilename2,    IO_Opts);
    Options::AddValueOption("-jmp",    FN,  jump,   "", settings.HasJumpDbStub,        settings.JumpDbStub,        IO_Opts);
    Options::AddValueOption("-mosaik", DIR, mosaik, "", settings.HasMosaikPath,        settings.MosaikPath,        IO_Opts);
    Options::AddValueOption("-out",    FN,  out,    "", settings.HasOutputFilename,    settings.OutputFilename,    IO_Opts);
    Options::AddValueOption("-ref",    FN,  ref,    "", settings.HasReferenceFilename, settings.ReferenceFilename, IO_Opts);
    Options::AddValueOption("-tmp",    DIR, tmp,    "", settings.HasScratchPath,       settings.ScratchPath,       IO_Opts, Defaults::ScratchPath);
    Options::AddOption("-keep",    keep,    settings.IsKeepGeneratedFiles, IO_Opts);
    Options::AddOption("-v",       verbose, settings.IsVerbose,            IO_Opts);
    Options::AddOption("-version", version, settings.IsVersionRequested,   IO_Opts);

    OptionGroup* PremoOpts = Options::CreateOptionGroup("Premo Bootstrapping Options");

    const string dfl("delta fragment length (fraction). Premo can stop when overall median fragment length changes by less than this amount after a new batch result");
    const string drl("delta read length (fraction). Premo can stop when overall median read length changes by less than this amount after a new batch result");
    const string n("# of pairs to align per batch");

    Options::AddValueOption("-delta-fl", "double", dfl, "", settings.HasDeltaFragmentLength, settings.DeltaFragmentLength, PremoOpts, Defaults::DeltaFragmentLength);
    Options::AddValueOption("-delta-rl", "double", drl, "", settings.HasDeltaReadLength,     settings.DeltaReadLength,     PremoOpts, Defaults::DeltaReadLength);
    Options::AddValueOption("-n",        "int",    n,   "", settings.HasBatchSize,           settings.BatchSize,           PremoOpts, Defaults::BatchSize);

    OptionGroup* MosaikOpts = Options::CreateOptionGroup("Mosaik Parameter-Generation Options");

    const string act("alignment candidate threshold. Generated MosaikAligner -act parameter will be ((ActSlope * ReadLength) + ActIntercept)");
    const string bwm("banded Smith-Waterman multiplier. Generated MosaikAligner -bw parameter will be (BwMultiplier * Mmp * ReadLength)");
    const string mhp("maximum hash positions. Used in premo batch runs, and included in generated parameter set");
    const string mmp("mismatch percent. Used in premo batch runs, and included in generated parameter set");
    const string st("sequencing technology: '454', 'helicos', 'illumina', 'illumina_long', 'sanger' or 'solid'. Required for premo batch runs, and included in generated parameter set");

    Options::AddValueOption("-act-intercept", "int",    act, "", settings.HasActIntercept, settings.ActIntercept, MosaikOpts, Defaults::ActIntercept);
    Options::AddValueOption("-act-slope",     "double", act, "", settings.HasActSlope,     settings.ActSlope,     MosaikOpts, Defaults::ActSlope);
    Options::AddValueOption("-bwm",           "int",    bwm, "", settings.HasBwMultiplier, settings.BwMultiplier, MosaikOpts, Defaults::BwMultiplier);
    Options::AddValueOption("-mhp",           "int",    mhp, "", settings.HasMhp,          settings.Mhp,          MosaikOpts, Defaults::Mhp);
    Options::AddValueOption("-mmp",           "double", mmp, "", settings.HasMmp,          settings.Mmp,          MosaikOpts, Defaults::Mmp);
    Options::AddValueOption("-st",            "string", st,  "", settings.HasSeqTech,      settings.SeqTech,      MosaikOpts);

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

    // initialize our Premo runner with cmdline settings
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
