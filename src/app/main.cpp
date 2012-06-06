// ***************************************************************************
// main.cpp (c) 2012 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 6 June 2012 (DB)
// ---------------------------------------------------------------------------
// Main entry point for the Premo app.
// ***************************************************************************

#include "premo_version.h"
#include "jsoncpp/json_value.h"
#include "bamtools/api/BamAlignment.h"
#include <iostream>
#include <string>
using namespace std;

static void printHelp(void) {



}

static void printVersion(void) {

    cout << endl 
         << "---------------------------" << endl
         << "Premo v" << PREMO_VERSION_MAJOR << "." << PREMO_VERSION_MINOR << "." << PREMO_VERSION_BUILD << endl
         << "(c) 2012 Derek Barnett, Gabor Marth" << endl
         << "Boston College, Biology Dept." << endl
         << "---------------------------" << endl
         << endl;
}

int main(int argc, char* argv[]) {

    Json::Value v;
    BamTools::BamAlignment a;

    cout << "Default RefID is -1? " << a.RefID << endl;

    printHelp();
    printVersion();

    return 0;
}
