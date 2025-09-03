



#include "matchingcommand.h"

MatchingCommand::MatchingCommand(const int argc, char **argv) : CommandParser(argc, argv) {
    
    options_key[OptionKeyword::Algorithm] = "-a";
    options_key[OptionKeyword::IndexType] = "-i";
    options_key[OptionKeyword::QueryGraphFile] = "-q";
    options_key[OptionKeyword::DataGraphFile] = "-d";
    options_key[OptionKeyword::ThreadCount] = "-n";
    options_key[OptionKeyword::DepthThreshold] = "-d0";
    options_key[OptionKeyword::WidthThreshold] = "-w0";
    options_key[OptionKeyword::Filter] = "-filter";
    options_key[OptionKeyword::Order] = "-order";
    options_key[OptionKeyword::Engine] = "-engine";
    options_key[OptionKeyword::MaxOutputEmbeddingNum] = "-num";
    options_key[OptionKeyword::SpectrumAnalysisTimeLimit] = "-time_limit";
    options_key[OptionKeyword::SpectrumAnalysisOrderNum] = "-order_num";
    options_key[OptionKeyword::DistributionFilePath] = "-dis_file";
    options_key[OptionKeyword::CSRFilePath] = "-csr";
    options_key[OptionKeyword::EnableSymmetry] = "-symmetry";
    options_key[OptionKeyword::Dataset] = "-dataset";
    options_key[OptionKeyword::QuerySize] = "-qsize";
    options_key[OptionKeyword::QueryNumber] = "-qnumber";
    options_key[OptionKeyword::QueryNumberL] = "-qnumberL";
    options_key[OptionKeyword::QueryProperty] = "-qprop";
    options_key[OptionKeyword::timeLimit] = "-time";
    options_key[OptionKeyword::outputF] = "-SF";
    options_key[OptionKeyword::FairT] = "-FairT";
    processOptions();
};

void MatchingCommand::processOptions() {
    options_value[OptionKeyword::outputF] = getCommandOption(options_key[OptionKeyword::outputF]);

    options_value[OptionKeyword::QueryGraphFile] = getCommandOption(options_key[OptionKeyword::QueryGraphFile]);;

    
    options_value[OptionKeyword::DataGraphFile] = getCommandOption(options_key[OptionKeyword::DataGraphFile]);

    
    options_value[OptionKeyword::Algorithm] = getCommandOption(options_key[OptionKeyword::Algorithm]);

    
    options_value[OptionKeyword::ThreadCount] = getCommandOption(options_key[OptionKeyword::ThreadCount]);

    
    options_value[OptionKeyword::DepthThreshold] = getCommandOption(options_key[OptionKeyword::DepthThreshold]);

    
    options_value[OptionKeyword::WidthThreshold] = getCommandOption(options_key[OptionKeyword::WidthThreshold]);

    
    options_value[OptionKeyword::IndexType] = getCommandOption(options_key[OptionKeyword::IndexType]);

    
    options_value[OptionKeyword::Filter] = getCommandOption(options_key[OptionKeyword::Filter]);

    
    options_value[OptionKeyword::Order] = getCommandOption(options_key[OptionKeyword::Order]);

    
    options_value[OptionKeyword::Engine] = getCommandOption(options_key[OptionKeyword::Engine]);

    
    options_value[OptionKeyword::MaxOutputEmbeddingNum] = getCommandOption(options_key[OptionKeyword::MaxOutputEmbeddingNum]);

    
    options_value[OptionKeyword::SpectrumAnalysisTimeLimit] = getCommandOption(options_key[OptionKeyword::SpectrumAnalysisTimeLimit]);

    
    options_value[OptionKeyword::SpectrumAnalysisOrderNum] = getCommandOption(options_key[OptionKeyword::SpectrumAnalysisOrderNum]);

    
    options_value[OptionKeyword::DistributionFilePath] = getCommandOption(options_key[OptionKeyword::DistributionFilePath]);

    
    options_value[OptionKeyword::CSRFilePath] = getCommandOption(options_key[OptionKeyword::CSRFilePath]);

    
    options_value[OptionKeyword::EnableSymmetry] = getCommandOption(options_key[OptionKeyword::EnableSymmetry]);

    options_value[OptionKeyword::timeLimit] = getCommandOption(options_key[OptionKeyword::timeLimit]);
    
    options_value[OptionKeyword::QueryNumber] = getCommandOption(options_key[OptionKeyword::QueryNumber]);

    options_value[OptionKeyword::QuerySize] = getCommandOption(options_key[OptionKeyword::QuerySize]);

    options_value[OptionKeyword::QueryProperty] = getCommandOption(options_key[OptionKeyword::QueryProperty]);
    options_value[OptionKeyword::FairT] = getCommandOption(options_key[OptionKeyword::FairT]);
    options_value[OptionKeyword::Dataset] = getCommandOption(options_key[OptionKeyword::Dataset]);
    options_value[OptionKeyword::QueryNumberL] = getCommandOption(options_key[OptionKeyword::QueryNumberL]);
}