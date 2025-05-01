#ifndef SUBGRAPHMATCHING_MATCHINGCOMMAND_H
#define SUBGRAPHMATCHING_MATCHINGCOMMAND_H

#include "utility/commandparser.h"
#include <map>
#include <iostream>
enum OptionKeyword {
    Algorithm = 0,
    QueryGraphFile = 1,
    DataGraphFile = 2,
    ThreadCount = 3,
    DepthThreshold = 4,
    WidthThreshold = 5,
    IndexType = 6,
    Filter = 7,
    Order = 8,
    ExportPlanPath = 9,
    MaxOutputEmbeddingNum = 10,
    SpectrumAnalysisTimeLimit = 11,
    SpectrumAnalysisOrderNum = 12,
    DistributionFilePath = 13,
    CSRFilePath = 14,
    ImportPlanPath = 15,
    InputOrder = 16,
    EnablePreprocessor = 17,
    OutputPath = 18
};

class MatchingCommand : public CommandParser{
private:
    std::map<OptionKeyword, std::string> options_key;
    std::map<OptionKeyword, std::string> options_value;

private:
    void processOptions();

public:
    MatchingCommand(int argc, char **argv);

    std::string getDataGraphFilePath() {
        return options_value[OptionKeyword::DataGraphFile];
    }

    std::string getQueryGraphFilePath() {
        return options_value[OptionKeyword::QueryGraphFile];
    }

    std::string getAlgorithm() {
        return options_value[OptionKeyword::Algorithm];
    }

    std::string getIndexType() {
        return options_value[OptionKeyword::IndexType] == "" ? "VertexCentric" : options_value[OptionKeyword::IndexType];
    }
    std::string getThreadCount() {
        return options_value[OptionKeyword::ThreadCount] == "" ? "1" : options_value[OptionKeyword::ThreadCount];
    }

    std::string getDepthThreshold() {
        return options_value[OptionKeyword::DepthThreshold] == "" ? "0" : options_value[OptionKeyword::DepthThreshold];
    }

    std::string getWidthThreshold() {
        return options_value[OptionKeyword::WidthThreshold] == "" ? "1" : options_value[OptionKeyword::WidthThreshold];
    }

    std::string getFilterType() {
        return options_value[OptionKeyword::Filter] == "" ? "CFL" : options_value[OptionKeyword::Filter];
    }

    std::string getOrderType() {
        return options_value[OptionKeyword::Order] == "" ? "nd" : options_value[OptionKeyword::Order];
    }

    std::string getExportPlanPath() {
        return options_value[OptionKeyword::ExportPlanPath] == "" ? "" : options_value[OptionKeyword::ExportPlanPath];
    }

    std::string getMaximumEmbeddingNum() {
        return options_value[OptionKeyword::MaxOutputEmbeddingNum] == "" ? "MAX" : options_value[OptionKeyword::MaxOutputEmbeddingNum];
    }

    std::string getTimeLimit() {
        return options_value[OptionKeyword::SpectrumAnalysisTimeLimit] == "" ? "60" : options_value[OptionKeyword::SpectrumAnalysisTimeLimit];
    }

    std::string getOrderNum() {
        return options_value[OptionKeyword::SpectrumAnalysisOrderNum] == "" ? "100" : options_value[OptionKeyword::SpectrumAnalysisOrderNum];
    }

    std::string getDistributionFilePath() {
        return options_value[OptionKeyword::DistributionFilePath] == "" ? "temp.distribution" : options_value[OptionKeyword::DistributionFilePath];
    }

    std::string getCSRFilePath() {
        return options_value[OptionKeyword::CSRFilePath] == "" ? "" : options_value[OptionKeyword::CSRFilePath];
    }

    std::string getImportPlanPath() {
        return options_value[OptionKeyword::ImportPlanPath] == "" ? "" : options_value[OptionKeyword::ImportPlanPath];
    }

    std::string getInputOrder() {
        return options_value[OptionKeyword::InputOrder] == "" ? "" : options_value[OptionKeyword::InputOrder];
    }

    std::string getPreprocessor() {
        return options_value[OptionKeyword::EnablePreprocessor] == "" ? "true" : options_value[OptionKeyword::EnablePreprocessor];
    }

    std::string getOutputPath() {
        return options_value[OptionKeyword::OutputPath] == "" ? "embeddings.bin" : options_value[OptionKeyword::OutputPath];
    }
};


#endif
