



#include <chrono>
#include <future>
#include <thread>
#include <fstream>

#include "matchingcommand.h"
#include "graph/graph.h"
#include "GenerateFilteringPlan.h"
#include "FilterVertices.h"
#include "BuildTable.h"
#include "GenerateQueryPlan.h"
#include "EvaluateQuery.h"

#include "utility/analyze_symmetry/analyze_symmetry.h"

#define NANOSECTOSEC(elapsed_time) ((elapsed_time)/(double)1000000000)
#define BYTESTOMB(memory_cost) ((memory_cost)/(double)(1024 * 1024))
//#define printingM
int parseLine(const char *line)
{
    // This assumes the format "VmSize: <value> kB"
    const char *p = line;
    // Move p to the first digit after "VmSize:"
    while (*p < '0' || *p > '9')
        p++;
    // Convert the number part of the line to an integer
    int value = atoi(p);
    return value;
}

int getValue1()
{
    FILE *file = fopen("/proc/self/status", "r");
    if (file == nullptr)
    {
        perror("fopen");
        return -1;
    }

    int result = -1;
    char line[256]; // Increase the buffer size to handle longer lines

    while (fgets(line, sizeof(line), file) != NULL)
    {
        // Print the line for debugging purposes
        // cout << "Read line: " << line;
        // Check if the line starts with "VmSize:"
        if (strncmp(line, "VmSize:", 7) == 0)
        {
            // cout << "Found VmSize line: " << line;
            result = parseLine(line);
            break;
        }
    }

    fclose(file);
    return result;
}
struct matching_algo_outputs
{   enumResult enumOutput;
    ui query_size;
    vector<set<ui>> candidate;
    ui candidate_count_sum_set;
    ui C_E;
    double total_time;
    double preprocessing_time;
    double enumeration_time;
    vector<ui> matching_order;
    string matching_order_string;
    ui *order_pointer = NULL;
    size_t call_count;
};


int main(int argc, char** argv) {
    MatchingCommand command(argc, argv);
    int MemSize=0;
    std::string input_query_graph_file = command.getQueryGraphFilePath();
    std::string input_data_graph_file = command.getDataGraphFilePath();
    std::string input_filter_type = command.getFilterType();
    std::string input_order_type = command.getOrderType();
    std::string input_engine_type = command.getEngineType();
    std::string input_max_embedding_num = command.getMaximumEmbeddingNum();
    std::string input_time_limit = command.getTimeLimit();
    std::string input_order_num = command.getOrderNum();
    std::string input_distribution_file_path = command.getDistributionFilePath();
    std::string input_csr_file_path = command.getCSRFilePath();
    std::string TimeL=command.getTime();
    std::string input_enable_symmetry = command.getEnableSymmetry();
    std::string QS = command.getQuerySize();
    std::string QP = command.getQueryProperty();
    std::string QN = command.getQueryNumber();
    std::string QNL = command.getQueryNumberL();
    int QNLR=stoi(QNL);
    int QNR=stoi(QN);
    //cout<<QNL<<"QNL"<<endl;
    std::string Dataset = command.getDatasetName();
    std::string QuerN = "../../../../../FairSM/dataset/";
    std::string QuerN1 = "../../../../../Pilos-Subgraph_Matching/dataset/";
    string Data_graph=QuerN+Dataset+"/data_graph/"+Dataset+".graph";
    string Data_graphQ=QuerN+Dataset+"/data_graph/"+Dataset+".graph";
    //Data_graph= "../../dataset/"+Dataset+"/data_graph/"+Dataset+".graph";
    input_data_graph_file=Data_graph;
    int FairT=stoi(command.getFairT());
    #ifdef printingM
    std::cout << "Command Line:" << std::endl;
    std::cout << "\tData Graph CSR: " << input_csr_file_path << std::endl;
    std::cout << "\tData Graph: " << input_data_graph_file << std::endl;
    std::cout << "\tQuery Graph: " << input_query_graph_file << std::endl;
    std::cout << "\tFilter Type: " << input_filter_type << std::endl;
    std::cout << "\tOrder Type: " << input_order_type << std::endl;
    std::cout << "\tEngine Type: " << input_engine_type << std::endl;
    std::cout << "\tOutput Limit: " << input_max_embedding_num << std::endl;
    std::cout << "\tTime Limit (seconds): " << input_time_limit << std::endl;
    std::cout << "\tOrder Num: " << input_order_num << std::endl;
    std::cout << "\tDistribution File Path: " << input_distribution_file_path << std::endl;

    std::cout << "\tEnable Symmetry: " << input_enable_symmetry << std::endl;

    std::cout << "--------------------------------------------------------------------" << std::endl;

    
    std::cout << "Load graphs..." << std::endl;
    #endif
    auto start = std::chrono::high_resolution_clock::now();

    
    

    Graph* data_graph = new Graph(true);

    if (input_csr_file_path.empty()) {
        data_graph->loadGraphFromFile(input_data_graph_file);
    }
    else {
        std::string degree_file_path = input_csr_file_path + "_deg.bin";
        std::string edge_file_path = input_csr_file_path + "_adj.bin";
        std::string label_file_path = input_csr_file_path + "_label.bin";
        data_graph->loadGraphFromFileCompressed(degree_file_path, edge_file_path, label_file_path);
    }

    auto end = std::chrono::high_resolution_clock::now();

    double load_graphs_time_in_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
    #ifdef printingM
    std::cout << "-----" << std::endl;
    std::cout << "Query Graph Meta Information" << std::endl;
    //query_graph->printGraphMetaData();
    std::cout << "-----" << std::endl;
    data_graph->printGraphMetaData();

    std::cout << "--------------------------------------------------------------------" << std::endl;
    #endif
    bool enable_symmetry = false;
    if (input_enable_symmetry == "1" || input_enable_symmetry == "true") {
        enable_symmetry = true;
    }
    for (QNR;QNR<=QNLR;QNR++){
    Graph* query_graph = new Graph(true);
    string Query_graph=QuerN1+Dataset+"/query_graph/query_"+QP+"_"+QS+"_"+to_string(QNR)+".graph";
    //Query_graph="../../dataset/"+Dataset+"/query_graph/query_"+QP+"_"+QS+"_"+to_string(QNR)+".graph";
    input_query_graph_file=Query_graph;
    query_graph->loadGraphFromFile(input_query_graph_file);
    query_graph->buildCoreTable();

    std::vector<std::set<std::pair<VertexID, VertexID>>> permutations;
    std::vector<std::pair<VertexID, VertexID>> constraints;
    std::unordered_map<VertexID, std::pair<std::set<VertexID>, std::set<VertexID>>> full_constraints;
    std::unordered_map<VertexID, std::pair<std::set<VertexID>, std::set<VertexID>>> ordered_constraints;
    if (enable_symmetry) {        
        std::unordered_map<VertexID, std::set<VertexID>> cosets = ANALYZE_SYMMETRY::analyze_symmetry(query_graph, permutations);
    
        ANALYZE_SYMMETRY::make_constraints(cosets, constraints);
        

        ANALYZE_SYMMETRY::make_full_constraints(constraints, full_constraints);
        
    }


    
    #ifdef printingM
    std::cout << "Start queries..." << std::endl;
    std::cout << "-----" << std::endl;
    std::cout << "Filter candidates..." << std::endl;
    #endif
    start = std::chrono::high_resolution_clock::now();

    ui** candidates = NULL;
    ui* candidates_count = NULL;
    ui* tso_order = NULL;
    TreeNode* tso_tree = NULL;
    ui* cfl_order = NULL;
    TreeNode* cfl_tree = NULL;
    ui* dpiso_order = NULL;
    TreeNode* dpiso_tree = NULL;
    TreeNode* veq_tree = NULL;
    ui* veq_order = NULL;
    TreeNode* ceci_tree = NULL;
    ui* ceci_order = NULL;
    catalog* storage = NULL;
    enumResult s;
    std::vector<std::unordered_map<VertexID, std::vector<VertexID >>> TE_Candidates;
    std::vector<std::vector<std::unordered_map<VertexID, std::vector<VertexID>>>> NTE_Candidates;
    ui *nodeID = new ui[data_graph->getVerticesCount()];
    std::fill(nodeID, nodeID + data_graph->getVerticesCount(), 0);
    if (input_filter_type == "LDF") {
        FilterVertices::LDFFilter(data_graph, query_graph, candidates, candidates_count);
    } else if (input_filter_type == "NLF") {
        FilterVertices::NLFFilter(data_graph, query_graph, candidates, candidates_count);
    } else if (input_filter_type == "GQL") {
        FilterVertices::GQLFilter(data_graph, query_graph, candidates, candidates_count);
    } else if (input_filter_type == "TSO") {
        FilterVertices::TSOFilter(data_graph, query_graph, candidates, candidates_count, tso_order, tso_tree);
    } else if (input_filter_type == "CFL") {
        FilterVertices::CFLFilter(data_graph, query_graph, candidates, candidates_count, cfl_order, cfl_tree);
    } else if (input_filter_type == "DPiso") {
        FilterVertices::DPisoFilter(data_graph, query_graph, candidates, candidates_count, dpiso_order, dpiso_tree);
    } else if (input_filter_type == "VEQ") {
        FilterVertices::VEQFilter(data_graph, query_graph, candidates, candidates_count, veq_order, veq_tree);
    } else if (input_filter_type == "CECI") {
        FilterVertices::CECIFilter(data_graph, query_graph, candidates, candidates_count, ceci_order, ceci_tree, TE_Candidates, NTE_Candidates);
    } else if (input_filter_type == "RM") {
        FilterVertices::RMFilter(data_graph, query_graph, candidates, candidates_count, storage);
    } else if (input_filter_type == "CaLiG") {
        FilterVertices::CaLiGFilter(data_graph, query_graph, candidates, candidates_count);
    }  else {
        std::cout << "The specified filter type '" << input_filter_type << "' is not supported." << std::endl;
        exit(-1);
    }
        if (getValue1() > MemSize)
        MemSize = getValue1();
    
    ui total_candidates_count = 0;
    for (ui i = 0; i < query_graph->getVerticesCount(); ++i) {
       total_candidates_count += candidates_count[i];
    }

    
    if (input_filter_type != "CECI")
        FilterVertices::sortCandidates(candidates, candidates_count, query_graph->getVerticesCount());

    end = std::chrono::high_resolution_clock::now();
    double filter_vertices_time_in_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

    
#ifdef OPTIMAL_CANDIDATES
    std::vector<ui> optimal_candidates_count;
    double avg_false_positive_ratio = FilterVertices::computeCandidatesFalsePositiveRatio(data_graph, query_graph, candidates,
                                                                                          candidates_count, optimal_candidates_count);
    FilterVertices::printCandidatesInfo(query_graph, candidates_count, optimal_candidates_count);
#endif
#ifdef printingM
    std::cout << "-----" << std::endl;
    std::cout << "Build indices..." << std::endl;

    start = std::chrono::high_resolution_clock::now();

#endif    
    Edges ***edge_matrix = NULL;
    if (input_filter_type != "CECI") {
        edge_matrix = new Edges **[query_graph->getVerticesCount()];
        for (ui i = 0; i < query_graph->getVerticesCount(); ++i) {
            edge_matrix[i] = new Edges *[query_graph->getVerticesCount()];
        }

        BuildTable::buildTables(data_graph, query_graph, candidates, candidates_count, edge_matrix);
    }
    end = std::chrono::high_resolution_clock::now();
    double build_table_time_in_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
    int sum = 0;
    int candidate_count_sum = accumulate(candidates_count, candidates_count + query_graph->getVerticesCount(), sum);
        unordered_set<ui> CandidateSet;
    for (int i = 0; i < query_graph->getVerticesCount(); i++){
        for (int j = 0; j < candidates_count[i]; j++)
        {
            CandidateSet.insert(candidates[i][j]);
        }
    }
    

    size_t memory_cost_in_bytes = 0;
    if (input_filter_type != "CECI") {
        memory_cost_in_bytes = BuildTable::computeMemoryCostInBytes(query_graph, candidates_count, edge_matrix);
        
    }
    else {
        memory_cost_in_bytes = BuildTable::computeMemoryCostInBytes(query_graph, candidates_count, ceci_order, ceci_tree,
                TE_Candidates, NTE_Candidates);
        
    }

#ifdef printingM
    std::cout << "-----" << std::endl;
    std::cout << "Generate a matching order..." << std::endl;
#endif
    start = std::chrono::high_resolution_clock::now();

    ui* matching_order = NULL;
    ui* pivots = NULL;
    ui** weight_array = NULL;

    size_t order_num = 0;
    sscanf(input_order_num.c_str(), "%zu", &order_num);

    std::vector<std::vector<ui>> spectrum;
    if (input_order_type == "QSI") {
        GenerateQueryPlan::generateQSIQueryPlan(data_graph, query_graph, edge_matrix, matching_order, pivots);
    } else if (input_order_type == "GQL") {
        GenerateQueryPlan::generateGQLQueryPlan(data_graph, query_graph, candidates_count, matching_order, pivots);
    }else if (input_order_type == "DSQL") {
        GenerateQueryPlan::generateDSQLQueryPlan(data_graph, query_graph, candidates_count, matching_order, pivots); 
    }else if (input_order_type == "TSO") {
        if (tso_tree == NULL) {
            GenerateFilteringPlan::generateTSOFilterPlan(data_graph, query_graph, tso_tree, tso_order);
        }
        GenerateQueryPlan::generateTSOQueryPlan(query_graph, edge_matrix, matching_order, pivots, tso_tree, tso_order);
    } else if (input_order_type == "CFL") {
        if (cfl_tree == NULL) {
            int level_count;
            ui* level_offset;
            GenerateFilteringPlan::generateCFLFilterPlan(data_graph, query_graph, cfl_tree, cfl_order, level_count, level_offset);
            delete[] level_offset;
        }
        GenerateQueryPlan::generateCFLQueryPlan(data_graph, query_graph, edge_matrix, matching_order, pivots, cfl_tree, cfl_order, candidates_count);
    } else if (input_order_type == "DPiso") {
        if (dpiso_tree == NULL) {
            GenerateFilteringPlan::generateDPisoFilterPlan(data_graph, query_graph, dpiso_tree, dpiso_order);
        }

        GenerateQueryPlan::generateDSPisoQueryPlan(query_graph, edge_matrix, matching_order, pivots, dpiso_tree, dpiso_order,
                                                    candidates_count, weight_array);
    } else if (input_order_type == "CECI") {
        GenerateQueryPlan::generateCECIQueryPlan(query_graph, ceci_tree, ceci_order, matching_order, pivots);
    } else if (input_order_type == "RI") {
        GenerateQueryPlan::generateRIQueryPlan(data_graph, query_graph, matching_order, pivots);
    } else if (input_order_type == "VF2PP") {
        GenerateQueryPlan::generateVF2PPQueryPlan(data_graph, query_graph, matching_order, pivots);
    } else if (input_order_type == "VF3") {
        GenerateQueryPlan::generateVF3QueryPlan(data_graph, query_graph, matching_order, pivots);
    } else if (input_order_type == "RM") {
        GenerateQueryPlan::generateRMQueryPlan(query_graph, matching_order, edge_matrix, pivots);
    }  else {
        std::cout << "The specified order type '" << input_order_type << "' is not supported." << std::endl;
    }

    end = std::chrono::high_resolution_clock::now();
    double generate_query_plan_time_in_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();


    GenerateQueryPlan::checkQueryPlanCorrectness(query_graph, matching_order, pivots);
    #ifdef printingM
    GenerateQueryPlan::printSimplifiedQueryPlan(query_graph, matching_order);
    #endif

        if (getValue1() > MemSize)
        MemSize = getValue1();
if (enable_symmetry) {
    ANALYZE_SYMMETRY::make_ordered_constraints(matching_order, query_graph->getVerticesCount(), full_constraints, ordered_constraints);
    
}
    #ifdef printingM
    std::cout << "-----" << std::endl;
    std::cout << "Enumerate..." << std::endl;
    #endif
    size_t output_limit = 0;
    size_t embedding_count = 0;
    if (input_max_embedding_num == "MAX") {
        output_limit = std::numeric_limits<size_t>::max();
    }
    else {
        sscanf(input_max_embedding_num.c_str(), "%zu", &output_limit);
    }


#if ENABLE_QFLITER == 1
    EvaluateQuery::qfliter_bsr_graph_ = BuildTable::qfliter_bsr_graph_;
#endif

    size_t call_count = 0;
    size_t time_limit = 0;
    size_t valid_vtx_count = 0;
    sscanf(input_time_limit.c_str(), "%zu", &time_limit); 
    int TimeL1= stoi(TimeL);
    start = std::chrono::high_resolution_clock::now();
    //cout<<input_engine_type<<endl;
    if (input_engine_type == "EXPLORE") {       
        embedding_count = EvaluateQuery::exploreGraph(data_graph, query_graph, edge_matrix, candidates,
                                                      candidates_count, matching_order, pivots, output_limit, call_count,
                                                      ordered_constraints);
    } else if (input_engine_type == "LFTJ") {  
        
        s = EvaluateQuery::LFTJ(data_graph, query_graph, nodeID,edge_matrix, candidates, candidates_count,
                                              matching_order, output_limit, call_count, valid_vtx_count,TimeL1,FairT,
                                              ordered_constraints);
                                              }
    else if (input_engine_type == "LFTJK") {  
        
        s = EvaluateQuery::LFTJKOPT(data_graph, query_graph, nodeID,edge_matrix, candidates, candidates_count,
                                              matching_order, output_limit, call_count, valid_vtx_count,TimeL1,FairT,
                                              ordered_constraints);
                                              }                                          
    else if (input_engine_type == "LFTJDLS") {  
        
        s = EvaluateQuery::LFTJDLS(data_graph, query_graph, nodeID,edge_matrix, candidates, candidates_count,
                                              matching_order, output_limit, call_count, valid_vtx_count,TimeL1,FairT,
                                              ordered_constraints);
    }
        else if (input_engine_type == "DLSBS") {  
        

        s = EvaluateQuery::DSQL(data_graph, query_graph, nodeID,edge_matrix, candidates, candidates_count,
                                              matching_order, output_limit, call_count, valid_vtx_count,TimeL1,FairT,
                                              ordered_constraints);
    }
    else if (input_engine_type == "DLSBSK") {  
        

        s = EvaluateQuery::DSQLTOPK(data_graph, query_graph, nodeID,edge_matrix, candidates, candidates_count,
                                              matching_order, output_limit, call_count, valid_vtx_count,TimeL1,FairT,
                                              ordered_constraints);
    }
    
    
        else if (input_engine_type == "LFSK") {  
        
        s = EvaluateQuery::LFTJDLSKOPT(data_graph, query_graph, nodeID,edge_matrix, candidates, candidates_count,
                                              matching_order, output_limit, call_count, valid_vtx_count,TimeL1,FairT,
                                              ordered_constraints);
    }
     else if (input_engine_type == "GQL") {    
        embedding_count = EvaluateQuery::exploreGraphQLStyle(data_graph, query_graph, candidates, candidates_count,
                                                             matching_order, output_limit, call_count,
                                                             ordered_constraints);
    } else if (input_engine_type == "QSI") {    
        embedding_count = EvaluateQuery::exploreQuickSIStyle(data_graph, query_graph, candidates, candidates_count,
                                                             matching_order, pivots, output_limit, call_count,
                                                             ordered_constraints);
    } else if (input_engine_type == "VF3") {    
        embedding_count = EvaluateQuery::exploreVF3Style(data_graph, query_graph, candidates, candidates_count,
                                                             matching_order, pivots, output_limit, call_count,
                                                             ordered_constraints);
    }
    else if (input_engine_type == "VEQ") {  
        
        //if (veq_tree == NULL) {
        //    GenerateFilteringPlan::generateDPisoFilterPlan(data_graph, query_graph, veq_tree, veq_order);
        //}
        
        s = EvaluateQuery::exploreVEQStyleOO(data_graph, query_graph,nodeID, veq_tree,veq_order,
                                                           edge_matrix, candidates, candidates_count,
                                                           output_limit, call_count,TimeL1,full_constraints);
        //s = EvaluateQuery::exploreVEQStyleOV(data_graph, query_graph, veq_tree,
        //                                                   edge_matrix, candidates, candidates_count,
        //                                                   output_limit, call_count,TimeL1,matching_order,full_constraints);
    }
        
    
        else if (input_engine_type == "VEQGM") {  
        //if (veq_tree == NULL) {
        //    GenerateFilteringPlan::generateDPisoFilterPlan(data_graph, query_graph, veq_tree, veq_order);
        //}
        s = EvaluateQuery::exploreVEQStyleOV(data_graph, query_graph, veq_tree,
                                                           edge_matrix, candidates, candidates_count,
                                                           output_limit, call_count,TimeL1,matching_order,ordered_constraints);
        //s = EvaluateQuery::exploreVEQStyle1(data_graph, query_graph, veq_tree,
          //                                                 edge_matrix, candidates, candidates_count,
          //                                                 output_limit, call_count,TimeL1, full_constraints);
    }
        else if (input_engine_type == "DV") {  
        //if (veq_tree == NULL) {
        //    GenerateFilteringPlan::generateDPisoFilterPlan(data_graph, query_graph, veq_tree, veq_order);        
        //}
        s = EvaluateQuery::DIVSM(data_graph, query_graph,nodeID,edge_matrix, candidates, candidates_count,matching_order,
                                                           output_limit, call_count,TimeL1,FairT, ordered_constraints);
        //s = EvaluateQuery::LFTJDIV1(data_graph, query_graph,edge_matrix, candidates, candidates_count,matching_order,
        //                                                   output_limit, call_count,TimeL1,FairT, ordered_constraints);
    }else if (input_engine_type == "DVS") {  
        //if (veq_tree == NULL) {
        //    GenerateFilteringPlan::generateDPisoFilterPlan(data_graph, query_graph, veq_tree, veq_order);        
        //}
        s = EvaluateQuery::DIVSMSQ(data_graph, query_graph,nodeID,edge_matrix, candidates, candidates_count,matching_order,
                                                           output_limit, call_count,TimeL1,FairT, ordered_constraints);
        //s = EvaluateQuery::LFTJDIV1(data_graph, query_graph,edge_matrix, candidates, candidates_count,matching_order,
        //                                                   output_limit, call_count,TimeL1,FairT, ordered_constraints);
    }

            else if (input_engine_type == "DVGM") {  
        //if (veq_tree == NULL) {
        //    GenerateFilteringPlan::generateDPisoFilterPlan(data_graph, query_graph, veq_tree, veq_order);
        
        //}

        s=EvaluateQuery::DIVTOPKOPT(data_graph, query_graph,nodeID,edge_matrix, candidates, candidates_count,matching_order,
                                                           output_limit, call_count,TimeL1,FairT, ordered_constraints);
        //s = EvaluateQuery::LFTJDIV1(data_graph, query_graph,edge_matrix, candidates, candidates_count,matching_order,
        //                                                   output_limit, call_count,TimeL1,FairT, ordered_constraints);
    }
            else if (input_engine_type == "DVGMSQ") {  
        //if (veq_tree == NULL) {
        //    GenerateFilteringPlan::generateDPisoFilterPlan(data_graph, query_graph, veq_tree, veq_order);
        
        //}

        s=EvaluateQuery::DIVTOPKOPTSQ(data_graph, query_graph,nodeID,edge_matrix, candidates, candidates_count,matching_order,
                                                           output_limit, call_count,TimeL1,FairT, ordered_constraints);
        //s = EvaluateQuery::LFTJDIV1(data_graph, query_graph,edge_matrix, candidates, candidates_count,matching_order,
        //                                                   output_limit, call_count,TimeL1,FairT, ordered_constraints);
    }
    else if (input_engine_type == "DPiso") { 
        embedding_count = EvaluateQuery::exploreDPisoStyle(data_graph, query_graph, dpiso_tree,
                                                           edge_matrix, candidates, candidates_count,
                                                           weight_array, dpiso_order, output_limit,
                                                           call_count,
                                                           full_constraints);
       s.embedding_cnt=      embedding_count;                                              
    }
    else if (input_engine_type == "RM") {
        s = EvaluateQuery::exploreRMStyle(query_graph, data_graph, storage, edge_matrix, candidates, candidates_count,
                                                        matching_order, output_limit, call_count,TimeL1);
    }
    else if (input_engine_type == "KSS") {
        s = EvaluateQuery::exploreKSSStyle(query_graph, data_graph, edge_matrix, candidates, candidates_count,
                                                        matching_order, output_limit, call_count,TimeL1);
    }
    else if (input_engine_type == "CECI") {   
        embedding_count = EvaluateQuery::exploreCECIStyle(data_graph, query_graph, ceci_tree, candidates, candidates_count, TE_Candidates,
                NTE_Candidates, ceci_order, output_limit, call_count,
                ordered_constraints);
    }
    else {
        std::cout << "The specified engine type '" << input_engine_type << "' is not supported." << std::endl;
        exit(-1);
    }

    end = std::chrono::high_resolution_clock::now();
    double enumeration_time_in_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

#ifdef DISTRIBUTION
    std::ofstream outfile (input_distribution_file_path , std::ofstream::binary);
    outfile.write((char*)EvaluateQuery::distribution_count_, sizeof(size_t) * data_graph->getVerticesCount());
    delete[] EvaluateQuery::distribution_count_;
#endif
    #ifdef printingM
    std::cout << "--------------------------------------------------------------------" << std::endl;
    std::cout << "Release memories..." << std::endl;
    #endif
    delete[] candidates_count;
    delete[] tso_order;
    delete[] tso_tree;
    delete[] cfl_order;
    delete[] cfl_tree;
    ordered_constraints.clear();
    delete[] dpiso_order;
    delete[] dpiso_tree;
    delete[] ceci_order;
    delete[] ceci_tree;
    delete[] matching_order;
    delete[] pivots;
    delete storage;
    
    for (ui i = 0; i < query_graph->getVerticesCount(); ++i) {
        delete[] candidates[i];
    }
    delete[] candidates;
    
    if (edge_matrix != NULL) {
        for (ui i = 0; i < query_graph->getVerticesCount(); ++i) {
            for (ui j = 0; j < query_graph->getVerticesCount(); ++j) {
                delete edge_matrix[i][j];
            }
            delete[] edge_matrix[i];
        }
        delete[] edge_matrix;
    }
    
    if (weight_array != NULL) {
        for (ui i = 0; i < query_graph->getVerticesCount(); ++i) {
            delete[] weight_array[i];
        }
        delete[] weight_array;
    }
    
    delete query_graph;
    
    //delete data_graph;
    double preprocessing_time_in_ns = filter_vertices_time_in_ns + build_table_time_in_ns + generate_query_plan_time_in_ns;
    double total_time_in_ns = preprocessing_time_in_ns + enumeration_time_in_ns;
    #ifdef printingM
    std::cout << "--------------------------------------------------------------------" << std::endl;
    
    
    printf("Load graphs time (seconds): %.4lf\n", NANOSECTOSEC(load_graphs_time_in_ns));
    printf("Filter vertices time (seconds): %.4lf\n", NANOSECTOSEC(filter_vertices_time_in_ns));
    printf("Build table time (seconds): %.4lf\n", NANOSECTOSEC(build_table_time_in_ns));
    printf("Generate query plan time (seconds): %.4lf\n", NANOSECTOSEC(generate_query_plan_time_in_ns));
    printf("Enumerate time (seconds): %.4lf\n", NANOSECTOSEC(enumeration_time_in_ns));
    printf("Preprocessing time (seconds): %.4lf\n", NANOSECTOSEC(preprocessing_time_in_ns));
    printf("Total time (seconds): %.4lf\n", NANOSECTOSEC(total_time_in_ns));
    printf("Memory cost (MB): %.4lf\n", BYTESTOMB(memory_cost_in_bytes));
    printf("#Embeddings: %zu\n", embedding_count);
    printf("Call Count: %zu\n", call_count);
    printf("Per Call Count Time (nanoseconds): %.4lf\n", enumeration_time_in_ns / (call_count == 0 ? 1 : call_count));
    std::cout << "End." << std::endl;

    #endif
    std::fstream output;
    output.open("/root/subgraph/test/output/survey_0408.csv", std::ios::out | std::ios::app);
    
    output << input_query_graph_file << ",";
    output << input_data_graph_file << ",";
    output << input_filter_type << ",";
    output << input_order_type << ",";
    output << input_engine_type << ",";
    output << NANOSECTOSEC(load_graphs_time_in_ns) << ",";
    output << NANOSECTOSEC(filter_vertices_time_in_ns) << ",";
    output << NANOSECTOSEC(build_table_time_in_ns) << ",";
    output << NANOSECTOSEC(generate_query_plan_time_in_ns) << ",";
    output << NANOSECTOSEC(enumeration_time_in_ns) << ",";
    output << NANOSECTOSEC(preprocessing_time_in_ns) << ",";
    output << NANOSECTOSEC(total_time_in_ns) << ",";
    output << embedding_count << ",";
    output << call_count << std::endl;
    double total_time=NANOSECTOSEC(total_time_in_ns);
    double preprocessing_time=NANOSECTOSEC(filter_vertices_time_in_ns)+NANOSECTOSEC(build_table_time_in_ns);
    double enumeration_time=NANOSECTOSEC(enumeration_time_in_ns);
    string StoreFile=command.getStoreFile();
    output.close();
    int row = 0;
    std::ostringstream oss;
    oss << QNR << " " << QS << " " << call_count << " " << s.embedding_cnt<<" "<<s.Can_embed<<" "<<CandidateSet.size();
    oss << " " << total_time << " " << candidate_count_sum; //<<LDF.first.enumOutput.embedding_cnt;
    oss << " " << preprocessing_time;
    oss << " " << enumeration_time<<" ";
    oss<<s.topk[0]<<" "<<s.topk[1]<<" "<<s.topk[2]<<" "<<s.topk[3]<<" "<<s.topk[4]<<" "<<s.topk[5]<<" "<<s.topk[6]<<" "<<s.topk[7];
    std::string var = oss.str();
    cout << var << endl;
    string file_path = "";
    file_path = "performance_experiment/" + StoreFile + "_" + input_engine_type + "_" + input_max_embedding_num + "_" + Dataset + "_" + QP +"_"+ QS +"_"+TimeL+ ".csv";
    std::ofstream myfile;
    myfile.open(file_path, std::ios_base::app);
    myfile << var << "\n";
    myfile.close();
    
}
return 0;
}

