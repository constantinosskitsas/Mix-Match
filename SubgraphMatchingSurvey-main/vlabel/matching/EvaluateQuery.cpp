
// #define ENABLE_FAILING_SET
#include "EvaluateQuery.h"
#include "utility/computesetintersection.h"
#include "utility/execution_tree/execution_tree_generator.h"
#include <vector>
#include <cstring>
#include <algorithm>
#include "utility/pretty_print.h"
#include <random>
#if ENABLE_QFLITER == 1
BSRGraph ***EvaluateQuery::qfliter_bsr_graph_;
int *EvaluateQuery::temp_bsr_base1_ = nullptr;
int *EvaluateQuery::temp_bsr_state1_ = nullptr;
int *EvaluateQuery::temp_bsr_base2_ = nullptr;
int *EvaluateQuery::temp_bsr_state2_ = nullptr;
#endif

#ifdef SPECTRUM
bool EvaluateQuery::exit_;
#endif

#ifdef DISTRIBUTION
size_t *EvaluateQuery::distribution_count_;
#endif

bool EvaluateQuery::nextCombination(std::vector<int>& indices, ui N, ui K) {
    int i = K - 1;
    while (i >= 0 && indices[i] == N - K + i) {
        --i;
    }

    // If no such index, we are done
    if (i < 0) return false;

    // Increment the current index
    ++indices[i];

    // Reset subsequent indices
    for (int j = i + 1; j < K; ++j) {
        indices[j] = indices[j - 1] + 1;
    }

    return true;
}
int parseLin2(const char *line)
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
int getValue3()
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
            result = parseLin2(line);
            break;
        }
    }

    fclose(file);
    return result;
}

bool isSorted(const vector<unsigned int> &vec)
{
    for (size_t i = 1; i < vec.size(); ++i)
    {
        if (vec[i] < vec[i - 1])
        {
            return false; // Found an element that is smaller than the previous one
        }
    }
    return true; // All elements are in order
}

void EvaluateQuery::calculateCell(const Graph *query_graph, Edges ***edge_matrix,
                                  ui *candidates_count, size_t **&candidatesHC, unordered_map<size_t, std::vector<VertexID>> *&idToValues, int count2, int i, int j)
{

    ui u_nbrs_count;
    bool add = true;
    const VertexID *u_nbrs = query_graph->getVertexNeighbors(i, u_nbrs_count);
    for (int k = 0; k < candidates_count[i]; k++)
    {
        add = true;
        if (k == j)
        {
            candidatesHC[i][k] = count2;
            idToValues[i][count2].push_back(k);
            continue;
        }
        if (candidatesHC[i][k] != 0)
            continue;
        for (int d = 0; d < u_nbrs_count; d++)
        {
            ui CVQ = u_nbrs[d];
            ui SP = edge_matrix[i][CVQ]->offset_[j];
            ui EP = edge_matrix[i][CVQ]->offset_[j + 1];
            ui SPC = edge_matrix[i][CVQ]->offset_[k];
            ui EPC = edge_matrix[i][CVQ]->offset_[k + 1];
            int Times = EP - SP;
            if ((EP - SP) == (EPC - SPC))
            {
                for (int ee = 0; ee <= Times; ee++)
                {
                    if (edge_matrix[i][CVQ]->edge_[SP + ee] != edge_matrix[i][CVQ]->edge_[SPC + ee])
                    {
                        add = false;
                        break;
                    }
                }
            }
            else
            {
                add = false;
            }
            if (add == false)
                break;
        }
        if (add == true)
        {
            candidatesHC[i][k] = count2;
            idToValues[i][count2].push_back(k);
        }
    }
}

void EvaluateQuery::calculateCellAC(const Graph *query_graph, Edges ***edge_matrix,
                                    ui *candidates_count, size_t **&candidatesHC, unordered_map<size_t, std::vector<VertexID>> *&idToValues, int count2, int i, int j)
{

    ui unbrs_count;
    bool equ = true;
    const ui *unbrs = query_graph->getVertexNeighbors(i, unbrs_count);

    candidatesHC[i][j] = count2;
    idToValues[i][count2].push_back(j);
    for (int k = 0; k < candidates_count[i]; k++)
    {
        bool equ = true;
        if (candidatesHC[i][k] != 0)
            continue;
        for (int u1 = 0; u1 < unbrs_count; u1++)
        {

            ui unbr = unbrs[u1];
            const Edges *edges = edge_matrix[i][unbr];
            if (edges->offset_[j + 1] - edges->offset_[j] != edges->offset_[k + 1] - edges->offset_[k])
            {
                equ = false;
                break;
            }

            for (ui u2 = 0; u2 < edges->offset_[j + 1] - edges->offset_[j]; u2++)
            {
                if (edges->edge_[u2 + edges->offset_[j]] != edges->edge_[u2 + edges->offset_[k]])
                {
                    equ = false;
                    break;
                }
            }
            if (equ == false)
                break;
        }

        if (equ == true)
        {
            candidatesHC[i][k] = count2;
            idToValues[i][count2].push_back(k);
        }
    }
}
void EvaluateQuery::calculateCellFN(const Graph *query_graph, Edges ***edge_matrix,
                                    ui *candidates_count, size_t **&candidatesHC, unordered_map<size_t, std::vector<VertexID>> *&idToValues, int count2, int i, int j, int ROQ[])
{

    ui u_nbrs_count;
    bool add = true;
    const VertexID *u_nbrs = query_graph->getVertexNeighbors(i, u_nbrs_count);
    ui PO = ROQ[i];
    for (int k = j; k < candidates_count[i]; k++)
    {
        add = true;

        if (k == j)
        {
            candidatesHC[i][k] = count2;
            idToValues[i][count2].push_back(k);
            continue;
        }
        if (candidatesHC[i][k] != 0)
            continue;
        for (int d = 0; d < u_nbrs_count; d++)
        {
            ui CVQ = u_nbrs[d];
            if (ROQ[CVQ] < PO)
                continue;
            ui SP = edge_matrix[i][CVQ]->offset_[j];
            ui EP = edge_matrix[i][CVQ]->offset_[j + 1];
            ui SPC = edge_matrix[i][CVQ]->offset_[k];
            ui EPC = edge_matrix[i][CVQ]->offset_[k + 1];
            int Times = EP - SP;
            if ((EP - SP) == (EPC - SPC))
            {
                for (int ee = 0; ee < Times; ee++)
                {
                    if (edge_matrix[i][CVQ]->edge_[SP + ee] != edge_matrix[i][CVQ]->edge_[SPC + ee])
                    {
                        add = false;
                        break;
                    }
                }
            }
            else
            {
                add = false;
            }
            if (add == false)
                break;
        }
        if (add == true)
        {
            candidatesHC[i][k] = count2;
            idToValues[i][count2].push_back(k);
        }
    }
}

void EvaluateQuery::rankSimple(ui u,ui *&nodeId, ui **candidates, ui **valid_candidate_idx, ui *idx, ui valid_idx, ui *idx_count)
{
    VertexID vN = candidates[u][valid_idx];

    if (nodeId[vN]==1)
    {

        ui ccIDX = idx[0] + 1;
        ui valid_idx1 = valid_candidate_idx[0][ccIDX];
        VertexID vC = candidates[u][valid_idx1];
        while (ccIDX < idx_count[0])
        {
            if (valid_idx1 != 10000000)
            {
                vC = candidates[u][valid_idx1];
                if (nodeId[vC] == 0)
                {
                    ui tempCC = valid_candidate_idx[0][ccIDX];
                    valid_candidate_idx[0][ccIDX] = valid_candidate_idx[0][idx[0]];
                    valid_candidate_idx[0][idx[0]] = tempCC;
                    break;
                }
            }
            ccIDX++;
            valid_idx1 = valid_candidate_idx[0][ccIDX];
        }
    }
}

void EvaluateQuery::rankLess(const Graph *query_graph, Edges ***edge_matrix, ui u, ui *&nodeId, ui **candidates, ui **valid_candidate_idx, ui *idx, ui valid_idx, ui *idx_count)
{
    vector<int> rankV;
    ui ccIDX = idx[0];

    int minVal = 100000;
    int minPos = ccIDX;
    int countmin = 0;
    while (ccIDX < idx_count[0])
    {
        ui valid_idx1 = valid_candidate_idx[0][ccIDX];
        if (valid_idx1 != 10000000)
        {
            VertexID vC = candidates[u][valid_idx1];
            if (nodeId[vC] == 0)
            {
                // add all elements to a vector. (not found yet)
                rankV.push_back(ccIDX);
            }
        }
        ccIDX++;
    }
    ui u_nbrs_count;
    for (int i = 0; i < rankV.size(); i++)
    {
        ccIDX = rankV[i];
        ui valid_idx1 = valid_candidate_idx[0][ccIDX];
        const VertexID *u_nbrs = query_graph->getVertexNeighbors(u, u_nbrs_count);
        countmin = 0;
        for (int d = 0; d < u_nbrs_count; d++)
        {
            ui CVQ = u_nbrs[d];
            ui SP = edge_matrix[u][CVQ]->offset_[valid_idx1];
            ui EP = edge_matrix[u][CVQ]->offset_[valid_idx1 + 1];
            int Times = EP - SP;
            for (int ee = 0; ee < Times; ee++)
            {
                ui vd = edge_matrix[u][CVQ]->edge_[SP + ee];
                ui vc = candidates[CVQ][vd];
                if (nodeId[vc]==1)
                    countmin++;
            }
        }
        if (countmin == 0)
        {
            minPos = ccIDX;
            break;
        }
        if (countmin < minVal)
        {
            minPos = ccIDX;
            minVal = countmin;
        }
    }
    ui tempCC = valid_candidate_idx[0][minPos];
    valid_candidate_idx[0][minPos] = valid_candidate_idx[0][idx[0]];
    valid_candidate_idx[0][idx[0]] = tempCC;
    rankV.clear();
    return;
}

void EvaluateQuery::rankMore(const Graph *query_graph, Edges ***edge_matrix, ui u, ui *&nodeId, ui **candidates, ui **valid_candidate_idx, ui *idx, ui valid_idx, ui *idx_count)
{
    vector<int> rankV;
    ui ccIDX = idx[0];
    int maxVal = 0;
    int maxPos = ccIDX;
    int countmax = 0;
    while (ccIDX < idx_count[0])
    {
        ui valid_idx1 = valid_candidate_idx[0][ccIDX];
        if (valid_idx1 != 10000000)
        {
            VertexID vC = candidates[u][valid_idx1];
            if (nodeId[vC]==0)
            {
                // add all elements to a vector. (not found yet)
                rankV.push_back(ccIDX);
            }
        }
        ccIDX++;
    }
    ui u_nbrs_count;
    for (int i = 0; i < rankV.size(); i++)
    {
        ccIDX = rankV[i];
        ui valid_idx1 = valid_candidate_idx[0][ccIDX];
        const VertexID *u_nbrs = query_graph->getVertexNeighbors(u, u_nbrs_count);
        countmax = 0;
        for (int d = 0; d < u_nbrs_count; d++)
        {
            ui CVQ = u_nbrs[d];
            ui SP = edge_matrix[u][CVQ]->offset_[valid_idx1];
            ui EP = edge_matrix[u][CVQ]->offset_[valid_idx1 + 1];
            int Times = EP - SP;
            for (int ee = 0; ee < Times; ee++)
            {
                ui vd = edge_matrix[u][CVQ]->edge_[SP + ee];
                ui vc = candidates[CVQ][vd];
                if (nodeId[vc]==0)
                    countmax++;
            }
        }
        if (countmax > maxVal)
        {
            maxPos = ccIDX;
            maxVal = countmax;
        }
    }
    ui tempCC = valid_candidate_idx[0][maxPos];
    valid_candidate_idx[0][maxPos] = valid_candidate_idx[0][idx[0]];
    valid_candidate_idx[0][idx[0]] = tempCC;
    return;
}

void EvaluateQuery::rankFinal(const Graph *query_graph, Edges ***edge_matrix, ui u, ui *&nodeId, ui **candidates, ui **valid_candidate_idx, ui *idx, ui valid_idx, ui *idx_count)
{
    vector<int> rankV;
    ui ccIDX = idx[0];

    int maxVal = 0; // FUCKKKSSSS
    int maxPos = ccIDX;
    int countmin = 0;
    int countmax = 0;
    while (ccIDX < idx_count[0])
    {
        ui valid_idx1 = valid_candidate_idx[0][ccIDX];
        if (valid_idx1 != 10000000)
        {
            VertexID vC = candidates[u][valid_idx1];
            if (nodeId[vC]==0)
            {
                // add all elements to a vector. (not found yet)
                rankV.push_back(ccIDX);
            }
        }
        ccIDX++;
    }
    ui u_nbrs_count;
    for (int i = 0; i < rankV.size(); i++)
    {
        ccIDX = rankV[i];
        ui valid_idx1 = valid_candidate_idx[0][ccIDX];
        const VertexID *u_nbrs = query_graph->getVertexNeighbors(u, u_nbrs_count);
        countmin = 0;
        countmax = 0;
        for (int d = 0; d < u_nbrs_count; d++)
        {
            ui CVQ = u_nbrs[d];
            ui SP = edge_matrix[u][CVQ]->offset_[valid_idx1];
            ui EP = edge_matrix[u][CVQ]->offset_[valid_idx1 + 1];
            int Times = EP - SP;
            for (int ee = 0; ee < Times; ee++)
            {
                ui vd = edge_matrix[u][CVQ]->edge_[SP + ee];
                ui vc = candidates[CVQ][vd];
                if (nodeId[vc]==1)
                    countmin++;
                else
                {
                    countmax++;
                }
            }
        }
        if (countmin == 0)
        {
            maxPos = ccIDX;
            break;
        }
        if (countmax > maxVal)
        {
            maxPos = ccIDX;
            maxVal = countmax;
        }
    }
    ui tempCC = valid_candidate_idx[0][maxPos];
    valid_candidate_idx[0][maxPos] = valid_candidate_idx[0][idx[0]];
    valid_candidate_idx[0][idx[0]] = tempCC;
    rankV.clear();
    return;
}

void EvaluateQuery::rankFinal1(const Graph *query_graph, Edges ***edge_matrix, ui u, ui *&nodeId, ui **candidates, ui **valid_candidate_idx, ui *idx, ui valid_idx, ui *idx_count)
{
    vector<int> rankV;
    ui ccIDX = idx[0];

    int maxVal = 0;
    int maxPos = ccIDX;
    int maxVal1 = 0;
    int maxPos1 = ccIDX;
    int countmin = 0;
    int countmax = 0;

    // get all the nodes that do not exists in solution set
    while (ccIDX < idx_count[0])
    {
        ui valid_idx1 = valid_candidate_idx[0][ccIDX];
        if (valid_idx1 != 10000000)
        {
            VertexID vC = candidates[u][valid_idx1];
            if (nodeId[vC]==0)
            {
                // add all elements to a vector. (not found yet)
                rankV.push_back(ccIDX);
            }
        }
        ccIDX++;
    }
    ui u_nbrs_count;
    for (int i = 0; i < rankV.size(); i++)
    {
        ccIDX = rankV[i];
        ui valid_idx1 = valid_candidate_idx[0][ccIDX];
        const VertexID *u_nbrs = query_graph->getVertexNeighbors(u, u_nbrs_count);
        countmin = 0;
        countmax = 0;
        // evaluate the neighboors of the node
        for (int d = 0; d < u_nbrs_count; d++)
        {
            ui CVQ = u_nbrs[d];
            ui SP = edge_matrix[u][CVQ]->offset_[valid_idx1];
            ui EP = edge_matrix[u][CVQ]->offset_[valid_idx1 + 1];
            int Times = EP - SP;
            for (int ee = 0; ee < Times; ee++)
            {
                ui vd = edge_matrix[u][CVQ]->edge_[SP + ee];
                ui vc = candidates[CVQ][vd];
                if (nodeId[vc]==1)
                    countmin++; // found duplicate
                else
                {
                    countmax++; // found unique
                }
            }
        }
        if (countmin == 0)
        { // if there is no duplicate
            if (countmax > maxVal1)
            {

                // check if among the only-uniques it has the more uniques
                maxPos1 = ccIDX;
                maxVal1 = countmax;
            }
        }
        if (countmax > maxVal)
        { // keep count the node with the more uniques
            maxPos = ccIDX;
            maxVal = countmax;
        }
    }
    if (maxVal1 != 0)
        maxPos = maxPos1; // we care about positions to change.
    ui tempCC = valid_candidate_idx[0][maxPos];
    valid_candidate_idx[0][maxPos] = valid_candidate_idx[0][idx[0]];
    valid_candidate_idx[0][idx[0]] = tempCC;
    rankV.clear();
    return;
}
void EvaluateQuery::rankRandom(const Graph *query_graph, Edges ***edge_matrix, ui u, ui *&nodeId, ui **candidates, ui **valid_candidate_idx, ui *idx, ui valid_idx, ui *idx_count)
{   
    ui ccIDX = idx[0];
    std::random_device rd;  // Obtain a random seed
    std::mt19937 gen(rd()); // Seed the generator with random_device
    std::uniform_int_distribution<> distr(ccIDX, idx_count[0]-1); // Define the range [0, N)
    
    // Generate a random number
    int maxVal = 0;
    int maxPos = ccIDX;
    int maxVal1 = 0;
    int maxPos1 = ccIDX;
    int countmin = 0;
    int countmax = 0;
    
    // get all the nodes that do not exists in solution set
    maxPos=distr(gen);
    if (maxPos>idx_count[0]-1 )
    cout<<"bigger"<<endl;
    if (maxPos<ccIDX)
    cout<< "smaller"<<endl;
    ui tempCC = valid_candidate_idx[0][maxPos];
    valid_candidate_idx[0][maxPos] = valid_candidate_idx[0][idx[0]];
    valid_candidate_idx[0][idx[0]] = tempCC;
    return;
}

std::function<bool(std::pair<std::pair<VertexID, ui>, ui>, std::pair<std::pair<VertexID, ui>, ui>)>
    EvaluateQuery::extendable_vertex_compare = [](std::pair<std::pair<VertexID, ui>, ui> l, std::pair<std::pair<VertexID, ui>, ui> r)
{
    if (l.first.second == 1 && r.first.second != 1)
    {
        return true;
    }
    else if (l.first.second != 1 && r.first.second == 1)
    {
        return false;
    }
    else
    {
        return l.second > r.second;
    }
};

void EvaluateQuery::generateBN(const Graph *query_graph, ui *order, ui *pivot, ui **&bn, ui *&bn_count)
{
    ui query_vertices_num = query_graph->getVerticesCount();
    bn_count = new ui[query_vertices_num];
    std::fill(bn_count, bn_count + query_vertices_num, 0);
    bn = new ui *[query_vertices_num];
    for (ui i = 0; i < query_vertices_num; ++i)
    {
        bn[i] = new ui[query_vertices_num];
    }

    std::vector<bool> visited_vertices(query_vertices_num, false);
    visited_vertices[order[0]] = true;
    for (ui i = 1; i < query_vertices_num; ++i)
    {
        VertexID vertex = order[i];

        ui nbrs_cnt;
        const ui *nbrs = query_graph->getVertexNeighbors(vertex, nbrs_cnt);
        for (ui j = 0; j < nbrs_cnt; ++j)
        {
            VertexID nbr = nbrs[j];

            if (visited_vertices[nbr] && nbr != pivot[i])
            {
                bn[i][bn_count[i]++] = nbr;
            }
        }

        visited_vertices[vertex] = true;
    }
}

void EvaluateQuery::generateBN(const Graph *query_graph, ui *order, ui **&bn, ui *&bn_count)
{
    ui query_vertices_num = query_graph->getVerticesCount();
    bn_count = new ui[query_vertices_num];
    std::fill(bn_count, bn_count + query_vertices_num, 0);
    bn = new ui *[query_vertices_num];
    for (ui i = 0; i < query_vertices_num; ++i)
    {
        bn[i] = new ui[query_vertices_num];
    }

    std::vector<bool> visited_vertices(query_vertices_num, false);
    visited_vertices[order[0]] = true;
    for (ui i = 1; i < query_vertices_num; ++i)
    {
        VertexID vertex = order[i];

        ui nbrs_cnt;
        const ui *nbrs = query_graph->getVertexNeighbors(vertex, nbrs_cnt);
        for (ui j = 0; j < nbrs_cnt; ++j)
        {
            VertexID nbr = nbrs[j];

            if (visited_vertices[nbr])
            {
                bn[i][bn_count[i]++] = nbr;
            }
        }

        visited_vertices[vertex] = true;
    }
}

void EvaluateQuery::generateBNLM(const Graph *query_graph, ui *order, ui **&bn, ui *&bn_count)
{
    ui query_vertices_num = query_graph->getVerticesCount();
    std::vector<bool> visited_vertices(query_vertices_num, false);
    visited_vertices[order[0]] = true;
    std::fill(bn_count, bn_count + query_vertices_num, 0);
    for (ui i = 1; i < query_vertices_num; ++i)
    {
        VertexID vertex = order[i];

        ui nbrs_cnt;
        const ui *nbrs = query_graph->getVertexNeighbors(vertex, nbrs_cnt);
        for (ui j = 0; j < nbrs_cnt; ++j)
        {
            VertexID nbr = nbrs[j];

            if (visited_vertices[nbr])
            {
                bn[i][bn_count[i]++] = nbr;
            }
        }

        visited_vertices[vertex] = true;
    }
}
size_t
EvaluateQuery::exploreGraph(const Graph *data_graph, const Graph *query_graph, Edges ***edge_matrix, ui **candidates,
                            ui *candidates_count, ui *order, ui *pivot, size_t output_limit_num, size_t &call_count,
                            const std::unordered_map<VertexID, std::pair<std::set<VertexID>, std::set<VertexID>>> &ordered_constraints)
{

    ui **bn;
    ui *bn_count;
    generateBN(query_graph, order, pivot, bn, bn_count);

    ui *idx;
    ui *idx_count;
    ui *embedding;
    ui *idx_embedding;
    ui *temp_buffer;
    ui **valid_candidate_idx;
    bool *visited_vertices;
    allocateBuffer(data_graph, query_graph, candidates_count, idx, idx_count, embedding, idx_embedding,
                   temp_buffer, valid_candidate_idx, visited_vertices);

    size_t embedding_cnt = 0;
    int cur_depth = 0;
    ui max_depth = query_graph->getVerticesCount();
    VertexID start_vertex = order[0];

    idx[cur_depth] = 0;
    idx_count[cur_depth] = candidates_count[start_vertex];

    for (ui i = 0; i < idx_count[cur_depth]; ++i)
    {
        valid_candidate_idx[cur_depth][i] = i;
    }

    while (true)
    {
        while (idx[cur_depth] < idx_count[cur_depth])
        {
            ui valid_idx = valid_candidate_idx[cur_depth][idx[cur_depth]];
            VertexID u = order[cur_depth];
            VertexID v = candidates[u][valid_idx];

            embedding[u] = v;
            idx_embedding[u] = valid_idx;
            visited_vertices[v] = true;
            idx[cur_depth] += 1;

            if (cur_depth == max_depth - 1)
            {
                embedding_cnt += 1;
                visited_vertices[v] = false;
                if (embedding_cnt >= output_limit_num)
                {
                    goto EXIT;
                }
            }
            else
            {
                call_count += 1;
                cur_depth += 1;
                idx[cur_depth] = 0;
                generateValidCandidateIndex(data_graph, cur_depth, embedding, idx_embedding, idx_count,
                                            valid_candidate_idx,
                                            edge_matrix, visited_vertices, bn, bn_count, order, pivot, candidates);
                pruneCandidatesIndexBySymmetryBreaking(cur_depth, embedding, order, candidates_count, candidates,
                                                       idx_count, valid_candidate_idx, ordered_constraints);
            }
        }

        cur_depth -= 1;
        if (cur_depth < 0)
            break;
        else
            visited_vertices[embedding[order[cur_depth]]] = false;
    }

EXIT:
    releaseBuffer(max_depth, idx, idx_count, embedding, idx_embedding, temp_buffer, valid_candidate_idx,
                  visited_vertices,
                  bn, bn_count);

    return embedding_cnt;
}
void EvaluateQuery::allocateBufferLM(const Graph *data_graph, const Graph *query_graph, ui *candidates_count, ui *&idx,
                                   ui *&idx_count, ui *&embedding, ui *&idx_embedding, ui *&temp_buffer,
                                   ui **&valid_candidate_idx, bool *&visited_vertices)
{
    ui query_vertices_num = query_graph->getVerticesCount();
    ui data_vertices_num = data_graph->getVerticesCount();
    std::fill(idx, idx + query_vertices_num, std::numeric_limits<ui>::max());
    std::fill(idx_count, idx_count + query_vertices_num, std::numeric_limits<ui>::max());
    std::fill(embedding, embedding + query_vertices_num, std::numeric_limits<ui>::max());
    std::fill(idx_embedding, idx_embedding + query_vertices_num, std::numeric_limits<ui>::max());
    //idx = new ui[query_vertices_num];
    //idx_count = new ui[query_vertices_num];
    //embedding = new ui[query_vertices_num];
    //idx_embedding = new ui[query_vertices_num];
    std::fill(visited_vertices, visited_vertices + data_vertices_num, false);
}
void EvaluateQuery::allocateBuffer(const Graph *data_graph, const Graph *query_graph, ui *candidates_count, ui *&idx,
                                   ui *&idx_count, ui *&embedding, ui *&idx_embedding, ui *&temp_buffer,
                                   ui **&valid_candidate_idx, bool *&visited_vertices)
{
    ui query_vertices_num = query_graph->getVerticesCount();
    ui data_vertices_num = data_graph->getVerticesCount();
    ui max_candidates_num = candidates_count[0];

    for (ui i = 1; i < query_vertices_num; ++i)
    {
        VertexID cur_vertex = i;
        ui cur_candidate_num = candidates_count[cur_vertex];

        if (cur_candidate_num > max_candidates_num)
        {
            max_candidates_num = cur_candidate_num;
        }
    }

    idx = new ui[query_vertices_num];
    idx_count = new ui[query_vertices_num];
    embedding = new ui[query_vertices_num];
    idx_embedding = new ui[query_vertices_num];
    visited_vertices = new bool[data_vertices_num];
    temp_buffer = new ui[max_candidates_num];
    valid_candidate_idx = new ui *[query_vertices_num];
    for (ui i = 0; i < query_vertices_num; ++i)
    {
        valid_candidate_idx[i] = new ui[max_candidates_num];
    }

    std::fill(visited_vertices, visited_vertices + data_vertices_num, false);
}

void EvaluateQuery::generateValidCandidateIndex(const Graph *data_graph, ui depth, ui *embedding, ui *idx_embedding,
                                                ui *idx_count, ui **valid_candidate_index, Edges ***edge_matrix,
                                                bool *visited_vertices, ui **bn, ui *bn_cnt, ui *order, ui *pivot,
                                                ui **candidates)
{
    VertexID u = order[depth];
    VertexID pivot_vertex = pivot[depth];
    ui idx_id = idx_embedding[pivot_vertex];
    Edges &edge = *edge_matrix[pivot_vertex][u];
    ui count = edge.offset_[idx_id + 1] - edge.offset_[idx_id];
    ui *candidate_idx = edge.edge_ + edge.offset_[idx_id];

    ui valid_candidate_index_count = 0;

    if (bn_cnt[depth] == 0)
    {
        for (ui i = 0; i < count; ++i)
        {
            ui temp_idx = candidate_idx[i];
            VertexID temp_v = candidates[u][temp_idx];

            if (!visited_vertices[temp_v])
                valid_candidate_index[depth][valid_candidate_index_count++] = temp_idx;
        }
    }
    else
    {
        for (ui i = 0; i < count; ++i)
        {
            ui temp_idx = candidate_idx[i];
            VertexID temp_v = candidates[u][temp_idx];

            if (!visited_vertices[temp_v])
            {
                bool valid = true;

                for (ui j = 0; j < bn_cnt[depth]; ++j)
                {
                    VertexID u_bn = bn[depth][j];
                    VertexID u_bn_v = embedding[u_bn];

                    if (!data_graph->checkEdgeExistence(temp_v, u_bn_v))
                    {
                        valid = false;
                        break;
                    }
                }

                if (valid)
                    valid_candidate_index[depth][valid_candidate_index_count++] = temp_idx;
            }
        }
    }

    idx_count[depth] = valid_candidate_index_count;
}

void EvaluateQuery::releaseBuffer(ui query_vertices_num, ui *idx, ui *idx_count, ui *embedding, ui *idx_embedding,
                                  ui *temp_buffer, ui **valid_candidate_idx, bool *visited_vertices, ui **bn,
                                  ui *bn_count)
{
    delete[] idx;
    delete[] idx_count;
    delete[] embedding;
    delete[] idx_embedding;
    delete[] visited_vertices;
    delete[] bn_count;
    delete[] temp_buffer;
    for (ui i = 0; i < query_vertices_num; ++i)
    {
        delete[] valid_candidate_idx[i];
        delete[] bn[i];
    }

    delete[] valid_candidate_idx;
    delete[] bn;
}

enumResult
EvaluateQuery::LFTJDLS(const Graph *data_graph, const Graph *query_graph,ui *&nodeId, Edges ***edge_matrix, ui **candidates,
                    ui *candidates_count, ui *order, size_t output_limit_num, size_t &call_count, size_t &valid_vtx_cnt, int TimeL, int FairT,
                    const std::unordered_map<VertexID, std::pair<std::set<VertexID>, std::set<VertexID>>> &ordered_constraints)
    {

    auto start = std::chrono::high_resolution_clock::now();
    std::vector<VertexID> priority_neighbors; // To store neighbors that are is_used[nbrs[i]]
    std::vector<VertexID> secondary_neighbors; // To store neighbors that are not is_used[nbrs[i]]
    enumResult s;
    ui **bn;
    ui *bn_count;
    generateBN(query_graph, order, bn, bn_count);
    int UNPM = 0;
    ui *idx;
    ui *idx_count;
    ui *embedding;
    ui *idx_embedding;
    ui *temp_buffer;
    ui **valid_candidate_idx;
    int i=0;
    bool *visited_vertices;
        int count1=0;
    int count2=0;
    int K_size2=0;
        ui tempstore=0;
    ui nbr_cnt;
    allocateBuffer(data_graph, query_graph, candidates_count, idx, idx_count, embedding, idx_embedding,
                   temp_buffer, valid_candidate_idx, visited_vertices);

    size_t embedding_cnt = 0;
    int cur_depth = 0;
    ui max_depth = query_graph->getVerticesCount();
    VertexID start_vertex = order[0];
    ui order_matrix[max_depth];
    bool is_used[max_depth];
    bool is_used_copy[max_depth];
    //unordered_set<ui> EmbSum;
    bool Match_BA[max_depth] = {false};
    idx[cur_depth] = 0;
    idx_count[cur_depth] = candidates_count[start_vertex];

    for (ui i = 0; i < idx_count[cur_depth]; ++i)
    {
        valid_candidate_idx[cur_depth][i] = i;
    }

    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> ancestors;
    computeAncestor(query_graph, bn, bn_count, order, ancestors);
    AncestorWithSymmetryBreaking(query_graph, ancestors, ordered_constraints);
    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> vec_failing_set(max_depth);
    std::unordered_map<VertexID, VertexID> reverse_embedding;
    reverse_embedding.reserve(MAXIMUM_QUERY_GRAPH_SIZE * 2);
    bool hasNext = true;
    ui K_size=1;
    TimeL = TimeL * 1000;
    vector<int> combination;
    double ens = 0;
    while (true)
    {
        while (idx[cur_depth] < idx_count[cur_depth])
        {
            Match_BA[cur_depth] = false;
            ui valid_idx = valid_candidate_idx[cur_depth][idx[cur_depth]];
            VertexID u = order[cur_depth];
            VertexID v = candidates[u][valid_idx];
            auto end = std::chrono::high_resolution_clock::now();
            ens = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
            if (ens > TimeL)
            { // 1000 1 sec
                goto EXIT;
            }
            if (visited_vertices[v])
            {
                idx[cur_depth] += 1;
                vec_failing_set[cur_depth] = ancestors[u];
                vec_failing_set[cur_depth] |= ancestors[reverse_embedding[v]];
                vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
                continue;
            }
            if(nodeId[v]==1){
                idx[cur_depth] += 1;
                //vec_failing_set[cur_depth] = ancestors[u];
                //vec_failing_set[cur_depth] |= ancestors[reverse_embedding[v]];
                
                if(cur_depth - 1>0){
                 vec_failing_set[cur_depth] |= ancestors[order[cur_depth]];       
                vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];

                }
                
                continue;
            }

            embedding[u] = v;
            idx_embedding[u] = valid_idx;
            visited_vertices[v] = true;
            idx[cur_depth] += 1;

            reverse_embedding[v] = u;
            if (cur_depth == max_depth - 1)
            {
                UNPM = 0;
                if(nodeId[embedding[order[cur_depth]]]==0){
                    UNPM++;
                    nodeId[embedding[order[cur_depth]]]=1;
                }
                embedding_cnt += 1;
                visited_vertices[v] = false;
                int ao = cur_depth; //-1;
                ui vqo = order[ao];
                nodeId[v]=1;
                while (ao > 0)
                {                
                    if(nodeId[embedding[vqo]]==0){
                    UNPM++;
                    nodeId[embedding[vqo]]=1;
                }
                //vec_failing_set[ao].set();
                //vec_failing_set[ao - 1] |= vec_failing_set[ao];
                //reverse_embedding.erase(embedding[embedding[vqo]]);
                    idx[ao] = idx_count[ao];
                    // Match_BA[ao] = true;
                    ao--;
                    vqo = order[ao];
                }
                nodeId[embedding[vqo]]=1;


                reverse_embedding.erase(embedding[u]);
                vec_failing_set[cur_depth].set();
                vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
                if (embedding_cnt >= output_limit_num)
                {
                    goto EXIT;
                }
            }
            else
            {
                call_count += 1;
                cur_depth += 1;

                idx[cur_depth] = 0;
                generateValidCandidateIndex(cur_depth, idx_embedding, idx_count, valid_candidate_idx, edge_matrix, bn,
                                            bn_count, order, temp_buffer);

                
                if (idx_count[cur_depth] == 0)
                {
                    vec_failing_set[cur_depth - 1] = ancestors[order[cur_depth]];
                }
                else
                {
                    vec_failing_set[cur_depth - 1].reset();
                }
            }
        }

        cur_depth -= 1;
        if (cur_depth < 0)
            break;
        else
        {
            VertexID u = order[cur_depth];
            reverse_embedding.erase(embedding[u]);
            if (cur_depth != 0)
            {
                if (!vec_failing_set[cur_depth].test(u))
                {
                    vec_failing_set[cur_depth - 1] = vec_failing_set[cur_depth];
                    idx[cur_depth] = idx_count[cur_depth];
                }
                else
                {
                    vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
                }
            }
            visited_vertices[embedding[u]] = false;

        }
    }
   
    hasNext = true;
     K_size=1;
        //for (int i=0;i<max_depth;i++)
        //cout<<order[i];
        //cout<<""<<endl;

    while (K_size<max_depth){
        vector<int> indices(K_size);        
     
        for (ui i = 0; i < K_size; ++i) {
        indices[i] = i;
        }
        hasNext = true;
        while (hasNext==true) {
            count1=0;
            count2=0;
            K_size2=0;
            tempstore=0;
            nbr_cnt=0;
            for (int i=0;i<max_depth;i++){
            is_used[i]=false;
            is_used_copy[i]=false;
            }
            combination.clear();
            
            for (ui i = 0; i < K_size; ++i) {
                combination.push_back(order[indices[i]]);
                is_used[order[indices[i]]] = true;
            }

            order_matrix[0] = combination[0];
            is_used_copy[combination[0]]=true;
            count1=1;
            count2=0;
            //count1++;
        /*
        while(count1<max_depth){
            const VertexID *nbrs = query_graph->getVertexNeighbors(order_matrix[count2], nbr_cnt);
            for (int i=0;i<nbr_cnt;i++){
                if(is_used_copy[nbrs[i]]==false){
                    is_used_copy[nbrs[i]]=true;
                    order_matrix[count1]=nbrs[i];
                    count1++;
                }
            }count2++;
        }*/
       
        while (count1 < max_depth) {
    const VertexID *nbrs = query_graph->getVertexNeighbors(order_matrix[count2], nbr_cnt);
    
        priority_neighbors.clear();
    secondary_neighbors.clear();
    for (int i = 0; i < nbr_cnt; i++) {
        if (is_used_copy[nbrs[i]] == false) {
            if (is_used[nbrs[i]]) {
                priority_neighbors.push_back(nbrs[i]); // Add to priority group
            } else {
                secondary_neighbors.push_back(nbrs[i]); // Add to secondary group
            }
        }
    }

    // Process priority neighbors first
    for (const auto &nbr : priority_neighbors) {
        if (count1 >= max_depth) break; // Ensure we don't exceed max_depth
        is_used_copy[nbr] = true;
        order_matrix[count1] = nbr;
        count1++;
    }

    // Process secondary neighbors next
    for (const auto &nbr : secondary_neighbors) {
        if (count1 >= max_depth) break; // Ensure we don't exceed max_depth
        is_used_copy[nbr] = true;
        order_matrix[count1] = nbr;
        count1++;
    }

    count2++;

}
        
        count1=max_depth;
        count2=0;
        i=1;
        while (i<count1){
            if(query_graph->getVertexDegree(order_matrix[i]==1)){
                tempstore=order_matrix[i];
                for (int j=i+1;j<max_depth;j++){
                    order_matrix[j-1]=order_matrix[j];
                }
                order_matrix[max_depth-1]=tempstore;
                count1--;
            }
            else{i++;}
        }   

        K_size2=1;
        i=0;
        while(is_used[order_matrix[i]]==true){
            K_size2++;
            i++;
        }
        //cout<<" #"<<K_size2<<endl;
        i=0;
        


        //generateBN(query_graph, order_matrix, bn, bn_count);
        generateBNLM(query_graph, order_matrix, bn, bn_count);
        computeAncestor(query_graph, bn, bn_count, order_matrix, ancestors);
        allocateBufferLM(data_graph, query_graph, candidates_count, idx, idx_count, embedding, idx_embedding,
                   temp_buffer, valid_candidate_idx, visited_vertices);
        //std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> vec_failing_set(max_depth);
        vec_failing_set.clear();
        vec_failing_set.resize(max_depth);
        start_vertex = order_matrix[0];
        cur_depth=0;
        idx[cur_depth] = 0;
        idx_count[cur_depth] = candidates_count[start_vertex];
        Match_BA[start_vertex] = false;
        for (ui i = 0; i < idx_count[cur_depth]; ++i)
        {
        valid_candidate_idx[cur_depth][i] = i;
        }
        while (true)
    {
        while (idx[cur_depth] < idx_count[cur_depth])
        {   
            
            ui valid_idx = valid_candidate_idx[cur_depth][idx[cur_depth]];
            VertexID u = order_matrix[cur_depth];
            VertexID v = candidates[u][valid_idx];
            auto end = std::chrono::high_resolution_clock::now();
            ens = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
            if (ens > TimeL)
            { // 1000 1 sec
                goto EXIT;
            }
            if(nodeId[v]==0&&is_used[u]){
                idx[cur_depth] += 1;
                vec_failing_set[cur_depth] = ancestors[u];
                vec_failing_set[cur_depth] |= ancestors[order_matrix[i]];
                if(cur_depth - 1>0)
                vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
                continue;
            }
            if(nodeId[v]==1&& (!is_used[u])){
                idx[cur_depth] += 1;
               vec_failing_set[cur_depth] = ancestors[u];
                vec_failing_set[cur_depth] |= ancestors[order_matrix[i]];
                if(cur_depth - 1>0)
                vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
                continue;
            }
            if (visited_vertices[v])
            {
                idx[cur_depth] += 1;
                vec_failing_set[cur_depth] = ancestors[u];
                vec_failing_set[cur_depth] |= ancestors[reverse_embedding[v]];
                vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
                continue;
            }



            embedding[u] = v;
            idx_embedding[u] = valid_idx;
            visited_vertices[v] = true;
            idx[cur_depth] += 1;

            reverse_embedding[v] = u;
            if (cur_depth == max_depth - 1)
            {   
                UNPM = 0;
                if(nodeId[embedding[order_matrix[cur_depth]]]==0){
                    UNPM++;
                    nodeId[embedding[order_matrix[cur_depth]]]=1;
                }
                embedding_cnt += 1;
                visited_vertices[v] = false;
                int ao = cur_depth; //-1;
                ui vqo = order_matrix[ao];

                while (ao > 0)
                {                
                    if(nodeId[embedding[vqo]]==0){
                    UNPM++;
                    nodeId[embedding[vqo]]=1;
                    //idx[ao] = idx_count[ao];
                }
                   if(ao>K_size2){
                    idx[ao] = idx_count[ao];

                    vec_failing_set[ao].set();
                    vec_failing_set[ao - 1] |= vec_failing_set[ao];
                    //vec_failing_set[ao - 1] = ancestors[order_matrix[ao]];
                    reverse_embedding.erase(embedding[ao]);
                    
                    
                   }
                    vec_failing_set[cur_depth].set();
                    vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
                    
                    Match_BA[vqo]=true;
                    ao--;
                    vqo = order_matrix[ao];
                }
                Match_BA[vqo]=true;
                nodeId[embedding[vqo]]=1;
                
                reverse_embedding.erase(embedding[u]);
                vec_failing_set[cur_depth].set();
                vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
                if (embedding_cnt >= output_limit_num)
                {
                    goto EXIT;
                }

            }
            else
            {
                call_count += 1;
                cur_depth += 1;   
                //Match_BA[order_matrix[cur_depth]]=false;          
                idx[cur_depth] = 0;
                generateValidCandidateIndex(cur_depth, idx_embedding, idx_count, valid_candidate_idx, edge_matrix, bn,
                                            bn_count, order_matrix, temp_buffer);
            
            //if(is_used[u]){
             //   while(nodeId[v]==0&&idx[cur_depth]<idx_count[cur_depth])
            //    idx[cur_depth] += 1;
            //}else{
           //     while(nodeId[v]==1&&idx[cur_depth]<idx_count[cur_depth])
            //    idx[cur_depth] += 1;
            //}
            //if(idx[cur_depth]==idx_count[cur_depth]){
            //    idx_count[cur_depth]=0;
            //}
            

                if (idx_count[cur_depth] == 0)
                {
                    vec_failing_set[cur_depth - 1] = ancestors[order_matrix[cur_depth]];
                }  
                else
                {   Match_BA[u]==false;
                    vec_failing_set[cur_depth - 1].reset();
                }
            }
        }
        cur_depth -= 1;
        if (cur_depth < 0)
            break;
        else
        {
            VertexID u = order_matrix[cur_depth];
            reverse_embedding.erase(embedding[u]);
            if (cur_depth != 0)
            {
                //if (!vec_failing_set[cur_depth].test(u)&&Match_BA[u]==false)
                if (!vec_failing_set[cur_depth].test(u))
                {
                    vec_failing_set[cur_depth - 1] = vec_failing_set[cur_depth];
                    //if(idx[cur_depth]!=idx_count[cur_depth])
                    //cout<<"in"<<idx[cur_depth]<<"-"<<idx_count[cur_depth];
                    idx[cur_depth] = idx_count[cur_depth];
                    
                }
                else
                {
                    vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
                }
            }
            visited_vertices[embedding[u]] = false;

        }
    }

        hasNext = nextCombination(indices, max_depth, K_size);
    }
        K_size++;
    }


EXIT:
    releaseBuffer(max_depth, idx, idx_count, embedding, idx_embedding, temp_buffer, valid_candidate_idx,
                  visited_vertices,
                  bn, bn_count);

        s.embedding_cnt = embedding_cnt;
        ui countSMU=0;
        for (int i=0;i<data_graph->getVerticesCount();i++){
            countSMU+=nodeId[i];
        }
        s.Can_embed = countSMU;
    // s.topk=greedysum;
        return s;
}

enumResult
EvaluateQuery::LFTJDLSKOPT(const Graph *data_graph, const Graph *query_graph,ui *&nodeId, Edges ***edge_matrix, ui **candidates,
                    ui *candidates_count, ui *order, size_t output_limit_num, size_t &call_count, size_t &valid_vtx_cnt, int TimeL, int FairT,
                    const std::unordered_map<VertexID, std::pair<std::set<VertexID>, std::set<VertexID>>> &ordered_constraints)
    {


    auto start = std::chrono::high_resolution_clock::now();
    enumResult s;
    ui **bn;
    ui *bn_count;
    generateBN(query_graph, order, bn, bn_count);
    int UNPM = 0;
    ui *idx;
    ui *idx_count;
    ui *embedding;
    ui *idx_embedding;
    ui *temp_buffer;
    ui **valid_candidate_idx;
    int i=0;
    bool *visited_vertices;
        int count1=0;
    int count2=0;
    int K_size2=0;
        ui tempstore=0;
    ui nbr_cnt;
    allocateBuffer(data_graph, query_graph, candidates_count, idx, idx_count, embedding, idx_embedding,
                   temp_buffer, valid_candidate_idx, visited_vertices);

    size_t embedding_cnt = 0;
    int cur_depth = 0;
    ui max_depth = query_graph->getVerticesCount();
    VertexID start_vertex = order[0];
    ui order_matrix[max_depth];
    bool is_used[max_depth];
    bool is_used_copy[max_depth];
    //unordered_set<ui> EmbSum;
    bool Match_BA[max_depth] = {false};
    idx[cur_depth] = 0;
    idx_count[cur_depth] = candidates_count[start_vertex];

    for (ui i = 0; i < idx_count[cur_depth]; ++i)
    {
        valid_candidate_idx[cur_depth][i] = i;
    }

    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> ancestors;
    computeAncestor(query_graph, bn, bn_count, order, ancestors);
    AncestorWithSymmetryBreaking(query_graph, ancestors, ordered_constraints);
    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> vec_failing_set(max_depth);
    std::unordered_map<VertexID, VertexID> reverse_embedding;
    reverse_embedding.reserve(MAXIMUM_QUERY_GRAPH_SIZE * 2);
    bool hasNext = true;
    ui K_size=1;
    TimeL = TimeL * 1000;
    vector<int> combination;
    double ens = 0;
    int qsiz = query_graph->getVerticesCount();
    int ksize = 10000;
    int SolPos[ksize+1][qsiz];
    std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, std::greater<>> pq;
    int tempPos[qsiz];
    pq.push({0,0});
    int heapcount = 0;
    while (true)
    {
        while (idx[cur_depth] < idx_count[cur_depth])
        {
            Match_BA[cur_depth] = false;
            ui valid_idx = valid_candidate_idx[cur_depth][idx[cur_depth]];
            VertexID u = order[cur_depth];
            VertexID v = candidates[u][valid_idx];
            auto end = std::chrono::high_resolution_clock::now();
            ens = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
            if (ens > TimeL)
            { // 1000 1 sec
                goto EXIT;
            }
            if (visited_vertices[v])
            {
                idx[cur_depth] += 1;
                vec_failing_set[cur_depth] = ancestors[u];
                vec_failing_set[cur_depth] |= ancestors[reverse_embedding[v]];
                vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
                continue;
            }
            if(nodeId[v]==1){
                idx[cur_depth] += 1;
                vec_failing_set[cur_depth] = ancestors[u];
                //vec_failing_set[cur_depth] |= ancestors[reverse_embedding[v]];
                vec_failing_set[cur_depth] |= ancestors[order[i]];
                if(cur_depth - 1>0)
                vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
                continue;
            }

            embedding[u] = v;
            idx_embedding[u] = valid_idx;
            visited_vertices[v] = true;
            idx[cur_depth] += 1;

            reverse_embedding[v] = u;
            if (cur_depth == max_depth - 1)
            {   
                UNPM = 0;
                if(nodeId[embedding[order[cur_depth]]]==0){
                    UNPM++;
                    nodeId[embedding[order[cur_depth]]]=1;
                }
                embedding_cnt += 1;
                visited_vertices[v] = false;
                int ao = cur_depth; //-1;
                tempPos[ao]=embedding[order[cur_depth]];
                ui vqo = order[ao];
                nodeId[v]=1;
                while (ao > 0)
                {                
                    if(nodeId[embedding[vqo]]==0){
                    UNPM++;
                    nodeId[embedding[vqo]]=1;
                }
                //vec_failing_set[ao].set();
                //vec_failing_set[ao - 1] |= vec_failing_set[ao];
                //reverse_embedding.erase(embedding[embedding[vqo]]);
                   tempPos[ao]=embedding[vqo];
                    idx[ao] = idx_count[ao];
                    // Match_BA[ao] = true;
                    ao--;
                    vqo = order[ao];
                }
                tempPos[ao]=embedding[vqo];
                if(nodeId[embedding[vqo]]==0){
                    UNPM++;
                    nodeId[embedding[vqo]]=1;
                    }

                if (UNPM > pq.top().first)
                            {                                       
                                int addPos=0;
                                if (heapcount >= ksize)
                                {
                                   auto [value, index] = pq.top();      
                                    pq.pop();
                                    pq.push({UNPM,index});
                                    addPos=index;
                                    
                                }else{
                                     pq.push({UNPM,heapcount});
                                     addPos=heapcount;            
                                }
                            for (int aa=0;aa<qsiz;aa++){
                                SolPos[addPos][aa]=tempPos[aa];
                                    }
                                heapcount++;
                            }

                reverse_embedding.erase(embedding[u]);
                vec_failing_set[cur_depth].set();
                vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
                if (embedding_cnt >= output_limit_num)
                {
                    goto EXIT;
                }
            }
            else
            {
                call_count += 1;
                cur_depth += 1;

                idx[cur_depth] = 0;
                generateValidCandidateIndex(cur_depth, idx_embedding, idx_count, valid_candidate_idx, edge_matrix, bn,
                                            bn_count, order, temp_buffer);

                
                if (idx_count[cur_depth] == 0)
                {
                    vec_failing_set[cur_depth - 1] = ancestors[order[cur_depth]];
                }
                else
                {
                    vec_failing_set[cur_depth - 1].reset();
                }
            }
        }

        cur_depth -= 1;
        if (cur_depth < 0)
            break;
        else
        {
            VertexID u = order[cur_depth];
            reverse_embedding.erase(embedding[u]);
            if (cur_depth != 0)
            {
                if (!vec_failing_set[cur_depth].test(u))
                {
                    vec_failing_set[cur_depth - 1] = vec_failing_set[cur_depth];
                    idx[cur_depth] = idx_count[cur_depth];
                }
                else
                {
                    vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
                }
            }
            visited_vertices[embedding[u]] = false;

        }
    }
   
    hasNext = true;
     K_size=1;
        //for (int i=0;i<max_depth;i++)
        //cout<<order[i];
        //cout<<""<<endl;

    while (K_size<max_depth){
        //cout<<K_size<<endl;
        vector<int> indices(K_size);        
     
        for (ui i = 0; i < K_size; ++i) {
        indices[i] = i;
        }
        hasNext = true;
        while (hasNext==true) {
            count1=0;
            count2=0;
            K_size2=0;
            tempstore=0;
            nbr_cnt=0;
            for (int i=0;i<max_depth;i++){
            is_used[i]=false;
            is_used_copy[i]=false;
            }
            combination.clear();
            
            for (ui i = 0; i < K_size; ++i) {
                combination.push_back(order[indices[i]]);
                is_used[order[indices[i]]] = true;
            }

            order_matrix[0] = combination[0];
            is_used_copy[combination[0]]=true;
            count1=1;
            count2=0;
            //count1++;
        /*
        while(count1<max_depth){
            const VertexID *nbrs = query_graph->getVertexNeighbors(order_matrix[count2], nbr_cnt);
            for (int i=0;i<nbr_cnt;i++){
                if(is_used_copy[nbrs[i]]==false){
                    is_used_copy[nbrs[i]]=true;
                    order_matrix[count1]=nbrs[i];
                    count1++;
                }
            }count2++;
        }*/
       
        while (count1 < max_depth) {
    const VertexID *nbrs = query_graph->getVertexNeighbors(order_matrix[count2], nbr_cnt);
    std::vector<VertexID> priority_neighbors; // To store neighbors that are is_used[nbrs[i]]
    std::vector<VertexID> secondary_neighbors; // To store neighbors that are not is_used[nbrs[i]]

    for (int i = 0; i < nbr_cnt; i++) {
        if (is_used_copy[nbrs[i]] == false) {
            if (is_used[nbrs[i]]) {
                priority_neighbors.push_back(nbrs[i]); // Add to priority group
            } else {
                secondary_neighbors.push_back(nbrs[i]); // Add to secondary group
            }
        }
    }

    // Process priority neighbors first
    for (const auto &nbr : priority_neighbors) {
        if (count1 >= max_depth) break; // Ensure we don't exceed max_depth
        is_used_copy[nbr] = true;
        order_matrix[count1] = nbr;
        count1++;
    }

    // Process secondary neighbors next
    for (const auto &nbr : secondary_neighbors) {
        if (count1 >= max_depth) break; // Ensure we don't exceed max_depth
        is_used_copy[nbr] = true;
        order_matrix[count1] = nbr;
        count1++;
    }

    count2++;
}
        
        count1=max_depth;
        count2=0;
        i=1;
        while (i<count1){
            if(query_graph->getVertexDegree(order_matrix[i]==1)){
                tempstore=order_matrix[i];
                for (int j=i+1;j<max_depth;j++){
                    order_matrix[j-1]=order_matrix[j];
                }
                order_matrix[max_depth-1]=tempstore;
                count1--;
            }
            else{i++;}
        }   

        K_size2=1;
        i=0;
        while(is_used[order_matrix[i]]==true){
            K_size2++;
            i++;
        }
        //cout<<" #"<<K_size2<<endl;
        i=0;
        


        generateBNLM(query_graph, order_matrix, bn, bn_count);
        computeAncestor(query_graph, bn, bn_count, order_matrix, ancestors);
        allocateBufferLM(data_graph, query_graph, candidates_count, idx, idx_count, embedding, idx_embedding,
                   temp_buffer, valid_candidate_idx, visited_vertices);
        std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> vec_failing_set(max_depth);
        start_vertex = order_matrix[0];
        cur_depth=0;
        idx[cur_depth] = 0;
        idx_count[cur_depth] = candidates_count[start_vertex];
        Match_BA[start_vertex] = false;
        for (ui i = 0; i < idx_count[cur_depth]; ++i)
        {
        valid_candidate_idx[cur_depth][i] = i;
        }
        while (true)
    {
        while (idx[cur_depth] < idx_count[cur_depth])
        {   
            
            ui valid_idx = valid_candidate_idx[cur_depth][idx[cur_depth]];
            VertexID u = order_matrix[cur_depth];
            VertexID v = candidates[u][valid_idx];
            auto end = std::chrono::high_resolution_clock::now();
            ens = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
            if (ens > TimeL)
            { // 1000 1 sec
                goto EXIT;
            }
            if(nodeId[v]==0&&is_used[u]){
                idx[cur_depth] += 1;
                vec_failing_set[cur_depth] = ancestors[u];
                vec_failing_set[cur_depth] |= ancestors[order_matrix[i]];
                if(cur_depth - 1>0)
                vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
                continue;
            }
            if(nodeId[v]==1&& (!is_used[u])){
                idx[cur_depth] += 1;
               vec_failing_set[cur_depth] = ancestors[u];
                vec_failing_set[cur_depth] |= ancestors[order_matrix[i]];
                if(cur_depth - 1>0)
                vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
                continue;
            }
            if (visited_vertices[v])
            {
                idx[cur_depth] += 1;
                vec_failing_set[cur_depth] = ancestors[u];
                vec_failing_set[cur_depth] |= ancestors[reverse_embedding[v]];
                vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
                continue;
            }



            embedding[u] = v;
            idx_embedding[u] = valid_idx;
            visited_vertices[v] = true;
            idx[cur_depth] += 1;

            reverse_embedding[v] = u;
            if (cur_depth == max_depth - 1)
            {   
                UNPM = 0;
                if(nodeId[embedding[order_matrix[cur_depth]]]==0){
                    UNPM++;
                    nodeId[embedding[order_matrix[cur_depth]]]=1;
                }
                embedding_cnt += 1;
                visited_vertices[v] = false;
                int ao = cur_depth; //-1;
                tempPos[ao]=embedding[order_matrix[cur_depth]];

                ui vqo = order_matrix[ao];

                while (ao > 0)
                {                
                    if(nodeId[embedding[vqo]]==0){
                    UNPM++;
                    nodeId[embedding[vqo]]=1;
                    //idx[ao] = idx_count[ao];
                }
                    tempPos[ao]=embedding[vqo];
                   if(ao>K_size2){
                    idx[ao] = idx_count[ao];

                    vec_failing_set[ao].set();
                    vec_failing_set[ao - 1] |= vec_failing_set[ao];
                    //vec_failing_set[ao - 1] = ancestors[order_matrix[ao]];
                    reverse_embedding.erase(embedding[ao]);
                    
                    
                   }

                    vec_failing_set[cur_depth].set();
                    vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
                    
                    Match_BA[vqo]=true;
                    ao--;
                    vqo = order_matrix[ao];
                }
                Match_BA[vqo]=true;
                if(nodeId[embedding[vqo]]==0){
                    UNPM++;
                    nodeId[embedding[vqo]]=1;
                    //idx[ao] = idx_count[ao];
                }
                tempPos[ao]=embedding[vqo];
                if (UNPM > pq.top().first)
                            {                                        
                                int addPos=0;
                                if (heapcount >= ksize)
                                {
                                   auto [value, index] = pq.top();      
                                    pq.pop();
                                    pq.push({UNPM,index});
                                    addPos=index;
                                    
                                }else{
                                     pq.push({UNPM,heapcount});
                                     addPos=heapcount;            
                                }
                            for (int aa=0;aa<qsiz;aa++){
                                SolPos[addPos][aa]=tempPos[aa];
                                    }
                                heapcount++;
                            }
                reverse_embedding.erase(embedding[u]);
                vec_failing_set[cur_depth].set();
                vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
                if (embedding_cnt >= output_limit_num)
                {
                    goto EXIT;
                }

            }
            else
            {
                call_count += 1;
                cur_depth += 1;   
                //Match_BA[order_matrix[cur_depth]]=false;          
                idx[cur_depth] = 0;
                generateValidCandidateIndex(cur_depth, idx_embedding, idx_count, valid_candidate_idx, edge_matrix, bn,
                                            bn_count, order_matrix, temp_buffer);
            
            //if(is_used[u]){
             //   while(nodeId[v]==0&&idx[cur_depth]<idx_count[cur_depth])
            //    idx[cur_depth] += 1;
            //}else{
           //     while(nodeId[v]==1&&idx[cur_depth]<idx_count[cur_depth])
            //    idx[cur_depth] += 1;
            //}
            //if(idx[cur_depth]==idx_count[cur_depth]){
            //    idx_count[cur_depth]=0;
            //}
            

                if (idx_count[cur_depth] == 0)
                {
                    vec_failing_set[cur_depth - 1] = ancestors[order_matrix[cur_depth]];
                }  
                else
                {   Match_BA[u]==false;
                    vec_failing_set[cur_depth - 1].reset();
                }
            }
        }
        cur_depth -= 1;
        if (cur_depth < 0)
            break;
        else
        {
            VertexID u = order_matrix[cur_depth];
            reverse_embedding.erase(embedding[u]);
            if (cur_depth != 0)
            {
                //if (!vec_failing_set[cur_depth].test(u)&&Match_BA[u]==false)
                if (!vec_failing_set[cur_depth].test(u))
                {
                    vec_failing_set[cur_depth - 1] = vec_failing_set[cur_depth];
                    //if(idx[cur_depth]!=idx_count[cur_depth])
                    //cout<<"in"<<idx[cur_depth]<<"-"<<idx_count[cur_depth];
                    idx[cur_depth] = idx_count[cur_depth];
                    
                }
                else
                {
                    vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
                }
            }
            visited_vertices[embedding[u]] = false;

        }
    }

        hasNext = nextCombination(indices, max_depth, K_size);
    }
        K_size++;
    }


EXIT:
    releaseBuffer(max_depth, idx, idx_count, embedding, idx_embedding, temp_buffer, valid_candidate_idx,
                  visited_vertices,
                  bn, bn_count);
        std::priority_queue<std::pair<int, int>> maxHeap;
    int topkcounter = 0;
    int greedysum = 0;
    unordered_set<ui> EmbGreedy;
    while (!pq.empty())
    {   if(pq.top().first!=0)
        maxHeap.push(pq.top());
        pq.pop();
    }
    while (!maxHeap.empty())
    {   auto [value, index]= maxHeap.top();
        for (int aa=0;aa<qsiz;aa++){
            EmbGreedy.insert(SolPos[index][aa]);
        }
        maxHeap.pop();
        topkcounter++;
        switch (topkcounter)
        {

        case VALUE_10:
            s.topk[0] = EmbGreedy.size();
            break;
        case VALUE_50:
            s.topk[1] = EmbGreedy.size();
            break;
        case VALUE_100:
            s.topk[2] = EmbGreedy.size();
            break;
        case VALUE_500:
            s.topk[3] = EmbGreedy.size();
            break;
        case VALUE_1000:
            s.topk[4] = EmbGreedy.size();
            break;
        case VALUE_2500:
            s.topk[5] = EmbGreedy.size();
            break;    
        case VALUE_5000:
            s.topk[6] = EmbGreedy.size();
            break;
        case VALUE_10000:
        s.topk[7] = EmbGreedy.size();
            break;    
        default:
            break;
        }
    }
    for (int i = 0; i < 8; i++)
    {
        if (s.topk[i] == 0)
            s.topk[i] = EmbGreedy.size();
    }

        s.embedding_cnt = embedding_cnt;
        ui countSMU=0;
        for (int i=0;i<data_graph->getVerticesCount();i++){
            countSMU+=nodeId[i];
        }
        s.Can_embed = countSMU;
    // s.topk=greedysum;
        return s;
}




enumResult
EvaluateQuery::LFTJDLSBASIC(const Graph *data_graph, const Graph *query_graph,ui *&nodeId, Edges ***edge_matrix, ui **candidates,
                    ui *candidates_count, ui *order, size_t output_limit_num, size_t &call_count, size_t &valid_vtx_cnt, int TimeL, int FairT,
                    const std::unordered_map<VertexID, std::pair<std::set<VertexID>, std::set<VertexID>>> &ordered_constraints)
    {

    auto start = std::chrono::high_resolution_clock::now();
    enumResult s;
    ui **bn;
    ui *bn_count;
    generateBN(query_graph, order, bn, bn_count);
    int UNPM = 0;
    ui *idx;
    ui *idx_count;
    ui *embedding;
    ui *idx_embedding;
    ui *temp_buffer;
    ui **valid_candidate_idx;
    int i=0;
    bool *visited_vertices;
        int count1=0;
    int count2=0;
    int K_size2=0;
        ui tempstore=0;
    ui nbr_cnt;
    allocateBuffer(data_graph, query_graph, candidates_count, idx, idx_count, embedding, idx_embedding,
                   temp_buffer, valid_candidate_idx, visited_vertices);

    size_t embedding_cnt = 0;
    int cur_depth = 0;
    ui max_depth = query_graph->getVerticesCount();
    VertexID start_vertex = order[0];
    ui order_matrix[max_depth];
    bool is_used[max_depth];
    bool is_used_copy[max_depth];
    //unordered_set<ui> EmbSum;
    bool Match_BA[max_depth] = {false};
    bool conflicts[max_depth][max_depth] = {false};
    bool conflictT[max_depth][max_depth] = {false};
        for (int i=0;i<max_depth;i++){
        for (int j=0;j<max_depth;j++)
        conflicts[i][j]=false;
    }
    for (int i=0;i<max_depth;i++){
    const VertexID *nbrs = query_graph->getVertexNeighbors(i, nbr_cnt);
            for (int j=0;j<nbr_cnt;j++){
                conflicts[i][nbrs[j]]=true;
            }
    }
    /*
    for (int i=0;i<max_depth;i++){
        for (int j=0;j<max_depth;j++)
        cout<< conflicts[i][j]<<" ";
    cout<<" "<<endl;
    }
    for (int i=0;i<max_depth;i++){
        cout<< order[i];
    }cout<<" "<<endl;*/
    idx[cur_depth] = 0;
    idx_count[cur_depth] = candidates_count[start_vertex];

    for (ui i = 0; i < idx_count[cur_depth]; ++i)
    {
        valid_candidate_idx[cur_depth][i] = i;
    }


    std::unordered_map<VertexID, VertexID> reverse_embedding;
    reverse_embedding.reserve(MAXIMUM_QUERY_GRAPH_SIZE * 2);
    bool hasNext = true;
    ui K_size=1;
    TimeL = TimeL * 1000;
    vector<int> combination;
    double ens = 0;
    while (true)
    {
        while (idx[cur_depth] < idx_count[cur_depth])
        {   
            ui valid_idx = valid_candidate_idx[cur_depth][idx[cur_depth]];
            VertexID u = order[cur_depth];
            VertexID v = candidates[u][valid_idx];
            auto end = std::chrono::high_resolution_clock::now();
            ens = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
            if (ens > TimeL)
            { // 1000 1 sec
                goto EXIT;
            }
            if (visited_vertices[v])
            {
                idx[cur_depth] += 1;
                conflictT[u][reverse_embedding[v]]=true;
                conflictT[reverse_embedding[v]][u]=true;
                continue;
            }
            if(nodeId[v]==1){
                idx[cur_depth] += 1;
                
                //Match_BA[cur_depth]=true;
                //for (int b=max_depth;b>=0;b--)
                //Match_BA[b]=true;
                //conflictT[u][reverse_embedding[v]]=true;
                continue;
            }
            Match_BA[cur_depth]=false;
            embedding[u] = v;
            idx_embedding[u] = valid_idx;
            visited_vertices[v] = true;
            idx[cur_depth] += 1;

            reverse_embedding[v] = u;
            if (cur_depth == max_depth - 1)
            {
                UNPM = 0;
                if(nodeId[embedding[order[cur_depth]]]==0){
                    UNPM++;
                    nodeId[embedding[order[cur_depth]]]=1;
                }
                embedding_cnt += 1;
                visited_vertices[v] = false;
                int ao = cur_depth; //-1;
                ui vqo = order[ao];
                nodeId[v]=1;
                while (ao > 0)
                {   
                    Match_BA[ao]=true;
                    if(nodeId[embedding[vqo]]==0){
                    UNPM++;
                    nodeId[embedding[vqo]]=1;
                }
                    idx[ao] = idx_count[ao];
                    ao--;
                    //reverse_embedding.erase(embedding[vqo]);
                    //visited_vertices[embedding[vqo]] = false;
                    vqo = order[ao];

                }
                Match_BA[ao]=true;
                nodeId[embedding[vqo]]=1;
                //reverse_embedding.erase(embedding[u]);
                //visited_vertices[embedding[vqo]] = false;
            }
            else
            {
                call_count += 1;
                cur_depth += 1;
                
                idx[cur_depth] = 0;
                generateValidCandidateIndex(cur_depth, idx_embedding, idx_count, valid_candidate_idx, edge_matrix, bn,
                                            bn_count, order, temp_buffer);
            //if (cur_depth!=0&&cur_depth!=max_depth-1)
            for (int i=0;i<max_depth;i++){
                conflictT[order[cur_depth]][i]=conflicts[order[cur_depth]][i];
            }
           // ui valid_idx1 = valid_candidate_idx[cur_depth][idx[cur_depth]];
           // while(nodeId[candidates[order[cur_depth]][valid_idx1]]==1&&idx[cur_depth]<idx_count[cur_depth]){
          ///  idx[cur_depth]++;
           // valid_idx1 = valid_candidate_idx[cur_depth][idx[cur_depth]];    }
            }
        }

        cur_depth -= 1;
        if (cur_depth < 0)
            break;
        else
        {
            VertexID u = order[cur_depth];
            reverse_embedding.erase(embedding[u]);
            if (cur_depth != 0&&Match_BA[cur_depth]==false&&cur_depth < max_depth-2)
            {
                //for (int k=cur_depth+1;k<max_depth-1;k++)
                //if(conflictT[order[cur_depth+1]][u]==false){
                if(conflictT[u][order[cur_depth+1]]==false){
                    //cout<<"in "<<u<<order[cur_depth+1]<<endl;
                    //idx[cur_depth] = idx_count[cur_depth];
                    //break;
                    ;
                }
            }
            visited_vertices[embedding[u]] = false;

        }
    }
   
    hasNext = true;
     K_size=1;
        //for (int i=0;i<max_depth;i++)
        //cout<<order[i];
        //cout<<""<<endl;

    while (K_size<max_depth){
        //cout<<K_size<<endl;
        vector<int> indices(K_size);        
     
        for (ui i = 0; i < K_size; ++i) {
        indices[i] = i;
        }
        hasNext = true;
        while (hasNext==true) {
            count1=0;
            count2=0;
            K_size2=0;
            tempstore=0;
            nbr_cnt=0;
            for (int i=0;i<max_depth;i++){
            is_used[i]=false;
            is_used_copy[i]=false;
            }
            combination.clear();
            
            for (ui i = 0; i < K_size; ++i) {
                combination.push_back(order[indices[i]]);
                is_used[order[indices[i]]] = true;
            }

            order_matrix[0] = combination[0];
            is_used_copy[combination[0]]=true;
            count1=1;
            count2=0;
            //count1++;
        while(count1<max_depth){
            const VertexID *nbrs = query_graph->getVertexNeighbors(order_matrix[count2], nbr_cnt);
            for (int i=0;i<nbr_cnt;i++){
                if(is_used_copy[nbrs[i]]==false){
                    is_used_copy[nbrs[i]]=true;
                    order_matrix[count1]=nbrs[i];
                    count1++;
                }
            }count2++;
        }

        
        count1=max_depth;
        count2=0;
        i=1;
        while (i<count1){
            if(query_graph->getVertexDegree(order_matrix[i])==1){
                tempstore=order_matrix[i];
                for (int j=i+1;j<max_depth;j++){
                    order_matrix[j-1]=order_matrix[j];
                }
                order_matrix[max_depth-1]=tempstore;
                count1--;
            }
            else{i++;}
        }   

        K_size2=1;
        i=0;
        while(is_used[order_matrix[i]]==true){
            K_size2++;
            i++;
        }
        //cout<<" #"<<K_size2<<endl;
        i=0;
        
        generateBNLM(query_graph, order_matrix, bn, bn_count);
        allocateBufferLM(data_graph, query_graph, candidates_count, idx, idx_count, embedding, idx_embedding,
                   temp_buffer, valid_candidate_idx, visited_vertices);

        std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> vec_failing_set(max_depth);
        start_vertex = order_matrix[0];
        cur_depth=0;
        idx[cur_depth] = 0;
        idx_count[cur_depth] = candidates_count[start_vertex];
        Match_BA[0] = false;
        for (ui i = 0; i < idx_count[cur_depth]; ++i)
        {
        valid_candidate_idx[cur_depth][i] = i;
        }
        while (true)
    {
        while (idx[cur_depth] < idx_count[cur_depth])
        {   
            
            ui valid_idx = valid_candidate_idx[cur_depth][idx[cur_depth]];
            VertexID u = order_matrix[cur_depth];
            VertexID v = candidates[u][valid_idx];
            auto end = std::chrono::high_resolution_clock::now();
            ens = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
            if (ens > TimeL)
            { // 1000 1 sec
                goto EXIT;
            }
            if (visited_vertices[v])
            {
                idx[cur_depth] += 1;
                conflictT[u][reverse_embedding[v]]=1;
                continue;
            }
            if(nodeId[v]==0&&is_used[u]){
                idx[cur_depth] += 1;
                conflictT[u][reverse_embedding[v]]=1;
                continue;
            }
            if(nodeId[v]==1&& (!is_used[u])){
                idx[cur_depth] += 1;
                conflictT[u][reverse_embedding[v]]=1;
                continue;
            }




            embedding[u] = v;
            idx_embedding[u] = valid_idx;
            visited_vertices[v] = true;
            idx[cur_depth] += 1;

            reverse_embedding[v] = u;
            
            if (cur_depth == max_depth - 1)
            {   
                UNPM = 0;
                if(nodeId[embedding[order_matrix[cur_depth]]]==0){
                    UNPM++;
                    nodeId[embedding[order_matrix[cur_depth]]]=1;
                }
                embedding_cnt += 1;
                visited_vertices[v] = false;
                int ao = cur_depth; //-1;
                ui vqo = order_matrix[ao];

                while (ao > 0)
                {                
                    if(nodeId[embedding[vqo]]==0){
                    UNPM++;
                    nodeId[embedding[vqo]]=1;
                    //idx[ao] = idx_count[ao];
                }
                   if(ao>=K_size2){
                    idx[ao] = idx_count[ao];            
                   }
                    
                    
                    Match_BA[ao]=true;
                    ao--;
                    vqo = order_matrix[ao];
                }
                Match_BA[ao]=true;
                nodeId[embedding[vqo]]=1;
                
                reverse_embedding.erase(embedding[u]);
                
                if (embedding_cnt >= output_limit_num)
                {
                    goto EXIT;
                }
            }
            else
            {
                call_count += 1;
                cur_depth += 1;   
                Match_BA[cur_depth]=false;          
                idx[cur_depth] = 0;
                generateValidCandidateIndex(cur_depth, idx_embedding, idx_count, valid_candidate_idx, edge_matrix, bn,
                                            bn_count, order_matrix, temp_buffer);
            

            if (cur_depth!=max_depth-1)
            for (int i=0;i<max_depth;i++){
                conflictT[order_matrix[cur_depth]][i]=conflicts[order_matrix[cur_depth]][i];
            }

            //if(is_used[u]){
            //    while(nodeId[v]==0&&idx[cur_depth]<idx_count[cur_depth])
            //    idx[cur_depth] += 1;
            //}else{
            //    while(nodeId[v]==1&&idx[cur_depth]<idx_count[cur_depth])
           //     idx[cur_depth] += 1;
           // }
            //if(idx[cur_depth]==idx_count[cur_depth]){
            //    idx_count[cur_depth]=0;
            //}
            
            
            }
        }
        cur_depth -= 1;
        if (cur_depth < 0)
            break;
        else
        {   
            VertexID u = order_matrix[cur_depth];
            reverse_embedding.erase(embedding[u]);

            if (cur_depth != 0&&Match_BA[cur_depth]==false&&cur_depth <  max_depth-2)
            {
                //for (int k=cur_depth+1;k<max_depth-1;k++)
                if(conflictT[order_matrix[cur_depth+1]][u]==false){
                    //idx[cur_depth] = idx_count[cur_depth];
                    ;
                    //break;
                }
            }
            visited_vertices[embedding[u]] = false;
           
        }
    }

        hasNext = nextCombination(indices, max_depth, K_size);
    }
        K_size++;
    }


EXIT:
    releaseBuffer(max_depth, idx, idx_count, embedding, idx_embedding, temp_buffer, valid_candidate_idx,
                  visited_vertices,
                  bn, bn_count);

        s.embedding_cnt = embedding_cnt;
        ui countSMU=0;
        for (int i=0;i<data_graph->getVerticesCount();i++){
            countSMU+=nodeId[i];
        }
        s.Can_embed = countSMU;
    // s.topk=greedysum;
        return s;
}






enumResult
EvaluateQuery::LFTJ(const Graph *data_graph, const Graph *query_graph,ui *&nodeId, Edges ***edge_matrix, ui **candidates,
                    ui *candidates_count, ui *order, size_t output_limit_num, size_t &call_count, size_t &valid_vtx_cnt, int TimeL, int FairT,
                    const std::unordered_map<VertexID, std::pair<std::set<VertexID>, std::set<VertexID>>> &ordered_constraints)
{
    
#ifdef TOPKGREEDY
    int qsiz = query_graph->getVerticesCount();
    int ksize = 1000;
    int SolPos[ksize+1][qsiz];
    std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, std::greater<>> pq;
    int tempPos[qsiz];
    pq.push({0,0});
    int heapcount = 0;
#endif
    auto start = std::chrono::high_resolution_clock::now();
    enumResult s;
    ui **bn;
    ui *bn_count;
    generateBN(query_graph, order, bn, bn_count);

    int UNPM = 0;
    ui *idx;
    ui *idx_count;
    ui *embedding;
    ui *idx_embedding;
    ui *temp_buffer;
    ui **valid_candidate_idx;
    bool *visited_vertices;
    allocateBuffer(data_graph, query_graph, candidates_count, idx, idx_count, embedding, idx_embedding,
                   temp_buffer, valid_candidate_idx, visited_vertices);

    size_t embedding_cnt = 0;
    int cur_depth = 0;
    ui max_depth = query_graph->getVerticesCount();
    VertexID start_vertex = order[0];
    //unordered_set<ui> EmbSum;
    bool Match_BA[max_depth] = {false};
    idx[cur_depth] = 0;
    idx_count[cur_depth] = candidates_count[start_vertex];

    for (ui i = 0; i < idx_count[cur_depth]; ++i)
    {
        valid_candidate_idx[cur_depth][i] = i;
    }

#ifdef ENABLE_FAILING_SET
    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> ancestors;
    computeAncestor(query_graph, bn, bn_count, order, ancestors);
    AncestorWithSymmetryBreaking(query_graph, ancestors, ordered_constraints);
    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> vec_failing_set(max_depth);
    std::unordered_map<VertexID, VertexID> reverse_embedding;
    reverse_embedding.reserve(MAXIMUM_QUERY_GRAPH_SIZE * 2);
#endif
    TimeL = TimeL * 1000;
    double ens = 0;
    while (true)
    {
        while (idx[cur_depth] < idx_count[cur_depth])
        {
            Match_BA[cur_depth] = false;
            ui valid_idx = valid_candidate_idx[cur_depth][idx[cur_depth]];

#ifdef ENABLE_DIVERSITY
            if (cur_depth == 0 && idx[cur_depth] != 0 && EmbSum.size() > 0 && FairT > 0)
            {
                if (FairT == 1)
                    rankSimple(order[0], nodeId, candidates, valid_candidate_idx, idx, valid_idx, idx_count);
                else if (FairT == 2) // already have experiments on that
                    rankLess(query_graph, edge_matrix, order[0], nodeId, candidates, valid_candidate_idx, idx, valid_idx, idx_count);
                else if (FairT == 3)
                    rankMore(query_graph, edge_matrix, order[0], nodeId, candidates, valid_candidate_idx, idx, valid_idx, idx_count);
                else if (FairT == 4)
                    rankFinal(query_graph, edge_matrix, order[0], nodeId, candidates, valid_candidate_idx, idx, valid_idx, idx_count);
                else if (FairT == 5)
                    rankFinal1(query_graph, edge_matrix, order[0], nodeId, candidates, valid_candidate_idx, idx, valid_idx, idx_count);

                valid_idx = valid_candidate_idx[cur_depth][idx[cur_depth]];
            }
            // from all the valid IDS
#endif

            VertexID u = order[cur_depth];
            VertexID v = candidates[u][valid_idx];
            auto end = std::chrono::high_resolution_clock::now();
            ens = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
            if (ens > TimeL)
            { // 1000 1 sec
                goto EXIT;
            }
            if (visited_vertices[v])
            {
                idx[cur_depth] += 1;
                vec_failing_set[cur_depth] = ancestors[u];
                vec_failing_set[cur_depth] |= ancestors[reverse_embedding[v]];
                vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
                continue;
            }

            embedding[u] = v;
            idx_embedding[u] = valid_idx;
            visited_vertices[v] = true;
            idx[cur_depth] += 1;

            reverse_embedding[v] = u;
            if (cur_depth == max_depth - 1)
            {
                UNPM = 0;
                if(nodeId[embedding[order[cur_depth]]]==0){
                    UNPM++;
                    nodeId[embedding[order[cur_depth]]]=1;
                }
                embedding_cnt += 1;
                visited_vertices[v] = false;
                int ao = cur_depth; //-1;
                ui vqo = order[ao];
                #ifdef TOPKGREEDY
                tempPos[ao]=embedding[order[cur_depth]];
                #endif
                while (ao >= 0 && Match_BA[ao] == false)
                {                
                    if(nodeId[embedding[vqo]]==0){
                    UNPM++;
                    nodeId[embedding[vqo]]=1;
                }
                #ifdef TOPKGREEDY
                tempPos[ao]=embedding[vqo];
                #endif
                    // Match_BA[ao] = true;
                    ao--;
                    vqo = order[ao];
                }
#ifdef TOPKGREEDY
                            if (UNPM > pq.top().first)
                            {                                        
                                int addPos=0;
                                if (heapcount >= ksize)
                                {
                                   auto [value, index] = pq.top();      
                                    pq.pop();
                                    pq.push({UNPM,index});
                                    addPos=index;
                                    
                                }else{
                                     pq.push({UNPM,heapcount});
                                     addPos=heapcount;            
                                }
                            for (int aa=0;aa<qsiz;aa++){
                                SolPos[addPos][aa]=tempPos[aa];
                                    }
                                heapcount++;
                            }
#endif

                reverse_embedding.erase(embedding[u]);
                vec_failing_set[cur_depth].set();
                vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
                if (embedding_cnt >= output_limit_num)
                {
                    goto EXIT;
                }
            }
            else
            {
                call_count += 1;
                cur_depth += 1;

                idx[cur_depth] = 0;
                generateValidCandidateIndex(cur_depth, idx_embedding, idx_count, valid_candidate_idx, edge_matrix, bn,
                                            bn_count, order, temp_buffer);

                pruneCandidatesIndexBySymmetryBreaking(cur_depth, embedding, order, candidates_count, candidates,
                                                       idx_count, valid_candidate_idx, ordered_constraints);

                if (idx_count[cur_depth] == 0)
                {
                    vec_failing_set[cur_depth - 1] = ancestors[order[cur_depth]];
                }
                else
                {
                    vec_failing_set[cur_depth - 1].reset();
                }
            }
        }

        cur_depth -= 1;
        if (cur_depth < 0)
            break;
        else
        {
            VertexID u = order[cur_depth];
            reverse_embedding.erase(embedding[u]);
            if (cur_depth != 0)
            {
                if (!vec_failing_set[cur_depth].test(u))
                {
                    vec_failing_set[cur_depth - 1] = vec_failing_set[cur_depth];
                    idx[cur_depth] = idx_count[cur_depth];
                }
                else
                {
                    vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
                }
            }
            visited_vertices[embedding[u]] = false;

        }
    }


EXIT:
    releaseBuffer(max_depth, idx, idx_count, embedding, idx_embedding, temp_buffer, valid_candidate_idx,
                  visited_vertices,
                  bn, bn_count);

#if ENABLE_QFLITER == 1
    delete[] temp_bsr_base1_;
    delete[] temp_bsr_base2_;
    delete[] temp_bsr_state1_;
    delete[] temp_bsr_state2_;

    for (ui i = 0; i < max_depth; ++i)
    {
        for (ui j = 0; j < max_depth; ++j)
        {
        }
        delete[] qfliter_bsr_graph_[i];
    }
    delete[] qfliter_bsr_graph_;
#endif
    /*
    #ifdef TOPKGREEDY
           int greedysum=0;
                   while(!pq.empty()){
                       greedysum+=pq.top();
                       pq.pop();
                   }
       #endif

   */
#ifdef TOPKGREEDY
    std::priority_queue<std::pair<int, int>> maxHeap;
    int topkcounter = 0;
    int greedysum = 0;
    unordered_set<ui> EmbGreedy;
    while (!pq.empty())
    {
        maxHeap.push(pq.top());
        pq.pop();
    }
    while (!maxHeap.empty())
    {   auto [value, index]= maxHeap.top();
        for (int aa=0;aa<qsiz;aa++){
            EmbGreedy.insert(SolPos[index][aa]);
        }
        maxHeap.pop();
        topkcounter++;
        switch (topkcounter)
        {
        case VALUE_5:
            s.topk[0] = EmbGreedy.size();
            break;
        case VALUE_10:
            s.topk[1] = EmbGreedy.size();
            break;
        case VALUE_50:
            s.topk[2] = EmbGreedy.size();
            break;
        case VALUE_100:
            s.topk[3] = EmbGreedy.size();
            break;
        case VALUE_250:
            s.topk[4] = EmbGreedy.size();
            break;
        case VALUE_500:
            s.topk[5] = EmbGreedy.size();
            break;
        case VALUE_750:
            s.topk[6] = EmbGreedy.size();
            break;
        case VALUE_1000:
            s.topk[7] = EmbGreedy.size();
            break;
        default:
            break;
        }
    }
    for (int i = 0; i < 8; i++)
    {
        if (s.topk[i] == 0)
            s.topk[i] = EmbGreedy.size();
    }

#endif
    s.embedding_cnt = embedding_cnt;
    ui countSMU=0;
    for (int i=0;i<data_graph->getVerticesCount();i++){
        countSMU+=nodeId[i];
    }
    s.Can_embed = countSMU;
    // s.topk=greedysum;
    return s;
}



enumResult
EvaluateQuery::LFTJKOPTSQ(const Graph *data_graph, const Graph *query_graph,ui *&nodeId, Edges ***edge_matrix, ui **candidates,
                    ui *candidates_count, ui *order, size_t output_limit_num, size_t &call_count, size_t &valid_vtx_cnt, int TimeL, int FairT,
                    const std::unordered_map<VertexID, std::pair<std::set<VertexID>, std::set<VertexID>>> &ordered_constraints)
{
    
    int qsiz = query_graph->getVerticesCount();
    int ksize = 10000;
    int SolPos[ksize+1][qsiz];
    std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, std::greater<>> pq;
    int tempPos[qsiz];
    pq.push({0,0});
    int heapcount = 0;
    auto start = std::chrono::high_resolution_clock::now();
    enumResult s;
    ui **bn;
    ui *bn_count;
    generateBN(query_graph, order, bn, bn_count);

    int UNPM = 0;
    ui *idx;
    ui *idx_count;
    ui *embedding;
    ui *idx_embedding;
    ui *temp_buffer;
    ui **valid_candidate_idx;
    bool *visited_vertices;
    allocateBuffer(data_graph, query_graph, candidates_count, idx, idx_count, embedding, idx_embedding,
                   temp_buffer, valid_candidate_idx, visited_vertices);

    size_t embedding_cnt = 0;
    int cur_depth = 0;
    ui max_depth = query_graph->getVerticesCount();
    VertexID start_vertex = order[0];
    //unordered_set<ui> EmbSum;
    bool Match_BA[max_depth] = {false};
    idx[cur_depth] = 0;
    idx_count[cur_depth] = candidates_count[start_vertex];

    for (ui i = 0; i < idx_count[cur_depth]; ++i)
    {
        valid_candidate_idx[cur_depth][i] = i;
    }

#ifdef ENABLE_FAILING_SET
    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> ancestors;
    computeAncestor(query_graph, bn, bn_count, order, ancestors);
    AncestorWithSymmetryBreaking(query_graph, ancestors, ordered_constraints);
    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> vec_failing_set(max_depth);
    std::unordered_map<VertexID, VertexID> reverse_embedding;
    reverse_embedding.reserve(MAXIMUM_QUERY_GRAPH_SIZE * 2);
#endif
    TimeL = TimeL * 1000;
    double ens = 0;
    while (true)
    {
        while (idx[cur_depth] < idx_count[cur_depth])
        {
            Match_BA[cur_depth] = false;
            ui valid_idx = valid_candidate_idx[cur_depth][idx[cur_depth]];

            VertexID u = order[cur_depth];
            VertexID v = candidates[u][valid_idx];
            auto end = std::chrono::high_resolution_clock::now();
            ens = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
            if (ens > TimeL)
            { // 1000 1 sec
                goto EXIT;
            }
            if (visited_vertices[v])
            {
                idx[cur_depth] += 1;
                vec_failing_set[cur_depth] = ancestors[u];
                vec_failing_set[cur_depth] |= ancestors[reverse_embedding[v]];
                vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
                continue;
            }

            embedding[u] = v;
            idx_embedding[u] = valid_idx;
            visited_vertices[v] = true;
            idx[cur_depth] += 1;

            reverse_embedding[v] = u;
            if (cur_depth == max_depth - 1)
            {
                UNPM = 0;
                if(nodeId[embedding[order[cur_depth]]]==0){
                    UNPM++;
                    nodeId[embedding[order[cur_depth]]]=1;
                }
                embedding_cnt += 1;
                visited_vertices[v] = false;
                int ao = cur_depth; //-1;
                ui vqo = order[ao];
                tempPos[ao]=embedding[order[cur_depth]];
                while (ao >= 0 && Match_BA[ao] == false)
                {                
                    if(nodeId[embedding[vqo]]==0){
                    UNPM++;
                    nodeId[embedding[vqo]]=1;
                }
                tempPos[ao]=embedding[vqo];
                    // Match_BA[ao] = true;
                    ao--;
                    vqo = order[ao];
                }
                            if (UNPM > pq.top().first)
                            {                                        
                                int addPos=0;
                                if (heapcount >= ksize)
                                {
                                   auto [value, index] = pq.top();      
                                    pq.pop();
                                    pq.push({UNPM,index});
                                    addPos=index;
                                    
                                }else{
                                     pq.push({UNPM,heapcount});
                                     addPos=heapcount;            
                                }
                            for (int aa=0;aa<qsiz;aa++){
                                SolPos[addPos][aa]=tempPos[aa];
                                    }
                                heapcount++;
                            }

                reverse_embedding.erase(embedding[u]);
                vec_failing_set[cur_depth].set();
                vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
                if (embedding_cnt >= output_limit_num)
                {
                    goto EXIT;
                }
            }
            else
            {
                call_count += 1;
                cur_depth += 1;

                idx[cur_depth] = 0;
                generateValidCandidateIndex(cur_depth, idx_embedding, idx_count, valid_candidate_idx, edge_matrix, bn,
                                            bn_count, order, temp_buffer);

                pruneCandidatesIndexBySymmetryBreaking(cur_depth, embedding, order, candidates_count, candidates,
                                                       idx_count, valid_candidate_idx, ordered_constraints);

                if (idx_count[cur_depth] == 0)
                {
                    vec_failing_set[cur_depth - 1] = ancestors[order[cur_depth]];
                }
                else
                {
                    vec_failing_set[cur_depth - 1].reset();
                }
            }
        }

        cur_depth -= 1;
        if (cur_depth < 0)
            break;
        else
        {
            VertexID u = order[cur_depth];
            reverse_embedding.erase(embedding[u]);
            if (cur_depth != 0)
            {
                if (!vec_failing_set[cur_depth].test(u))
                {
                    vec_failing_set[cur_depth - 1] = vec_failing_set[cur_depth];
                    idx[cur_depth] = idx_count[cur_depth];
                }
                else
                {
                    vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
                }
            }
            visited_vertices[embedding[u]] = false;

        }
    }


EXIT:
    releaseBuffer(max_depth, idx, idx_count, embedding, idx_embedding, temp_buffer, valid_candidate_idx,
                  visited_vertices,
                  bn, bn_count);

#if ENABLE_QFLITER == 1
    delete[] temp_bsr_base1_;
    delete[] temp_bsr_base2_;
    delete[] temp_bsr_state1_;
    delete[] temp_bsr_state2_;

    for (ui i = 0; i < max_depth; ++i)
    {
        for (ui j = 0; j < max_depth; ++j)
        {
        }
        delete[] qfliter_bsr_graph_[i];
    }
    delete[] qfliter_bsr_graph_;
#endif

    std::priority_queue<std::pair<int, int>> maxHeap;
    int topkcounter = 0;
    int greedysum = 0;
    unordered_set<ui> EmbGreedy;
    while (!pq.empty())
    {
        maxHeap.push(pq.top());
        pq.pop();
    }
    while (!maxHeap.empty())
    {   auto [value, index]= maxHeap.top();
        for (int aa=0;aa<qsiz;aa++){
            EmbGreedy.insert(SolPos[index][aa]);
        }
        maxHeap.pop();
        topkcounter++;
 switch (topkcounter)
        {

        case VALUE_10:
            s.topk[0] = EmbGreedy.size();
            break;
        case VALUE_50:
            s.topk[1] = EmbGreedy.size();
            break;
        case VALUE_100:
            s.topk[2] = EmbGreedy.size();
            break;
        case VALUE_500:
            s.topk[3] = EmbGreedy.size();
            break;
        case VALUE_1000:
            s.topk[4] = EmbGreedy.size();
            break;
        case VALUE_2500:
            s.topk[5] = EmbGreedy.size();
            break;    
        case VALUE_5000:
            s.topk[6] = EmbGreedy.size();
            break;
        case VALUE_10000:
        s.topk[7] = EmbGreedy.size();
            break;    
        default:
            break;
        }
    }
    for (int i = 0; i < 8; i++)
    {
        if (s.topk[i] == 0)
            s.topk[i] = EmbGreedy.size();
    }

    s.embedding_cnt = embedding_cnt;
    ui countSMU=0;
    for (int i=0;i<data_graph->getVerticesCount();i++){
        countSMU+=nodeId[i];
    }
    s.Can_embed = countSMU;
    // s.topk=greedysum;
    return s;
}
enumResult
EvaluateQuery::LFTJKOPT(const Graph *data_graph, const Graph *query_graph,ui *&nodeId, Edges ***edge_matrix, ui **candidates,
                    ui *candidates_count, ui *order, size_t output_limit_num, size_t &call_count, size_t &valid_vtx_cnt, int TimeL, int FairT,
                    const std::unordered_map<VertexID, std::pair<std::set<VertexID>, std::set<VertexID>>> &ordered_constraints)
{
    
    int qsiz = query_graph->getVerticesCount();
    int ksize = 10000;
    int SolPos[ksize+1][qsiz];
    std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, std::greater<>> pq;
    int tempPos[qsiz];
    pq.push({0,0});
    int heapcount = 0;
    auto start = std::chrono::high_resolution_clock::now();
    enumResult s;
    ui **bn;
    ui *bn_count;
    generateBN(query_graph, order, bn, bn_count);

    int UNPM = 0;
    ui *idx;
    ui *idx_count;
    ui *embedding;
    ui *idx_embedding;
    ui *temp_buffer;
    ui **valid_candidate_idx;
    bool *visited_vertices;
    allocateBuffer(data_graph, query_graph, candidates_count, idx, idx_count, embedding, idx_embedding,
                   temp_buffer, valid_candidate_idx, visited_vertices);

    size_t embedding_cnt = 0;
    int cur_depth = 0;
    ui max_depth = query_graph->getVerticesCount();
    VertexID start_vertex = order[0];
    //unordered_set<ui> EmbSum;
    bool Match_BA[max_depth] = {false};
    idx[cur_depth] = 0;
    idx_count[cur_depth] = candidates_count[start_vertex];

    for (ui i = 0; i < idx_count[cur_depth]; ++i)
    {
        valid_candidate_idx[cur_depth][i] = i;
    }

#ifdef ENABLE_FAILING_SET
    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> ancestors;
    computeAncestor(query_graph, bn, bn_count, order, ancestors);
    AncestorWithSymmetryBreaking(query_graph, ancestors, ordered_constraints);
    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> vec_failing_set(max_depth);
    std::unordered_map<VertexID, VertexID> reverse_embedding;
    reverse_embedding.reserve(MAXIMUM_QUERY_GRAPH_SIZE * 2);
#endif
    TimeL = TimeL * 1000;
    double ens = 0;
    while (true)
    {
        while (idx[cur_depth] < idx_count[cur_depth])
        {
            Match_BA[cur_depth] = false;
            ui valid_idx = valid_candidate_idx[cur_depth][idx[cur_depth]];

            VertexID u = order[cur_depth];
            VertexID v = candidates[u][valid_idx];
            auto end = std::chrono::high_resolution_clock::now();
            ens = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
            if (ens > TimeL)
            { // 1000 1 sec
                goto EXIT;
            }
            if (visited_vertices[v])
            {
                idx[cur_depth] += 1;
                vec_failing_set[cur_depth] = ancestors[u];
                vec_failing_set[cur_depth] |= ancestors[reverse_embedding[v]];
                vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
                continue;
            }

            embedding[u] = v;
            idx_embedding[u] = valid_idx;
            visited_vertices[v] = true;
            idx[cur_depth] += 1;

            reverse_embedding[v] = u;
            if (cur_depth == max_depth - 1)
            {
                UNPM = 0;
                if(nodeId[embedding[order[cur_depth]]]==0){
                    UNPM++;
                    nodeId[embedding[order[cur_depth]]]=1;
                }
                embedding_cnt += 1;
                visited_vertices[v] = false;
                int ao = cur_depth; //-1;
                ui vqo = order[ao];
                tempPos[ao]=embedding[order[cur_depth]];
                while (ao >= 0 && Match_BA[ao] == false)
                {                
                    if(nodeId[embedding[vqo]]==0){
                    UNPM++;
                    nodeId[embedding[vqo]]=1;
                }
                tempPos[ao]=embedding[vqo];
                    // Match_BA[ao] = true;
                    ao--;
                    vqo = order[ao];
                }
                            if (UNPM > pq.top().first)
                            {                                        
                                int addPos=0;
                                if (heapcount >= ksize)
                                {
                                   auto [value, index] = pq.top();      
                                    pq.pop();
                                    pq.push({UNPM,index});
                                    addPos=index;
                                    
                                }else{
                                     pq.push({UNPM,heapcount});
                                     addPos=heapcount;            
                                }
                            for (int aa=0;aa<qsiz;aa++){
                                SolPos[addPos][aa]=tempPos[aa];
                                    }
                                heapcount++;
                            }

                reverse_embedding.erase(embedding[u]);
                vec_failing_set[cur_depth].set();
                vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
                if (embedding_cnt >= output_limit_num)
                {
                    goto EXIT;
                }
            }
            else
            {
                call_count += 1;
                cur_depth += 1;

                idx[cur_depth] = 0;
                generateValidCandidateIndex(cur_depth, idx_embedding, idx_count, valid_candidate_idx, edge_matrix, bn,
                                            bn_count, order, temp_buffer);

                pruneCandidatesIndexBySymmetryBreaking(cur_depth, embedding, order, candidates_count, candidates,
                                                       idx_count, valid_candidate_idx, ordered_constraints);

                if (idx_count[cur_depth] == 0)
                {
                    vec_failing_set[cur_depth - 1] = ancestors[order[cur_depth]];
                }
                else
                {
                    vec_failing_set[cur_depth - 1].reset();
                }
            }
        }

        cur_depth -= 1;
        if (cur_depth < 0)
            break;
        else
        {
            VertexID u = order[cur_depth];
            reverse_embedding.erase(embedding[u]);
            if (cur_depth != 0)
            {
                if (!vec_failing_set[cur_depth].test(u))
                {
                    vec_failing_set[cur_depth - 1] = vec_failing_set[cur_depth];
                    idx[cur_depth] = idx_count[cur_depth];
                }
                else
                {
                    vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
                }
            }
            visited_vertices[embedding[u]] = false;

        }
    }


EXIT:
    releaseBuffer(max_depth, idx, idx_count, embedding, idx_embedding, temp_buffer, valid_candidate_idx,
                  visited_vertices,
                  bn, bn_count);

#if ENABLE_QFLITER == 1
    delete[] temp_bsr_base1_;
    delete[] temp_bsr_base2_;
    delete[] temp_bsr_state1_;
    delete[] temp_bsr_state2_;

    for (ui i = 0; i < max_depth; ++i)
    {
        for (ui j = 0; j < max_depth; ++j)
        {
        }
        delete[] qfliter_bsr_graph_[i];
    }
    delete[] qfliter_bsr_graph_;
#endif

    std::priority_queue<std::pair<int, int>> maxHeap;
    int topkcounter = 0;
    int greedysum = 0;
    unordered_set<ui> EmbGreedy;
    while (!pq.empty())
    {   if(pq.top().first!=0)
        maxHeap.push(pq.top());
        pq.pop();
    }
    while (!maxHeap.empty())
    {   auto [value, index]= maxHeap.top();
        for (int aa=0;aa<qsiz;aa++){
            EmbGreedy.insert(SolPos[index][aa]);
        }
        maxHeap.pop();
        topkcounter++;
 switch (topkcounter)
        {

        case VALUE_10:
            s.topk[0] = EmbGreedy.size();
            break;
        case VALUE_50:
            s.topk[1] = EmbGreedy.size();
            break;
        case VALUE_100:
            s.topk[2] = EmbGreedy.size();
            break;
        case VALUE_500:
            s.topk[3] = EmbGreedy.size();
            break;
        case VALUE_1000:
            s.topk[4] = EmbGreedy.size();
            break;
        case VALUE_2500:
            s.topk[5] = EmbGreedy.size();
            break;    
        case VALUE_5000:
            s.topk[6] = EmbGreedy.size();
            break;
        case VALUE_10000:
        s.topk[7] = EmbGreedy.size();
            break;    
        default:
            break;
        }
    }
    for (int i = 0; i < 8; i++)
    {
        if (s.topk[i] == 0)
            s.topk[i] = EmbGreedy.size();
    }

    s.embedding_cnt = embedding_cnt;
    ui countSMU=0;
    for (int i=0;i<data_graph->getVerticesCount();i++){
        countSMU+=nodeId[i];
    }
    s.Can_embed = countSMU;
    // s.topk=greedysum;
    return s;
}


void EvaluateQuery::generateValidCandidateIndex(ui depth, ui *idx_embedding, ui *idx_count, ui **valid_candidate_index,
                                                Edges ***edge_matrix, ui **bn, ui *bn_cnt, ui *order,
                                                ui *&temp_buffer)
{   

    VertexID u = order[depth];
    VertexID previous_bn = bn[depth][0];
    ui previous_index_id = idx_embedding[previous_bn];
    ui valid_candidates_count = 0;
    
#if ENABLE_QFLITER == 1
    BSRGraph &bsr_graph = *qfliter_bsr_graph_[previous_bn][u];
    BSRSet &bsr_set = bsr_graph.bsrs[previous_index_id];

    if (bsr_set.size_ != 0)
    {
        offline_bsr_trans_uint(bsr_set.base_, bsr_set.states_, bsr_set.size_,
                               (int *)valid_candidate_index[depth]);

        valid_candidates_count = bsr_set.size_;
    }

    if (bn_cnt[depth] > 0)
    {
        if (temp_bsr_base1_ == nullptr)
        {
            temp_bsr_base1_ = new int[1024 * 1024];
        }
        if (temp_bsr_state1_ == nullptr)
        {
            temp_bsr_state1_ = new int[1024 * 1024];
        }
        if (temp_bsr_base2_ == nullptr)
        {
            temp_bsr_base2_ = new int[1024 * 1024];
        }
        if (temp_bsr_state2_ == nullptr)
        {
            temp_bsr_state2_ = new int[1024 * 1024];
        }
        int *res_base_ = temp_bsr_base1_;
        int *res_state_ = temp_bsr_state1_;
        int *input_base_ = temp_bsr_base2_;
        int *input_state_ = temp_bsr_state2_;

        memcpy(input_base_, bsr_set.base_, sizeof(int) * bsr_set.size_);
        memcpy(input_state_, bsr_set.states_, sizeof(int) * bsr_set.size_);

        for (ui i = 1; i < bn_cnt[depth]; ++i)
        {
            VertexID current_bn = bn[depth][i];
            ui current_index_id = idx_embedding[current_bn];
            BSRGraph &cur_bsr_graph = *qfliter_bsr_graph_[current_bn][u];
            BSRSet &cur_bsr_set = cur_bsr_graph.bsrs[current_index_id];

            if (valid_candidates_count == 0 || cur_bsr_set.size_ == 0)
            {
                valid_candidates_count = 0;
                break;
            }

            valid_candidates_count = intersect_qfilter_bsr_b4_v2(cur_bsr_set.base_, cur_bsr_set.states_,
                                                                 cur_bsr_set.size_,
                                                                 input_base_, input_state_, valid_candidates_count,
                                                                 res_base_, res_state_);

            swap(res_base_, input_base_);
            swap(res_state_, input_state_);
        }

        if (valid_candidates_count != 0)
        {
            valid_candidates_count = offline_bsr_trans_uint(input_base_, input_state_, valid_candidates_count,
                                                            (int *)valid_candidate_index[depth]);
        }
    }

    idx_count[depth] = valid_candidates_count;

#ifdef YCHE_DEBUG
    Edges &previous_edge = *edge_matrix[previous_bn][u];

    auto gt_valid_candidates_count =
        previous_edge.offset_[previous_index_id + 1] - previous_edge.offset_[previous_index_id];
    ui *previous_candidates = previous_edge.edge_ + previous_edge.offset_[previous_index_id];
    ui *gt_valid_candidate_index = new ui[1024 * 1024];
    memcpy(gt_valid_candidate_index, previous_candidates, gt_valid_candidates_count * sizeof(ui));
    ui temp_count;
    for (ui i = 1; i < bn_cnt[depth]; ++i)
    {
        VertexID current_bn = bn[depth][i];
        Edges &current_edge = *edge_matrix[current_bn][u];
        ui current_index_id = idx_embedding[current_bn];

        ui current_candidates_count =
            current_edge.offset_[current_index_id + 1] - current_edge.offset_[current_index_id];
        ui *current_candidates = current_edge.edge_ + current_edge.offset_[current_index_id];

        ComputeSetIntersection::ComputeCandidates(current_candidates, current_candidates_count,
                                                  gt_valid_candidate_index, gt_valid_candidates_count, temp_buffer,
                                                  temp_count);

        std::swap(temp_buffer, gt_valid_candidate_index);
        gt_valid_candidates_count = temp_count;
    }
    assert(valid_candidates_count == gt_valid_candidates_count);
    
    cout << "Ret, Level:" << bn_cnt[depth] << ", BSR:"
         << pretty_print_array(valid_candidate_index[depth], valid_candidates_count)
         << "; GT: " << pretty_print_array(gt_valid_candidate_index, gt_valid_candidates_count) << "\n";

    for (auto i = 0; i < valid_candidates_count; i++)
    {
        assert(gt_valid_candidate_index[i] == valid_candidate_index[depth][i]);
    }
    delete[] gt_valid_candidate_index;
#endif
#else
    Edges &previous_edge = *edge_matrix[previous_bn][u];

    valid_candidates_count = previous_edge.offset_[previous_index_id + 1] - previous_edge.offset_[previous_index_id];
    ui *previous_candidates = previous_edge.edge_ + previous_edge.offset_[previous_index_id];

    memcpy(valid_candidate_index[depth], previous_candidates, valid_candidates_count * sizeof(ui));

    ui temp_count;
    for (ui i = 1; i < bn_cnt[depth]; ++i)
    {
        VertexID current_bn = bn[depth][i];
        Edges &current_edge = *edge_matrix[current_bn][u];
        ui current_index_id = idx_embedding[current_bn];

        ui current_candidates_count = current_edge.offset_[current_index_id + 1] - current_edge.offset_[current_index_id];
        ui *current_candidates = current_edge.edge_ + current_edge.offset_[current_index_id];
        ComputeSetIntersection::ComputeCandidates(current_candidates, current_candidates_count, valid_candidate_index[depth], valid_candidates_count,
                                                  temp_buffer, temp_count);
        std::swap(temp_buffer, valid_candidate_index[depth]);
        valid_candidates_count = temp_count;
    }
    idx_count[depth] = valid_candidates_count;
#endif
}

size_t EvaluateQuery::exploreGraphQLStyle(const Graph *data_graph, const Graph *query_graph, ui **candidates,
                                          ui *candidates_count, ui *order,
                                          size_t output_limit_num, size_t &call_count,
                                          const std::unordered_map<VertexID, std::pair<std::set<VertexID>, std::set<VertexID>>> &ordered_constraints)
{
    size_t embedding_cnt = 0;
    int cur_depth = 0;
    ui max_depth = query_graph->getVerticesCount();
    VertexID start_vertex = order[0];

    ui **bn;
    ui *bn_count;

    bn = new ui *[max_depth];
    for (ui i = 0; i < max_depth; ++i)
    {
        bn[i] = new ui[max_depth];
    }

    bn_count = new ui[max_depth];
    std::fill(bn_count, bn_count + max_depth, 0);

    std::vector<bool> visited_query_vertices(max_depth, false);
    visited_query_vertices[start_vertex] = true;
    for (ui i = 1; i < max_depth; ++i)
    {
        VertexID cur_vertex = order[i];
        ui nbr_cnt;
        const VertexID *nbrs = query_graph->getVertexNeighbors(cur_vertex, nbr_cnt);

        for (ui j = 0; j < nbr_cnt; ++j)
        {
            VertexID nbr = nbrs[j];

            if (visited_query_vertices[nbr])
            {
                bn[i][bn_count[i]++] = nbr;
            }
        }

        visited_query_vertices[cur_vertex] = true;
    }

    ui *idx;
    ui *idx_count;
    ui *embedding;
    VertexID **valid_candidate;
    bool *visited_vertices;

    idx = new ui[max_depth];
    idx_count = new ui[max_depth];
    embedding = new ui[max_depth];
    visited_vertices = new bool[data_graph->getVerticesCount()];
    std::fill(visited_vertices, visited_vertices + data_graph->getVerticesCount(), false);
    valid_candidate = new ui *[max_depth];

    for (ui i = 0; i < max_depth; ++i)
    {
        VertexID cur_vertex = order[i];
        ui max_candidate_count = candidates_count[cur_vertex];
        valid_candidate[i] = new VertexID[max_candidate_count];
    }

    idx[cur_depth] = 0;
    idx_count[cur_depth] = candidates_count[start_vertex];
    std::copy(candidates[start_vertex], candidates[start_vertex] + candidates_count[start_vertex],
              valid_candidate[cur_depth]);

    while (true)
    {
        while (idx[cur_depth] < idx_count[cur_depth])
        {
            VertexID u = order[cur_depth];
            VertexID v = valid_candidate[cur_depth][idx[cur_depth]];
            embedding[u] = v;
            visited_vertices[v] = true;
            idx[cur_depth] += 1;

            if (cur_depth == max_depth - 1)
            {
                embedding_cnt += 1;
                visited_vertices[v] = false;
                if (embedding_cnt >= output_limit_num)
                {
                    goto EXIT;
                }
            }
            else
            {
                call_count += 1;
                cur_depth += 1;
                idx[cur_depth] = 0;
                generateValidCandidates(data_graph, cur_depth, embedding, idx_count, valid_candidate,
                                        visited_vertices, bn, bn_count, order, candidates, candidates_count);
                pruneCandidatesBySymmetryBreaking(cur_depth, embedding, order,
                                                  idx_count, valid_candidate, ordered_constraints);
            }
        }

        cur_depth -= 1;
        if (cur_depth < 0)
            break;
        else
            visited_vertices[embedding[order[cur_depth]]] = false;
    }

EXIT:
    delete[] bn_count;
    delete[] idx;
    delete[] idx_count;
    delete[] embedding;
    delete[] visited_vertices;
    for (ui i = 0; i < max_depth; ++i)
    {
        delete[] bn[i];
        delete[] valid_candidate[i];
    }

    delete[] bn;
    delete[] valid_candidate;

    return embedding_cnt;
}

void EvaluateQuery::generateValidCandidates(const Graph *data_graph, ui depth, ui *embedding, ui *idx_count,
                                            ui **valid_candidate, bool *visited_vertices, ui **bn, ui *bn_cnt,
                                            ui *order, ui **candidates, ui *candidates_count)
{
    VertexID u = order[depth];

    idx_count[depth] = 0;

    for (ui i = 0; i < candidates_count[u]; ++i)
    {
        VertexID v = candidates[u][i];

        if (!visited_vertices[v])
        {
            bool valid = true;

            for (ui j = 0; j < bn_cnt[depth]; ++j)
            {
                VertexID u_nbr = bn[depth][j];
                VertexID u_nbr_v = embedding[u_nbr];

                if (!data_graph->checkEdgeExistence(v, u_nbr_v))
                {
                    valid = false;
                    break;
                }
            }

            if (valid)
            {
                valid_candidate[depth][idx_count[depth]++] = v;
            }
        }
    }
}

size_t EvaluateQuery::exploreQuickSIStyle(const Graph *data_graph, const Graph *query_graph, ui **candidates,
                                          ui *candidates_count, ui *order,
                                          ui *pivot, size_t output_limit_num, size_t &call_count,
                                          const std::unordered_map<VertexID, std::pair<std::set<VertexID>, std::set<VertexID>>> &ordered_constraints)
{
    size_t embedding_cnt = 0;
    int cur_depth = 0;
    ui max_depth = query_graph->getVerticesCount();
    VertexID start_vertex = order[0];

    ui **bn;
    ui *bn_count;
    generateBN(query_graph, order, pivot, bn, bn_count);

    ui *idx;
    ui *idx_count;
    ui *embedding;
    VertexID **valid_candidate;
    bool *visited_vertices;

    idx = new ui[max_depth];
    idx_count = new ui[max_depth];
    embedding = new ui[max_depth];
    visited_vertices = new bool[data_graph->getVerticesCount()];
    std::fill(visited_vertices, visited_vertices + data_graph->getVerticesCount(), false);
    valid_candidate = new ui *[max_depth];

    ui max_candidate_count = data_graph->getGraphMaxLabelFrequency();
    for (ui i = 0; i < max_depth; ++i)
    {
        valid_candidate[i] = new VertexID[max_candidate_count];
    }

    idx[cur_depth] = 0;
    idx_count[cur_depth] = candidates_count[start_vertex];
    std::copy(candidates[start_vertex], candidates[start_vertex] + candidates_count[start_vertex],
              valid_candidate[cur_depth]);

    while (true)
    {
        while (idx[cur_depth] < idx_count[cur_depth])
        {
            VertexID u = order[cur_depth];
            VertexID v = valid_candidate[cur_depth][idx[cur_depth]];
            embedding[u] = v;
            visited_vertices[v] = true;
            idx[cur_depth] += 1;

            if (cur_depth == max_depth - 1)
            {
                embedding_cnt += 1;
                visited_vertices[v] = false;
                if (embedding_cnt >= output_limit_num)
                {
                    goto EXIT;
                }
            }
            else
            {
                call_count += 1;
                cur_depth += 1;
                idx[cur_depth] = 0;
                generateValidCandidates(query_graph, data_graph, cur_depth, embedding, idx_count, valid_candidate,
                                        visited_vertices, bn, bn_count, order, pivot);
                pruneCandidatesBySymmetryBreaking(cur_depth, embedding, order,
                                                  idx_count, valid_candidate, ordered_constraints);
            }
        }

        cur_depth -= 1;
        if (cur_depth < 0)
            break;
        else
            visited_vertices[embedding[order[cur_depth]]] = false;
    }

EXIT:
    delete[] bn_count;
    delete[] idx;
    delete[] idx_count;
    delete[] embedding;
    delete[] visited_vertices;
    for (ui i = 0; i < max_depth; ++i)
    {
        delete[] bn[i];
        delete[] valid_candidate[i];
    }

    delete[] bn;
    delete[] valid_candidate;

    return embedding_cnt;
}

void EvaluateQuery::generateValidCandidates(const Graph *query_graph, const Graph *data_graph, ui depth, ui *embedding,
                                            ui *idx_count, ui **valid_candidate, bool *visited_vertices, ui **bn,
                                            ui *bn_cnt,
                                            ui *order, ui *pivot)
{
    VertexID u = order[depth];
    LabelID u_label = query_graph->getVertexLabel(u);
    ui u_degree = query_graph->getVertexDegree(u);

    idx_count[depth] = 0;

    VertexID p = embedding[pivot[depth]];
    ui nbr_cnt;
    const VertexID *nbrs = data_graph->getVertexNeighbors(p, nbr_cnt);

    for (ui i = 0; i < nbr_cnt; ++i)
    {
        VertexID v = nbrs[i];

        if (!visited_vertices[v] && u_label == data_graph->getVertexLabel(v) &&
            u_degree <= data_graph->getVertexDegree(v))
        {
            bool valid = true;

            for (ui j = 0; j < bn_cnt[depth]; ++j)
            {
                VertexID u_nbr = bn[depth][j];
                VertexID u_nbr_v = embedding[u_nbr];

                if (!data_graph->checkEdgeExistence(v, u_nbr_v))
                {
                    valid = false;
                    break;
                }
            }

            if (valid)
            {
                valid_candidate[depth][idx_count[depth]++] = v;
            }
        }
    }
}

size_t EvaluateQuery::exploreDPisoStyle(const Graph *data_graph, const Graph *query_graph, TreeNode *tree,
                                        Edges ***edge_matrix, ui **candidates, ui *candidates_count,
                                        ui **weight_array, ui *order, size_t output_limit_num,
                                        size_t &call_count,
                                        const std::unordered_map<VertexID, std::pair<std::set<VertexID>, std::set<VertexID>>> &full_constraints)
{
    ui max_depth = query_graph->getVerticesCount();

    ui *extendable = new ui[max_depth];
    for (ui i = 0; i < max_depth; ++i)
    {
        extendable[i] = tree[i].bn_count_;
    }

    ui **bn;
    ui *bn_count;
    generateBN(query_graph, order, bn, bn_count);

    ui *idx;
    ui *idx_count;
    ui *embedding;
    ui *idx_embedding;
    ui *temp_buffer;
    ui **valid_candidate_idx;
    bool *visited_vertices;
    allocateBuffer(data_graph, query_graph, candidates_count, idx, idx_count, embedding, idx_embedding,
                   temp_buffer, valid_candidate_idx, visited_vertices);

    bool *visited_u = new bool[max_depth];
    bool *violated_symmetry_u = new bool[max_depth];
    std::fill(visited_u, visited_u + max_depth, false);

    size_t embedding_cnt = 0;
    int cur_depth = 0;

#ifdef ENABLE_FAILING_SET
    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> ancestors;
    computeAncestor(query_graph, tree, order, ancestors);
    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> vec_failing_set(max_depth);
    std::unordered_map<VertexID, VertexID> reverse_embedding;
    reverse_embedding.reserve(MAXIMUM_QUERY_GRAPH_SIZE * 2);
#endif

    VertexID start_vertex = order[0];
    std::vector<dpiso_min_pq> vec_rank_queue;

    for (ui i = 0; i < candidates_count[start_vertex]; ++i)
    {
        VertexID v = candidates[start_vertex][i];
        embedding[start_vertex] = v;
        idx_embedding[start_vertex] = i;
        visited_vertices[v] = true;
        visited_u[start_vertex] = true;

#ifdef ENABLE_FAILING_SET
        reverse_embedding[v] = start_vertex;
#endif
        vec_rank_queue.emplace_back(dpiso_min_pq(extendable_vertex_compare));
        updateExtendableVertex(idx_embedding, idx_count, valid_candidate_idx, edge_matrix, temp_buffer, weight_array,
                               tree, start_vertex, extendable,
                               vec_rank_queue, query_graph);

        VertexID u = vec_rank_queue.back().top().first.first;
        vec_rank_queue.back().pop();

#ifdef ENABLE_FAILING_SET
        if (idx_count[u] == 0)
        {
            vec_failing_set[cur_depth] = ancestors[u];
        }
        else
        {
            vec_failing_set[cur_depth].reset();
        }
#endif

        call_count += 1;
        cur_depth += 1;
        order[cur_depth] = u;
        idx[u] = 0;
        while (cur_depth > 0)
        {
            while (idx[u] < idx_count[u])
            {
                ui valid_idx = valid_candidate_idx[u][idx[u]];
                v = candidates[u][valid_idx];

                bool is_voilate_symmetry = false;

                if (!full_constraints.empty())
                {
                    is_voilate_symmetry = checkSingeVertexBySymmetryBreaking(cur_depth, embedding, order, visited_u, v,
                                                                             full_constraints,
                                                                             violated_symmetry_u,
                                                                             max_depth);
                }
                if (visited_vertices[v] || is_voilate_symmetry)
                {
                    idx[u] += 1;
#ifdef ENABLE_FAILING_SET
                    vec_failing_set[cur_depth] = ancestors[u];
                    if (visited_vertices[v])
                    {
                        vec_failing_set[cur_depth] |= ancestors[reverse_embedding[v]];
                    }
                    if (is_voilate_symmetry)
                    {
                        for (ui i = 0; i < cur_depth; ++i)
                        {
                            if (violated_symmetry_u[order[i]])
                            {
                                vec_failing_set[cur_depth] |= ancestors[order[i]];
                            }
                        }
                    }
                    vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
#endif
                    continue;
                }

                embedding[u] = v;
                idx_embedding[u] = valid_idx;
                visited_vertices[v] = true;
                visited_u[u] = true;
                idx[u] += 1;

#ifdef ENABLE_FAILING_SET
                reverse_embedding[v] = u;
#endif

                if (cur_depth == max_depth - 1)
                {
                    embedding_cnt += 1;
                    visited_vertices[v] = false;
                    visited_u[u] = false;
#ifdef ENABLE_FAILING_SET
                    reverse_embedding.erase(embedding[u]);
                    vec_failing_set[cur_depth].set();
                    vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];

#endif
                    if (embedding_cnt >= output_limit_num)
                    {
                        goto EXIT;
                    }
                }
                else
                {
                    call_count += 1;
                    cur_depth += 1;
                    vec_rank_queue.emplace_back(vec_rank_queue.back());
                    updateExtendableVertex(idx_embedding, idx_count, valid_candidate_idx, edge_matrix, temp_buffer,
                                           weight_array, tree, u, extendable,
                                           vec_rank_queue, query_graph);

                    u = vec_rank_queue.back().top().first.first;
                    vec_rank_queue.back().pop();
                    idx[u] = 0;
                    order[cur_depth] = u;

#ifdef ENABLE_FAILING_SET
                    if (idx_count[u] == 0)
                    {

                        vec_failing_set[cur_depth - 1] = ancestors[u];
                    }
                    else
                    {
                        vec_failing_set[cur_depth - 1].reset();
                    }
#endif
                }
            }

            cur_depth -= 1;
            vec_rank_queue.pop_back();
            u = order[cur_depth];
            visited_vertices[embedding[u]] = false;
            visited_u[u] = false;
            restoreExtendableVertex(tree, u, extendable);
#ifdef ENABLE_FAILING_SET
            reverse_embedding.erase(embedding[u]);
            if (cur_depth != 0)
            {
                if (!vec_failing_set[cur_depth].test(u))
                {
                    vec_failing_set[cur_depth - 1] = vec_failing_set[cur_depth];

                    idx[u] = idx_count[u];
                }
                else
                {
                    vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
                }
            }
#endif
        }
    }

EXIT:
    releaseBuffer(max_depth, idx, idx_count, embedding, idx_embedding, temp_buffer, valid_candidate_idx,
                  visited_vertices,
                  bn, bn_count);
    delete[] visited_u;

    return embedding_cnt;
}

void EvaluateQuery::updateExtendableVertex(ui *idx_embedding, ui *idx_count, ui **valid_candidate_index,
                                           Edges ***edge_matrix, ui *&temp_buffer, ui **weight_array,
                                           TreeNode *tree, VertexID mapped_vertex, ui *extendable,
                                           std::vector<dpiso_min_pq> &vec_rank_queue, const Graph *query_graph)
{
    TreeNode &node = tree[mapped_vertex];
    for (ui i = 0; i < node.fn_count_; ++i)
    {
        VertexID u = node.fn_[i];
        extendable[u] -= 1;
        if (extendable[u] == 0)
        {
            generateValidCandidateIndex(u, idx_embedding, idx_count, valid_candidate_index[u], edge_matrix, tree[u].bn_,
                                        tree[u].bn_count_, temp_buffer);

            ui weight = 0;
            for (ui j = 0; j < idx_count[u]; ++j)
            {
                ui idx = valid_candidate_index[u][j];
                weight += weight_array[u][idx];
            }
            vec_rank_queue.back().emplace(std::make_pair(std::make_pair(u, query_graph->getVertexDegree(u)), weight));
        }
    }
}

void EvaluateQuery::restoreExtendableVertex(TreeNode *tree, VertexID unmapped_vertex, ui *extendable)
{
    TreeNode &node = tree[unmapped_vertex];
    for (ui i = 0; i < node.fn_count_; ++i)
    {
        VertexID u = node.fn_[i];
        extendable[u] += 1;
    }
}

void EvaluateQuery::generateValidCandidateIndex(VertexID u, ui *idx_embedding, ui *idx_count, ui *&valid_candidate_index,
                                                Edges ***edge_matrix, ui *bn, ui bn_cnt, ui *&temp_buffer)
{  
    VertexID previous_bn = bn[0];
    Edges &previous_edge = *edge_matrix[previous_bn][u];
    ui previous_index_id = idx_embedding[previous_bn];
    
    ui previous_candidates_count =
        previous_edge.offset_[previous_index_id + 1] - previous_edge.offset_[previous_index_id];
    ui *previous_candidates = previous_edge.edge_ + previous_edge.offset_[previous_index_id];

    ui valid_candidates_count = 0;
    for (ui i = 0; i < previous_candidates_count; ++i)
    {
        valid_candidate_index[valid_candidates_count++] = previous_candidates[i];
    }
    ui temp_count;
    for (ui i = 1; i < bn_cnt; ++i)
    {
        VertexID current_bn = bn[i];
        Edges &current_edge = *edge_matrix[current_bn][u];
        ui current_index_id = idx_embedding[current_bn];

        ui current_candidates_count =
            current_edge.offset_[current_index_id + 1] - current_edge.offset_[current_index_id];
        ui *current_candidates = current_edge.edge_ + current_edge.offset_[current_index_id];

        ComputeSetIntersection::ComputeCandidates(current_candidates, current_candidates_count, valid_candidate_index,
                                                  valid_candidates_count,
                                                  temp_buffer, temp_count);

        std::swap(temp_buffer, valid_candidate_index);
        valid_candidates_count = temp_count;
    }
    idx_count[u] = valid_candidates_count;
}

void EvaluateQuery::computeAncestor(const Graph *query_graph, TreeNode *tree, VertexID *order,
                                    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> &ancestors)
{
    ui query_vertices_num = query_graph->getVerticesCount();
    ancestors.resize(query_vertices_num);

    for (ui i = 0; i < query_vertices_num; ++i)
    {
        VertexID u = order[i];
        TreeNode &u_node = tree[u];
        ancestors[u].set(u);
        for (ui j = 0; j < u_node.bn_count_; ++j)
        {
            VertexID u_bn = u_node.bn_[j];
            ancestors[u] |= ancestors[u_bn];
        }
    }
}



size_t EvaluateQuery::exploreDPisoRecursiveStyle(const Graph *data_graph, const Graph *query_graph, TreeNode *tree,
                                                 Edges ***edge_matrix, ui **candidates, ui *candidates_count,
                                                 ui **weight_array, ui *order, size_t output_limit_num,
                                                 size_t &call_count,
                                                 const std::unordered_map<VertexID, std::pair<std::set<VertexID>, std::set<VertexID>>> &order_constraints)
{
    ui max_depth = query_graph->getVerticesCount();

    ui *extendable = new ui[max_depth];
    for (ui i = 0; i < max_depth; ++i)
    {
        extendable[i] = tree[i].bn_count_;
    }

    ui **bn;
    ui *bn_count;
    generateBN(query_graph, order, bn, bn_count);

    ui *idx;
    ui *idx_count;
    ui *embedding;
    ui *idx_embedding;
    ui *temp_buffer;
    ui **valid_candidate_idx;
    bool *visited_vertices;
    allocateBuffer(data_graph, query_graph, candidates_count, idx, idx_count, embedding, idx_embedding,
                   temp_buffer, valid_candidate_idx, visited_vertices);

    size_t embedding_cnt = 0;

    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> ancestors;
    computeAncestor(query_graph, tree, order, ancestors);

    std::unordered_map<VertexID, VertexID> reverse_embedding;
    reverse_embedding.reserve(MAXIMUM_QUERY_GRAPH_SIZE * 2);
    VertexID start_vertex = order[0];

    for (ui i = 0; i < candidates_count[start_vertex]; ++i)
    {
        VertexID v = candidates[start_vertex][i];
        embedding[start_vertex] = v;
        idx_embedding[start_vertex] = i;
        visited_vertices[v] = true;
        reverse_embedding[v] = start_vertex;
        call_count += 1;

        exploreDPisoBacktrack(max_depth, 1, start_vertex, tree, idx_embedding, embedding, reverse_embedding,
                              visited_vertices, idx_count, valid_candidate_idx, edge_matrix,
                              ancestors, dpiso_min_pq(extendable_vertex_compare), weight_array, temp_buffer, extendable,
                              candidates, embedding_cnt,
                              call_count, nullptr);

        visited_vertices[v] = false;
        reverse_embedding.erase(v);
    }

    releaseBuffer(max_depth, idx, idx_count, embedding, idx_embedding, temp_buffer, valid_candidate_idx,
                  visited_vertices,
                  bn, bn_count);

    return embedding_cnt;
}

std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>
EvaluateQuery::exploreDPisoBacktrack(ui max_depth, ui depth, VertexID mapped_vertex, TreeNode *tree, ui *idx_embedding,
                                     ui *embedding, std::unordered_map<VertexID, VertexID> &reverse_embedding,
                                     bool *visited_vertices, ui *idx_count, ui **valid_candidate_index,
                                     Edges ***edge_matrix,
                                     std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> &ancestors,
                                     dpiso_min_pq rank_queue, ui **weight_array, ui *&temp_buffer, ui *extendable,
                                     ui **candidates, size_t &embedding_count, size_t &call_count,
                                     const Graph *query_graph)
{

    TreeNode &node = tree[mapped_vertex];
    for (ui i = 0; i < node.fn_count_; ++i)
    {
        VertexID u = node.fn_[i];
        extendable[u] -= 1;
        if (extendable[u] == 0)
        {
            generateValidCandidateIndex(u, idx_embedding, idx_count, valid_candidate_index[u], edge_matrix, tree[u].bn_,
                                        tree[u].bn_count_, temp_buffer);

            ui weight = 0;
            for (ui j = 0; j < idx_count[u]; ++j)
            {
                ui idx = valid_candidate_index[u][j];
                weight += weight_array[u][idx];
            }
            rank_queue.emplace(std::make_pair(std::make_pair(u, query_graph->getVertexDegree(u)), weight));
        }
    }

    VertexID u = rank_queue.top().first.first;
    rank_queue.pop();

    if (idx_count[u] == 0)
    {
        restoreExtendableVertex(tree, mapped_vertex, extendable);
        return ancestors[u];
    }
    else
    {
        std::bitset<MAXIMUM_QUERY_GRAPH_SIZE> current_fs;
        std::bitset<MAXIMUM_QUERY_GRAPH_SIZE> child_fs;

        for (ui i = 0; i < idx_count[u]; ++i)
        {
            ui valid_index = valid_candidate_index[u][i];
            VertexID v = candidates[u][valid_index];

            if (!visited_vertices[v])
            {
                embedding[u] = v;
                idx_embedding[u] = valid_index;
                visited_vertices[v] = true;
                reverse_embedding[v] = u;
                if (depth != max_depth - 1)
                {
                    call_count += 1;
                    child_fs = exploreDPisoBacktrack(max_depth, depth + 1, u, tree, idx_embedding, embedding,
                                                     reverse_embedding, visited_vertices, idx_count,
                                                     valid_candidate_index, edge_matrix,
                                                     ancestors, rank_queue, weight_array, temp_buffer, extendable,
                                                     candidates, embedding_count,
                                                     call_count, query_graph);
                }
                else
                {
                    embedding_count += 1;
                    child_fs.set();
                }
                visited_vertices[v] = false;
                reverse_embedding.erase(v);

                if (!child_fs.test(u))
                {
                    current_fs = child_fs;
                    break;
                }
            }
            else
            {
                child_fs.reset();
                child_fs |= ancestors[u];
                child_fs |= ancestors[reverse_embedding[v]];
            }

            current_fs |= child_fs;
        }

        restoreExtendableVertex(tree, mapped_vertex, extendable);
        return current_fs;
    }
}

size_t
EvaluateQuery::exploreCECIStyle(const Graph *data_graph, const Graph *query_graph, TreeNode *tree, ui **candidates,
                                ui *candidates_count,
                                std::vector<std::unordered_map<VertexID, std::vector<VertexID>>> &TE_Candidates,
                                std::vector<std::vector<std::unordered_map<VertexID, std::vector<VertexID>>>> &NTE_Candidates,
                                ui *order, size_t &output_limit_num, size_t &call_count,
                                const std::unordered_map<VertexID, std::pair<std::set<VertexID>, std::set<VertexID>>> &ordered_constraints)
{

    ui max_depth = query_graph->getVerticesCount();
    ui data_vertices_count = data_graph->getVerticesCount();
    ui max_valid_candidates_count = 0;
    for (ui i = 0; i < max_depth; ++i)
    {
        if (candidates_count[i] > max_valid_candidates_count)
        {
            max_valid_candidates_count = candidates_count[i];
        }
    }

    ui *idx = new ui[max_depth];
    ui *idx_count = new ui[max_depth];
    ui *embedding = new ui[max_depth];
    ui *temp_buffer = new ui[max_valid_candidates_count];
    ui **valid_candidates = new ui *[max_depth];
    for (ui i = 0; i < max_depth; ++i)
    {
        valid_candidates[i] = new ui[max_valid_candidates_count];
    }
    bool *visited_vertices = new bool[data_vertices_count];
    std::fill(visited_vertices, visited_vertices + data_vertices_count, false);

    size_t embedding_cnt = 0;
    int cur_depth = 0;
    VertexID start_vertex = order[0];

    idx[cur_depth] = 0;
    idx_count[cur_depth] = candidates_count[start_vertex];

    for (ui i = 0; i < idx_count[cur_depth]; ++i)
    {
        valid_candidates[cur_depth][i] = candidates[start_vertex][i];
    }

#ifdef ENABLE_FAILING_SET
    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> ancestors;
    computeAncestor(query_graph, order, ancestors);
    AncestorWithSymmetryBreaking(query_graph, ancestors, ordered_constraints);
    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> vec_failing_set(max_depth);
    std::unordered_map<VertexID, VertexID> reverse_embedding;
    reverse_embedding.reserve(MAXIMUM_QUERY_GRAPH_SIZE * 2);
#endif

    while (true)
    {
        while (idx[cur_depth] < idx_count[cur_depth])
        {
            VertexID u = order[cur_depth];
            VertexID v = valid_candidates[cur_depth][idx[cur_depth]];
            idx[cur_depth] += 1;

            if (visited_vertices[v])
            {
#ifdef ENABLE_FAILING_SET
                vec_failing_set[cur_depth] = ancestors[u];
                vec_failing_set[cur_depth] |= ancestors[reverse_embedding[v]];
                vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
#endif
                continue;
            }

            embedding[u] = v;
            visited_vertices[v] = true;

#ifdef ENABLE_FAILING_SET
            reverse_embedding[v] = u;
#endif
            if (cur_depth == max_depth - 1)
            {
                embedding_cnt += 1;
                visited_vertices[v] = false;
#ifdef ENABLE_FAILING_SET
                reverse_embedding.erase(embedding[u]);
                vec_failing_set[cur_depth].set();
                vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
#endif
                if (embedding_cnt >= output_limit_num)
                {
                    goto EXIT;
                }
            }
            else
            {
                call_count += 1;
                cur_depth += 1;
                idx[cur_depth] = 0;
                generateValidCandidates(cur_depth, embedding, idx_count, valid_candidates, order, temp_buffer, tree,
                                        TE_Candidates,
                                        NTE_Candidates);
                pruneCandidatesBySymmetryBreaking(cur_depth, embedding, order,
                                                  idx_count, valid_candidates, ordered_constraints);
#ifdef ENABLE_FAILING_SET
                if (idx_count[cur_depth] == 0)
                {
                    vec_failing_set[cur_depth - 1] = ancestors[order[cur_depth]];
                }
                else
                {
                    vec_failing_set[cur_depth - 1].reset();
                }
#endif
            }
        }

        cur_depth -= 1;
        if (cur_depth < 0)
            break;
        else
        {
            VertexID u = order[cur_depth];
#ifdef ENABLE_FAILING_SET
            reverse_embedding.erase(embedding[u]);
            if (cur_depth != 0)
            {
                if (!vec_failing_set[cur_depth].test(u))
                {
                    vec_failing_set[cur_depth - 1] = vec_failing_set[cur_depth];
                    idx[cur_depth] = idx_count[cur_depth];
                }
                else
                {
                    vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
                }
            }
#endif
            visited_vertices[embedding[u]] = false;
        }
    }

EXIT:
    delete[] idx;
    delete[] idx_count;
    delete[] embedding;
    delete[] temp_buffer;
    delete[] visited_vertices;
    for (ui i = 0; i < max_depth; ++i)
    {
        delete[] valid_candidates[i];
    }
    delete[] valid_candidates;

    return embedding_cnt;
}

void EvaluateQuery::generateValidCandidates(ui depth, ui *embedding, ui *idx_count, ui **valid_candidates, ui *order,
                                            ui *&temp_buffer, TreeNode *tree,
                                            std::vector<std::unordered_map<VertexID, std::vector<VertexID>>> &TE_Candidates,
                                            std::vector<std::vector<std::unordered_map<VertexID, std::vector<VertexID>>>> &NTE_Candidates)
{

    VertexID u = order[depth];
    idx_count[depth] = 0;
    ui valid_candidates_count = 0;
    {
        VertexID u_p = tree[u].parent_;
        VertexID v_p = embedding[u_p];

        auto iter = TE_Candidates[u].find(v_p);
        if (iter == TE_Candidates[u].end() || iter->second.empty())
        {
            return;
        }

        valid_candidates_count = iter->second.size();
        VertexID *v_p_nbrs = iter->second.data();

        for (ui i = 0; i < valid_candidates_count; ++i)
        {
            valid_candidates[depth][i] = v_p_nbrs[i];
        }
    }
    ui temp_count;
    for (ui i = 0; i < tree[u].bn_count_; ++i)
    {
        VertexID u_p = tree[u].bn_[i];
        VertexID v_p = embedding[u_p];

        auto iter = NTE_Candidates[u][u_p].find(v_p);
        if (iter == NTE_Candidates[u][u_p].end() || iter->second.empty())
        {
            return;
        }

        ui current_candidates_count = iter->second.size();
        ui *current_candidates = iter->second.data();

        ComputeSetIntersection::ComputeCandidates(current_candidates, current_candidates_count,
                                                  valid_candidates[depth], valid_candidates_count,
                                                  temp_buffer, temp_count);

        std::swap(temp_buffer, valid_candidates[depth]);
        valid_candidates_count = temp_count;
    }

    idx_count[depth] = valid_candidates_count;
}

void EvaluateQuery::computeAncestor(const Graph *query_graph, ui **bn, ui *bn_cnt, VertexID *order,
                                    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> &ancestors)
{
    ui query_vertices_num = query_graph->getVerticesCount();
    ancestors.resize(query_vertices_num);

    for (ui i = 0; i < query_vertices_num; ++i)
    {
        VertexID u = order[i];
        ancestors[u].set(u);
        for (ui j = 0; j < bn_cnt[i]; ++j)
        {
            VertexID u_bn = bn[i][j];
            ancestors[u] |= ancestors[u_bn];
        }
    }
}

void EvaluateQuery::computeAncestor(const Graph *query_graph, VertexID *order,
                                    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> &ancestors)
{
    ui query_vertices_num = query_graph->getVerticesCount();
    ancestors.resize(query_vertices_num);

    for (ui i = 0; i < query_vertices_num; ++i)
    {
        VertexID u = order[i];
        ancestors[u].set(u);
        for (ui j = 0; j < i; ++j)
        {
            VertexID u_bn = order[j];
            if (query_graph->checkEdgeExistence(u, u_bn))
            {
                ancestors[u] |= ancestors[u_bn];
            }
        }
    }
}

size_t EvaluateQuery::exploreVF3Style(const Graph *data_graph, const Graph *query_graph, ui **candidates,
                                      ui *candidates_count, ui *order,
                                      ui *pivot, size_t output_limit_num, size_t &call_count,
                                      const std::unordered_map<VertexID, std::pair<std::set<VertexID>, std::set<VertexID>>> &ordered_constraints)
{
    ui query_vertices_num = query_graph->getVerticesCount();
    ui query_labels_num = query_graph->getLabelsCount();
    ui query_max_label_fre = query_graph->getGraphMaxLabelFrequency();
    ui data_vertices_num = data_graph->getVerticesCount();
    ui data_labels_num = data_graph->getLabelsCount();
    ui data_max_label_fre = data_graph->getGraphMaxLabelFrequency();

    ui ***query_feasibility = NULL;
    ui **query_feasibility_count = NULL;
    query_feasibility = new ui **[query_vertices_num + 1];
    query_feasibility_count = new ui *[query_vertices_num + 1];
    for (ui i = 0; i <= query_vertices_num; i++)
    {
        query_feasibility[i] = new ui *[query_labels_num];
        query_feasibility_count[i] = new ui[query_labels_num];
        memset(query_feasibility_count[i], 0, query_labels_num * sizeof(ui));
        for (ui j = 0; j < query_labels_num; j++)
        {
            query_feasibility[i][j] = new ui[query_max_label_fre];
        }
    }

    generateFeasibility(query_graph, order, query_feasibility, query_feasibility_count);

    ui ***data_feasibility = NULL;
    ui **data_feasibility_count = NULL;
    data_feasibility = new ui **[query_vertices_num + 1];
    data_feasibility_count = new ui *[query_vertices_num + 1];
    for (ui i = 0; i <= query_vertices_num; i++)
    {
        data_feasibility[i] = new ui *[data_labels_num];
        data_feasibility_count[i] = new ui[data_labels_num];
        memset(data_feasibility_count[i], 0, data_labels_num * sizeof(ui));
        for (ui j = 0; j < data_labels_num; j++)
        {
            data_feasibility[i][j] = new ui[data_max_label_fre];
        }
    }

    ui max_depth = query_graph->getVerticesCount();
    ui *idx = new ui[max_depth];
    ui *idx_count = new ui[max_depth];
    ui *embedding = new ui[max_depth];
    ui embedding_cnt = 0;
    ui depth = 0;
    ui **valid_candidates = new ui *[max_depth];

    for (ui i = 0; i < max_depth; ++i)
    {
        valid_candidates[i] = new ui[data_max_label_fre];
    }
    bool *visited_vertices = new bool[data_vertices_num];
    memset(visited_vertices, 0, data_vertices_num * sizeof(bool));
    exploreVF3Backtrack(data_graph, query_graph, order, pivot, output_limit_num, call_count, embedding_cnt,
                        depth, max_depth, idx, idx_count, embedding, valid_candidates, visited_vertices,
                        query_feasibility, query_feasibility_count, data_feasibility, data_feasibility_count, ordered_constraints);

    for (ui i = 0; i <= query_vertices_num; i++)
    {
        for (ui j = 0; j < query_labels_num; j++)
        {
            delete[] query_feasibility[i][j];
        }
        delete[] query_feasibility[i];
        delete[] query_feasibility_count[i];
    }
    delete[] query_feasibility;
    delete[] query_feasibility_count;
    for (ui i = 0; i <= query_vertices_num; i++)
    {
        for (ui j = 0; j < data_labels_num; j++)
        {
            delete[] data_feasibility[i][j];
        }
        delete[] data_feasibility[i];
        delete[] data_feasibility_count[i];
    }
    delete[] data_feasibility;
    delete[] data_feasibility_count;
    delete[] idx;
    delete[] idx_count;
    delete[] embedding;
    for (ui i = 0; i < max_depth; i++)
    {
        delete[] valid_candidates[i];
    }
    delete[] valid_candidates;
    delete[] visited_vertices;
    return embedding_cnt;
}

void EvaluateQuery::generateFeasibility(const Graph *graph, ui v, bool *matched, ui level, ui ***feasibility, ui **feasibility_count)
{
    if (level < (ui)1)
    {

        exit(1);
    }

    ui labels_num = graph->getLabelsCount();
    for (ui j = 0; j < labels_num; j++)
    {
        feasibility_count[level][j] = 0;

        for (ui k = 0; k < feasibility_count[level - 1][j]; k++)
        {
            if (feasibility[level - 1][j][k] != v)
            {
                feasibility[level][j][feasibility_count[level][j]++] = feasibility[level - 1][j][k];
            }
        }
    }

    ui nbrs_cnt;
    const ui *nbrs = graph->getVertexNeighbors(v, nbrs_cnt);
    for (ui j = 0; j < nbrs_cnt; j++)
    {
        if (!matched[nbrs[j]])
        {
            ui label = graph->getVertexLabel(nbrs[j]);
            bool f_add = true;
            for (ui k = 0; k < feasibility_count[level][label]; k++)
            {
                if (nbrs[j] == feasibility[level][label][k])
                {
                    f_add = false;
                    break;
                }
            }
            if (f_add == true)
            {
                ui count = feasibility_count[level][label]++;
                feasibility[level][label][count] = nbrs[j];
            }
        }
    }
}

void EvaluateQuery::generateFeasibility(const Graph *query_graph, const ui *order, ui ***feasibility, ui **feasibility_count)
{
    ui vertices_num = query_graph->getVerticesCount();

    bool *query_matched = new bool[vertices_num];
    memset(query_matched, 0, vertices_num * sizeof(bool));

    for (ui i = 1; i < vertices_num; i++)
    {

        ui u = order[i - 1];
        query_matched[u] = true;
        generateFeasibility(query_graph, u, query_matched, i, feasibility, feasibility_count);
    }
    delete[] query_matched;
}

bool EvaluateQuery::exploreVF3Backtrack(const Graph *data_graph, const Graph *query_graph,
                                        ui *order, ui *pivot, size_t output_limit_num, size_t &call_count,
                                        ui &embedding_cnt, ui depth, ui max_depth, ui *idx, ui *idx_count,
                                        ui *embedding, ui **valid_candidates, bool *visited_vertices,
                                        ui ***query_feasibility, ui **query_feasibility_count,
                                        ui ***data_feasibility, ui **data_feasibility_count,
                                        const std::unordered_map<VertexID, std::pair<std::set<VertexID>, std::set<VertexID>>> &ordered_constraints)
{
    if (depth == max_depth)
    {

        embedding_cnt++;
        if (embedding_cnt >= output_limit_num)
        {
            for (ui i = 0; i < depth; i++)
            {
                idx[i] = idx_count[i];
            }
            return true;
        }
        return true;
    }
    ui u = order[depth];

    const ui *candidates;
    ui candidate_count;
    idx[depth] = 0;
    idx_count[depth] = 0;

    if (pivot[depth] == (ui)-1)
    {

        ui query_label = query_graph->getVertexLabel(u);
        candidates = data_graph->getVerticesByLabel(query_label, candidate_count);
        for (ui i = 0; i < candidate_count; i++)
        {
            if (!visited_vertices[candidates[i]])
            {
                valid_candidates[depth][idx_count[depth]++] = candidates[i];
            }
        }
    }
    else
    {

        if (!visited_vertices[embedding[pivot[depth]]])
        {

            exit(1);
        }

        ui v = embedding[pivot[depth]];
        candidates = data_graph->getVertexNeighbors(v, candidate_count);
        ui query_label = query_graph->getVertexLabel(u);
        for (ui i = 0; i < candidate_count; i++)
        {
            ui data_label = data_graph->getVertexLabel(candidates[i]);
            if (!visited_vertices[candidates[i]] && data_label == query_label)
            {
                valid_candidates[depth][idx_count[depth]++] = candidates[i];
            }
        }
    }

    pruneCandidatesBySymmetryBreaking(depth, embedding, order,
                                      idx_count, valid_candidates, ordered_constraints);

    bool result = false;
    while (idx[depth] < idx_count[depth])
    {

        if (isFeasibility(data_graph, query_graph, depth, order[depth], valid_candidates[depth][idx[depth]],
                          embedding, order, query_feasibility, query_feasibility_count,
                          data_feasibility, data_feasibility_count) == true)
        {
            embedding[order[depth]] = valid_candidates[depth][idx[depth]];
            visited_vertices[embedding[order[depth]]] = true;
            call_count++;
            if (exploreVF3Backtrack(data_graph, query_graph, order, pivot, output_limit_num, call_count, embedding_cnt,
                                    depth + 1, max_depth, idx, idx_count, embedding, valid_candidates, visited_vertices,
                                    query_feasibility, query_feasibility_count, data_feasibility, data_feasibility_count, ordered_constraints) == true)
                result = true;
            visited_vertices[embedding[order[depth]]] = false;
        }
        idx[depth]++;
    }
    return result;
}

bool EvaluateQuery::isFeasibility(const Graph *data_graph, const Graph *query_graph, ui depth, ui cur_u, ui cur_v,
                                  ui *embedding, ui *order, ui ***query_feasibility, ui **query_feasibility_count,
                                  ui ***data_feasibility, ui **data_feasibility_count)
{
    ui data_vertices_num = data_graph->getVerticesCount();
    ui query_vertices_num = query_graph->getVerticesCount();
    ui data_labels_num = data_graph->getLabelsCount();
    ui query_labels_num = query_graph->getLabelsCount();
    ui labels_num = query_labels_num;
    if (data_labels_num < query_labels_num)
    {

        exit(1);
    }

    bool *data_matched = new bool[data_vertices_num];
    bool *query_matched = new bool[query_vertices_num];
    memset(data_matched, false, data_vertices_num * sizeof(bool));
    memset(query_matched, false, query_vertices_num * sizeof(bool));
    for (ui i = 0; i < depth; i++)
    {
        data_matched[embedding[order[i]]] = true;
        query_matched[order[i]] = true;
    }

    generateFeasibility(data_graph, cur_v, data_matched, depth + 1, data_feasibility, data_feasibility_count);

    bool *data_fmatched = new bool[data_vertices_num];
    bool *query_fmatched = new bool[query_vertices_num];
    memset(data_fmatched, false, data_vertices_num * sizeof(bool));
    memset(query_fmatched, false, query_vertices_num * sizeof(bool));
    for (ui i = 0; i < query_labels_num; i++)
    {
        for (ui j = 0; j < query_feasibility_count[depth + 1][i]; j++)
        {
            query_fmatched[query_feasibility[depth + 1][i][j]] = true;
        }
    }
    for (ui i = 0; i < data_labels_num; i++)
    {
        for (ui j = 0; j < data_feasibility_count[depth + 1][i]; j++)
        {
            data_fmatched[data_feasibility[depth + 1][i][j]] = true;
        }
    }

    ui unbrs_count = 0;
    const ui *unbrs = query_graph->getVertexNeighbors(cur_u, unbrs_count);
    ui vnbrs_count = 0;
    const ui *vnbrs = data_graph->getVertexNeighbors(cur_v, vnbrs_count);
    if (vnbrs < unbrs)
        return false;

    ui *unbrs_flabeled_num = new ui[query_labels_num];
    ui *vnbrs_flabeled_num = new ui[data_labels_num];
    memset(unbrs_flabeled_num, 0, query_labels_num * sizeof(ui));
    memset(vnbrs_flabeled_num, 0, data_labels_num * sizeof(ui));

    ui *unbrs_labeled_num = new ui[query_labels_num];
    ui *vnbrs_labeled_num = new ui[data_labels_num];
    memset(unbrs_labeled_num, 0, query_labels_num * sizeof(ui));
    memset(vnbrs_labeled_num, 0, data_labels_num * sizeof(ui));

    for (ui i = 0; i < unbrs_count; i++)
    {
        ui unbr = unbrs[i];
        ui unbr_label = query_graph->getVertexLabel(unbr);

        if (query_matched[unbr] == true)
        {
            if (!data_graph->checkEdgeExistence(cur_v, embedding[unbr]))
                return false;
        }
        else if (query_fmatched[unbr] == true)
        {
            unbrs_flabeled_num[unbr_label]++;
        }
        else
        {
            unbrs_labeled_num[unbr_label]++;
        }
    }

    for (ui i = 0; i < vnbrs_count; i++)
    {
        ui vnbr = vnbrs[i];
        ui vnbr_label = data_graph->getVertexLabel(vnbr);

        if (data_matched[vnbr] == true)
        {
            true;
        }
        else if (data_fmatched[vnbr] == true)
        {
            vnbrs_flabeled_num[vnbr_label]++;
        }
        else
        {
            vnbrs_labeled_num[vnbr_label]++;
        }
    }

    for (ui i = 0; i < labels_num; i++)
    {
        if (vnbrs_labeled_num[i] < unbrs_labeled_num[i])
        {

            return false;
        }
        if (vnbrs_flabeled_num[i] < unbrs_flabeled_num[i])
        {

            return false;
        }
    }

    return true;
}

enumResult EvaluateQuery::exploreVEQStyle1(const Graph *data_graph, const Graph *query_graph, TreeNode *tree,
                                           Edges ***edge_matrix, ui **candidates, ui *candidates_count,
                                           size_t output_limit_num, size_t &call_count, int TimeL,
                                           const std::unordered_map<VertexID, std::pair<std::set<VertexID>, std::set<VertexID>>> &full_constraints)
{

    auto start = std::chrono::high_resolution_clock::now();
    ui query_vertices_num = query_graph->getVerticesCount();
    ui data_vertices_num = data_graph->getVerticesCount();
    ui **nec = new ui *[query_vertices_num];
    memset(nec, 0, sizeof(ui *) * query_vertices_num);
    computeNEC(query_graph, nec);
    ui vec_count = 0;
    enumResult s;

    int depth = 0;
    ui max_depth = query_graph->getVerticesCount();
    ui *embedding = new ui[query_vertices_num];
    ui *order = new ui[max_depth];
    ui embedding_cnt = 0;
    ui *extendable = new ui[query_vertices_num];
    for (ui i = 0; i < query_vertices_num; ++i)
    {
        extendable[i] = tree[i].bn_count_;
    }
    ui *idx = new ui[max_depth];
    ui *idx_count = new ui[max_depth];
    ui **valid_cans = new ui *[query_vertices_num];
    for (ui i = 0; i < query_vertices_num; i++)
    {
        valid_cans[i] = new ui[candidates_count[i]];
    }
    ui *valid_cans_count = new ui[query_vertices_num];
    bool *visited_u = new bool[query_vertices_num];
    memset(visited_u, 0, query_vertices_num * sizeof(bool));
    bool *visited_vertices = new bool[data_vertices_num];
    memset(visited_vertices, 0, data_vertices_num * sizeof(bool));
    bool *violated_symmetry_u = new bool[max_depth];

    priority_queue<int, vector<int>, greater<int>> pq;
    pq.push(0);
    int heapcount = 1;
    int ksize = 1000;
    bool Match_BA_t[max_depth] = {false};
    unordered_set<ui> EmbSum_t;
    int UNPM = 0;

    TimeL = TimeL * 1000;
    double ens = 0;

    ui **TM = new ui *[query_vertices_num];
    for (ui i = 0; i < query_vertices_num; i++)
    {
        TM[i] = new ui[candidates_count[i] + 1];
    }
    memset(TM[0], 0, sizeof(ui) * (candidates_count[0] + 1));
    // TM has size TM[qsiz][candi[i]]+1 element.
    // TM initialized the first row with zeros?
    std::vector<std::vector<ui>> vec_index(query_vertices_num);
    for (ui i = 0; i < query_vertices_num; i++)
    {
        vec_index[i].resize(candidates_count[i]);
        std::fill(vec_index[i].begin(), vec_index[i].end(), (ui)-1);
    }
    // vec_index is initialized with largest number ui-1
    std::vector<std::vector<ui>> vec_set;
    // computeNEC(query_graph, edge_matrix, candidates_count, candidates, vec_index, vec_set);

    std::vector<std::vector<ui>> pi_m_index(query_vertices_num);
    for (ui i = 0; i < query_vertices_num; i++)
    {
        pi_m_index[i].resize(candidates_count[i], (ui)-1);
    }
    std::vector<std::vector<ui>> pi_m;
    ui *pi_m_count = new ui[max_depth];
    pi_m_count[0] = 0;

    std::vector<std::vector<ui>> dm(query_vertices_num);
    std::unordered_map<VertexID, VertexID> reverse_embedding;
    reverse_embedding.reserve(query_vertices_num);

    ui start_vertex = (ui)-1;
    for (ui i = 0; i < query_vertices_num; i++)
    {
        if (extendable[i] == 0)
        {
            start_vertex = i;
            break;
        }
    }
    assert(start_vertex != (ui)-1);

    order[depth] = start_vertex;
    visited_u[start_vertex] = true;

    idx[depth] = 0;
    valid_cans_count[start_vertex] = candidates_count[start_vertex];
    idx_count[depth] = valid_cans_count[start_vertex];
    for (ui i = 0; i < candidates_count[start_vertex]; i++)
    {
        valid_cans[start_vertex][i] = candidates[start_vertex][i];
    }

    std::fill(pi_m_index[start_vertex].begin(), pi_m_index[start_vertex].end(), (ui)-1);

    bool Match_BA[max_depth] = {false};
    unordered_set<ui> EmbSum;

    while (true)
    {
        while (idx[depth] < idx_count[depth])
        {
            Match_BA[depth] = false;
            Match_BA_t[depth] = false;
            auto end = std::chrono::high_resolution_clock::now();
            ens = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
            if (ens > TimeL)
            { // 1000 1 sec
                goto EXIT;
            }
#ifdef ENABLE_DYNAMIC_CANS

            // if (idx[order[depth]] != valid_cans[order[depth]].begin())
            if (idx[order[depth]] != valid_cans[order[depth]][0])
            {
                RestoreValidCans(query_graph, data_graph, visited_u, order[depth], embedding[order[depth]], valid_cans);
            }
#endif
            VertexID u = order[depth];
            VertexID v = valid_cans[u][idx[depth]];
            ui v_index1 = 0;
            for (; v_index1 < candidates_count[u]; v_index1++)
            {
                if (candidates[u][v_index1] == v)
                    break;
            }
            if (vec_index[u][v_index1] == (ui)-1)
            {
                computeNEC_Adaptive(query_graph, edge_matrix, candidates_count, candidates, vec_index, vec_set, u, v_index1, vec_count);
                vec_count++;
            }

            bool is_voilate_symmetry = false;

            if (!full_constraints.empty())
            {
                is_voilate_symmetry = checkSingeVertexBySymmetryBreaking(depth, embedding, order, visited_u, v,
                                                                         full_constraints,
                                                                         violated_symmetry_u,
                                                                         max_depth);
            }

            if (visited_vertices[v] || is_voilate_symmetry)
            {

#ifdef ENABLE_EQUIVALENT_SET
                TM[u][idx[depth]] = 0;
                // keep only common in conflix- like pmnegative
                if (visited_vertices[v])
                {
                    VertexID con_u = reverse_embedding[v];

                    ui con_v_index = 0, v_index = 0;
                    for (; con_v_index < valid_cans_count[con_u]; con_v_index++)
                    {
                        if (valid_cans[con_u][con_v_index] == v)
                            break;
                    }
                    for (; v_index < candidates_count[u]; v_index++)
                    {
                        if (candidates[u][v_index] == v)
                            break;
                    }

                    assert(con_v_index < valid_cans_count[con_u]);
                    assert(v_index < candidates_count[u]);

                    auto &con_uv_idx = pi_m_index[con_u][con_v_index];

                    for (ui i = 0; i < pi_m[con_uv_idx].size(); i++)
                    {
                        bool f_in = false;
                        for (ui j = 0; j < vec_set[vec_index[u][v_index]].size(); j++)
                        {
                            if (vec_set[vec_index[u][v_index]][j] == pi_m[con_uv_idx][i])
                            {
                                f_in = true;
                                break;
                            }
                        }
                        if (f_in == false)
                        {
                            pi_m[con_uv_idx][i] = pi_m[con_uv_idx][pi_m[con_uv_idx].size() - 1];
                            pi_m[con_uv_idx].pop_back();
                            i--;
                        }
                    }
                }
                if (is_voilate_symmetry)
                {
                    for (ui i = 0; i < depth; ++i)
                    {

                        if (violated_symmetry_u[order[i]])
                        {
                            VertexID sym_u = order[i];
                            VertexID sym_v = embedding[sym_u];
                            ui sym_v_index = 0, v_index = 0;
                            for (; sym_v_index < valid_cans_count[sym_u]; sym_v_index++)
                            {
                                if (valid_cans[sym_u][sym_v_index] == sym_v)
                                    break;
                            }
                            for (; v_index < candidates_count[u]; v_index++)
                            {
                                if (candidates[u][v_index] == sym_v)
                                    break;
                            }

                            auto &sym_uv_idx = pi_m_index[sym_u][sym_v_index];

                            for (ui i = 0; i < pi_m[sym_uv_idx].size(); i++)
                            {
                                bool f_in = false;
                                if ((vec_index[u][v_index]) == (ui)-1)
                                {
                                    computeNEC_Adaptive(query_graph, edge_matrix, candidates_count, candidates, vec_index, vec_set, u, v_index, vec_count);
                                    vec_count++;
                                }
                                for (ui j = 0; j < vec_set[vec_index[u][v_index]].size(); j++)
                                {
                                    if (vec_set[vec_index[u][v_index]][j] == pi_m[sym_uv_idx][i])
                                    {
                                        f_in = true;
                                        break;
                                    }
                                }
                                if (f_in == false)
                                {
                                    pi_m[sym_uv_idx][i] = pi_m[sym_uv_idx][pi_m[sym_uv_idx].size() - 1];
                                    pi_m[sym_uv_idx].pop_back();
                                    i--;
                                }
                            }
                        }
                    }
                }

#endif

                idx[depth]++;
                continue;
            }

#ifdef ENABLE_EQUIVALENT_SET
            if (pi_m_index[u][idx[depth]] != (ui)-1)
            {
                VertexID equ_v = pi_m[pi_m_index[u][idx[depth]]][0];
                ui equ_v_index = 0;

                for (; equ_v_index < idx_count[depth]; equ_v_index++)
                {
                    if (equ_v == valid_cans[u][equ_v_index])
                        break;
                }
                assert(equ_v_index < idx_count[depth]);

                embedding_cnt += TM[u][equ_v_index];
                TM[u][candidates_count[u]] += TM[u][equ_v_index];
                idx[depth]++;

                continue;
            }
            reverse_embedding[v] = u;
#endif
            embedding[u] = v;
            visited_vertices[v] = true;
            ui cur_idx = idx[depth]++;
#ifdef ENABLE_EQUIVALENT_SET
            pi_m_index[u][cur_idx] = pi_m_count[depth]++;
            ui v_index = 0;
            for (; v_index < candidates_count[u]; v_index++)
            {
                if (candidates[u][v_index] == v)
                    break;
            }
            if ((vec_index[u][v_index]) == (ui)-1)
            {
                computeNEC_Adaptive(query_graph, edge_matrix, candidates_count, candidates, vec_index, vec_set, u, v_index, vec_count);
                vec_count++;
            }
            auto &pi = vec_set[vec_index[u][v_index]];
            pi_m.push_back(pi);
            assert(pi_m.size() == pi_m_count[depth]);
            dm[u].clear();
            for (ui i = 0; i < depth; i++)
            {
                bool va_in_pi = false, sec_empty = true;
                ui va_index = 0;
                VertexID ua = order[i];
                VertexID va = embedding[ua];
                for (; va_index < candidates_count[ua]; va_index++)
                {
                    if (candidates[ua][va_index] == va)
                        break;
                }
                assert(va_index < candidates_count[ua]);
                if ((vec_index[ua][va_index]) == (ui)-1)
                {
                    computeNEC_Adaptive(query_graph, edge_matrix, candidates_count, candidates, vec_index, vec_set, ua, v_index, vec_count);
                    vec_count++;
                }
                auto &pi_a = vec_set[vec_index[ua][va_index]];
                for (ui j = 0; j < pi.size(); j++)
                {
                    if (pi[j] == embedding[order[i]])
                    {
                        va_in_pi = true;
                        break;
                    }
                    for (ui k = 0; k < pi_a.size(); k++)
                    {
                        if (pi_a[k] == pi[j])
                        {
                            sec_empty = false;
                            break;
                        }
                    }
                    if (sec_empty == false)
                        break;
                }
                if (va_in_pi == false && sec_empty == false)
                {
                    ui dm_size = dm[ua].size();
                    for (ui j = 0; j < pi.size(); j++)
                    {
                        bool f_add = true;
                        for (ui k = 0; k < dm_size; k++)
                        {
                            if (dm[ua][k] == pi[j])
                            {
                                f_add = false;
                                break;
                            }
                        }
                        if (f_add == true)
                        {
                            dm[ua].push_back(pi[j]);
                        }
                    }
                }
            }
#endif
            if (depth == max_depth - 1)
            {
                embedding_cnt += 1;
                visited_vertices[v] = false;
                visited_u[u] = false;
                // EmbSum.insert(embedding[order[depth]]);
                int ao = depth; //-1;
                ui vqo = order[ao];
                int UNPMtemp = 0;
                int UNPM2 = 0;
                UNPM = 0;
                while (ao >= 0)
                {
                    UNPMtemp = EmbSum.insert(embedding[vqo]).second;

                    if (Match_BA_t[ao] == false)
                    {
                        if (pi_m_index[vqo][idx[ao] - 1] != (ui)-1)
                        {
                            for (auto eq_node : pi_m[pi_m_index[vqo][idx[ao] - 1]])
                            {
                                EmbSum_t.insert(eq_node);
                            }
                        }
                        Match_BA_t[ao] = true;
                    }

                    if (UNPMtemp == 0 && Match_BA[ao] == false)
                    {
                        Match_BA[ao] = true;
                        if (pi_m_index[vqo][idx[ao] - 1] != (ui)-1)
                        {
                            for (auto eq_node : pi_m[pi_m_index[vqo][idx[ao] - 1]])
                            {
                                UNPMtemp += EmbSum.insert(eq_node).second;
                                if (UNPMtemp == 1)
                                    break;
                            }
                        }
                        if (UNPMtemp == 1)
                        {
                            Match_BA[ao] = false;
                        }
                    }

                    UNPM += UNPMtemp;
                    ao--;
                    vqo = order[ao];
                    UNPMtemp = 0;
                }
                if (UNPM > pq.top())
                {
                    if (heapcount >= ksize)
                    {
                        pq.pop();
                    }
                    pq.push(UNPM);
                    heapcount++;
                }
                if (embedding_cnt >= output_limit_num)
                {
                    goto EXIT;
                }
#ifdef ENABLE_EQUIVALENT_SET
                // for (auto eq_node : pi_m[pi_m_index[u][idx[depth]]]) {

                reverse_embedding.erase(v);
                TM[u][cur_idx] = 1;
                TM[u][candidates_count[u]]++;
                auto &uv_idx = pi_m_index[u][cur_idx];
                for (ui i = 0; i < pi_m[uv_idx].size(); i++)
                {
                    bool f_del = false;
                    for (ui j = 0; j < dm[u].size(); j++)
                    {
                        if (pi_m[uv_idx][i] == dm[u][j])
                        {
                            f_del = true;
                            break;
                        }
                    }
                    if (f_del == true)
                    {
                        pi_m[uv_idx][i] = pi_m[uv_idx][pi_m[uv_idx].size() - 1];
                        pi_m[uv_idx].pop_back();
                        i--;
                    }
                }
                for (ui i = 1; i < pi_m[uv_idx].size(); i++)
                {
                    ui v_equ_index = 0;
                    for (; v_equ_index < valid_cans_count[u]; v_equ_index++)
                    {
                        if (valid_cans[u][v_equ_index] == pi_m[uv_idx][i])
                            break;
                    }
                    assert(v_equ_index < valid_cans_count[u]);
                    pi_m_index[u][v_equ_index] = uv_idx;
                }

                if (pi_m[uv_idx].size() > 1)
                {
                    auto pi_v_idx = std::find(pi_m[uv_idx].begin(), pi_m[uv_idx].end(), v);
                    assert(pi_v_idx != pi_m[uv_idx].end());
                    ui tmp = *pi_v_idx;
                    *pi_v_idx = pi_m[uv_idx][0];
                    pi_m[uv_idx][0] = tmp;
                }
#endif
            }
            else
            {
                call_count += 1;
                depth += 1;
                order[depth] = generateNextU(data_graph, query_graph, candidates, candidates_count, valid_cans,
                                             valid_cans_count, extendable, nec, depth, embedding,
                                             edge_matrix, visited_vertices, visited_u, order, tree);

                if (order[depth] == (ui)-1)
                {
                    break;
                }
                else
                {
                    visited_u[order[depth]] = true;
                    idx[depth] = 0;
                    idx_count[depth] = valid_cans_count[order[depth]];
                }
#ifdef ENABLE_EQUIVALENT_SET
                memset(TM[order[depth]], 0, sizeof(ui) * (candidates_count[order[depth]] + 1));
                std::fill(pi_m_index[order[depth]].begin(), pi_m_index[order[depth]].end(), (ui)-1);
                pi_m_count[depth] = pi_m_count[depth - 1];
#endif
            }
        }

        depth -= 1;

        if (depth < 0)
            break;
        VertexID u = order[depth];
        ui cur_idx = idx[depth] - 1;
        visited_vertices[embedding[u]] = false;
        restoreExtendableVertex(tree, u, extendable);
        if (order[depth + 1] != (ui)-1)
        {
            VertexID last_u = order[depth + 1];
            visited_u[last_u] = false;
            if (nec[last_u] != NULL)
                (*(nec[last_u]))++;
#ifdef ENABLE_DYNAMIC_CANS
            RestoreValidCans(query_graph, data_graph, visited_u, last_u, last_v, valid_cans);
#endif
        }
#ifdef ENABLE_EQUIVALENT_SET
        if (order[depth + 1] != (ui)-1)
        {
            TM[u][cur_idx] = TM[order[depth + 1]][candidates_count[order[depth + 1]]];
        }
        else
        {
            TM[u][cur_idx] = 0;
        }
        TM[u][candidates_count[u]] += TM[u][cur_idx];
        auto &uv_idx = pi_m_index[u][cur_idx];
        if (TM[u][cur_idx] != 0)
        {

            for (ui i = 0; i < pi_m[uv_idx].size(); i++)
            {
                bool f_del = false;
                for (ui j = 0; j < dm[u].size(); j++)
                {
                    if (pi_m[uv_idx][i] == dm[u][j])
                    {
                        f_del = true;
                        break;
                    }
                }
                if (f_del == true)
                {
                    pi_m[uv_idx][i] = pi_m[uv_idx][pi_m[uv_idx].size() - 1];
                    pi_m[uv_idx].pop_back();
                    i--;
                }
            }

            for (ui i = 0; i < pi_m[uv_idx].size(); i++)
            {
            }

            if (pi_m[uv_idx].size() > 1)
            {
                auto pi_v_idx = std::find(pi_m[uv_idx].begin(), pi_m[uv_idx].end(), embedding[u]);
                assert(pi_v_idx != pi_m[uv_idx].end());
                ui tmp = *pi_v_idx;
                *pi_v_idx = pi_m[uv_idx][0];
                pi_m[uv_idx][0] = tmp;
            }
        }
        for (ui i = 1; i < pi_m[uv_idx].size(); i++)
        {
            ui v_equ_index = 0;
            for (; v_equ_index < valid_cans_count[u]; v_equ_index++)
            {
                if (valid_cans[u][v_equ_index] == pi_m[uv_idx][i])
                    break;
            }
            assert(v_equ_index < valid_cans_count[u]);

            pi_m_index[u][v_equ_index] = uv_idx;
        }
        pi_m.resize(pi_m_count[depth]);
#endif
    }

EXIT:

    for (ui i = 0; i < query_vertices_num; i++)
    {
    }
    delete[] nec;
    delete[] embedding;
    delete[] visited_u;
    delete[] visited_vertices;
    delete[] order;
    delete[] extendable;
    delete[] idx;
    delete[] idx_count;
    for (ui i = 0; i < query_vertices_num; i++)
    {
        delete[] valid_cans[i];
    }
    delete[] valid_cans;
    delete[] valid_cans_count;
#ifdef ENABLE_EQUIVALENT_SET
    delete[] TM;
#endif

    s.embedding_cnt = embedding_cnt;
    s.Can_embed = EmbSum.size();
    int tempC = 0;
    for (auto it = EmbSum_t.begin(); it != EmbSum_t.end(); ++it)
    {
        tempC += EmbSum.insert(*it).second; // Inserting elements into set2
    }
    EmbSum_t.clear();
    while (tempC > 0 && pq.top() < 1)
    {
        if (heapcount >= ksize)
        {
            pq.pop();
        }
        pq.push(1);
        heapcount++;
        tempC--;
    }
#ifdef TOPKGREEDY
    std::priority_queue<int> maxHeap;
    int topkcounter = 0;
    int greedysum = 0;
    while (!pq.empty())
    {
        maxHeap.push(pq.top());
        pq.pop();
    }
    while (!maxHeap.empty())
    {
        greedysum += maxHeap.top();
        maxHeap.pop();
        topkcounter++;
        switch (topkcounter)
        {
        case VALUE_5:
            s.topk[0] = greedysum;
            break;
        case VALUE_10:
            s.topk[1] = greedysum;
            break;
        case VALUE_50:
            s.topk[2] = greedysum;
            break;
        case VALUE_100:
            s.topk[3] = greedysum;
            break;
        case VALUE_250:
            s.topk[4] = greedysum;
            break;
        case VALUE_500:
            s.topk[5] = greedysum;
            break;
        case VALUE_750:
            s.topk[6] = greedysum;
            break;
        case VALUE_1000:
            s.topk[7] = greedysum;
            break;
        default:
            break;
        }
    }
    for (int i = 0; i < 8; i++)
    {
        if (s.topk[i] == 0)
            s.topk[i] = greedysum;
    }
#endif
    return s;
}

enumResult EvaluateQuery::exploreVEQStyle(const Graph *data_graph, const Graph *query_graph, TreeNode *tree,
                                          Edges ***edge_matrix, ui **candidates, ui *candidates_count,
                                          size_t output_limit_num, size_t &call_count, int TimeL, ui *order1,
                                          const std::unordered_map<VertexID, std::pair<std::set<VertexID>, std::set<VertexID>>> &full_constraints)
{

    auto start = std::chrono::high_resolution_clock::now();
    ui query_vertices_num = query_graph->getVerticesCount();
    ui data_vertices_num = data_graph->getVerticesCount();
    ui **nec = new ui *[query_vertices_num];
    memset(nec, 0, sizeof(ui *) * query_vertices_num);
    computeNEC(query_graph, nec);
    ui vec_count = 0;
    enumResult s;

    int depth = 0;
    ui max_depth = query_graph->getVerticesCount();
    ui *embedding = new ui[query_vertices_num];
    ui *order = new ui[max_depth];
    ui embedding_cnt = 0;
    ui *extendable = new ui[query_vertices_num];
    for (ui i = 0; i < query_vertices_num; ++i)
    {
        extendable[i] = tree[i].bn_count_;
    }
    ui *idx = new ui[max_depth];
    ui *idx_count = new ui[max_depth];
    ui **valid_cans = new ui *[query_vertices_num];
    for (ui i = 0; i < query_vertices_num; i++)
    {
        valid_cans[i] = new ui[candidates_count[i]];
    }
    ui *valid_cans_count = new ui[query_vertices_num];
    bool *visited_u = new bool[query_vertices_num];
    memset(visited_u, 0, query_vertices_num * sizeof(bool));
    bool *visited_vertices = new bool[data_vertices_num];
    memset(visited_vertices, 0, data_vertices_num * sizeof(bool));
    bool *violated_symmetry_u = new bool[max_depth];

#ifdef TOPKGREEDY
    priority_queue<int, vector<int>, greater<int>> pq;
    pq.push(0);
    int heapcount = 1;
    int ksize = 1000;
#endif
    int UNPM = 0;

    TimeL = TimeL * 1000;
    double ens = 0;
#ifdef ENABLE_EQUIVALENT_SET
    ui **TM = new ui *[query_vertices_num];
    for (ui i = 0; i < query_vertices_num; i++)
    {
        TM[i] = new ui[candidates_count[i] + 1];
    }
    memset(TM[0], 0, sizeof(ui) * (candidates_count[0] + 1));
    // TM has size TM[qsiz][candi[i]]+1 element.
    // TM initialized the first row with zeros?
    std::vector<std::vector<ui>> vec_index(query_vertices_num);
    for (ui i = 0; i < query_vertices_num; i++)
    {
        vec_index[i].resize(candidates_count[i]);
        std::fill(vec_index[i].begin(), vec_index[i].end(), (ui)-1);
    }
    // vec_index is initialized with largest number ui-1
    std::vector<std::vector<ui>> vec_set;
    computeNEC(query_graph, edge_matrix, candidates_count, candidates, vec_index, vec_set);

    std::vector<std::vector<ui>> pi_m_index(query_vertices_num);
    for (ui i = 0; i < query_vertices_num; i++)
    {
        pi_m_index[i].resize(candidates_count[i], (ui)-1);
    }
    std::vector<std::vector<ui>> pi_m;
    ui *pi_m_count = new ui[max_depth];
    pi_m_count[0] = 0;

    std::vector<std::vector<ui>> dm(query_vertices_num);
    std::unordered_map<VertexID, VertexID> reverse_embedding;
    reverse_embedding.reserve(query_vertices_num);
#endif

    ui start_vertex = (ui)-1;
    for (ui i = 0; i < query_vertices_num; i++)
    {
        if (extendable[i] == 0)
        {
            start_vertex = i;
            break;
        }
    }
    assert(start_vertex != (ui)-1);

    order[depth] = start_vertex;
    visited_u[start_vertex] = true;

    idx[depth] = 0;
    valid_cans_count[start_vertex] = candidates_count[start_vertex];
    idx_count[depth] = valid_cans_count[start_vertex];
    for (ui i = 0; i < candidates_count[start_vertex]; i++)
    {
        valid_cans[start_vertex][i] = candidates[start_vertex][i];
    }
#ifdef ENABLE_EQUIVALENT_SET
    std::fill(pi_m_index[start_vertex].begin(), pi_m_index[start_vertex].end(), (ui)-1);

#endif
    bool Match_BA[max_depth] = {false};
    unordered_set<ui> EmbSum;

    while (true)
    {
        while (idx[depth] < idx_count[depth])
        {
            Match_BA[depth] = false;
            auto end = std::chrono::high_resolution_clock::now();
            ens = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
            if (ens > TimeL)
            { // 1000 1 sec
                goto EXIT;
            }
#ifdef ENABLE_DYNAMIC_CANS

            // if (idx[order[depth]] != valid_cans[order[depth]].begin())
            if (idx[order[depth]] != valid_cans[order[depth]][0])
            {
                RestoreValidCans(query_graph, data_graph, visited_u, order[depth], embedding[order[depth]], valid_cans);
            }
#endif
            VertexID u = order[depth];
            VertexID v = valid_cans[u][idx[depth]];
            // ui v_index1=0;
            //   for (; v_index1 < candidates_count[u]; v_index1++)
            //   {
            //       if (candidates[u][v_index1] == v)
            //           break;
            //   }
            //  if (vec_index[u][v_index1] == (ui)-1){
            //      computeNEC_Adaptive(query_graph, edge_matrix, candidates_count, candidates, vec_index, vec_set,u,v_index1,vec_count);
            //     vec_count ++;
            // }
            bool is_voilate_symmetry = false;

            if (!full_constraints.empty())
            {
                is_voilate_symmetry = checkSingeVertexBySymmetryBreaking(depth, embedding, order, visited_u, v,
                                                                         full_constraints,
                                                                         violated_symmetry_u,
                                                                         max_depth);
            }

            if (visited_vertices[v] || is_voilate_symmetry)
            {

#ifdef ENABLE_EQUIVALENT_SET
                TM[u][idx[depth]] = 0;
                // keep only common in conflix- like pmnegative
                if (visited_vertices[v])
                {
                    VertexID con_u = reverse_embedding[v];

                    ui con_v_index = 0, v_index = 0;
                    for (; con_v_index < valid_cans_count[con_u]; con_v_index++)
                    {
                        if (valid_cans[con_u][con_v_index] == v)
                            break;
                    }
                    for (; v_index < candidates_count[u]; v_index++)
                    {
                        if (candidates[u][v_index] == v)
                            break;
                    }

                    assert(con_v_index < valid_cans_count[con_u]);
                    assert(v_index < candidates_count[u]);

                    auto &con_uv_idx = pi_m_index[con_u][con_v_index];

                    for (ui i = 0; i < pi_m[con_uv_idx].size(); i++)
                    {
                        bool f_in = false;
                        // if ((vec_index[u][v_index])== (ui)-1){
                        //         computeNEC_Adaptive(query_graph, edge_matrix, candidates_count, candidates, vec_index, vec_set,u,v_index,vec_count);
                        //         vec_count++;
                        //        }
                        for (ui j = 0; j < vec_set[vec_index[u][v_index]].size(); j++)
                        {
                            if (vec_set[vec_index[u][v_index]][j] == pi_m[con_uv_idx][i])
                            {
                                f_in = true;
                                break;
                            }
                        }
                        if (f_in == false)
                        {
                            pi_m[con_uv_idx][i] = pi_m[con_uv_idx][pi_m[con_uv_idx].size() - 1];
                            pi_m[con_uv_idx].pop_back();
                            i--;
                        }
                    }
                }
                if (is_voilate_symmetry)
                {
                    for (ui i = 0; i < depth; ++i)
                    {

                        if (violated_symmetry_u[order[i]])
                        {
                            VertexID sym_u = order[i];
                            VertexID sym_v = embedding[sym_u];
                            ui sym_v_index = 0, v_index = 0;
                            for (; sym_v_index < valid_cans_count[sym_u]; sym_v_index++)
                            {
                                if (valid_cans[sym_u][sym_v_index] == sym_v)
                                    break;
                            }
                            for (; v_index < candidates_count[u]; v_index++)
                            {
                                if (candidates[u][v_index] == sym_v)
                                    break;
                            }

                            auto &sym_uv_idx = pi_m_index[sym_u][sym_v_index];

                            for (ui i = 0; i < pi_m[sym_uv_idx].size(); i++)
                            {
                                bool f_in = false;
                                // if ((vec_index[u][v_index])== (ui)-1){
                                //     computeNEC_Adaptive(query_graph, edge_matrix, candidates_count, candidates, vec_index, vec_set,u,v_index,vec_count);
                                // vec_count++;
                                // }
                                for (ui j = 0; j < vec_set[vec_index[u][v_index]].size(); j++)
                                {
                                    if (vec_set[vec_index[u][v_index]][j] == pi_m[sym_uv_idx][i])
                                    {
                                        f_in = true;
                                        break;
                                    }
                                }
                                if (f_in == false)
                                {
                                    pi_m[sym_uv_idx][i] = pi_m[sym_uv_idx][pi_m[sym_uv_idx].size() - 1];
                                    pi_m[sym_uv_idx].pop_back();
                                    i--;
                                }
                            }
                        }
                    }
                }

#endif

                idx[depth]++;
                continue;
            }

#ifdef ENABLE_EQUIVALENT_SET

            if (pi_m_index[u][idx[depth]] != (ui)-1)
            {
                VertexID equ_v = pi_m[pi_m_index[u][idx[depth]]][0];
                ui equ_v_index = 0;

                for (; equ_v_index < idx_count[depth]; equ_v_index++)
                {
                    if (equ_v == valid_cans[u][equ_v_index])
                        break;
                }
                assert(equ_v_index < idx_count[depth]);
                embedding_cnt += TM[u][equ_v_index];
                // EmbSum.insert(equ_v);

                TM[u][candidates_count[u]] += TM[u][equ_v_index];
                if (TM[u][equ_v_index] > 0)
                    EmbSum.insert(v);
                // for (auto eq_node : pi_m[pi_m_index[u][idx[depth]]][]) {
                //  EmbSum.insert(eq_node);

                //    }
                idx[depth]++;
                continue;
            }
            reverse_embedding[v] = u;
#endif
            embedding[u] = v;
            visited_vertices[v] = true;
            ui cur_idx = idx[depth]++;
#ifdef ENABLE_EQUIVALENT_SET
            pi_m_index[u][cur_idx] = pi_m_count[depth]++;
            ui v_index = 0;
            for (; v_index < candidates_count[u]; v_index++)
            {
                if (candidates[u][v_index] == v)
                    break;
            }
            // if ((vec_index[u][v_index])== (ui)-1){
            //     computeNEC_Adaptive(query_graph, edge_matrix, candidates_count, candidates, vec_index, vec_set,u,v_index,vec_count);
            //    vec_count++;
            //    }
            auto &pi = vec_set[vec_index[u][v_index]];
            pi_m.push_back(pi);
            assert(pi_m.size() == pi_m_count[depth]);
            dm[u].clear();
            for (ui i = 0; i < depth; i++)
            {
                bool va_in_pi = false, sec_empty = true;
                ui va_index = 0;
                VertexID ua = order[i];
                VertexID va = embedding[ua];
                for (; va_index < candidates_count[ua]; va_index++)
                {
                    if (candidates[ua][va_index] == va)
                        break;
                }
                assert(va_index < candidates_count[ua]);
                // if ((vec_index[u][v_index])== (ui)-1){
                //    computeNEC_Adaptive(query_graph, edge_matrix, candidates_count, candidates, vec_index, vec_set,ua,va_index,vec_count);
                //    vec_count++;
                //    }
                auto &pi_a = vec_set[vec_index[ua][va_index]];
                for (ui j = 0; j < pi.size(); j++)
                {
                    if (pi[j] == embedding[order[i]])
                    {
                        va_in_pi = true;
                        break;
                    }
                    for (ui k = 0; k < pi_a.size(); k++)
                    {
                        if (pi_a[k] == pi[j])
                        {
                            sec_empty = false;
                            break;
                        }
                    }
                    if (sec_empty == false)
                        break;
                }
                if (va_in_pi == false && sec_empty == false)
                {
                    ui dm_size = dm[ua].size();
                    for (ui j = 0; j < pi.size(); j++)
                    {
                        bool f_add = true;
                        for (ui k = 0; k < dm_size; k++)
                        {
                            if (dm[ua][k] == pi[j])
                            {
                                f_add = false;
                                break;
                            }
                        }
                        if (f_add == true)
                        {
                            dm[ua].push_back(pi[j]);
                        }
                    }
                }
            }
#endif
            if (depth == max_depth - 1)
            {
                embedding_cnt += 1;
                visited_vertices[v] = false;
                visited_u[u] = false;
                // EmbSum.insert(embedding[order[depth]]);
                int ao = depth; //-1;
                ui vqo = order[ao];
                int UNPMtemp = 0;
                int UNPM2 = 0;
                UNPM = 0;
                while (ao >= 0)
                {
                    UNPMtemp = EmbSum.insert(embedding[vqo]).second;
// UNPMtemp=EmbSum.insert(v).second;
// Match_BA[ao] = true;
// pi_m[pi_m_index[u][idx[depth]]]
#ifdef ENABLE_EQUIVALENT_SET
                    if (pi_m_index[vqo][idx[ao] - 1] != (ui)-1)
                    {
                        for (auto eq_node : pi_m[pi_m_index[vqo][idx[ao] - 1]])
                        {
                            UNPMtemp += EmbSum.insert(eq_node).second;
                        }
                    }
#endif
                    if (UNPMtemp > 1)
                    {
                        UNPM2 += (UNPMtemp - 1);
                        UNPMtemp = 1;
                    }
                    UNPM += UNPMtemp;
                    ao--;
                    vqo = order[ao];
                }
#ifdef TOPKGREEDY
                if (UNPM > pq.top())
                {
                    if (heapcount >= ksize)
                    {
                        pq.pop();
                    }
                    pq.push(UNPM);
                    heapcount++;
                }
                while (UNPM2 > 0 && pq.top() < 1)
                {
                    if (heapcount >= ksize)
                    {
                        pq.pop();
                    }
                    UNPM2--;
                    pq.push(1);
                    heapcount++;
                }
#endif

                if (embedding_cnt >= output_limit_num)
                {
                    goto EXIT;
                }
#ifdef ENABLE_EQUIVALENT_SET
                // for (auto eq_node : pi_m[pi_m_index[u][idx[depth]]]) {
                EmbSum.insert(v);
                reverse_embedding.erase(v);
                TM[u][cur_idx] = 1;
                TM[u][candidates_count[u]]++;

                auto &uv_idx = pi_m_index[u][cur_idx];
                for (ui i = 0; i < pi_m[uv_idx].size(); i++)
                {
                    bool f_del = false;
                    for (ui j = 0; j < dm[u].size(); j++)
                    {
                        if (pi_m[uv_idx][i] == dm[u][j])
                        {
                            f_del = true;
                            break;
                        }
                    }
                    if (f_del == true)
                    {
                        pi_m[uv_idx][i] = pi_m[uv_idx][pi_m[uv_idx].size() - 1];
                        pi_m[uv_idx].pop_back();
                        i--;
                    }
                }
                for (ui i = 1; i < pi_m[uv_idx].size(); i++)
                {
                    ui v_equ_index = 0;
                    for (; v_equ_index < valid_cans_count[u]; v_equ_index++)
                    {
                        if (valid_cans[u][v_equ_index] == pi_m[uv_idx][i])
                            break;
                    }
                    assert(v_equ_index < valid_cans_count[u]);
                    pi_m_index[u][v_equ_index] = uv_idx;
                }

                if (pi_m[uv_idx].size() > 1)
                {
                    auto pi_v_idx = std::find(pi_m[uv_idx].begin(), pi_m[uv_idx].end(), v);
                    assert(pi_v_idx != pi_m[uv_idx].end());
                    ui tmp = *pi_v_idx;
                    *pi_v_idx = pi_m[uv_idx][0];
                    pi_m[uv_idx][0] = tmp;
                }
#endif
            }
            else
            {
                call_count += 1;
                depth += 1;

                order[depth] = generateNextU(data_graph, query_graph, candidates, candidates_count, valid_cans,
                                             valid_cans_count, extendable, nec, depth, embedding,
                                             edge_matrix, visited_vertices, visited_u, order, tree);

                if (order[depth] == (ui)-1)
                {
                    break;
                }
                else
                {
                    visited_u[order[depth]] = true;
                    idx[depth] = 0;
                    idx_count[depth] = valid_cans_count[order[depth]];
                }
#ifdef ENABLE_EQUIVALENT_SET
                memset(TM[order[depth]], 0, sizeof(ui) * (candidates_count[order[depth]] + 1));
                std::fill(pi_m_index[order[depth]].begin(), pi_m_index[order[depth]].end(), (ui)-1);
                pi_m_count[depth] = pi_m_count[depth - 1];
#endif
            }
        }

        depth -= 1;

        if (depth < 0)
            break;
        VertexID u = order[depth];
        ui cur_idx = idx[depth] - 1;
        visited_vertices[embedding[u]] = false;
        restoreExtendableVertex(tree, u, extendable);
        if (order[depth + 1] != (ui)-1)
        {
            VertexID last_u = order[depth + 1];
            visited_u[last_u] = false;
            if (nec[last_u] != NULL)
                (*(nec[last_u]))++;
#ifdef ENABLE_DYNAMIC_CANS
            RestoreValidCans(query_graph, data_graph, visited_u, last_u, last_v, valid_cans);
#endif
        }
#ifdef ENABLE_EQUIVALENT_SET
        if (order[depth + 1] != (ui)-1)
        {
            TM[u][cur_idx] = TM[order[depth + 1]][candidates_count[order[depth + 1]]];
        }
        else
        {
            TM[u][cur_idx] = 0;
        }
        TM[u][candidates_count[u]] += TM[u][cur_idx];
        // add here
        auto &uv_idx = pi_m_index[u][cur_idx];
        if (TM[u][cur_idx] != 0)
        {

            for (ui i = 0; i < pi_m[uv_idx].size(); i++)
            {
                bool f_del = false;
                for (ui j = 0; j < dm[u].size(); j++)
                {
                    if (pi_m[uv_idx][i] == dm[u][j])
                    {
                        f_del = true;
                        break;
                    }
                }
                if (f_del == true)
                {
                    pi_m[uv_idx][i] = pi_m[uv_idx][pi_m[uv_idx].size() - 1];
                    pi_m[uv_idx].pop_back();
                    i--;
                }
            }

            for (ui i = 0; i < pi_m[uv_idx].size(); i++)
            {
            }

            if (pi_m[uv_idx].size() > 1)
            {
                auto pi_v_idx = std::find(pi_m[uv_idx].begin(), pi_m[uv_idx].end(), embedding[u]);
                assert(pi_v_idx != pi_m[uv_idx].end());
                ui tmp = *pi_v_idx;
                *pi_v_idx = pi_m[uv_idx][0];
                pi_m[uv_idx][0] = tmp;
            }
        }
        for (ui i = 1; i < pi_m[uv_idx].size(); i++)
        {
            ui v_equ_index = 0;
            for (; v_equ_index < valid_cans_count[u]; v_equ_index++)
            {
                if (valid_cans[u][v_equ_index] == pi_m[uv_idx][i])
                    break;
            }
            assert(v_equ_index < valid_cans_count[u]);

            pi_m_index[u][v_equ_index] = uv_idx;
        }
        pi_m.resize(pi_m_count[depth]);
#endif
    }

EXIT:

    for (ui i = 0; i < query_vertices_num; i++)
    {
    }
    delete[] nec;
    delete[] embedding;
    delete[] visited_u;
    delete[] visited_vertices;
    delete[] order;
    delete[] extendable;
    delete[] idx;
    delete[] idx_count;
    for (ui i = 0; i < query_vertices_num; i++)
    {
        delete[] valid_cans[i];
    }
    delete[] valid_cans;
    delete[] valid_cans_count;
#ifdef ENABLE_EQUIVALENT_SET
    delete[] TM;
#endif
    s.embedding_cnt = embedding_cnt;
    s.Can_embed = EmbSum.size();
#ifdef TOPKGREEDY
    std::priority_queue<int> maxHeap;
    int topkcounter = 0;
    int greedysum = 0;
    while (!pq.empty())
    {
        maxHeap.push(pq.top());
        pq.pop();
    }
    cout << maxHeap.top() << endl;
    while (!maxHeap.empty())
    {
        greedysum += maxHeap.top();
        maxHeap.pop();
        topkcounter++;
        switch (topkcounter)
        {
        case VALUE_5:
            s.topk[0] = greedysum;
            break;
        case VALUE_10:
            s.topk[1] = greedysum;
            break;
        case VALUE_50:
            s.topk[2] = greedysum;
            break;
        case VALUE_100:
            s.topk[3] = greedysum;
            break;
        case VALUE_250:
            s.topk[4] = greedysum;
            break;
        case VALUE_500:
            s.topk[5] = greedysum;
            break;
        case VALUE_750:
            s.topk[6] = greedysum;
            break;
        case VALUE_1000:
            s.topk[7] = greedysum;
            break;
        default:
            break;
        }
    }
    for (int i = 0; i < 8; i++)
    {
        if (s.topk[i] == 0)
            s.topk[i] = greedysum;
    }
#endif
    return s;
}

enumResult
EvaluateQuery::LFTJDIV1(const Graph *data_graph, const Graph *query_graph, Edges ***edge_matrix, ui **candidates,
                        ui *candidates_count, ui *order, size_t output_limit_num, size_t &call_count, int TimeL, int FairT, const std::unordered_map<VertexID, std::pair<std::set<VertexID>, std::set<VertexID>>> &ordered_constraints)
{
    auto start = std::chrono::high_resolution_clock::now();
    enumResult s;
    // Generate bn.

    ui **bn;
    ui *bn_count;
    generateBN(query_graph, order, bn, bn_count);
    unordered_set<ui> EmbSum;
    bool NotifyALL = true;
    // Allocate the memory buffer.
    ui *idx;
    ui *idx_count;
    ui *embedding;
    ui *idx_embedding;
    ui *temp_buffer;
    ui **valid_candidate_idx;
    bool *visited_vertices;
    int max_depth = query_graph->getVerticesCount();
    priority_queue<int, vector<int>, greater<int>> pq;
    pq.push(0);
    int heapcount = 1;
    int ksize = 1000;
    bool Match_BA_t[max_depth] = {false};
    unordered_set<ui> EmbSum_t;
    int UNPM = 0;
    allocateBuffer(data_graph, query_graph, candidates_count, idx, idx_count, embedding, idx_embedding,
                   temp_buffer, valid_candidate_idx, visited_vertices);
    vector<ui> VN[query_graph->getVerticesCount()];
    vector<ui> VNTemp;
    int qsiz = query_graph->getVerticesCount();

    int RQ[query_graph->getVerticesCount()];
    for (int i = 0; i < query_graph->getVerticesCount(); i++)
    {
        RQ[order[i]] = i;
    }
    bool *embd_vertices = new bool[data_graph->getVerticesCount()];
    std::fill(embd_vertices, embd_vertices + data_graph->getVerticesCount(), false);
    size_t embedding_cnt = 0;
    int cur_depth = 0;

    VertexID start_vertex = order[0];
    bool Match_BA[max_depth] = {false};

    idx[cur_depth] = 0;                                    // set the starting position to 0 for the current depth
    idx_count[cur_depth] = candidates_count[start_vertex]; // set the size for current depth
    // all the nodes are possible candidates

    int count2 = 1;
    size_t **candidatesHC = NULL;
    unordered_map<size_t, vector<ui>> *idToValues;
    idToValues = new unordered_map<size_t, vector<ui>>[query_graph->getVerticesCount()];
    candidatesHC = new size_t *[qsiz];

    for (ui i = 0; i < qsiz; ++i)
    {
        candidatesHC[i] = new size_t[candidates_count[i]];
        memset(candidatesHC[i], 0, sizeof(size_t) * candidates_count[i]);
    }
    calculateCell(query_graph, edge_matrix, candidates_count, candidatesHC, idToValues, count2, start_vertex, 0);
    count2++;
    VN[start_vertex] = idToValues[start_vertex][candidatesHC[start_vertex][0]];

    for (ui i = 0; i < idx_count[cur_depth]; ++i)
    {
        valid_candidate_idx[cur_depth][i] = i;
    }

    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> ancestors;
    computeAncestor(query_graph, bn, bn_count, order, ancestors);
    AncestorWithSymmetryBreaking(query_graph, ancestors, ordered_constraints);
    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> vec_failing_set(max_depth);
    std::unordered_map<VertexID, VertexID> reverse_embedding;
    reverse_embedding.reserve(MAXIMUM_QUERY_GRAPH_SIZE * 2);

    int k;

    VertexID utemp; // u is current depth qid
    VertexID vtemp;
    int counter = 0;
    double ens = 0;
    int ordCand[idx_count[0]] = {100000};
    TimeL = TimeL * 1000;
    while (true)
    {
        // current position smaller than size
        while (idx[cur_depth] < idx_count[cur_depth])
        {
            Match_BA[cur_depth] = false;
            Match_BA_t[cur_depth] = false;
            auto end = std::chrono::high_resolution_clock::now();
            ens = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
            if (ens > TimeL)
            { // 1000 1 sec
                goto EXIT;
            }
            ui valid_idx = valid_candidate_idx[cur_depth][idx[cur_depth]]; // valid id
            if (valid_idx == 10000000)
            {
                idx[cur_depth] += 1;
                continue;
            }
            
            if (cur_depth == 0 && idx[cur_depth] != 0 && EmbSum.size() > 0 && FairT > 0)
            {

                int tempC = 0;
                for (auto it = EmbSum_t.begin(); it != EmbSum_t.end(); ++it)
                {
                    tempC += EmbSum.insert(*it).second; // Inserting elements into set2
                }
                EmbSum_t.clear();
                while (tempC > 0 && pq.top() < 1)
                {
                    if (heapcount >= ksize)
                    {
                        pq.pop();
                    }
                    pq.push(1);
                    heapcount++;
                    tempC--;
                }
/*
                if (FairT == 1)
                    rankSimple(order[0], nodeId, candidates, valid_candidate_idx, idx, valid_idx, idx_count);
                else if (FairT == 2) // already have experiments on that
                    rankLess(query_graph, edge_matrix, order[0], nodeId, candidates, valid_candidate_idx, idx, valid_idx, idx_count);
                else if (FairT == 3)
                    rankMore(query_graph, edge_matrix, order[0], nodeId, candidates, valid_candidate_idx, idx, valid_idx, idx_count);
                else if (FairT == 4)
                    rankFinal(query_graph, edge_matrix, order[0], nodeId, candidates, valid_candidate_idx, idx, valid_idx, idx_count);
                else if (FairT == 5)
                    rankFinal1(query_graph, edge_matrix, order[0], nodeId, candidates, valid_candidate_idx, idx, valid_idx, idx_count);
                valid_idx = valid_candidate_idx[cur_depth][idx[cur_depth]];
            */}
// from all the valid IDS
            VertexID u = order[cur_depth]; // u is current depth qid

            VertexID v = candidates[u][valid_idx]; // v is FCS[qid][vid] valid_idx-> position in CS

            size_t candidate = candidatesHC[u][valid_idx];
            if (candidate == 0 && cur_depth == 0)
            {
                calculateCell(query_graph, edge_matrix, candidates_count, candidatesHC, idToValues, count2, u, valid_idx);
                count2++;
            }

            else if (candidate == 0 && cur_depth != max_depth - 1)
            { // calculateCell(query_graph, edge_matrix, candidates_count, candidatesHC, idToValues, count2, u, valid_idx);

                calculateCellFN(query_graph, edge_matrix, candidates_count, candidatesHC, idToValues, count2, u, valid_idx, RQ);
                count2++;
            }
            candidate = candidatesHC[u][valid_idx];
            VN[u].clear();
            if (cur_depth != max_depth - 1) // change here
            {
                VN[u] = idToValues[u][candidatesHC[u][valid_idx]];
                if (cur_depth != 0)
                {
                    int ia = idx[cur_depth];
                    int io = 0;
                    VNTemp.clear();
                    while (ia < idx_count[cur_depth] && io < VN[u].size())
                    {
                        ui valid_idxN = valid_candidate_idx[cur_depth][ia];
                        if (valid_idxN == 10000000)
                        {
                            ia += 1;
                            continue;
                        }
                        if (VN[u][io] == valid_idxN)
                        {
                            VNTemp.push_back(VN[u][io]);
                            ia++;
                            io++;
                        }

                        else if (VN[u][io] > valid_idxN)
                            ia++;
                        else
                            io++;
                    }
                    VN[u] = VNTemp;
                }
            }
            // Match_BA[u] = false; // set embedings found to empty for start

            if (visited_vertices[v])
            {
#ifdef ENABLE_EQUIVALENT_SET
                VNTemp.clear();
                ui utemp = reverse_embedding[v];

                if (cur_depth != max_depth - 1)
                {
                    if (VN[u].size() > 0 && VN[utemp].size() > 0)
                    {
                        int ia = 0;
                        int io = 0;

                        while (ia < VN[u].size() && io < VN[utemp].size())
                        {
                            if (candidates[u][VN[u][ia]] == candidates[utemp][VN[utemp][io]])
                            {
                                VNTemp.push_back(VN[utemp][io]);
                                ia++;
                                io++;
                            }
                            else if (candidates[u][VN[u][ia]] < candidates[utemp][VN[utemp][io]])
                                ia++;
                            else
                                io++;
                        }
                    }
                    VN[utemp] = VNTemp;
                    VNTemp.clear();
                    // break;
                }
                else
                { // enumerate the results if it has equivalent nodes.
                    if (VN[utemp].size() > 1)
                    {
                        VNTemp.clear();
                        for (int ai = 0; ai < VN[utemp].size(); ai++)
                        {
                            if (!visited_vertices[candidates[utemp][VN[utemp][ai]]])
                                VNTemp.push_back(VN[utemp][ai]);
                        }
                        if (VNTemp.size() >= 1)
                        {
                            // heuristic there is a match so we count.
                            //  we can also exclude by just using VTEMP
                            // for utemp to count the correct number of matchings

                            int ao = cur_depth; //-1;
                            ui vqo = order[ao];
                            int UNPMtemp = 0;
                            int UNPM2 = 0;
                            UNPM = 0;
                            while (ao >= 0)
                            {

                                // EmbSum.insert(embedding[vqo]);
                                UNPMtemp = EmbSum.insert(embedding[vqo]).second;

                                if (Match_BA_t[ao] == false)
                                {
                                    if (VN[vqo].size() > 0)
                                    {
                                        for (int dd = 0; dd < VN[vqo].size(); dd++)
                                        {
                                            if (VN[vqo][dd] != 10000000)
                                                EmbSum_t.insert(candidates[vqo][VN[vqo][dd]]);
                                        }
                                    }
                                    Match_BA_t[ao] = true;
                                }
                                if (UNPMtemp == 0)
                                {
                                    Match_BA[ao] = true;
                                    if (VN[vqo].size() > 0)
                                    {
                                        for (int dd = 0; dd < VN[vqo].size(); dd++)
                                        {
                                            if (VN[vqo][dd] != 10000000)
                                                UNPMtemp += EmbSum.insert(candidates[vqo][VN[vqo][dd]]).second;
                                            if (UNPMtemp == 1)
                                                break;
                                        }
                                    }
                                    if (UNPMtemp == 1)
                                        Match_BA[ao] = false;
                                }

                                ao--;
                                vqo = order[ao];
                                UNPM += UNPMtemp;
                            }
                            if (UNPM > pq.top())
                            {
                                if (heapcount >= ksize)
                                {
                                    pq.pop();
                                }
                                pq.push(UNPM);
                                heapcount++;
                            }
#ifdef ENABLE_FAILING_SET
                            // reverse_embedding.erase(embedding[u]);
                            vec_failing_set[cur_depth].set();
                            vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
#endif
                        }
                        else
                        { // update equivalent class of VN[utemp]
                            int ia = 0;
                            int io = 0;
                            VNTemp.clear();
                            while (ia < idx_count[cur_depth] && io < VN[utemp].size())
                            {

                                ui valid_idxN = valid_candidate_idx[cur_depth][ia];
                                if (valid_idxN == 10000000)
                                {
                                    ia += 1;
                                    continue;
                                }
                                if (VN[utemp][io] == valid_idxN)
                                {
                                    // VNTemp.push_back(VN[u][io]);
                                    VNTemp.push_back(VN[utemp][io]);
                                    ia++;
                                    io++;
                                }

                                else if (VN[utemp][io] > valid_idxN)
                                    ia++;
                                else
                                    io++;
                                // ia++;
                            }
                            VN[utemp] = VNTemp;
                            VNTemp.clear();
                        }
                    }
                }
                // }
                //}
                VN[u].clear();
#endif
                idx[cur_depth] += 1;
#ifdef ENABLE_FAILING_SET
                vec_failing_set[cur_depth] = ancestors[u];
                vec_failing_set[cur_depth] |= ancestors[reverse_embedding[v]];
                vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
#endif
                continue;
            }
            embd_vertices[v] = true;
            counter++;
            embedding[u] = v; // set the embeding of depth qid to v
            // valid_idx = valid_candidate_idx[cur_depth][idx[cur_depth]];
            idx_embedding[u] = valid_idx; // reverse the id of the embdeing to the position
            visited_vertices[v] = true;
            idx[cur_depth] += 1; // next element
            while (cur_depth < idx_count[cur_depth] && valid_candidate_idx[cur_depth][idx[cur_depth]] == 10000000)
            {
                idx[cur_depth] += 1;
            }
            reverse_embedding[v] = u;

            if (cur_depth == max_depth - 1)
            {
                embedding_cnt += 1;
                visited_vertices[v] = false;
                embd_vertices[v] = false;
                counter--;
                // EmbSum.insert(embedding[order[cur_depth]]);

                int ao = cur_depth; //-1;
                ui vqo = order[ao];
                int UNPMtemp = 0;
                int UNPM2 = 0;
                UNPM = 0;
                while (ao >= 0)
                {

                    // EmbSum.insert(embedding[vqo]);
                    UNPMtemp = EmbSum.insert(embedding[vqo]).second;

                    if (Match_BA_t[ao] == false)
                    {
                        if (VN[vqo].size() > 0)
                        {
                            for (int dd = 0; dd < VN[vqo].size(); dd++)
                            {
                                if (VN[vqo][dd] != 10000000)
                                    EmbSum_t.insert(candidates[vqo][VN[vqo][dd]]);
                            }
                        }
                        Match_BA_t[ao] = true;
                    }
                    if (UNPMtemp == 0)
                    {
                        Match_BA[ao] = true;
                        if (VN[vqo].size() > 0)
                        {
                            for (int dd = 0; dd < VN[vqo].size(); dd++)
                            {
                                if (VN[vqo][dd] != 10000000)
                                    UNPMtemp += EmbSum.insert(candidates[vqo][VN[vqo][dd]]).second;
                                if (UNPMtemp == 1)
                                    break;
                            }
                        }
                        if (UNPMtemp == 1)
                            Match_BA[ao] = false;
                    }

                    ao--;
                    vqo = order[ao];
                    UNPM += UNPMtemp;
                }
                if (UNPM > pq.top())
                {
                    if (heapcount >= ksize)
                    {
                        pq.pop();
                    }
                    pq.push(UNPM);
                    heapcount++;
                }

                reverse_embedding.erase(embedding[u]);
                vec_failing_set[cur_depth].set();
                vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];

                if (embedding_cnt >= output_limit_num)
                {
                    goto EXIT;
                }
            }
            else
            {
                call_count += 1;
                cur_depth += 1;     // next depth
                idx[cur_depth] = 0; // set the element to 0
                generateValidCandidateIndex(cur_depth, idx_embedding, idx_count, valid_candidate_idx, edge_matrix, bn,
                                            bn_count, order, temp_buffer, candidates, candidates_count, query_graph);
                //    generateValidCandidateIndex(cur_depth, idx_embedding, idx_count, valid_candidate_idx, edge_matrix, bn,
                //                            bn_count, order, temp_buffer);
                pruneCandidatesIndexBySymmetryBreaking(cur_depth, embedding, order, candidates_count, candidates,
                                                       idx_count, valid_candidate_idx, ordered_constraints);

                // get the candidates for next depth
                if (idx_count[cur_depth] == 0)
                {
                    vec_failing_set[cur_depth - 1] = ancestors[order[cur_depth]];
                }
                else
                {
                    vec_failing_set[cur_depth - 1].reset();
                }
            }
        }

        cur_depth -= 1;
        int ao = cur_depth;

        if (cur_depth == max_depth - 1)
        {
            continue;
        }
        if (cur_depth < 0)
            break;
        else
        {
            ui vq = order[cur_depth];
        // if (Match_BA[vq] == false && idx[cur_depth] < idx_count[cur_depth])
        //{ positive and negative removal now. dont need to count positive symmetric
        if (VN[vq].size() > 1)
        { // something to remove
            for (int dd = 0; dd < VN[vq].size(); dd++)
            {
                for (int kk = (idx[cur_depth]); kk < idx_count[cur_depth]; kk++)
                {
                    ui valid_idxn = valid_candidate_idx[cur_depth][kk];
                    if (valid_idxn == 10000000)
                        continue;
                    if (valid_idxn == VN[vq][dd])
                    {
                        valid_candidate_idx[cur_depth][kk] = 10000000;
                        break;
                    }
                }
            }
        }VN[vq].clear();


            VertexID u = order[cur_depth];
            reverse_embedding.erase(embedding[u]);
            if (cur_depth != 0)
            {
                if (!vec_failing_set[cur_depth].test(u))
                {
                    vec_failing_set[cur_depth - 1] = vec_failing_set[cur_depth];
                    idx[cur_depth] = idx_count[cur_depth];
                }
                else
                {
                    vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
                }
            }
            visited_vertices[embedding[u]] = false;
            embd_vertices[embedding[u]] = false;
            counter--;
        }

    }

EXIT:
    releaseBuffer(max_depth, idx, idx_count, embedding, idx_embedding, temp_buffer, valid_candidate_idx,
                  visited_vertices,
                  bn, bn_count);
    int tempC = 0;
    for (auto it = EmbSum_t.begin(); it != EmbSum_t.end(); ++it)
    {
        tempC += EmbSum.insert(*it).second; // Inserting elements into set2
    }
    EmbSum_t.clear();
    while (tempC > 0 && pq.top() < 1)
    {
        if (heapcount >= ksize)
        {
            pq.pop();
        }
        pq.push(1);
        heapcount++;
        tempC--;
    }

#ifdef TOPKGREEDY
    std::priority_queue<int> maxHeap;
    int topkcounter = 0;
    int greedysum = 0;
    while (!pq.empty())
    {
        maxHeap.push(pq.top());
        pq.pop();
    }
    while (!maxHeap.empty())
    {
        greedysum += maxHeap.top();
        maxHeap.pop();
        topkcounter++;
        switch (topkcounter)
        {
        case VALUE_5:
            s.topk[0] = greedysum;
            break;
        case VALUE_10:
            s.topk[1] = greedysum;
            break;
        case VALUE_50:
            s.topk[2] = greedysum;
            break;
        case VALUE_100:
            s.topk[3] = greedysum;
            break;
        case VALUE_250:
            s.topk[4] = greedysum;
            break;
        case VALUE_500:
            s.topk[5] = greedysum;
            break;
        case VALUE_750:
            s.topk[6] = greedysum;
            break;
        case VALUE_1000:
            s.topk[7] = greedysum;
            break;
        default:
            break;
        }
    }
    for (int i = 0; i < 8; i++)
    {
        if (s.topk[i] == 0)
            s.topk[i] = greedysum;
    }
#endif
    int true_cand_sum = 0;
    s.Can_embed = EmbSum.size();
    EmbSum.clear();
    s.embedding_cnt = embedding_cnt;
    s.candidate_true_count_sum = true_cand_sum;
    // s.topk=greedysum;
    return s;
}

enumResult
EvaluateQuery::LFTJDIV(const Graph *data_graph, const Graph *query_graph, Edges ***edge_matrix, ui **candidates,
                       ui *candidates_count, ui *order, size_t output_limit_num, size_t &call_count, int TimeL, int FairT, const std::unordered_map<VertexID, std::pair<std::set<VertexID>, std::set<VertexID>>> &ordered_constraints)
{
    //assumming that creating matrix can be done before we do not count the time.

    auto start = std::chrono::high_resolution_clock::now();
    enumResult s;
    // Generate bn.

    ui **bn;
    ui *bn_count;
    generateBN(query_graph, order, bn, bn_count);
    unordered_set<ui> EmbSum;
    bool NotifyALL = true;
    // Allocate the memory buffer.
    ui *idx;
    ui *idx_count;
    ui *embedding;
    ui *idx_embedding;
    ui *temp_buffer;
    ui **valid_candidate_idx;
    bool *visited_vertices;
#ifdef TOPKGREEDY
    priority_queue<int, vector<int>, greater<int>> pq;
    pq.push(0);
    int heapcount = 1;
    int ksize = 1000;
#endif
    int UNPM = 0;
    allocateBuffer(data_graph, query_graph, candidates_count, idx, idx_count, embedding, idx_embedding,
                   temp_buffer, valid_candidate_idx, visited_vertices);
    vector<ui> VN[query_graph->getVerticesCount()];
    vector<ui> VNTemp;
    int qsiz = query_graph->getVerticesCount();

    int RQ[query_graph->getVerticesCount()];
    for (int i = 0; i < query_graph->getVerticesCount(); i++)
    {
        RQ[order[i]] = i; // depth RQ[u]=
    }
    bool *embd_vertices = new bool[data_graph->getVerticesCount()];
    std::fill(embd_vertices, embd_vertices + data_graph->getVerticesCount(), false);
    size_t embedding_cnt = 0;
    int cur_depth = 0;
    int max_depth = query_graph->getVerticesCount();
    VertexID start_vertex = order[0];
    bool Match_BA[max_depth] = {false};
    idx[cur_depth] = 0;                                    // set the starting position to 0 for the current depth
    idx_count[cur_depth] = candidates_count[start_vertex]; // set the size for current depth
    // all the nodes are possible candidates

    // end1 = std::chrono::high_resolution_clock::now();
    int count2 = 1;
#ifdef ENABLE_EQUIVALENT_SET
    size_t **candidatesHC = NULL;
    unordered_map<size_t, vector<ui>> *idToValues;
    idToValues = new unordered_map<size_t, vector<ui>>[query_graph->getVerticesCount()];
    candidatesHC = new size_t *[qsiz];

    for (ui i = 0; i < qsiz; ++i)
    {
        candidatesHC[i] = new size_t[candidates_count[i]];
        memset(candidatesHC[i], 0, sizeof(size_t) * candidates_count[i]);
    }
    calculateCell(query_graph, edge_matrix, candidates_count, candidatesHC, idToValues, count2, start_vertex, 0);
    // calculateCellFN(query_graph, edge_matrix, candidates_count, candidatesHC, idToValues, count2, start_vertex, 0,RQ);
    count2++;
    VN[start_vertex] = idToValues[start_vertex][candidatesHC[start_vertex][0]];
#endif

    for (ui i = 0; i < idx_count[cur_depth]; ++i)
    {
        valid_candidate_idx[cur_depth][i] = i;
    }

#ifdef ENABLE_FAILING_SET
    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> ancestors;
    computeAncestor(query_graph, bn, bn_count, order, ancestors);
    AncestorWithSymmetryBreaking(query_graph, ancestors, ordered_constraints);
    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> vec_failing_set(max_depth);
    std::unordered_map<VertexID, VertexID> reverse_embedding;
    reverse_embedding.reserve(MAXIMUM_QUERY_GRAPH_SIZE * 2);
    bool del = false;
#endif
    int k;
    VertexID utemp; // u is current depth qid
    VertexID vtemp;
    int counter = 0;
    double ens = 0;
    int ordCand[idx_count[0]] = {100000};
    TimeL = TimeL * 1000;
    while (true)
    {
        // current position smaller than size
        while (idx[cur_depth] < idx_count[cur_depth])
        {
            Match_BA[cur_depth] = false;
            auto end = std::chrono::high_resolution_clock::now();
            ens = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
            if (ens > TimeL)
            { // 1000 1 sec
                goto EXIT;
            }
            ui valid_idx = valid_candidate_idx[cur_depth][idx[cur_depth]]; // valid id
            if (valid_idx == 10000000)                                     //||
            {
                idx[cur_depth] += 1;
                continue;
            }
/*
#ifdef ENABLE_DIVERSITY
            if (cur_depth == 0 && idx[cur_depth] != 0 && EmbSum.size() > 0 && FairT > 0)
            {
                if (FairT == 1)
                    rankSimple(order[0], nodeId, candidates, valid_candidate_idx, idx, valid_idx, idx_count);
                else if (FairT == 2) // already have experiments on that
                    rankLess(query_graph, edge_matrix, order[0], nodeId, candidates, valid_candidate_idx, idx, valid_idx, idx_count);
                else if (FairT == 3)
                    rankMore(query_graph, edge_matrix, order[0], nodeId, candidates, valid_candidate_idx, idx, valid_idx, idx_count);
                else if (FairT == 4)
                    rankFinal(query_graph, edge_matrix, order[0], nodeId, candidates, valid_candidate_idx, idx, valid_idx, idx_count);
                else if (FairT == 5)
                    rankFinal1(query_graph, edge_matrix, order[0], nodeId, candidates, valid_candidate_idx, idx, valid_idx, idx_count);

                valid_idx = valid_candidate_idx[cur_depth][idx[cur_depth]];
            }*/
            // from all the valid IDS
//#endif
            VertexID u = order[cur_depth]; // u is current depth qid

            VertexID v = candidates[u][valid_idx]; // v is FCS[qid][vid] valid_idx-> position in CS

#ifdef ENABLE_EQUIVALENT_SET
            size_t candidate = candidatesHC[u][valid_idx];
            if (candidate == 0 && cur_depth == 0)
            {
                calculateCell(query_graph, edge_matrix, candidates_count, candidatesHC, idToValues, count2, u, valid_idx);
                // calculateCellFN(query_graph, edge_matrix, candidates_count, candidatesHC, idToValues, count2, u, valid_idx,RQ);
                count2++;
            }
            // else if (candidate == 0 &&cur_depth != max_depth - 1)
            else if (candidate == 0 && cur_depth != max_depth - 1)
            { // calculateCell(query_graph, edge_matrix, candidates_count, candidatesHC, idToValues, count2, u, valid_idx);

                calculateCellFN(query_graph, edge_matrix, candidates_count, candidatesHC, idToValues, count2, u, valid_idx, RQ);
                count2++;
            }
            // candidate = candidatesHC[u][valid_idx];
            VN[u].clear();
            if (cur_depth != max_depth - 1) // change here
            {
                // VN[u] = idToValues[u][candidate];
                VN[u] = idToValues[u][candidatesHC[u][valid_idx]];

                if (cur_depth != 0)
                {
                    int ia = idx[cur_depth];
                    int io = 0;
                    VNTemp.clear();
                    while (ia < idx_count[cur_depth] && io < VN[u].size())
                    {

                        ui valid_idxN = valid_candidate_idx[cur_depth][ia];
                        if (valid_idxN == 10000000)
                        {
                            ia += 1;
                            continue;
                        }
                        if (VN[u][io] == valid_idxN)
                        {
                            // VNTemp.push_back(VN[u][io]);
                            VNTemp.push_back(VN[u][io]);
                            ia++;
                            io++;
                        }

                        else if (VN[u][io] > valid_idxN)
                            ia++;
                        else
                            io++;
                        // ia++;
                    }
                    VN[u] = VNTemp;
                }
            }
#endif
            // set embedings found to empty for start

            if (visited_vertices[v])
            { // if conflict on the last level then equivalence are match.
#ifdef ENABLE_EQUIVALENT_SET
                VNTemp.clear();
                ui utemp = reverse_embedding[v];


                if (cur_depth != max_depth - 1)
                {
                    if (VN[u].size() > 0 && VN[utemp].size() > 0)
                    {
                        int ia = 0;
                        int io = 0;

                        while (ia < VN[u].size() && io < VN[utemp].size())
                        {
                            if (candidates[u][VN[u][ia]] == candidates[utemp][VN[utemp][io]])
                            {
                                VNTemp.push_back(VN[utemp][io]);
                                ia++;
                                io++;
                            }
                            else if (candidates[u][VN[u][ia]] < candidates[utemp][VN[utemp][io]])
                                ia++;
                            else
                                io++;
                        }
                    }
                    VN[utemp] = VNTemp;
                    VNTemp.clear();
                    // break;
                }
                else
                { // enumerate the results if it has equivalent nodes.
                    if (VN[utemp].size() > 1)
                    {
                        VNTemp.clear();
                        for (int ai = 0; ai < VN[utemp].size(); ai++)
                        {
                            if (!visited_vertices[candidates[utemp][VN[utemp][ai]]])
                                VNTemp.push_back(VN[utemp][ai]);
                        }
                        if (VNTemp.size() >= 1)
                        {
                            // heuristic there is a match so we count.
                            //  we can also exclude by just using VTEMP
                            // for utemp to count the correct number of matchings

                            int UNPMtemp = 0;
                            UNPMtemp = EmbSum.insert(v).second;
                            int ao = cur_depth - 1; //-1;
                            ui vqo = order[ao];
                            int UNPM2 = 0;
                            UNPM = 0; //
                            while (ao >= 0)
                            {
                                Match_BA[ao] = true;
                                // EmbSum.insert(embedding[vqo]);
                                UNPMtemp = EmbSum.insert(embedding[vqo]).second;

#ifdef ENABLE_EQUIVALENT_SET
                                //&&ao!=max_depth-1
                                if (VN[vqo].size() > 0)
                                {
                                    for (int dd = 0; dd < VN[vqo].size(); dd++)
                                    {
                                        if (VN[vqo][dd] != 10000000)
                                        {
                                            UNPMtemp += EmbSum.insert(candidates[vqo][VN[vqo][dd]]).second;
                                        }
                                    }
                                }
                                if (UNPMtemp > 1)
                                {
                                    UNPM2 += (UNPMtemp - 1);
                                    UNPMtemp = 1;
                                }
#endif
                                ao--;
                                vqo = order[ao];
                                UNPM += UNPMtemp;
                            }
#ifdef TOPKGREEDY
                            if (UNPM > pq.top())
                            {
                                if (heapcount >= ksize)
                                {
                                    pq.pop();
                                }
                                pq.push(UNPM);
                                heapcount++;
                            }
                            embedding_cnt += UNPM2 - 1;
                            while (UNPM2 > 0 && pq.top() < 1)
                            {
                                if (heapcount >= ksize)
                                {
                                    pq.pop();
                                }
                                pq.push(1);
                                heapcount++;
                                UNPM2--;
                            }
#endif
#ifdef ENABLE_FAILING_SET
                            // reverse_embedding.erase(embedding[u]);
                            vec_failing_set[cur_depth].set();
                            vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
#endif
                        }
                        else
                        { // update equivalent class of VN[utemp]
                            int ia = 0;
                            int io = 0;
                            VNTemp.clear();
                            while (ia < idx_count[cur_depth] && io < VN[utemp].size())
                            {

                                ui valid_idxN = valid_candidate_idx[cur_depth][ia];
                                if (valid_idxN == 10000000)
                                {
                                    ia += 1;
                                    continue;
                                }
                                if (VN[utemp][io] == valid_idxN)
                                {
                                    // VNTemp.push_back(VN[u][io]);
                                    VNTemp.push_back(VN[utemp][io]);
                                    ia++;
                                    io++;
                                }

                                else if (VN[utemp][io] > valid_idxN)
                                    ia++;
                                else
                                    io++;
                                // ia++;
                            }
                            VN[utemp] = VNTemp;
                            VNTemp.clear();
                        }
                    }
                }
                // }
                //}
                VN[u].clear();
#endif
                idx[cur_depth] += 1;
#ifdef ENABLE_FAILING_SET
                vec_failing_set[cur_depth] = ancestors[u];
                vec_failing_set[cur_depth] |= ancestors[reverse_embedding[v]];
                vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
#endif
                continue;
            }
            embd_vertices[v] = true;
            counter++;
            embedding[u] = v;             // set the embeding of depth qid to v
            idx_embedding[u] = valid_idx; // reverse the id of the embdeing to the position
            visited_vertices[v] = true;
            idx[cur_depth] += 1; // next element

            while (cur_depth < idx_count[cur_depth] && valid_candidate_idx[cur_depth][idx[cur_depth]] == 10000000)
            {
                idx[cur_depth] += 1;
            }

#ifdef ENABLE_FAILING_SET
            reverse_embedding[v] = u;
#endif

            if (cur_depth == max_depth - 1)
            {
                embedding_cnt += 1;
                visited_vertices[v] = false;
                embd_vertices[v] = false;
                counter--;
                int UNPMtemp = 0;
                // UNPMtemp=EmbSum.insert(embedding[order[cur_depth]]).second;

                int ao = cur_depth; //-1;
                ui vqo = order[ao];
                int UNPM2 = 0;
                UNPM = 0; //
                while (ao >= 0)
                {
                    Match_BA[ao] = true;
                    // EmbSum.insert(embedding[vqo]);
                    UNPMtemp = EmbSum.insert(embedding[vqo]).second;

#ifdef ENABLE_EQUIVALENT_SET
                    //&&ao!=max_depth-1
                    if (VN[vqo].size() > 0)
                    {
                        for (int dd = 0; dd < VN[vqo].size(); dd++)
                        {
                            if (VN[vqo][dd] != 10000000)
                            {
                                UNPMtemp += EmbSum.insert(candidates[vqo][VN[vqo][dd]]).second;
                            }
                        }
                    }
                    if (UNPMtemp > 1)
                    {
                        UNPM2 += (UNPMtemp - 1);
                        UNPMtemp = 1;
                    }
#endif
                    // Match_BA[ao] = true;
                    ao--;
                    vqo = order[ao];
                    UNPM += UNPMtemp;
                }
#ifdef TOPKGREEDY
                if (UNPM > pq.top())
                {
                    if (heapcount >= ksize)
                    {
                        pq.pop();
                    }
                    pq.push(UNPM);
                    heapcount++;
                }
                embedding_cnt += UNPM2;
                while (UNPM2 > 0 && pq.top() < 1)
                {
                    if (heapcount >= ksize)
                    {
                        pq.pop();
                    }
                    pq.push(1);
                    heapcount++;
                    UNPM2--;
                }
#endif

#ifdef ENABLE_FAILING_SET
                reverse_embedding.erase(embedding[u]);
                vec_failing_set[cur_depth].set();
                vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
#endif

                if (embedding_cnt >= output_limit_num)
                {
                    goto EXIT;
                }
            }
            else
            {
                call_count += 1;
                cur_depth += 1;     // next depth
                idx[cur_depth] = 0; // set the element to 0
                // generateValidCandidateIndex(cur_depth, idx_embedding, idx_count, valid_candidate_idx, edge_matrix, bn,
                //                             bn_count, order, temp_buffer, candidates, candidates_count, query_graph);
                generateValidCandidateIndex(cur_depth, idx_embedding, idx_count, valid_candidate_idx, edge_matrix, bn,
                                            bn_count, order, temp_buffer);
                pruneCandidatesIndexBySymmetryBreaking(cur_depth, embedding, order, candidates_count, candidates,
                                                       idx_count, valid_candidate_idx, ordered_constraints);

                // get the candidates for next depth
#ifdef ENABLE_FAILING_SET
                if (idx_count[cur_depth] == 0)
                {
                    vec_failing_set[cur_depth - 1] = ancestors[order[cur_depth]];
                }
                else
                {
                    vec_failing_set[cur_depth - 1].reset();
                }
#endif
            }
        }

        cur_depth -= 1;

        if (cur_depth == max_depth - 1)
        {
            continue;
        }

        if (cur_depth < 0)
            break;
        else
        {

#ifdef ENABLE_EQUIVALENT_SET

            // if(idx[cur_depth]>idx_count[cur_depth])
            // continue;
            // if(valid_candidate_idx[cur_depth][idx[cur_depth]]==10000000)
            // continue;

            // if (valid_candidate_idx[cur_depth][idx[cur_depth]]==10000000)  {
            //     idx[cur_depth]++;
            //     continue;
            // }

#endif
            ui vq = order[cur_depth];
            if (VN[vq].size() > 1)
            { // something to remove
                // if (idx[cur_depth]==0){
                //     cout<<"my logic is wrong"<<endl;
                // }
                int add = 0;
                for (int dd = 0; dd < VN[vq].size(); dd++)
                {

                    //  if (VN[vq][dd]==idx[cur_depth]-1){
                    //  cout<<"oups"<<endl;
                    //  continue;
                    //}
                    for (int kk = (idx[cur_depth]); kk < idx_count[cur_depth]; kk++)
                    {
                        ui valid_idxn = valid_candidate_idx[cur_depth][kk];
                        if (valid_idxn == 10000000)
                        {

                            continue;
                        }

                        if (valid_idxn == VN[vq][dd])
                        {
                            valid_candidate_idx[cur_depth][kk] = 10000000;
                            break;
                        }
                    }
                }
            }
            VN[vq].clear();

            VertexID u = order[cur_depth];
#ifdef ENABLE_FAILING_SET
            reverse_embedding.erase(embedding[u]);
            if (cur_depth != 0)
            {
                if (!vec_failing_set[cur_depth].test(u))
                {
                    vec_failing_set[cur_depth - 1] = vec_failing_set[cur_depth];
                    idx[cur_depth] = idx_count[cur_depth];
                }
                else
                {
                    vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
                }
            }
#endif
            visited_vertices[embedding[u]] = false;
            embd_vertices[embedding[u]] = false;
            counter--;
            // if (idx[cur_depth] == idx_count[cur_depth])
            // continue;

            // if (Match_BA[vq] == false && idx[cur_depth] < idx_count[cur_depth])
            //{ positive and negative removal now. dont need to count positive symmetric
        }
    }

EXIT:
    releaseBuffer(max_depth, idx, idx_count, embedding, idx_embedding, temp_buffer, valid_candidate_idx,
                  visited_vertices,
                  bn, bn_count);

#ifdef TOPKGREEDY
    std::priority_queue<int> maxHeap;
    int topkcounter = 0;
    int greedysum = 0;
    while (!pq.empty())
    {
        maxHeap.push(pq.top());
        pq.pop();
    }
    cout << maxHeap.top() << endl;
    while (!maxHeap.empty())
    {
        greedysum += maxHeap.top();
        maxHeap.pop();
        topkcounter++;
        switch (topkcounter)
        {
        case VALUE_5:
            s.topk[0] = greedysum;
            cout << greedysum << endl;
            break;
        case VALUE_10:
            s.topk[1] = greedysum;
            break;
        case VALUE_50:
            s.topk[2] = greedysum;
            break;
        case VALUE_100:
            s.topk[3] = greedysum;
            break;
        case VALUE_250:
            s.topk[4] = greedysum;
            break;
        case VALUE_500:
            s.topk[5] = greedysum;
            break;
        case VALUE_750:
            s.topk[6] = greedysum;
            break;
        case VALUE_1000:
            s.topk[7] = greedysum;
            break;
        default:
            break;
        }
    }
    for (int i = 0; i < 8; i++)
    {
        if (s.topk[i] == 0)
            s.topk[i] = greedysum;
    }
#endif

    int true_cand_sum = 0;
    s.Can_embed = EmbSum.size();
    EmbSum.clear();
    s.embedding_cnt = embedding_cnt;
    s.candidate_true_count_sum = true_cand_sum;
    return s;
}

enumResult
EvaluateQuery::DIVSM(const Graph *data_graph, const Graph *query_graph,ui *&nodeId, Edges ***edge_matrix, ui **candidates,
                       ui *candidates_count, ui *order, size_t output_limit_num, size_t &call_count, int TimeL, int FairT, const std::unordered_map<VertexID, std::pair<std::set<VertexID>, std::set<VertexID>>> &ordered_constraints)
{
    //assumming that creating matrix can be done before we do not count the time.
    int qsiz = query_graph->getVerticesCount();
    int ksize = 1000;
    int SolPos[ksize+1][qsiz];
    std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, std::greater<>> pq;
    int tempPos[qsiz];
    
    auto start = std::chrono::high_resolution_clock::now();

    enumResult s;
    // Generate bn.

    ui **bn;
    ui *bn_count;
    generateBN(query_graph, order, bn, bn_count);
    //unordered_set<ui> EmbSum;
    bool NotifyALL = true;
    // Allocate the memory buffer.
    ui *idx;
    ui *idx_count;
    ui *embedding;
    ui *idx_embedding;
    ui *temp_buffer;
    ui **valid_candidate_idx;
    bool *visited_vertices;
    int addPos=0;


    int UNPM = 0;
    allocateBuffer(data_graph, query_graph, candidates_count, idx, idx_count, embedding, idx_embedding,
                   temp_buffer, valid_candidate_idx, visited_vertices);
    vector<ui> VN[query_graph->getVerticesCount()];
    vector<ui> VNTemp;
    

    int RQ[query_graph->getVerticesCount()];
    for (int i = 0; i < query_graph->getVerticesCount(); i++)
    {
        RQ[order[i]] = i; // depth RQ[u]=
    }
    bool *embd_vertices = new bool[data_graph->getVerticesCount()];
    std::fill(embd_vertices, embd_vertices + data_graph->getVerticesCount(), false);
    size_t embedding_cnt = 0;
    int cur_depth = 0;
    int max_depth = query_graph->getVerticesCount();
    VertexID start_vertex = order[0];
    bool Match_BA[max_depth] = {false};
    idx[cur_depth] = 0;                                    // set the starting position to 0 for the current depth
    idx_count[cur_depth] = candidates_count[start_vertex]; // set the size for current depth
    // all the nodes are possible candidates

    // end1 = std::chrono::high_resolution_clock::now();
    int count2 = 1;

    size_t **candidatesHC = NULL;
    unordered_map<size_t, vector<ui>> *idToValues;
    idToValues = new unordered_map<size_t, vector<ui>>[query_graph->getVerticesCount()];
    candidatesHC = new size_t *[qsiz];

    for (ui i = 0; i < qsiz; ++i)
    {
        candidatesHC[i] = new size_t[candidates_count[i]];
        memset(candidatesHC[i], 0, sizeof(size_t) * candidates_count[i]);
    }
    calculateCell(query_graph, edge_matrix, candidates_count, candidatesHC, idToValues, count2, start_vertex, 0);
    // calculateCellFN(query_graph, edge_matrix, candidates_count, candidatesHC, idToValues, count2, start_vertex, 0,RQ);
    count2++;
    VN[start_vertex] = idToValues[start_vertex][candidatesHC[start_vertex][0]];


    for (ui i = 0; i < idx_count[cur_depth]; ++i)
    {
        valid_candidate_idx[cur_depth][i] = i;
    }


    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> ancestors;
    computeAncestor(query_graph, bn, bn_count, order, ancestors);
    AncestorWithSymmetryBreaking(query_graph, ancestors, ordered_constraints);
    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> vec_failing_set(max_depth);
    std::unordered_map<VertexID, VertexID> reverse_embedding;
    reverse_embedding.reserve(MAXIMUM_QUERY_GRAPH_SIZE * 2);
    bool del = false;

    int k;
    VertexID utemp; // u is current depth qid
    VertexID vtemp;
    int counter = 0;
    double ens = 0;
    int ordCand[idx_count[0]] = {100000};
    TimeL = TimeL * 1000;


      



    while (true)
    {
        // current position smaller than size
        while (idx[cur_depth] < idx_count[cur_depth])
        {
            Match_BA[cur_depth] = false;
            auto end = std::chrono::high_resolution_clock::now();
            ens = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
            if (ens > TimeL)
            { // 1000 1 sec
                goto EXIT;
            }
            ui valid_idx = valid_candidate_idx[cur_depth][idx[cur_depth]]; // valid id
            if (valid_idx == 10000000)                                     //||
            {
                idx[cur_depth] += 1;
                continue;
            }
            
            if (cur_depth == 0 && idx[cur_depth] != 0  && FairT > 0)
            {
                if (FairT == 1)
                    rankSimple(order[0], nodeId, candidates, valid_candidate_idx, idx, valid_idx, idx_count);
                else if (FairT == 2) 
                    rankLess(query_graph, edge_matrix, order[0], nodeId, candidates, valid_candidate_idx, idx, valid_idx, idx_count);
                else if (FairT == 3)
                    rankMore(query_graph, edge_matrix, order[0], nodeId, candidates, valid_candidate_idx, idx, valid_idx, idx_count);
                else if (FairT == 4)
                    rankFinal(query_graph, edge_matrix, order[0], nodeId, candidates, valid_candidate_idx, idx, valid_idx, idx_count);
                else if (FairT == 5)
                    rankFinal1(query_graph, edge_matrix, order[0], nodeId, candidates, valid_candidate_idx, idx, valid_idx, idx_count);
                else if (FairT == 6)
                    rankRandom(query_graph, edge_matrix, order[0], nodeId, candidates, valid_candidate_idx, idx, valid_idx, idx_count);
                valid_idx = valid_candidate_idx[cur_depth][idx[cur_depth]];
            if (valid_idx == 10000000)                                     //||
            {
                idx[cur_depth] += 1;
                continue;
            }
            }
            // from all the valid IDS

            VertexID u = order[cur_depth]; // u is current depth qid

            VertexID v = candidates[u][valid_idx]; // v is FCS[qid][vid] valid_idx-> position in CS

            size_t candidate = candidatesHC[u][valid_idx];
            if (candidate == 0 && cur_depth == 0)
            {
                calculateCell(query_graph, edge_matrix, candidates_count, candidatesHC, idToValues, count2, u, valid_idx);
                //calculateCell(query_graph, edge_matrix, candidates_count, candidatesHC, idToValues, count2, u, valid_idx,RQ);
                count2++;
            }
            // else if (candidate == 0 &&cur_depth != max_depth - 1)
            else if (candidate == 0 && cur_depth != max_depth - 1)
            {
                //calculateCell(query_graph, edge_matrix, candidates_count, candidatesHC, idToValues, count2, u, valid_idx);
                //calculateCell(query_graph, edge_matrix, candidates_count, candidatesHC, idToValues, count2, u, valid_idx);
                calculateCellFN(query_graph, edge_matrix, candidates_count, candidatesHC, idToValues, count2, u, valid_idx, RQ);
                count2++;
            }
            // candidate = candidatesHC[u][valid_idx];
            VN[u].clear();
            if (cur_depth != max_depth - 1) // change here
            {
                // VN[u] = idToValues[u][candidate];
                VN[u] = idToValues[u][candidatesHC[u][valid_idx]];

                if (cur_depth != 0)
                {
                    int ia = idx[cur_depth];
                    int io = 0;
                    VNTemp.clear();
                    while (ia < idx_count[cur_depth] && io < VN[u].size())
                    {

                        ui valid_idxN = valid_candidate_idx[cur_depth][ia];
                        if (valid_idxN == 10000000)
                        {
                            ia += 1;
                            continue;
                        }
                        if (VN[u][io] == valid_idxN)
                        {
                            VNTemp.push_back(VN[u][io]);
                            ia++;
                            io++;
                        }

                        else if (VN[u][io] > valid_idxN)
                            ia++;
                        else
                            io++;
                    }
                    VN[u] = VNTemp;
                }
            }

            // set embedings found to empty for start

            if (visited_vertices[v])
            { // if conflict on the last level then equivalence are match.

                VNTemp.clear();
                ui utemp = reverse_embedding[v];


                if (cur_depth != max_depth - 1)
                {
                    if (VN[u].size() > 0 && VN[utemp].size() > 0)
                    {
                        int ia = 0;
                        int io = 0;

                        while (ia < VN[u].size() && io < VN[utemp].size())
                        {
                            if (candidates[u][VN[u][ia]] == candidates[utemp][VN[utemp][io]])
                            {
                                VNTemp.push_back(VN[utemp][io]);
                                ia++;
                                io++;
                            }
                            else if (candidates[u][VN[u][ia]] < candidates[utemp][VN[utemp][io]])
                                ia++;
                            else
                                io++;
                        }
                    }
                    VN[utemp] = VNTemp;
                    VNTemp.clear();
                    // break;
                }
                else
                { // enumerate the results if it has equivalent nodes.
                    if (VN[utemp].size() > 1)
                    {
                        VNTemp.clear();
                        for (int ai = 0; ai < VN[utemp].size(); ai++)
                        {
                            if (!visited_vertices[candidates[utemp][VN[utemp][ai]]])
                                VNTemp.push_back(VN[utemp][ai]);
                        }
                                
                                    
                                
                        if (VNTemp.size() >= 1)
                        {
                            // heuristic there is a match so we count.
                            //  we can also exclude by just using VTEMP
                            // for utemp to count the correct number of matchings
                            //temp_sol[utemp].push_back(candidates[utemp][VNTemp[0]]);
                            int UNPMtemp = 0;
                            //UNPMtemp = EmbSum.insert(v).second;
                            nodeId[v]=1;
                            Match_BA[0] = true;
                            int ao = cur_depth - 1; //-1;
                            ui vqo = order[ao];
                            UNPM = 0; //
                            while (ao >= 0)
                            {  
                                Match_BA[ao] = true;
                                // EmbSum.insert(embedding[vqo]);
                                if(embedding[vqo]!=v){   
                                    nodeId[embedding[vqo]]=1;                                 
                                }else
                                {
                                    nodeId[candidates[vqo][VNTemp[0]]]=1; 
                                }
                                //&&ao!=max_depth-1
                                if (VN[vqo].size() > 0)
                                {
                                    for (int dd = 0; dd < VN[vqo].size(); dd++)
                                    {
                                        if (VN[vqo][dd] != 10000000)
                                        {     
                                            nodeId[candidates[vqo][VN[vqo][dd]]]=1;                                                                                                                                  
                                        }
                                    }
                                }

                                ao--;
                                vqo = order[ao];
                            }
                            // reverse_embedding.erase(embedding[u]);
                            vec_failing_set[cur_depth].set();
                            vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];

                        }
                        else
                        { // update equivalent class of VN[utemp]
                            int ia = 0;
                            int io = 0;
                            VNTemp.clear();
                            while (ia < idx_count[cur_depth] && io < VN[utemp].size())
                            {

                                ui valid_idxN = valid_candidate_idx[cur_depth][ia];
                                if (valid_idxN == 10000000)
                                {
                                    ia += 1;
                                    continue;
                                }
                                if (VN[utemp][io] == valid_idxN)
                                {
                                    // VNTemp.push_back(VN[u][io]);
                                    VNTemp.push_back(VN[utemp][io]);
                                    ia++;
                                    io++;
                                }

                                else if (VN[utemp][io] > valid_idxN)
                                    ia++;
                                else
                                    io++;
                                // ia++;
                            }
                            VN[utemp] = VNTemp;
                            VNTemp.clear();
                        }
                    }
                }
                // }
                //}
                VN[u].clear();

                idx[cur_depth] += 1;
                vec_failing_set[cur_depth] = ancestors[u];
                vec_failing_set[cur_depth] |= ancestors[reverse_embedding[v]];
                vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];

                continue;
            }
            embd_vertices[v] = true;
            counter++;
            embedding[u] = v;             // set the embeding of depth qid to v
            idx_embedding[u] = valid_idx; // reverse the id of the embdeing to the position
            visited_vertices[v] = true;
            idx[cur_depth] += 1; // next element

            while (cur_depth < idx_count[cur_depth] && valid_candidate_idx[cur_depth][idx[cur_depth]] == 10000000)
            {
                idx[cur_depth] += 1;
            }

            reverse_embedding[v] = u;

            if (cur_depth == max_depth - 1)
            {
                embedding_cnt += 1;
                visited_vertices[v] = false;
                embd_vertices[v] = false;
                counter--;
                int UNPMtemp = 0;
                // UNPMtemp=EmbSum.insert(embedding[order[cur_depth]]).second;

                int ao = cur_depth; //-1;
                ui vqo = order[ao];
                int UNPM2 = 0;
                UNPM = 0; //
                while (ao >= 0)
                {
                    Match_BA[ao] = true;
                    // EmbSum.insert(embedding[vqo]);
                    //UNPMtemp = EmbSum.insert(embedding[vqo]).second;
                     nodeId[embedding[vqo]]=1;  
                    //&&ao!=max_depth-1
                    if (VN[vqo].size() > 0)
                    {
                        for (int dd = 0; dd < VN[vqo].size(); dd++)
                        {
                            if (VN[vqo][dd] != 10000000)
                            {  
                                nodeId[candidates[vqo][VN[vqo][dd]]]=1;                                
                            }
                        }
                    }

                    // Match_BA[ao] = true;
                    ao--;
                    vqo = order[ao];
                }

                reverse_embedding.erase(embedding[u]);
                vec_failing_set[cur_depth].set();
                vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];


                if (embedding_cnt >= output_limit_num)
                {
                    goto EXIT;
                }
            }
            else
            {
                call_count += 1;
                cur_depth += 1;     // next depth
                idx[cur_depth] = 0; // set the element to 0
                // generateValidCandidateIndex(cur_depth, idx_embedding, idx_count, valid_candidate_idx, edge_matrix, bn,
                //                             bn_count, order, temp_buffer, candidates, candidates_count, query_graph);
                generateValidCandidateIndex(cur_depth, idx_embedding, idx_count, valid_candidate_idx, edge_matrix, bn,
                                            bn_count, order, temp_buffer);
                
                pruneCandidatesIndexBySymmetryBreakingC(cur_depth, embedding, order, candidates_count, candidates,
                                                       idx_count, valid_candidate_idx, ordered_constraints);

                // get the candidates for next depth
                if (idx_count[cur_depth] == 0)
                {
                    vec_failing_set[cur_depth - 1] = ancestors[order[cur_depth]];
                }
                else
                {
                    vec_failing_set[cur_depth - 1].reset();
                }
            }
        }

        cur_depth -= 1;

        if (cur_depth == max_depth - 1)
        {
            continue;
        }

        if (cur_depth < 0)
            break;
        else
        {
            ui vq = order[cur_depth];
            if (VN[vq].size() > 1)
            { 
                int add = 0;
                for (int dd = 0; dd < VN[vq].size(); dd++)
                {
                    for (int kk = (idx[cur_depth]); kk < idx_count[cur_depth]; kk++)
                    {
                        ui valid_idxn = valid_candidate_idx[cur_depth][kk];
                        if (valid_idxn == 10000000)
                        {
                            continue;
                        }

                        if (valid_idxn == VN[vq][dd])
                        {
                            valid_candidate_idx[cur_depth][kk] = 10000000;
                            break;
                        }
                    }
                }
            }
            VN[vq].clear();

            VertexID u = order[cur_depth];
            reverse_embedding.erase(embedding[u]);
            if (cur_depth != 0)
            {
                if (!vec_failing_set[cur_depth].test(u))
                {
                    vec_failing_set[cur_depth - 1] = vec_failing_set[cur_depth];
                    idx[cur_depth] = idx_count[cur_depth];
                }
                else
                {
                    vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
                }
            }

            visited_vertices[embedding[u]] = false;
            embd_vertices[embedding[u]] = false;
            counter--;
        }
    }

    EXIT:
    releaseBuffer(max_depth, idx, idx_count, embedding, idx_embedding, temp_buffer, valid_candidate_idx,
                  visited_vertices,
                  bn, bn_count);

    std::priority_queue<std::pair<int, int>> maxHeap;
    int topkcounter = 0;
    int greedysum = 0;
    int true_cand_sum = 0;
    int countSol=0;
    for (int i=0;i<data_graph->getVerticesCount();i++){
        countSol+=nodeId[i];
    }
    s.Can_embed = countSol;
    s.embedding_cnt = embedding_cnt;
    s.candidate_true_count_sum = true_cand_sum;
    return s;
}

enumResult
EvaluateQuery::DIVSMSQ(const Graph *data_graph, const Graph *query_graph,ui *&nodeId, Edges ***edge_matrix, ui **candidates,
                       ui *candidates_count, ui *order, size_t output_limit_num, size_t &call_count, int TimeL, int FairT, const std::unordered_map<VertexID, std::pair<std::set<VertexID>, std::set<VertexID>>> &ordered_constraints)
{
    //assumming that creating matrix can be done before we do not count the time.
    int qsiz = query_graph->getVerticesCount();
    int ksize = 1000;
    int SolPos[ksize+1][qsiz];
    std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, std::greater<>> pq;
    int tempPos[qsiz];
    
    auto start = std::chrono::high_resolution_clock::now();

    enumResult s;
    // Generate bn.

    ui **bn;
    ui *bn_count;
    generateBN(query_graph, order, bn, bn_count);
    //unordered_set<ui> EmbSum;
    bool NotifyALL = true;
    // Allocate the memory buffer.
    ui *idx;
    ui *idx_count;
    ui *embedding;
    ui *idx_embedding;
    ui *temp_buffer;
    ui **valid_candidate_idx;
    bool *visited_vertices;
    int addPos=0;


    int UNPM = 0;
    allocateBuffer(data_graph, query_graph, candidates_count, idx, idx_count, embedding, idx_embedding,
                   temp_buffer, valid_candidate_idx, visited_vertices);
    vector<ui> VN[query_graph->getVerticesCount()];
    vector<ui> VNTemp;
    

    int RQ[query_graph->getVerticesCount()];
    for (int i = 0; i < query_graph->getVerticesCount(); i++)
    {
        RQ[order[i]] = i; // depth RQ[u]=
    }
    bool *embd_vertices = new bool[data_graph->getVerticesCount()];
    std::fill(embd_vertices, embd_vertices + data_graph->getVerticesCount(), false);
    size_t embedding_cnt = 0;
    int cur_depth = 0;
    int max_depth = query_graph->getVerticesCount();
    VertexID start_vertex = order[0];
    bool Match_BA[max_depth] = {false};
    idx[cur_depth] = 0;                                    // set the starting position to 0 for the current depth
    idx_count[cur_depth] = candidates_count[start_vertex]; // set the size for current depth
    // all the nodes are possible candidates

    // end1 = std::chrono::high_resolution_clock::now();
    int count2 = 1;

    size_t **candidatesHC = NULL;
    unordered_map<size_t, vector<ui>> *idToValues;
    idToValues = new unordered_map<size_t, vector<ui>>[query_graph->getVerticesCount()];
    candidatesHC = new size_t *[qsiz];

    for (ui i = 0; i < qsiz; ++i)
    {
        candidatesHC[i] = new size_t[candidates_count[i]];
        memset(candidatesHC[i], 0, sizeof(size_t) * candidates_count[i]);
    }
    calculateCell(query_graph, edge_matrix, candidates_count, candidatesHC, idToValues, count2, start_vertex, 0);
    // calculateCellFN(query_graph, edge_matrix, candidates_count, candidatesHC, idToValues, count2, start_vertex, 0,RQ);
    count2++;
    VN[start_vertex] = idToValues[start_vertex][candidatesHC[start_vertex][0]];


    for (ui i = 0; i < idx_count[cur_depth]; ++i)
    {
        valid_candidate_idx[cur_depth][i] = i;
    }

    
    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> ancestors;
    computeAncestor(query_graph, bn, bn_count, order, ancestors);
    AncestorWithSymmetryBreaking(query_graph, ancestors, ordered_constraints);
    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> vec_failing_set(max_depth);
    std::unordered_map<VertexID, VertexID> reverse_embedding;
    reverse_embedding.reserve(MAXIMUM_QUERY_GRAPH_SIZE * 2);
    bool del = false;

    int k;
    VertexID utemp; // u is current depth qid
    VertexID vtemp;
    int counter = 0;
    double ens = 0;
    int ordCand[idx_count[0]] = {100000};
    TimeL = TimeL * 1000;
    
     while (true)
    {
        while (idx[cur_depth] < idx_count[cur_depth])
        {
            Match_BA[cur_depth] = false;
            ui valid_idx = valid_candidate_idx[cur_depth][idx[cur_depth]];
            VertexID u = order[cur_depth];
            VertexID v = candidates[u][valid_idx];
            auto end = std::chrono::high_resolution_clock::now();
            ens = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
            if (ens > TimeL)
            { // 1000 1 sec
                goto EXIT;
            }
            if (visited_vertices[v])
            {
                idx[cur_depth] += 1;
                vec_failing_set[cur_depth] = ancestors[u];
                vec_failing_set[cur_depth] |= ancestors[reverse_embedding[v]];
                vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
                continue;
            }
            if(nodeId[v]==1){
                idx[cur_depth] += 1;
                vec_failing_set[cur_depth] = ancestors[u];
                //vec_failing_set[cur_depth] |= ancestors[reverse_embedding[v]];
                //vec_failing_set[cur_depth] |= ancestors[order[i]];
                if(cur_depth - 1>0)
                vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
                continue;
            }

            embedding[u] = v;
            idx_embedding[u] = valid_idx;
            visited_vertices[v] = true;
            idx[cur_depth] += 1;

            reverse_embedding[v] = u;
            if (cur_depth == max_depth - 1)
            {
                UNPM = 0;
                if(nodeId[embedding[order[cur_depth]]]==0){
                    UNPM++;
                    nodeId[embedding[order[cur_depth]]]=1;
                }
                embedding_cnt += 1;
                visited_vertices[v] = false;
                int ao = cur_depth; //-1;
                ui vqo = order[ao];
                nodeId[v]=1;
                while (ao > 0)
                {                
                    if(nodeId[embedding[vqo]]==0){
                    UNPM++;
                    nodeId[embedding[vqo]]=1;
                }   if(ao>1)
                    idx[ao] = idx_count[ao];
                    ao--;
                    vqo = order[ao];
                }
                nodeId[embedding[vqo]]=1;


                reverse_embedding.erase(embedding[u]);
                vec_failing_set[cur_depth].set();
                vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
                if (embedding_cnt >= output_limit_num)
                {
                    goto EXIT;
                }
            }
            else
            {
                call_count += 1;
                cur_depth += 1;

                idx[cur_depth] = 0;
                generateValidCandidateIndex(cur_depth, idx_embedding, idx_count, valid_candidate_idx, edge_matrix, bn,
                                            bn_count, order, temp_buffer);

                
                if (idx_count[cur_depth] == 0)
                {
                    vec_failing_set[cur_depth - 1] = ancestors[order[cur_depth]];
                }
                else
                {
                    vec_failing_set[cur_depth - 1].reset();
                }
            }
        }

        cur_depth -= 1;
        if (cur_depth < 0)
            break;
        else
        {
            VertexID u = order[cur_depth];
            reverse_embedding.erase(embedding[u]);
            if (cur_depth != 0)
            {
                if (!vec_failing_set[cur_depth].test(u))
                {
                    vec_failing_set[cur_depth - 1] = vec_failing_set[cur_depth];
                    idx[cur_depth] = idx_count[cur_depth];
                }
                else
                {
                    vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
                }
            }
            visited_vertices[embedding[u]] = false;

        }
    }
    
    
    allocateBufferLM(data_graph, query_graph, candidates_count, idx, idx_count, embedding, idx_embedding,
                   temp_buffer, valid_candidate_idx, visited_vertices);
    ancestors.clear();
    computeAncestor(query_graph, bn, bn_count, order, ancestors);
    AncestorWithSymmetryBreaking(query_graph, ancestors, ordered_constraints);
     vec_failing_set.clear();
     reverse_embedding.clear();
    cur_depth = 0;
    idx[cur_depth] = 0;       
    start_vertex = order[0];                             // set the starting position to 0 for the current depth
    idx_count[cur_depth] = candidates_count[start_vertex]; // set the size for current depth
    while (true)
    {
        // current position smaller than size
        while (idx[cur_depth] < idx_count[cur_depth])
        {
            Match_BA[cur_depth] = false;
            auto end = std::chrono::high_resolution_clock::now();
            ens = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
            if (ens > TimeL)
            { // 1000 1 sec
                goto EXIT;
            }
            ui valid_idx = valid_candidate_idx[cur_depth][idx[cur_depth]]; // valid id
            if (valid_idx == 10000000)                                     //||
            {
                idx[cur_depth] += 1;
                continue;
            }
            
            if (cur_depth == 0 && idx[cur_depth] != 0  && FairT > 0)
            {
                if (FairT == 1)
                    rankSimple(order[0], nodeId, candidates, valid_candidate_idx, idx, valid_idx, idx_count);
                else if (FairT == 2) 
                    rankLess(query_graph, edge_matrix, order[0], nodeId, candidates, valid_candidate_idx, idx, valid_idx, idx_count);
                else if (FairT == 3)
                    rankMore(query_graph, edge_matrix, order[0], nodeId, candidates, valid_candidate_idx, idx, valid_idx, idx_count);
                else if (FairT == 4)
                    rankFinal(query_graph, edge_matrix, order[0], nodeId, candidates, valid_candidate_idx, idx, valid_idx, idx_count);
                else if (FairT == 5)
                    rankFinal1(query_graph, edge_matrix, order[0], nodeId, candidates, valid_candidate_idx, idx, valid_idx, idx_count);
                else if (FairT == 6)
                    rankRandom(query_graph, edge_matrix, order[0], nodeId, candidates, valid_candidate_idx, idx, valid_idx, idx_count);
                valid_idx = valid_candidate_idx[cur_depth][idx[cur_depth]];
            if (valid_idx == 10000000)                                     //||
            {
                idx[cur_depth] += 1;
                continue;
            }
            }
            // from all the valid IDS

            VertexID u = order[cur_depth]; // u is current depth qid

            VertexID v = candidates[u][valid_idx]; // v is FCS[qid][vid] valid_idx-> position in CS

            size_t candidate = candidatesHC[u][valid_idx];
            if (candidate == 0 && cur_depth == 0)
            {
                calculateCell(query_graph, edge_matrix, candidates_count, candidatesHC, idToValues, count2, u, valid_idx);
                // calculateCellFN(query_graph, edge_matrix, candidates_count, candidatesHC, idToValues, count2, u, valid_idx,RQ);
                count2++;
            }
            // else if (candidate == 0 &&cur_depth != max_depth - 1)
            else if (candidate == 0 && cur_depth != max_depth - 1)
            { // calculateCell(query_graph, edge_matrix, candidates_count, candidatesHC, idToValues, count2, u, valid_idx);

                calculateCellFN(query_graph, edge_matrix, candidates_count, candidatesHC, idToValues, count2, u, valid_idx, RQ);
                count2++;
            }
            // candidate = candidatesHC[u][valid_idx];
            VN[u].clear();
            if (cur_depth != max_depth - 1) // change here
            {
                // VN[u] = idToValues[u][candidate];
                VN[u] = idToValues[u][candidatesHC[u][valid_idx]];

                if (cur_depth != 0)
                {
                    int ia = idx[cur_depth];
                    int io = 0;
                    VNTemp.clear();
                    while (ia < idx_count[cur_depth] && io < VN[u].size())
                    {

                        ui valid_idxN = valid_candidate_idx[cur_depth][ia];
                        if (valid_idxN == 10000000)
                        {
                            ia += 1;
                            continue;
                        }
                        if (VN[u][io] == valid_idxN)
                        {
                            VNTemp.push_back(VN[u][io]);
                            ia++;
                            io++;
                        }

                        else if (VN[u][io] > valid_idxN)
                            ia++;
                        else
                            io++;
                    }
                    VN[u] = VNTemp;
                }
            }

            // set embedings found to empty for start

            if (visited_vertices[v])
            { // if conflict on the last level then equivalence are match.

                VNTemp.clear();
                ui utemp = reverse_embedding[v];


                if (cur_depth != max_depth - 1)
                {
                    if (VN[u].size() > 0 && VN[utemp].size() > 0)
                    {
                        int ia = 0;
                        int io = 0;

                        while (ia < VN[u].size() && io < VN[utemp].size())
                        {
                            if (candidates[u][VN[u][ia]] == candidates[utemp][VN[utemp][io]])
                            {
                                VNTemp.push_back(VN[utemp][io]);
                                ia++;
                                io++;
                            }
                            else if (candidates[u][VN[u][ia]] < candidates[utemp][VN[utemp][io]])
                                ia++;
                            else
                                io++;
                        }
                    }
                    VN[utemp] = VNTemp;
                    VNTemp.clear();
                    // break;
                }
                else
                { // enumerate the results if it has equivalent nodes.
                    if (VN[utemp].size() > 1)
                    {
                        VNTemp.clear();
                        for (int ai = 0; ai < VN[utemp].size(); ai++)
                        {
                            if (!visited_vertices[candidates[utemp][VN[utemp][ai]]])
                                VNTemp.push_back(VN[utemp][ai]);
                        }
                                
                                    
                                
                        if (VNTemp.size() >= 1)
                        {
                            // heuristic there is a match so we count.
                            //  we can also exclude by just using VTEMP
                            // for utemp to count the correct number of matchings
                            //temp_sol[utemp].push_back(candidates[utemp][VNTemp[0]]);
                            int UNPMtemp = 0;
                            //UNPMtemp = EmbSum.insert(v).second;
                            nodeId[v]=1;
                            int ao = cur_depth - 1; //-1;
                            ui vqo = order[ao];
                            UNPM = 0; //
                            while (ao >= 0)
                            {  
                                Match_BA[ao] = true;
                                // EmbSum.insert(embedding[vqo]);
                                if(embedding[vqo]!=v){   
                                    nodeId[embedding[vqo]]=1;                                 
                                }else
                                {
                                    nodeId[candidates[vqo][VNTemp[0]]]=1; 
                                }
                                //&&ao!=max_depth-1
                                if (VN[vqo].size() > 0)
                                {
                                    for (int dd = 0; dd < VN[vqo].size(); dd++)
                                    {
                                        if (VN[vqo][dd] != 10000000)
                                        {     
                                            nodeId[candidates[vqo][VN[vqo][dd]]]=1;                                                                                                                                  
                                        }
                                    }
                                }

                                ao--;
                                vqo = order[ao];
                            }
                            // reverse_embedding.erase(embedding[u]);
                            vec_failing_set[cur_depth].set();
                            vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];

                        }
                        else
                        { // update equivalent class of VN[utemp]
                            int ia = 0;
                            int io = 0;
                            VNTemp.clear();
                            while (ia < idx_count[cur_depth] && io < VN[utemp].size())
                            {

                                ui valid_idxN = valid_candidate_idx[cur_depth][ia];
                                if (valid_idxN == 10000000)
                                {
                                    ia += 1;
                                    continue;
                                }
                                if (VN[utemp][io] == valid_idxN)
                                {
                                    // VNTemp.push_back(VN[u][io]);
                                    VNTemp.push_back(VN[utemp][io]);
                                    ia++;
                                    io++;
                                }

                                else if (VN[utemp][io] > valid_idxN)
                                    ia++;
                                else
                                    io++;
                                // ia++;
                            }
                            VN[utemp] = VNTemp;
                            VNTemp.clear();
                        }
                    }
                }
                // }
                //}
                VN[u].clear();

                idx[cur_depth] += 1;
                vec_failing_set[cur_depth] = ancestors[u];
                vec_failing_set[cur_depth] |= ancestors[reverse_embedding[v]];
                vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];

                continue;
            }
            embd_vertices[v] = true;
            counter++;
            embedding[u] = v;             // set the embeding of depth qid to v
            idx_embedding[u] = valid_idx; // reverse the id of the embdeing to the position
            visited_vertices[v] = true;
            idx[cur_depth] += 1; // next element

            while (cur_depth < idx_count[cur_depth] && valid_candidate_idx[cur_depth][idx[cur_depth]] == 10000000)
            {
                idx[cur_depth] += 1;
            }

            reverse_embedding[v] = u;

            if (cur_depth == max_depth - 1)
            {
                embedding_cnt += 1;
                visited_vertices[v] = false;
                embd_vertices[v] = false;
                counter--;
                int UNPMtemp = 0;
                // UNPMtemp=EmbSum.insert(embedding[order[cur_depth]]).second;

                int ao = cur_depth; //-1;
                ui vqo = order[ao];
                int UNPM2 = 0;
                UNPM = 0; //
                while (ao >= 0)
                {
                    Match_BA[ao] = true;
                    // EmbSum.insert(embedding[vqo]);
                    //UNPMtemp = EmbSum.insert(embedding[vqo]).second;
                     nodeId[embedding[vqo]]=1;  
                    //&&ao!=max_depth-1
                    if (VN[vqo].size() > 0)
                    {
                        for (int dd = 0; dd < VN[vqo].size(); dd++)
                        {
                            if (VN[vqo][dd] != 10000000)
                            {  
                                nodeId[candidates[vqo][VN[vqo][dd]]]=1;                                
                            }
                        }
                    }

                    // Match_BA[ao] = true;
                    ao--;
                    vqo = order[ao];
                }

                reverse_embedding.erase(embedding[u]);
                vec_failing_set[cur_depth].set();
                vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];


                if (embedding_cnt >= output_limit_num)
                {
                    goto EXIT;
                }
            }
            else
            {
                call_count += 1;
                cur_depth += 1;     // next depth
                idx[cur_depth] = 0; // set the element to 0
                // generateValidCandidateIndex(cur_depth, idx_embedding, idx_count, valid_candidate_idx, edge_matrix, bn,
                //                             bn_count, order, temp_buffer, candidates, candidates_count, query_graph);
                generateValidCandidateIndex(cur_depth, idx_embedding, idx_count, valid_candidate_idx, edge_matrix, bn,
                                            bn_count, order, temp_buffer);
                pruneCandidatesIndexBySymmetryBreakingC(cur_depth, embedding, order, candidates_count, candidates,
                                                       idx_count, valid_candidate_idx, ordered_constraints);

                // get the candidates for next depth
                if (idx_count[cur_depth] == 0)
                {
                    vec_failing_set[cur_depth - 1] = ancestors[order[cur_depth]];
                }
                else
                {
                    vec_failing_set[cur_depth - 1].reset();
                }
            }
        }

        cur_depth -= 1;

        if (cur_depth == max_depth - 1)
        {
            continue;
        }

        if (cur_depth < 0)
            break;
        else
        {
            ui vq = order[cur_depth];
            if (VN[vq].size() > 1)
            { 
                int add = 0;
                for (int dd = 0; dd < VN[vq].size(); dd++)
                {
                    for (int kk = (idx[cur_depth]); kk < idx_count[cur_depth]; kk++)
                    {
                        ui valid_idxn = valid_candidate_idx[cur_depth][kk];
                        if (valid_idxn == 10000000)
                        {
                            continue;
                        }

                        if (valid_idxn == VN[vq][dd])
                        {
                            valid_candidate_idx[cur_depth][kk] = 10000000;
                            break;
                        }
                    }
                }
            }
            VN[vq].clear();

            VertexID u = order[cur_depth];
            reverse_embedding.erase(embedding[u]);
            if (cur_depth != 0)
            {
                if (!vec_failing_set[cur_depth].test(u))
                {
                    vec_failing_set[cur_depth - 1] = vec_failing_set[cur_depth];
                    idx[cur_depth] = idx_count[cur_depth];
                }
                else
                {
                    vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
                }
            }

            visited_vertices[embedding[u]] = false;
            embd_vertices[embedding[u]] = false;
            counter--;
        }
    }

    EXIT:
    releaseBuffer(max_depth, idx, idx_count, embedding, idx_embedding, temp_buffer, valid_candidate_idx,
                  visited_vertices,
                  bn, bn_count);

    std::priority_queue<std::pair<int, int>> maxHeap;
    int topkcounter = 0;
    int greedysum = 0;
    int true_cand_sum = 0;
    int countSol=0;
    for (int i=0;i<data_graph->getVerticesCount();i++){
        countSol+=nodeId[i];
    }
    s.Can_embed = countSol;
    s.embedding_cnt = embedding_cnt;
    s.candidate_true_count_sum = true_cand_sum;
    return s;
}



enumResult
EvaluateQuery::DIVTOPK(const Graph *data_graph, const Graph *query_graph,ui *&nodeId, Edges ***edge_matrix, ui **candidates,
                       ui *candidates_count, ui *order, size_t output_limit_num, size_t &call_count, int TimeL, int FairT, const std::unordered_map<VertexID, std::pair<std::set<VertexID>, std::set<VertexID>>> &ordered_constraints)
{
    //assumming that creating matrix can be done before we do not count the time.
    int qsiz = query_graph->getVerticesCount();
    int ksize = 1000;
    int SolPos[ksize+1][qsiz];
    std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, std::greater<>> pq;
    int tempPos[qsiz];
    
    auto start = std::chrono::high_resolution_clock::now();
    pq.push({0,0});
    int heapcount = 0;
    enumResult s;
    std::vector<std::vector<int>> temp_sol(qsiz);
    // Generate bn.

    ui **bn;
    ui *bn_count;
    generateBN(query_graph, order, bn, bn_count);
    //unordered_set<ui> EmbSum;
    bool NotifyALL = true;
    // Allocate the memory buffer.
    ui *idx;
    ui *idx_count;
    ui *embedding;
    ui *idx_embedding;
    ui *temp_buffer;
    ui **valid_candidate_idx;
    bool *visited_vertices;
    int addPos=0;


    int UNPM = 0;
    allocateBuffer(data_graph, query_graph, candidates_count, idx, idx_count, embedding, idx_embedding,
                   temp_buffer, valid_candidate_idx, visited_vertices);
    vector<ui> VN[query_graph->getVerticesCount()];
    vector<ui> VNTemp;
    

    int RQ[query_graph->getVerticesCount()];
    for (int i = 0; i < query_graph->getVerticesCount(); i++)
    {
        RQ[order[i]] = i; // depth RQ[u]=
    }
    bool *embd_vertices = new bool[data_graph->getVerticesCount()];
    std::fill(embd_vertices, embd_vertices + data_graph->getVerticesCount(), false);
    size_t embedding_cnt = 0;
    int cur_depth = 0;
    int max_depth = query_graph->getVerticesCount();
    VertexID start_vertex = order[0];
    bool Match_BA[max_depth] = {false};
    idx[cur_depth] = 0;                                    // set the starting position to 0 for the current depth
    idx_count[cur_depth] = candidates_count[start_vertex]; // set the size for current depth
    // all the nodes are possible candidates

    // end1 = std::chrono::high_resolution_clock::now();
    int count2 = 1;

    size_t **candidatesHC = NULL;
    unordered_map<size_t, vector<ui>> *idToValues;
    idToValues = new unordered_map<size_t, vector<ui>>[query_graph->getVerticesCount()];
    candidatesHC = new size_t *[qsiz];

    for (ui i = 0; i < qsiz; ++i)
    {
        candidatesHC[i] = new size_t[candidates_count[i]];
        memset(candidatesHC[i], 0, sizeof(size_t) * candidates_count[i]);
    }
    calculateCell(query_graph, edge_matrix, candidates_count, candidatesHC, idToValues, count2, start_vertex, 0);
    // calculateCellFN(query_graph, edge_matrix, candidates_count, candidatesHC, idToValues, count2, start_vertex, 0,RQ);
    count2++;
    VN[start_vertex] = idToValues[start_vertex][candidatesHC[start_vertex][0]];


    for (ui i = 0; i < idx_count[cur_depth]; ++i)
    {
        valid_candidate_idx[cur_depth][i] = i;
    }


    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> ancestors;
    computeAncestor(query_graph, bn, bn_count, order, ancestors);
    AncestorWithSymmetryBreaking(query_graph, ancestors, ordered_constraints);
    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> vec_failing_set(max_depth);
    std::unordered_map<VertexID, VertexID> reverse_embedding;
    reverse_embedding.reserve(MAXIMUM_QUERY_GRAPH_SIZE * 2);
    bool del = false;

    int k;
    VertexID utemp; // u is current depth qid
    VertexID vtemp;
    int counter = 0;
    double ens = 0;
    int ordCand[idx_count[0]] = {100000};
    TimeL = TimeL * 1000;
    while (true)
    {
        // current position smaller than size
        while (idx[cur_depth] < idx_count[cur_depth])
        {
            Match_BA[cur_depth] = false;
            auto end = std::chrono::high_resolution_clock::now();
            ens = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
            if (ens > TimeL)
            { // 1000 1 sec
                goto EXIT;
            }
            ui valid_idx = valid_candidate_idx[cur_depth][idx[cur_depth]]; // valid id
            if (valid_idx == 10000000)                                     //||
            {
                idx[cur_depth] += 1;
                continue;
            }
            /*
            if (cur_depth == 0 && idx[cur_depth] != 0 && EmbSum.size() > 0 && FairT > 0)
            {
                if (FairT == 1)
                    rankSimple(order[0], EmbSum, candidates, valid_candidate_idx, idx, valid_idx, idx_count);
                else if (FairT == 2) // already have experiments on that
                    rankLess(query_graph, edge_matrix, order[0], EmbSum, candidates, valid_candidate_idx, idx, valid_idx, idx_count);
                else if (FairT == 3)
                    rankMore(query_graph, edge_matrix, order[0], EmbSum, candidates, valid_candidate_idx, idx, valid_idx, idx_count);
                else if (FairT == 4)
                    rankFinal(query_graph, edge_matrix, order[0], EmbSum, candidates, valid_candidate_idx, idx, valid_idx, idx_count);
                else if (FairT == 5)
                    rankFinal1(query_graph, edge_matrix, order[0], EmbSum, candidates, valid_candidate_idx, idx, valid_idx, idx_count);

                valid_idx = valid_candidate_idx[cur_depth][idx[cur_depth]];
            }*/
            // from all the valid IDS

            VertexID u = order[cur_depth]; // u is current depth qid

            VertexID v = candidates[u][valid_idx]; // v is FCS[qid][vid] valid_idx-> position in CS

            size_t candidate = candidatesHC[u][valid_idx];
            if (candidate == 0 && cur_depth == 0)
            {
                calculateCell(query_graph, edge_matrix, candidates_count, candidatesHC, idToValues, count2, u, valid_idx);
                // calculateCellFN(query_graph, edge_matrix, candidates_count, candidatesHC, idToValues, count2, u, valid_idx,RQ);
                count2++;
            }
            // else if (candidate == 0 &&cur_depth != max_depth - 1)
            else if (candidate == 0 && cur_depth != max_depth - 1)
            { // calculateCell(query_graph, edge_matrix, candidates_count, candidatesHC, idToValues, count2, u, valid_idx);

                calculateCellFN(query_graph, edge_matrix, candidates_count, candidatesHC, idToValues, count2, u, valid_idx, RQ);
                count2++;
            }
            // candidate = candidatesHC[u][valid_idx];
            VN[u].clear();
            if (cur_depth != max_depth - 1) // change here
            {
                // VN[u] = idToValues[u][candidate];
                VN[u] = idToValues[u][candidatesHC[u][valid_idx]];

                if (cur_depth != 0)
                {
                    int ia = idx[cur_depth];
                    int io = 0;
                    VNTemp.clear();
                    while (ia < idx_count[cur_depth] && io < VN[u].size())
                    {

                        ui valid_idxN = valid_candidate_idx[cur_depth][ia];
                        if (valid_idxN == 10000000)
                        {
                            ia += 1;
                            continue;
                        }
                        if (VN[u][io] == valid_idxN)
                        {
                            VNTemp.push_back(VN[u][io]);
                            ia++;
                            io++;
                        }

                        else if (VN[u][io] > valid_idxN)
                            ia++;
                        else
                            io++;
                    }
                    VN[u] = VNTemp;
                }
            }

            // set embedings found to empty for start

            if (visited_vertices[v])
            { // if conflict on the last level then equivalence are match.

                VNTemp.clear();
                ui utemp = reverse_embedding[v];


                if (cur_depth != max_depth - 1)
                {
                    if (VN[u].size() > 0 && VN[utemp].size() > 0)
                    {
                        int ia = 0;
                        int io = 0;

                        while (ia < VN[u].size() && io < VN[utemp].size())
                        {
                            if (candidates[u][VN[u][ia]] == candidates[utemp][VN[utemp][io]])
                            {
                                VNTemp.push_back(VN[utemp][io]);
                                ia++;
                                io++;
                            }
                            else if (candidates[u][VN[u][ia]] < candidates[utemp][VN[utemp][io]])
                                ia++;
                            else
                                io++;
                        }
                    }
                    VN[utemp] = VNTemp;
                    VNTemp.clear();
                    // break;
                }
                else
                { // enumerate the results if it has equivalent nodes.
                    if (VN[utemp].size() > 1)
                    {
                        VNTemp.clear();
                        for (int ai = 0; ai < VN[utemp].size(); ai++)
                        {
                            if (!visited_vertices[candidates[utemp][VN[utemp][ai]]])
                                VNTemp.push_back(VN[utemp][ai]);
                        }
                                
                                    
                                
                        if (VNTemp.size() >= 1)
                        {
                            // heuristic there is a match so we count.
                            //  we can also exclude by just using VTEMP
                            // for utemp to count the correct number of matchings
                            //temp_sol[utemp].push_back(candidates[utemp][VNTemp[0]]);
                            int UNPMtemp = 0;
                            //UNPMtemp = EmbSum.insert(v).second;
                            temp_sol[cur_depth].push_back(v);
                            int ao = cur_depth - 1; //-1;
                            ui vqo = order[ao];
                            int UNPM2 = 0;
                            UNPM = 0; //
                            while (ao >= 0)
                            {  
                                Match_BA[ao] = true;
                                // EmbSum.insert(embedding[vqo]);
                                if(embedding[vqo]!=v){                                
                                temp_sol[ao].push_back(embedding[vqo]);
                                }else
                                {
                                    temp_sol[ao].push_back(candidates[vqo][VNTemp[0]]);
                                }
                                //&&ao!=max_depth-1
                                if (VN[vqo].size() > 0)
                                {
                                    for (int dd = 0; dd < VN[vqo].size(); dd++)
                                    {
                                        if (VN[vqo][dd] != 10000000)
                                        {                                                   
                                            if (nodeId[candidates[vqo][VN[vqo][dd]]]==0)                                            
                                            temp_sol[ao].push_back(candidates[vqo][VN[vqo][dd]]);
                                        }
                                    }
                                }

                                ao--;
                                vqo = order[ao];
                                //UNPM += UNPMtemp;
                            }
         
                            
                            for (int aa=0;aa<qsiz;aa++){
                                if (nodeId[temp_sol[aa][0]]==0){
                                    nodeId[temp_sol[aa][0]]=1;
                                    UNPM ++;
                                }
                            }

                            if (UNPM > pq.top().first)
                            {                                        
                                
                                if (heapcount >= ksize)
                                {
                                   auto [value, index] = pq.top();      
                                    pq.pop();
                                    pq.push({UNPM,index});
                                    addPos=index;
                                    
                                }else{
                                     pq.push({UNPM,heapcount});
                                     addPos=heapcount;            
                                }
                            for (int aa=0;aa<qsiz;aa++){
                                SolPos[addPos][aa]=temp_sol[aa][0];
                                //tempPos[aa]=temp_sol[aa][0];
                                    }
                                heapcount++;
                            }
                            //embedding_cnt += UNPM2 - 1;
                            int aa=0;
                            bool go=true;
                            int ii=1;
                            
                            while (aa<qsiz && pq.top().first < 1)
                            {   
                                if (temp_sol[aa].size()==1||temp_sol[aa].size()==ii){
                                aa++;
                                ii=1;
                                continue;
                            }
                            if(nodeId[temp_sol[aa][ii]]==0)
                                {
                                  nodeId[temp_sol[aa][ii]]=1;        
                                if (heapcount >= ksize)
                                {
                                   auto [value, index] = pq.top();      
                                    pq.pop();
                                    pq.push({1,index});
                                    addPos=index;                                    
                                }else{
                                     pq.push({1,heapcount});  
                                    addPos=heapcount;                               
                                }
                                    for (int ia=0;ia<qsiz;ia++){
                                        if(ia==aa){
                                          SolPos[addPos][aa]=temp_sol[aa][ii];  
                                        }
                                        SolPos[addPos][aa]=tempPos[aa];
                                    }
                                heapcount++;      
                            }
                            ii++;
                            }
                            for (auto& inner_vec : temp_sol) {
                                inner_vec.clear();  // Clears each inner vector
                            }

                            // reverse_embedding.erase(embedding[u]);
                            vec_failing_set[cur_depth].set();
                            vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];

                        }
                        else
                        { // update equivalent class of VN[utemp]
                            int ia = 0;
                            int io = 0;
                            VNTemp.clear();
                            while (ia < idx_count[cur_depth] && io < VN[utemp].size())
                            {

                                ui valid_idxN = valid_candidate_idx[cur_depth][ia];
                                if (valid_idxN == 10000000)
                                {
                                    ia += 1;
                                    continue;
                                }
                                if (VN[utemp][io] == valid_idxN)
                                {
                                    // VNTemp.push_back(VN[u][io]);
                                    VNTemp.push_back(VN[utemp][io]);
                                    ia++;
                                    io++;
                                }

                                else if (VN[utemp][io] > valid_idxN)
                                    ia++;
                                else
                                    io++;
                                // ia++;
                            }
                            VN[utemp] = VNTemp;
                            VNTemp.clear();
                        }
                    }
                }
                // }
                //}
                VN[u].clear();

                idx[cur_depth] += 1;
                vec_failing_set[cur_depth] = ancestors[u];
                vec_failing_set[cur_depth] |= ancestors[reverse_embedding[v]];
                vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];

                continue;
            }
            embd_vertices[v] = true;
            counter++;
            embedding[u] = v;             // set the embeding of depth qid to v
            idx_embedding[u] = valid_idx; // reverse the id of the embdeing to the position
            visited_vertices[v] = true;
            idx[cur_depth] += 1; // next element

            while (cur_depth < idx_count[cur_depth] && valid_candidate_idx[cur_depth][idx[cur_depth]] == 10000000)
            {
                idx[cur_depth] += 1;
            }

            reverse_embedding[v] = u;

            if (cur_depth == max_depth - 1)
            {
                embedding_cnt += 1;
                visited_vertices[v] = false;
                embd_vertices[v] = false;
                counter--;
                int UNPMtemp = 0;
                // UNPMtemp=EmbSum.insert(embedding[order[cur_depth]]).second;

                int ao = cur_depth; //-1;
                ui vqo = order[ao];
                int UNPM2 = 0;
                UNPM = 0; //
                while (ao >= 0)
                {
                    Match_BA[ao] = true;
                    // EmbSum.insert(embedding[vqo]);
                    //UNPMtemp = EmbSum.insert(embedding[vqo]).second;
                     temp_sol[ao].push_back(embedding[vqo]);

                    //&&ao!=max_depth-1
                    if (VN[vqo].size() > 0)
                    {
                        for (int dd = 0; dd < VN[vqo].size(); dd++)
                        {
                            if (VN[vqo][dd] != 10000000)
                            {  if (nodeId[candidates[vqo][VN[vqo][dd]]]==0)        
                                    temp_sol[ao].push_back(candidates[vqo][VN[vqo][dd]]);
                            }
                        }
                    }

                    // Match_BA[ao] = true;
                    ao--;
                    vqo = order[ao];
                }

                    for (int aa=0;aa<qsiz;aa++){
                                if (nodeId[temp_sol[aa][0]]==0){
                                    nodeId[temp_sol[aa][0]]=1;
                                    UNPM ++;
                                }
                            }

                            if (UNPM > pq.top().first)
                            {                                        
                                int addPos=0;
                                if (heapcount >= ksize)
                                {
                                   auto [value, index] = pq.top();      
                                    pq.pop();
                                    pq.push({UNPM,index});
                                    addPos=index;
                                    
                                }else{
                                     pq.push({UNPM,heapcount});
                                     addPos=heapcount;            
                                }
                            for (int aa=0;aa<qsiz;aa++){
                                SolPos[addPos][aa]=temp_sol[aa][0];
                                tempPos[aa]=temp_sol[aa][0];
                                    }
                                heapcount++;
                            }
                            //embedding_cnt += UNPM2 - 1;
                            int aa=0;
                            bool go=true;
                            int ii=1;
                            
                            while (aa<qsiz && pq.top().first < 1)
                            {   
                                if (temp_sol[aa].size()==1||temp_sol[aa].size()==ii){
                                aa++;
                                ii=1;
                                continue;
                            }
                            if(nodeId[temp_sol[aa][ii]]==0){        
                                nodeId[temp_sol[aa][ii]]=1;
                                if (heapcount >= ksize)
                                {
                                   auto [value, index] = pq.top();      
                                    pq.pop();
                                    pq.push({1,index});
                                    addPos=index;                                    
                                }else{
                                     pq.push({1,heapcount});  
                                    addPos=heapcount;                               
                                }
                                    for (int ia=0;ia<qsiz;ia++){
                                        if(ia==aa){
                                          SolPos[addPos][aa]=temp_sol[aa][ii];  
                                        }
                                        SolPos[addPos][aa]=tempPos[aa];
                                    }
                                heapcount++;      
                            }
                            ii++;
                            }
                            for (auto& inner_vec : temp_sol) {
                                inner_vec.clear();  // Clears each inner vector
                            }

                reverse_embedding.erase(embedding[u]);
                vec_failing_set[cur_depth].set();
                vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];


                if (embedding_cnt >= output_limit_num)
                {
                    goto EXIT;
                }
            }
            else
            {
                call_count += 1;
                cur_depth += 1;     // next depth
                idx[cur_depth] = 0; // set the element to 0
                // generateValidCandidateIndex(cur_depth, idx_embedding, idx_count, valid_candidate_idx, edge_matrix, bn,
                //                             bn_count, order, temp_buffer, candidates, candidates_count, query_graph);
                generateValidCandidateIndex(cur_depth, idx_embedding, idx_count, valid_candidate_idx, edge_matrix, bn,
                                            bn_count, order, temp_buffer);
                pruneCandidatesIndexBySymmetryBreakingC(cur_depth, embedding, order, candidates_count, candidates,
                                                       idx_count, valid_candidate_idx, ordered_constraints);

                // get the candidates for next depth
                if (idx_count[cur_depth] == 0)
                {
                    vec_failing_set[cur_depth - 1] = ancestors[order[cur_depth]];
                }
                else
                {
                    vec_failing_set[cur_depth - 1].reset();
                }
            }
        }

        cur_depth -= 1;

        if (cur_depth == max_depth - 1)
        {
            continue;
        }

        if (cur_depth < 0)
            break;
        else
        {


            // if(idx[cur_depth]>idx_count[cur_depth])
            // continue;
            // if(valid_candidate_idx[cur_depth][idx[cur_depth]]==10000000)
            // continue;

            // if (valid_candidate_idx[cur_depth][idx[cur_depth]]==10000000)  {
            //     idx[cur_depth]++;
            //     continue;
            // }

            ui vq = order[cur_depth];
            if (VN[vq].size() > 1)
            { // something to remove
                // if (idx[cur_depth]==0){
                //     cout<<"my logic is wrong"<<endl;
                // }
                int add = 0;
                for (int dd = 0; dd < VN[vq].size(); dd++)
                {

                    //  if (VN[vq][dd]==idx[cur_depth]-1){
                    //  cout<<"oups"<<endl;
                    //  continue;
                    //}
                    for (int kk = (idx[cur_depth]); kk < idx_count[cur_depth]; kk++)
                    {
                        ui valid_idxn = valid_candidate_idx[cur_depth][kk];
                        if (valid_idxn == 10000000)
                        {

                            continue;
                        }

                        if (valid_idxn == VN[vq][dd])
                        {
                            valid_candidate_idx[cur_depth][kk] = 10000000;
                            break;
                        }
                    }
                }
            }
            VN[vq].clear();

            VertexID u = order[cur_depth];
            reverse_embedding.erase(embedding[u]);
            if (cur_depth != 0)
            {
                if (!vec_failing_set[cur_depth].test(u))
                {
                    vec_failing_set[cur_depth - 1] = vec_failing_set[cur_depth];
                    idx[cur_depth] = idx_count[cur_depth];
                }
                else
                {
                    vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
                }
            }

            visited_vertices[embedding[u]] = false;
            embd_vertices[embedding[u]] = false;
            counter--;
            // if (idx[cur_depth] == idx_count[cur_depth])
            // continue;

            // if (Match_BA[vq] == false && idx[cur_depth] < idx_count[cur_depth])
            //{ positive and negative removal now. dont need to count positive symmetric
        }
    }

    EXIT:
    releaseBuffer(max_depth, idx, idx_count, embedding, idx_embedding, temp_buffer, valid_candidate_idx,
                  visited_vertices,
                  bn, bn_count);

    std::priority_queue<std::pair<int, int>> maxHeap;
    int topkcounter = 0;
    int greedysum = 0;
    unordered_set<ui> EmbGreedy;
    while (!pq.empty())
    {
        maxHeap.push(pq.top());
        pq.pop();
    }
    while (!maxHeap.empty())
    {   auto [value, index]= maxHeap.top();
        for (int aa=0;aa<qsiz;aa++){
            EmbGreedy.insert(SolPos[index][aa]);
        }
        maxHeap.pop();
        topkcounter++;
        switch (topkcounter)
        {
        case VALUE_5:
            s.topk[0] = EmbGreedy.size();
            break;
        case VALUE_10:
            s.topk[1] = EmbGreedy.size();
            break;
        case VALUE_50:
            s.topk[2] = EmbGreedy.size();
            break;
        case VALUE_100:
            s.topk[3] = EmbGreedy.size();
            break;
        case VALUE_250:
            s.topk[4] = EmbGreedy.size();
            break;
        case VALUE_500:
            s.topk[5] = EmbGreedy.size();
            break;
        case VALUE_750:
            s.topk[6] = EmbGreedy.size();
            break;
        case VALUE_1000:
            s.topk[7] = EmbGreedy.size();
            break;
        default:
            break;
        }
    }
    for (int i = 0; i < 8; i++)
    {
        if (s.topk[i] == 0)
            s.topk[i] = EmbGreedy.size();
    }
    int true_cand_sum = 0;
        int countSol=0;
    for (int i=0;i<data_graph->getVerticesCount();i++){
        countSol+=nodeId[i];
    }
    s.Can_embed = countSol;
    s.embedding_cnt = embedding_cnt;
    s.candidate_true_count_sum = true_cand_sum;
    return s;
}



enumResult
EvaluateQuery::DIVTOPKOPT(const Graph *data_graph, const Graph *query_graph,ui *&nodeId, Edges ***edge_matrix, ui **candidates,
                       ui *candidates_count, ui *order, size_t output_limit_num, size_t &call_count, int TimeL, int FairT, const std::unordered_map<VertexID, std::pair<std::set<VertexID>, std::set<VertexID>>> &ordered_constraints)
{
    //assumming that creating matrix can be done before we do not count the time.
    int qsiz = query_graph->getVerticesCount();
    int ksize = 10000;
    int SolPos[ksize+1][qsiz];
    std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, std::greater<>> pq;
    int tempPos[qsiz];
    int tempPos1[qsiz];
    auto start = std::chrono::high_resolution_clock::now();
    pq.push({0,0});
    int heapcount = 0;
    enumResult s;
    std::vector<std::vector<int>> temp_sol(qsiz);
    // Generate bn.
    int LastMatc=-1;
    ui **bn;
    ui *bn_count;
    generateBN(query_graph, order, bn, bn_count);
    //unordered_set<ui> EmbSum;
    bool NotifyALL = true;
    // Allocate the memory buffer.
    ui *idx;
    ui *idx_count;
    ui *embedding;
    ui *idx_embedding;
    ui *temp_buffer;
    ui **valid_candidate_idx;
    bool *visited_vertices;
    int addPos=0;


    int UNPM = 0;
    allocateBuffer(data_graph, query_graph, candidates_count, idx, idx_count, embedding, idx_embedding,
                   temp_buffer, valid_candidate_idx, visited_vertices);
    vector<ui> VN[query_graph->getVerticesCount()];
    vector<ui> VNTemp;
    

    int RQ[query_graph->getVerticesCount()];
    for (int i = 0; i < query_graph->getVerticesCount(); i++)
    {
        RQ[order[i]] = i; // depth RQ[u]=
    }
    bool *embd_vertices = new bool[data_graph->getVerticesCount()];
    std::fill(embd_vertices, embd_vertices + data_graph->getVerticesCount(), false);
    size_t embedding_cnt = 0;
    int cur_depth = 0;
    int max_depth = query_graph->getVerticesCount();
    VertexID start_vertex = order[0];
    bool Match_BA[max_depth] = {false};
    idx[cur_depth] = 0;                                    // set the starting position to 0 for the current depth
    idx_count[cur_depth] = candidates_count[start_vertex]; // set the size for current depth
    // all the nodes are possible candidates

    // end1 = std::chrono::high_resolution_clock::now();
    int count2 = 1;
    size_t **candidatesHC = NULL;
    unordered_map<size_t, vector<ui>> *idToValues;
    idToValues = new unordered_map<size_t, vector<ui>>[query_graph->getVerticesCount()];
    candidatesHC = new size_t *[qsiz];
    for (ui i = 0; i < qsiz; ++i)
    {
        candidatesHC[i] = new size_t[candidates_count[i]];
        memset(candidatesHC[i], 0, sizeof(size_t) * candidates_count[i]);
    }
    calculateCell(query_graph, edge_matrix, candidates_count, candidatesHC, idToValues, count2, start_vertex, 0);
    // calculateCellFN(query_graph, edge_matrix, candidates_count, candidatesHC, idToValues, count2, start_vertex, 0,RQ);
    count2++;
    VN[start_vertex] = idToValues[start_vertex][candidatesHC[start_vertex][0]];

    for (ui i = 0; i < idx_count[cur_depth]; ++i)
    {
        valid_candidate_idx[cur_depth][i] = i;
    }

    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> ancestors;
    computeAncestor(query_graph, bn, bn_count, order, ancestors);
    AncestorWithSymmetryBreaking(query_graph, ancestors, ordered_constraints);
    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> vec_failing_set(max_depth);
    std::unordered_map<VertexID, VertexID> reverse_embedding;
    reverse_embedding.reserve(MAXIMUM_QUERY_GRAPH_SIZE * 2);
    bool del = false;
    int k;
    VertexID utemp; // u is current depth qid
    VertexID vtemp;
    int counter = 0;
    double ens = 0;
    int ordCand[idx_count[0]] = {100000};
    TimeL = TimeL * 1000;
    while (true)
    {
        // current position smaller than size
        while (idx[cur_depth] < idx_count[cur_depth])
        {
            Match_BA[cur_depth] = false;
            auto end = std::chrono::high_resolution_clock::now();
            ens = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
            if (ens > TimeL)
            { // 1000 1 sec
                goto EXIT;
            }
            ui valid_idx = valid_candidate_idx[cur_depth][idx[cur_depth]]; // valid id
            if (valid_idx == 10000000)                                     //||
            {
                idx[cur_depth] += 1;
                continue;
            }
            /*
            if (cur_depth == 0 && idx[cur_depth] != 0 && EmbSum.size() > 0 && FairT > 0)
            {
                if (FairT == 1)
                    rankSimple(order[0], EmbSum, candidates, valid_candidate_idx, idx, valid_idx, idx_count);
                else if (FairT == 2) // already have experiments on that
                    rankLess(query_graph, edge_matrix, order[0], EmbSum, candidates, valid_candidate_idx, idx, valid_idx, idx_count);
                else if (FairT == 3)
                    rankMore(query_graph, edge_matrix, order[0], EmbSum, candidates, valid_candidate_idx, idx, valid_idx, idx_count);
                else if (FairT == 4)
                    rankFinal(query_graph, edge_matrix, order[0], EmbSum, candidates, valid_candidate_idx, idx, valid_idx, idx_count);
                else if (FairT == 5)
                    rankFinal1(query_graph, edge_matrix, order[0], EmbSum, candidates, valid_candidate_idx, idx, valid_idx, idx_count);

                valid_idx = valid_candidate_idx[cur_depth][idx[cur_depth]];
            }*/
            // from all the valid IDS
            if (cur_depth == 0 && idx[cur_depth] != 0  && FairT > 0)
            {
                if (FairT == 1)
                    rankSimple(order[0], nodeId, candidates, valid_candidate_idx, idx, valid_idx, idx_count);
                else if (FairT == 2) 
                    rankLess(query_graph, edge_matrix, order[0], nodeId, candidates, valid_candidate_idx, idx, valid_idx, idx_count);
                else if (FairT == 3)
                    rankMore(query_graph, edge_matrix, order[0], nodeId, candidates, valid_candidate_idx, idx, valid_idx, idx_count);
                else if (FairT == 4)
                    rankFinal(query_graph, edge_matrix, order[0], nodeId, candidates, valid_candidate_idx, idx, valid_idx, idx_count);
                else if (FairT == 5)
                    rankFinal1(query_graph, edge_matrix, order[0], nodeId, candidates, valid_candidate_idx, idx, valid_idx, idx_count);
                else if (FairT == 6)
                    rankRandom(query_graph, edge_matrix, order[0], nodeId, candidates, valid_candidate_idx, idx, valid_idx, idx_count);
                valid_idx = valid_candidate_idx[cur_depth][idx[cur_depth]];
            if (valid_idx == 10000000)                                     //||
            {
                idx[cur_depth] += 1;
                continue;
            }
            }
            VertexID u = order[cur_depth]; // u is current depth qid
            VertexID v = candidates[u][valid_idx]; // v is FCS[qid][vid] valid_idx-> position in CS

            size_t candidate = candidatesHC[u][valid_idx];
            if (candidate == 0 && cur_depth == 0)
            {
                calculateCell(query_graph, edge_matrix, candidates_count, candidatesHC, idToValues, count2, u, valid_idx);
                // calculateCellFN(query_graph, edge_matrix, candidates_count, candidatesHC, idToValues, count2, u, valid_idx,RQ);
                count2++;
            }
            // else if (candidate == 0 &&cur_depth != max_depth - 1)
            else if (candidate == 0 && cur_depth != max_depth - 1)
            { // calculateCell(query_graph, edge_matrix, candidates_count, candidatesHC, idToValues, count2, u, valid_idx);

                calculateCellFN(query_graph, edge_matrix, candidates_count, candidatesHC, idToValues, count2, u, valid_idx, RQ);
                count2++;
            }
            // candidate = candidatesHC[u][valid_idx];
            VN[u].clear();
            if (cur_depth != max_depth - 1) // change here
            {
                // VN[u] = idToValues[u][candidate];
                VN[u] = idToValues[u][candidatesHC[u][valid_idx]];

                if (cur_depth != 0)
                {
                    int ia = idx[cur_depth];
                    int io = 0;
                    VNTemp.clear();
                    while (ia < idx_count[cur_depth] && io < VN[u].size())
                    {

                        ui valid_idxN = valid_candidate_idx[cur_depth][ia];
                        if (valid_idxN == 10000000)
                        {
                            ia += 1;
                            continue;
                        }
                        if (VN[u][io] == valid_idxN)
                        {
                            VNTemp.push_back(VN[u][io]);
                            ia++;
                            io++;
                        }

                        else if (VN[u][io] > valid_idxN)
                            ia++;
                        else
                            io++;
                    }
                    VN[u] = VNTemp;
                }
            }

            // set embedings found to empty for start

            if (visited_vertices[v])
            { // if conflict on the last level then equivalence are match.
                VNTemp.clear();
                ui utemp = reverse_embedding[v];


                if (cur_depth != max_depth - 1)
                {
                    if (VN[u].size() > 0 && VN[utemp].size() > 0)
                    {
                        int ia = 0;
                        int io = 0;

                        while (ia < VN[u].size() && io < VN[utemp].size())
                        {
                            if (candidates[u][VN[u][ia]] == candidates[utemp][VN[utemp][io]])
                            {
                                VNTemp.push_back(VN[utemp][io]);
                                ia++;
                                io++;
                            }
                            else if (candidates[u][VN[u][ia]] < candidates[utemp][VN[utemp][io]])
                                ia++;
                            else
                                io++;
                        }
                    }
                    VN[utemp] = VNTemp;
                    VNTemp.clear();
                    // break;
                }
                else
                { // enumerate the results if it has equivalent nodes.
                    if (VN[utemp].size() > 1)
                    {
                        VNTemp.clear();
                        for (int ai = 0; ai < VN[utemp].size(); ai++)
                        {
                            if (!visited_vertices[candidates[utemp][VN[utemp][ai]]])
                                VNTemp.push_back(VN[utemp][ai]);
                        }
                        if (VNTemp.size() >= 1)
                        {
                            // heuristic there is a match so we count.
                            //  we can also exclude by just using VTEMP
                            // for utemp to count the correct number of matchings

                            int UNPMtemp = 0;
                            //UNPMtemp = EmbSum.insert(v).second;
                            temp_sol[cur_depth].push_back(v);
                             Match_BA[cur_depth] = true;
                            int ao = cur_depth-1 ; //-1;
                            ui vqo = order[ao];
                            int UNPM2 = 0;
                            UNPM = 0; //
                            while (ao >= 0)
                            {  
                                Match_BA[ao] = true;
                                // EmbSum.insert(embedding[vqo]);
                                if(embedding[vqo]!=v){                                
                                temp_sol[ao].push_back(embedding[vqo]);
                                }else
                                {
                                    temp_sol[ao].push_back(candidates[vqo][VNTemp[0]]);
                                }
                                //&&ao!=max_depth-1
                                if (VN[vqo].size() > 0)
                                {
                                    for (int dd = 0; dd < VN[vqo].size(); dd++)
                                    {
                                        if (VN[vqo][dd] != 10000000)
                                        {       
                                            if (nodeId[candidates[vqo][VN[vqo][dd]]]==0)                                                                                            
                                                temp_sol[ao].push_back(candidates[vqo][VN[vqo][dd]]);                                        
                                        }
                                    }
                                }
                                ao--;
                                vqo = order[ao];
                                //UNPM += UNPMtemp;
                            }
         
                            
                            for (int aa=0;aa<qsiz;aa++){
                                if (nodeId[temp_sol[aa][0]]==0){
                                    nodeId[temp_sol[aa][0]]=1;
                                    UNPM ++;
                                }
                                        
                                else if(temp_sol[aa].size()==1) 
                                    continue;  
                                else{
                                   int ja=temp_sol[aa].size()-1;
                                   //tempPos[aa]=temp_sol[aa][0];
                                   while(ja>0){
                                    
                                    if (nodeId[temp_sol[aa][ja]]==0){
                                        nodeId[temp_sol[aa][ja]]=1;
                                        temp_sol[aa][0]=temp_sol[aa][ja];
                                        temp_sol[aa].pop_back();
                                        UNPM ++;
                                        //tempPos[aa]=temp_sol[aa][0];
                                        break;
                                    }else{
                                        temp_sol[aa].pop_back();
                                        ja--;
                                    }

                                   }
                                }
                            }

                            if (UNPM > pq.top().first&&UNPM!=0)
                            {                                        
                                
                                if (heapcount >= ksize)
                                {
                                   auto [value, index] = pq.top();      
                                    pq.pop();
                                    pq.push({UNPM,index});
                                    addPos=index;
                                    
                                }else{
                                     pq.push({UNPM,heapcount});
                                     addPos=heapcount;            
                                }
                            for (int aa=0;aa<qsiz;aa++){
                                SolPos[addPos][aa]=temp_sol[aa][0];
                                //tempPos[aa]=temp_sol[aa][0];
                                    }
                                heapcount++;
                            }
                            //embedding_cnt += UNPM2 - 1;
                            int aa=0;
                            bool go=true;
                            int ii=1;
                            //reverse greedy
                            bool OG=true;
                            if(idx[cur_depth] == idx_count[cur_depth]-1)
                            while (OG)
                            {  
                                 OG=false;
                                UNPM=0;
                                for (int aa=0;aa<qsiz;aa++){
                                if (temp_sol[aa].size()==1){
                                continue;
                                    }
                                int ja=temp_sol[aa].size()-1;
                                while(ja>0){
                                    
                                    if (nodeId[temp_sol[aa][ja]]==0){
                                        nodeId[temp_sol[aa][ja]]=1;
                                        temp_sol[aa][0]=temp_sol[aa][ja];
                                        temp_sol[aa].pop_back();
                                        //tempPos[aa]=temp_sol[aa][0];
                                        UNPM++;
                                        break;
                                    }else{
                                        temp_sol[aa].pop_back();
                                        ja--;
                                    }
                                   }
                            }
                            if (UNPM > pq.top().first&&UNPM!=0)
                            {                                        
                                OG=true;
                                if (heapcount >= ksize)
                                {
                                   auto [value, index] = pq.top();      
                                    pq.pop();
                                    pq.push({UNPM,index});
                                    addPos=index;
                                    
                                }else{
                                     pq.push({UNPM,heapcount});
                                     addPos=heapcount;            
                                }
                            for (int aa=0;aa<qsiz;aa++){
                                SolPos[addPos][aa]=temp_sol[aa][0];
                                //tempPos[aa]=temp_sol[aa][0];
                                    }
                                heapcount++;
                            }
                            }
                            
                            for (auto& inner_vec : temp_sol) {
                                inner_vec.clear();  // Clears each inner vector
                            }

                            // reverse_embedding.erase(embedding[u]);
                            vec_failing_set[cur_depth].set();
                            vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
                        }
                        else
                        { // update equivalent class of VN[utemp]
                            int ia = 0;
                            int io = 0;
                            VNTemp.clear();
                            while (ia < idx_count[cur_depth] && io < VN[utemp].size())
                            {

                                ui valid_idxN = valid_candidate_idx[cur_depth][ia];
                                if (valid_idxN == 10000000)
                                {
                                    ia += 1;
                                    continue;
                                }
                                if (VN[utemp][io] == valid_idxN)
                                {
                                    // VNTemp.push_back(VN[u][io]);
                                    VNTemp.push_back(VN[utemp][io]);
                                    ia++;
                                    io++;
                                }

                                else if (VN[utemp][io] > valid_idxN)
                                    ia++;
                                else
                                    io++;
                                // ia++;
                            }
                            VN[utemp] = VNTemp;
                            VNTemp.clear();
                        }
                    }
                }
                // }
                //}
                VN[u].clear();

                idx[cur_depth] += 1;

                vec_failing_set[cur_depth] = ancestors[u];
                vec_failing_set[cur_depth] |= ancestors[reverse_embedding[v]];
                vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];

                continue;
            }
            embd_vertices[v] = true;
            counter++;
            embedding[u] = v;             // set the embeding of depth qid to v
            idx_embedding[u] = valid_idx; // reverse the id of the embdeing to the position
            visited_vertices[v] = true;
            idx[cur_depth] += 1; // next element

            while (cur_depth < idx_count[cur_depth] && valid_candidate_idx[cur_depth][idx[cur_depth]] == 10000000)
            {
                idx[cur_depth] += 1;
            }

            reverse_embedding[v] = u;


            if (cur_depth == max_depth - 1)
            {
                embedding_cnt += 1;
                visited_vertices[v] = false;
                embd_vertices[v] = false;
                counter--;
                int UNPMtemp = 0;
                // UNPMtemp=EmbSum.insert(embedding[order[cur_depth]]).second;

                int ao = cur_depth; //-1;
                ui vqo = order[ao];
                int UNPM2 = 0;
                UNPM = 0; //
                while (ao >= 0)
                            {  
                                Match_BA[ao] = true;
                                // EmbSum.insert(embedding[vqo]);                              
                                temp_sol[ao].push_back(embedding[vqo]);
                                //&&ao!=max_depth-1
                                if (VN[vqo].size() > 0)
                                {
                                    for (int dd = 0; dd < VN[vqo].size(); dd++)
                                    {
                                        if (VN[vqo][dd] != 10000000)
                                        {       
                                            if (nodeId[candidates[vqo][VN[vqo][dd]]]==0)                                            
                                                temp_sol[ao].push_back(candidates[vqo][VN[vqo][dd]]);                                        
                                        }
                                    }
                                }

                                ao--;
                                vqo = order[ao];
                                //UNPM += UNPMtemp;
                            }
         
                            
                            for (int aa=0;aa<qsiz;aa++){
                                if(nodeId[temp_sol[aa][0]]==0){                                    
                                    nodeId[temp_sol[aa][0]]=1;
                                    UNPM ++;
                                    tempPos[aa]=temp_sol[aa][0];
                                }
                                        
                                else if(temp_sol[aa].size()==1) 
                                    continue;  
                                else{
                                   int ja=temp_sol[aa].size()-1;
                                   while(ja>0){
                                    
                                    if (nodeId[temp_sol[aa][ja]]==0){
                                        nodeId[temp_sol[aa][ja]]=1;
                                        temp_sol[aa][0]=temp_sol[aa][ja];
                                        temp_sol[aa].pop_back();
                                        UNPM ++;
                                        break;
                                    }else{
                                        temp_sol[aa].pop_back();
                                        ja--;
                                    }

                                   }
                                }
                            }

                            if (UNPM > pq.top().first&&UNPM!=0)
                            {                                        
                                
                                if (heapcount >= ksize)
                                {
                                   auto [value, index] = pq.top();      
                                    pq.pop();
                                    pq.push({UNPM,index});
                                    addPos=index;
                                    
                                }else{
                                     pq.push({UNPM,heapcount});
                                     addPos=heapcount;            
                                }
                            for (int aa=0;aa<qsiz;aa++){
                                SolPos[addPos][aa]=temp_sol[aa][0];
                                //tempPos[aa]=temp_sol[aa][0];
                                    }
                                heapcount++;
                            }
                            //embedding_cnt += UNPM2 - 1;
                            int aa=0;
                            bool go=true;
                            int ii=1;
                            //reverse greedy
                            bool OG=true;
                            //if (true)
                            if(idx[cur_depth] == idx_count[cur_depth])
                            while (OG)
                            {   OG=false;
                                UNPM=0;
                                for (int aa=0;aa<qsiz;aa++){
                                if (temp_sol[aa].size()==1){
                                continue;
                                    }
                                int ja=temp_sol[aa].size()-1;
                                while(ja>0){
                                    
                                    if (nodeId[temp_sol[aa][ja]]==0){
                                        nodeId[temp_sol[aa][ja]]=1;
                                        temp_sol[aa][0]=temp_sol[aa][ja];
                                        temp_sol[aa].pop_back();
                                        UNPM++;
                                        break;
                                    }else{
                                        temp_sol[aa].pop_back();
                                        ja--;
                                    }
                                   }
                            }
                            if (UNPM > pq.top().first&&UNPM!=0)
                            {                                        
                                OG=true;
                                if (heapcount >= ksize)
                                {
                                   auto [value, index] = pq.top();      
                                    pq.pop();
                                    pq.push({UNPM,index});
                                    addPos=index;
                                    
                                }else{
                                     pq.push({UNPM,heapcount});
                                     addPos=heapcount;            
                                }
                            for (int aa=0;aa<qsiz;aa++){
                                SolPos[addPos][aa]=temp_sol[aa][0];
                                //tempPos[aa]=temp_sol[aa][0];
                                    }
                                heapcount++;
                            }

                            }
                            for (auto& inner_vec : temp_sol) {
                                inner_vec.clear();  // Clears each inner vector
                            }

                reverse_embedding.erase(embedding[u]);
                vec_failing_set[cur_depth].set();
                vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];


                if (embedding_cnt >= output_limit_num)
                {
                    goto EXIT;
                }
            }
            else
            {
                call_count += 1;
                cur_depth += 1;     // next depth
                idx[cur_depth] = 0; // set the element to 0
                // generateValidCandidateIndex(cur_depth, idx_embedding, idx_count, valid_candidate_idx, edge_matrix, bn,
                //                             bn_count, order, temp_buffer, candidates, candidates_count, query_graph);
                generateValidCandidateIndex(cur_depth, idx_embedding, idx_count, valid_candidate_idx, edge_matrix, bn,
                                            bn_count, order, temp_buffer);
                pruneCandidatesIndexBySymmetryBreakingC(cur_depth, embedding, order, candidates_count, candidates,
                                                       idx_count, valid_candidate_idx, ordered_constraints);

                // get the candidates for next depth

                if (idx_count[cur_depth] == 0)
                {
                    vec_failing_set[cur_depth - 1] = ancestors[order[cur_depth]];
                }
                else
                {
                    vec_failing_set[cur_depth - 1].reset();
                }

            }
        }

        cur_depth -= 1;
        LastMatc=-1;
        if (cur_depth == max_depth - 1)
        {
            continue;
        }

        if (cur_depth < 0)
            break;
        else
        {

            ui vq = order[cur_depth];
            if (VN[vq].size() > 1)
            { 
                int add = 0;
                for (int dd = 0; dd < VN[vq].size(); dd++)
                {

                    for (int kk = (idx[cur_depth]); kk < idx_count[cur_depth]; kk++)
                    {
                        ui valid_idxn = valid_candidate_idx[cur_depth][kk];
                        if (valid_idxn == 10000000)
                        {

                            continue;
                        }

                        if (valid_idxn == VN[vq][dd])
                        {
                            valid_candidate_idx[cur_depth][kk] = 10000000;
                            if(Match_BA[cur_depth]==true &&nodeId[candidates[order[cur_depth]][valid_idxn]]==0){
                                nodeId[candidates[order[cur_depth]][valid_idxn]]=1;
                                for (int mm=0;mm<max_depth;mm++){
                                    if (mm<cur_depth)
                                    temp_sol[mm][0]=embedding[order[mm]];
                                    else
                                    temp_sol[mm][0]=candidates[order[cur_depth]][valid_idxn];
                                }
                             if (1 > pq.top().first)
                            {                                        
                                
                                if (heapcount >= ksize)
                                {
                                   auto [value, index] = pq.top();      
                                    pq.pop();
                                    pq.push({1,index});
                                    addPos=index;                                    
                                }else{
                                     pq.push({1,heapcount});
                                     addPos=heapcount;            
                                }
                            for (int aa=0;aa<qsiz;aa++){
                                SolPos[addPos][aa]=temp_sol[aa][0];
                                //tempPos[aa]=temp_sol[aa][0];
                                    }
                                heapcount++;
                            }       
                            }
                            break;
                        }
                    }
                }
            }
            VN[vq].clear();

            VertexID u = order[cur_depth];

            reverse_embedding.erase(embedding[u]);
            if (cur_depth != 0)
            {
                if (!vec_failing_set[cur_depth].test(u))
                {
                    vec_failing_set[cur_depth - 1] = vec_failing_set[cur_depth];
                    idx[cur_depth] = idx_count[cur_depth];
                }
                else
                {
                    vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
                }
            }

            visited_vertices[embedding[u]] = false;
            embd_vertices[embedding[u]] = false;
            counter--;
            // if (idx[cur_depth] == idx_count[cur_depth])
            // continue;

            // if (Match_BA[vq] == false && idx[cur_depth] < idx_count[cur_depth])
            //{ positive and negative removal now. dont need to count positive symmetric
        }
    }

    EXIT:
    releaseBuffer(max_depth, idx, idx_count, embedding, idx_embedding, temp_buffer, valid_candidate_idx,
                  visited_vertices,
                  bn, bn_count);


    std::priority_queue<std::pair<int, int>> maxHeap;
    int topkcounter = 0;
    int greedysum = 0;
    unordered_set<ui> EmbGreedy;
    while (!pq.empty())
    {   if(pq.top().first!=0)
        maxHeap.push(pq.top());
        pq.pop();
    }


    while (!maxHeap.empty())
    {   auto [value, index]= maxHeap.top();
        for (int aa=0;aa<qsiz;aa++){
            EmbGreedy.insert(SolPos[index][aa]);
        }
        maxHeap.pop();
        topkcounter++;
        switch (topkcounter)
        {

        case VALUE_10:
            s.topk[0] = EmbGreedy.size();
            break;
        case VALUE_50:
            s.topk[1] = EmbGreedy.size();
            break;
        case VALUE_100:
            s.topk[2] = EmbGreedy.size();
            break;
        case VALUE_500:
            s.topk[3] = EmbGreedy.size();
            break;
        case VALUE_1000:
            s.topk[4] = EmbGreedy.size();
            break;
        case VALUE_2500:
            s.topk[5] = EmbGreedy.size();
            break;    
        case VALUE_5000:
            s.topk[6] = EmbGreedy.size();
            break;
        case VALUE_10000:
            s.topk[7] = EmbGreedy.size();
            break;    
        default:
            break;
        }
    }
    for (int i = 0; i < 8; i++)
    {
        if (s.topk[i] == 0)
            s.topk[i] = EmbGreedy.size();
    }
    int true_cand_sum = 0;
    int countSol=0;
    for (int i=0;i<data_graph->getVerticesCount();i++){
        countSol+=nodeId[i];
    }
    s.Can_embed = countSol;
    //EmbSum.clear();
    s.embedding_cnt = embedding_cnt;
    s.candidate_true_count_sum = true_cand_sum;
    return s;
}

enumResult
EvaluateQuery::DIVTOPKOPTSQ(const Graph *data_graph, const Graph *query_graph,ui *&nodeId, Edges ***edge_matrix, ui **candidates,
                       ui *candidates_count, ui *order, size_t output_limit_num, size_t &call_count, int TimeL, int FairT, const std::unordered_map<VertexID, std::pair<std::set<VertexID>, std::set<VertexID>>> &ordered_constraints)
{
    //assumming that creating matrix can be done before we do not count the time.
    int qsiz = query_graph->getVerticesCount();
    int ksize = 10000;
    int SolPos[ksize+1][qsiz];
    std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, std::greater<>> pq;
    int tempPos[qsiz];
    int tempPos1[qsiz];
    auto start = std::chrono::high_resolution_clock::now();
    pq.push({0,0});
    int heapcount = 0;
    enumResult s;
    std::vector<std::vector<int>> temp_sol(qsiz);
    // Generate bn.
    int LastMatc=-1;
    ui **bn;
    ui *bn_count;
    generateBN(query_graph, order, bn, bn_count);
    //unordered_set<ui> EmbSum;
    bool NotifyALL = true;
    // Allocate the memory buffer.
    ui *idx;
    ui *idx_count;
    ui *embedding;
    ui *idx_embedding;
    ui *temp_buffer;
    ui **valid_candidate_idx;
    bool *visited_vertices;
    int addPos=0;


    int UNPM = 0;
    allocateBuffer(data_graph, query_graph, candidates_count, idx, idx_count, embedding, idx_embedding,
                   temp_buffer, valid_candidate_idx, visited_vertices);
    vector<ui> VN[query_graph->getVerticesCount()];
    vector<ui> VNTemp;
    

    int RQ[query_graph->getVerticesCount()];
    for (int i = 0; i < query_graph->getVerticesCount(); i++)
    {
        RQ[order[i]] = i; // depth RQ[u]=
    }
    bool *embd_vertices = new bool[data_graph->getVerticesCount()];
    std::fill(embd_vertices, embd_vertices + data_graph->getVerticesCount(), false);
    size_t embedding_cnt = 0;
    int cur_depth = 0;
    int max_depth = query_graph->getVerticesCount();
    VertexID start_vertex = order[0];
    bool Match_BA[max_depth] = {false};
    idx[cur_depth] = 0;                                    // set the starting position to 0 for the current depth
    idx_count[cur_depth] = candidates_count[start_vertex]; // set the size for current depth
    // all the nodes are possible candidates

    // end1 = std::chrono::high_resolution_clock::now();
    int count2 = 1;
    size_t **candidatesHC = NULL;
    unordered_map<size_t, vector<ui>> *idToValues;
    idToValues = new unordered_map<size_t, vector<ui>>[query_graph->getVerticesCount()];
    candidatesHC = new size_t *[qsiz];
    for (ui i = 0; i < qsiz; ++i)
    {
        candidatesHC[i] = new size_t[candidates_count[i]];
        memset(candidatesHC[i], 0, sizeof(size_t) * candidates_count[i]);
    }
    calculateCell(query_graph, edge_matrix, candidates_count, candidatesHC, idToValues, count2, start_vertex, 0);
    // calculateCellFN(query_graph, edge_matrix, candidates_count, candidatesHC, idToValues, count2, start_vertex, 0,RQ);
    count2++;
    VN[start_vertex] = idToValues[start_vertex][candidatesHC[start_vertex][0]];

    for (ui i = 0; i < idx_count[cur_depth]; ++i)
    {
        valid_candidate_idx[cur_depth][i] = i;
    }

    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> ancestors;
    computeAncestor(query_graph, bn, bn_count, order, ancestors);
    AncestorWithSymmetryBreaking(query_graph, ancestors, ordered_constraints);
    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> vec_failing_set(max_depth);
    std::unordered_map<VertexID, VertexID> reverse_embedding;
    reverse_embedding.reserve(MAXIMUM_QUERY_GRAPH_SIZE * 2);
    bool del = false;
    int k;
    VertexID utemp; // u is current depth qid
    VertexID vtemp;
    int counter = 0;
    double ens = 0;
    int ordCand[idx_count[0]] = {100000};
    TimeL = TimeL * 1000;

    while (true)
    {
        while (idx[cur_depth] < idx_count[cur_depth])
        {
            Match_BA[cur_depth] = false;
            ui valid_idx = valid_candidate_idx[cur_depth][idx[cur_depth]];
            VertexID u = order[cur_depth];
            VertexID v = candidates[u][valid_idx];
            auto end = std::chrono::high_resolution_clock::now();
            ens = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
            if (ens > TimeL)
            { // 1000 1 sec
                goto EXIT;
            }
            if (visited_vertices[v])
            {
                idx[cur_depth] += 1;
                vec_failing_set[cur_depth] = ancestors[u];
                vec_failing_set[cur_depth] |= ancestors[reverse_embedding[v]];
                vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
                continue;
            }
            if(nodeId[v]==1){
                idx[cur_depth] += 1;
                vec_failing_set[cur_depth] = ancestors[u];
                //vec_failing_set[cur_depth] |= ancestors[reverse_embedding[v]];
                //vec_failing_set[cur_depth] |= ancestors[order[i]];
                if(cur_depth - 1>0)
                vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
                continue;
            }

            embedding[u] = v;
            idx_embedding[u] = valid_idx;
            visited_vertices[v] = true;
            idx[cur_depth] += 1;

            reverse_embedding[v] = u;
            if (cur_depth == max_depth - 1)
            {   
                UNPM = 0;
                if(nodeId[embedding[order[cur_depth]]]==0){
                    UNPM++;
                    nodeId[embedding[order[cur_depth]]]=1;
                }
                embedding_cnt += 1;
                visited_vertices[v] = false;
                int ao = cur_depth; //-1;
                tempPos[ao]=embedding[order[cur_depth]];
                ui vqo = order[ao];
                nodeId[v]=1;
                while (ao > 0)
                {                
                    if(nodeId[embedding[vqo]]==0){
                    UNPM++;
                    nodeId[embedding[vqo]]=1;
                }
                //vec_failing_set[ao].set();
                //vec_failing_set[ao - 1] |= vec_failing_set[ao];
                //reverse_embedding.erase(embedding[embedding[vqo]]);
                   tempPos[ao]=embedding[vqo];
                    idx[ao] = idx_count[ao];
                    // Match_BA[ao] = true;
                    ao--;
                    vqo = order[ao];
                }
                tempPos[ao]=embedding[vqo];
                if(nodeId[embedding[vqo]]==0){
                    UNPM++;
                    nodeId[embedding[vqo]]=1;
                    }

                if (UNPM > pq.top().first)
                            {                                       
                                int addPos=0;
                                if (heapcount >= ksize)
                                {
                                   auto [value, index] = pq.top();      
                                    pq.pop();
                                    pq.push({UNPM,index});
                                    addPos=index;
                                    
                                }else{
                                     pq.push({UNPM,heapcount});
                                     addPos=heapcount;            
                                }
                            for (int aa=0;aa<qsiz;aa++){
                                SolPos[addPos][aa]=tempPos[aa];
                                    }
                                heapcount++;
                            }

                reverse_embedding.erase(embedding[u]);
                vec_failing_set[cur_depth].set();
                vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
                if (embedding_cnt >= output_limit_num)
                {
                    goto EXIT;
                }
            }
            else
            {
                call_count += 1;
                cur_depth += 1;

                idx[cur_depth] = 0;
                generateValidCandidateIndex(cur_depth, idx_embedding, idx_count, valid_candidate_idx, edge_matrix, bn,
                                            bn_count, order, temp_buffer);

                
                if (idx_count[cur_depth] == 0)
                {
                    vec_failing_set[cur_depth - 1] = ancestors[order[cur_depth]];
                }
                else
                {
                    vec_failing_set[cur_depth - 1].reset();
                }
            }
        }

        cur_depth -= 1;
        if (cur_depth < 0)
            break;
        else
        {
            VertexID u = order[cur_depth];
            reverse_embedding.erase(embedding[u]);
            if (cur_depth != 0)
            {
                if (!vec_failing_set[cur_depth].test(u))
                {
                    vec_failing_set[cur_depth - 1] = vec_failing_set[cur_depth];
                    idx[cur_depth] = idx_count[cur_depth];
                }
                else
                {
                    vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
                }
            }
            visited_vertices[embedding[u]] = false;

        }
    }




    allocateBufferLM(data_graph, query_graph, candidates_count, idx, idx_count, embedding, idx_embedding,
                   temp_buffer, valid_candidate_idx, visited_vertices);
    ancestors.clear();
    computeAncestor(query_graph, bn, bn_count, order, ancestors);
    AncestorWithSymmetryBreaking(query_graph, ancestors, ordered_constraints);
     vec_failing_set.clear();
     reverse_embedding.clear();
    cur_depth = 0;
    idx[cur_depth] = 0;       
    start_vertex = order[0];                             // set the starting position to 0 for the current depth
    idx_count[cur_depth] = candidates_count[start_vertex]; // set the size for current depth



    while (true)
    {
        // current position smaller than size
        while (idx[cur_depth] < idx_count[cur_depth])
        {
            Match_BA[cur_depth] = false;
            auto end = std::chrono::high_resolution_clock::now();
            ens = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
            if (ens > TimeL)
            { // 1000 1 sec
                goto EXIT;
            }
            ui valid_idx = valid_candidate_idx[cur_depth][idx[cur_depth]]; // valid id
            if (valid_idx == 10000000)                                     //||
            {
                idx[cur_depth] += 1;
                continue;
            }
            /*
            if (cur_depth == 0 && idx[cur_depth] != 0 && EmbSum.size() > 0 && FairT > 0)
            {
                if (FairT == 1)
                    rankSimple(order[0], EmbSum, candidates, valid_candidate_idx, idx, valid_idx, idx_count);
                else if (FairT == 2) // already have experiments on that
                    rankLess(query_graph, edge_matrix, order[0], EmbSum, candidates, valid_candidate_idx, idx, valid_idx, idx_count);
                else if (FairT == 3)
                    rankMore(query_graph, edge_matrix, order[0], EmbSum, candidates, valid_candidate_idx, idx, valid_idx, idx_count);
                else if (FairT == 4)
                    rankFinal(query_graph, edge_matrix, order[0], EmbSum, candidates, valid_candidate_idx, idx, valid_idx, idx_count);
                else if (FairT == 5)
                    rankFinal1(query_graph, edge_matrix, order[0], EmbSum, candidates, valid_candidate_idx, idx, valid_idx, idx_count);

                valid_idx = valid_candidate_idx[cur_depth][idx[cur_depth]];
            }*/
            // from all the valid IDS
            if (cur_depth == 0 && idx[cur_depth] != 0  && FairT > 0)
            {
                if (FairT == 1)
                    rankSimple(order[0], nodeId, candidates, valid_candidate_idx, idx, valid_idx, idx_count);
                else if (FairT == 2) 
                    rankLess(query_graph, edge_matrix, order[0], nodeId, candidates, valid_candidate_idx, idx, valid_idx, idx_count);
                else if (FairT == 3)
                    rankMore(query_graph, edge_matrix, order[0], nodeId, candidates, valid_candidate_idx, idx, valid_idx, idx_count);
                else if (FairT == 4)
                    rankFinal(query_graph, edge_matrix, order[0], nodeId, candidates, valid_candidate_idx, idx, valid_idx, idx_count);
                else if (FairT == 5)
                    rankFinal1(query_graph, edge_matrix, order[0], nodeId, candidates, valid_candidate_idx, idx, valid_idx, idx_count);
                else if (FairT == 6)
                    rankRandom(query_graph, edge_matrix, order[0], nodeId, candidates, valid_candidate_idx, idx, valid_idx, idx_count);
                valid_idx = valid_candidate_idx[cur_depth][idx[cur_depth]];
            if (valid_idx == 10000000)                                     //||
            {
                idx[cur_depth] += 1;
                continue;
            }
            }
            VertexID u = order[cur_depth]; // u is current depth qid
            VertexID v = candidates[u][valid_idx]; // v is FCS[qid][vid] valid_idx-> position in CS

            size_t candidate = candidatesHC[u][valid_idx];
            if (candidate == 0 && cur_depth == 0)
            {
                calculateCell(query_graph, edge_matrix, candidates_count, candidatesHC, idToValues, count2, u, valid_idx);
                // calculateCellFN(query_graph, edge_matrix, candidates_count, candidatesHC, idToValues, count2, u, valid_idx,RQ);
                count2++;
            }
            // else if (candidate == 0 &&cur_depth != max_depth - 1)
            else if (candidate == 0 && cur_depth != max_depth - 1)
            { // calculateCell(query_graph, edge_matrix, candidates_count, candidatesHC, idToValues, count2, u, valid_idx);

                calculateCellFN(query_graph, edge_matrix, candidates_count, candidatesHC, idToValues, count2, u, valid_idx, RQ);
                count2++;
            }
            // candidate = candidatesHC[u][valid_idx];
            VN[u].clear();
            if (cur_depth != max_depth - 1) // change here
            {
                // VN[u] = idToValues[u][candidate];
                VN[u] = idToValues[u][candidatesHC[u][valid_idx]];

                if (cur_depth != 0)
                {
                    int ia = idx[cur_depth];
                    int io = 0;
                    VNTemp.clear();
                    while (ia < idx_count[cur_depth] && io < VN[u].size())
                    {

                        ui valid_idxN = valid_candidate_idx[cur_depth][ia];
                        if (valid_idxN == 10000000)
                        {
                            ia += 1;
                            continue;
                        }
                        if (VN[u][io] == valid_idxN)
                        {
                            VNTemp.push_back(VN[u][io]);
                            ia++;
                            io++;
                        }

                        else if (VN[u][io] > valid_idxN)
                            ia++;
                        else
                            io++;
                    }
                    VN[u] = VNTemp;
                }
            }

            // set embedings found to empty for start

            if (visited_vertices[v])
            { // if conflict on the last level then equivalence are match.
                VNTemp.clear();
                ui utemp = reverse_embedding[v];


                if (cur_depth != max_depth - 1)
                {
                    if (VN[u].size() > 0 && VN[utemp].size() > 0)
                    {
                        int ia = 0;
                        int io = 0;

                        while (ia < VN[u].size() && io < VN[utemp].size())
                        {
                            if (candidates[u][VN[u][ia]] == candidates[utemp][VN[utemp][io]])
                            {
                                VNTemp.push_back(VN[utemp][io]);
                                ia++;
                                io++;
                            }
                            else if (candidates[u][VN[u][ia]] < candidates[utemp][VN[utemp][io]])
                                ia++;
                            else
                                io++;
                        }
                    }
                    VN[utemp] = VNTemp;
                    VNTemp.clear();
                    // break;
                }
                else
                { // enumerate the results if it has equivalent nodes.
                    if (VN[utemp].size() > 1)
                    {
                        VNTemp.clear();
                        for (int ai = 0; ai < VN[utemp].size(); ai++)
                        {
                            if (!visited_vertices[candidates[utemp][VN[utemp][ai]]])
                                VNTemp.push_back(VN[utemp][ai]);
                        }
                        if (VNTemp.size() >= 1)
                        {
                            // heuristic there is a match so we count.
                            //  we can also exclude by just using VTEMP
                            // for utemp to count the correct number of matchings

                            int UNPMtemp = 0;
                            //UNPMtemp = EmbSum.insert(v).second;
                            temp_sol[cur_depth].push_back(v);
                            int ao = cur_depth - 1; //-1;
                            ui vqo = order[ao];
                            int UNPM2 = 0;
                            UNPM = 0; //
                            while (ao >= 0)
                            {  
                                Match_BA[ao] = true;
                                // EmbSum.insert(embedding[vqo]);
                                if(embedding[vqo]!=v){                                
                                temp_sol[ao].push_back(embedding[vqo]);
                                }else
                                {
                                    temp_sol[ao].push_back(candidates[vqo][VNTemp[0]]);
                                }
                                //&&ao!=max_depth-1
                                if (VN[vqo].size() > 0)
                                {
                                    for (int dd = 0; dd < VN[vqo].size(); dd++)
                                    {
                                        if (VN[vqo][dd] != 10000000)
                                        {       
                                            if (nodeId[candidates[vqo][VN[vqo][dd]]]==0)                                                                                            
                                                temp_sol[ao].push_back(candidates[vqo][VN[vqo][dd]]);                                        
                                        }
                                    }
                                }
                                ao--;
                                vqo = order[ao];
                                //UNPM += UNPMtemp;
                            }
         
                            
                            for (int aa=0;aa<qsiz;aa++){
                                if (nodeId[temp_sol[aa][0]]==0){
                                    nodeId[temp_sol[aa][0]]=1;
                                    UNPM ++;
                                Match_BA[aa]=true;
                                }
                                        
                                else if(temp_sol[aa].size()==1) 
                                    continue;  
                                else{
                                   int ja=temp_sol[aa].size()-1;
                                   //tempPos[aa]=temp_sol[aa][0];
                                   while(ja>0){
                                    
                                    if (nodeId[temp_sol[aa][ja]]==0){
                                        nodeId[temp_sol[aa][ja]]=1;
                                        temp_sol[aa][0]=temp_sol[aa][ja];
                                        temp_sol[aa].pop_back();
                                        UNPM ++;
                                        //tempPos[aa]=temp_sol[aa][0];
                                        break;
                                    }else{
                                        temp_sol[aa].pop_back();
                                        ja--;
                                    }

                                   }
                                }
                            }

                            if (UNPM > pq.top().first&&UNPM!=0)
                            {                                        
                                
                                if (heapcount >= ksize)
                                {
                                   auto [value, index] = pq.top();      
                                    pq.pop();
                                    pq.push({UNPM,index});
                                    addPos=index;
                                    
                                }else{
                                     pq.push({UNPM,heapcount});
                                     addPos=heapcount;            
                                }
                            for (int aa=0;aa<qsiz;aa++){
                                SolPos[addPos][aa]=temp_sol[aa][0];
                                //tempPos[aa]=temp_sol[aa][0];
                                    }
                                heapcount++;
                            }
                            //embedding_cnt += UNPM2 - 1;
                            int aa=0;
                            bool go=true;
                            int ii=1;
                            //reverse greedy
                            bool OG=true;
                            if(idx[cur_depth] == idx_count[cur_depth]-1)
                            while (OG)
                            {  
                                 OG=false;
                                UNPM=0;
                                for (int aa=0;aa<qsiz;aa++){
                                if (temp_sol[aa].size()==1){
                                continue;
                                    }
                                int ja=temp_sol[aa].size()-1;
                                while(ja>0){
                                    
                                    if (nodeId[temp_sol[aa][ja]]==0){
                                        nodeId[temp_sol[aa][ja]]=1;
                                        temp_sol[aa][0]=temp_sol[aa][ja];
                                        temp_sol[aa].pop_back();
                                        //tempPos[aa]=temp_sol[aa][0];
                                        UNPM++;
                                        break;
                                    }else{
                                        temp_sol[aa].pop_back();
                                        ja--;
                                    }
                                   }
                            }
                            if (UNPM > pq.top().first&&UNPM!=0)
                            {                                        
                                OG=true;
                                if (heapcount >= ksize)
                                {
                                   auto [value, index] = pq.top();      
                                    pq.pop();
                                    pq.push({UNPM,index});
                                    addPos=index;
                                    
                                }else{
                                     pq.push({UNPM,heapcount});
                                     addPos=heapcount;            
                                }
                            for (int aa=0;aa<qsiz;aa++){
                                SolPos[addPos][aa]=temp_sol[aa][0];
                                //tempPos[aa]=temp_sol[aa][0];
                                    }
                                heapcount++;
                            }
                            }
                            
                            for (auto& inner_vec : temp_sol) {
                                inner_vec.clear();  // Clears each inner vector
                            }

                            // reverse_embedding.erase(embedding[u]);
                            vec_failing_set[cur_depth].set();
                            vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
                        }
                        else
                        { // update equivalent class of VN[utemp]
                            int ia = 0;
                            int io = 0;
                            VNTemp.clear();
                            while (ia < idx_count[cur_depth] && io < VN[utemp].size())
                            {

                                ui valid_idxN = valid_candidate_idx[cur_depth][ia];
                                if (valid_idxN == 10000000)
                                {
                                    ia += 1;
                                    continue;
                                }
                                if (VN[utemp][io] == valid_idxN)
                                {
                                    // VNTemp.push_back(VN[u][io]);
                                    VNTemp.push_back(VN[utemp][io]);
                                    ia++;
                                    io++;
                                }

                                else if (VN[utemp][io] > valid_idxN)
                                    ia++;
                                else
                                    io++;
                                // ia++;
                            }
                            VN[utemp] = VNTemp;
                            VNTemp.clear();
                        }
                    }
                }
                // }
                //}
                VN[u].clear();

                idx[cur_depth] += 1;

                vec_failing_set[cur_depth] = ancestors[u];
                vec_failing_set[cur_depth] |= ancestors[reverse_embedding[v]];
                vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];

                continue;
            }
            embd_vertices[v] = true;
            counter++;
            embedding[u] = v;             // set the embeding of depth qid to v
            idx_embedding[u] = valid_idx; // reverse the id of the embdeing to the position
            visited_vertices[v] = true;
            idx[cur_depth] += 1; // next element

            while (cur_depth < idx_count[cur_depth] && valid_candidate_idx[cur_depth][idx[cur_depth]] == 10000000)
            {
                idx[cur_depth] += 1;
            }

            reverse_embedding[v] = u;


            if (cur_depth == max_depth - 1)
            {
                embedding_cnt += 1;
                visited_vertices[v] = false;
                embd_vertices[v] = false;
                counter--;
                int UNPMtemp = 0;
                // UNPMtemp=EmbSum.insert(embedding[order[cur_depth]]).second;

                int ao = cur_depth; //-1;
                ui vqo = order[ao];
                int UNPM2 = 0;
                UNPM = 0; //
                while (ao >= 0)
                            {  
                                Match_BA[ao] = true;
                                // EmbSum.insert(embedding[vqo]);                              
                                temp_sol[ao].push_back(embedding[vqo]);
                                //&&ao!=max_depth-1
                                if (VN[vqo].size() > 0)
                                {
                                    for (int dd = 0; dd < VN[vqo].size(); dd++)
                                    {
                                        if (VN[vqo][dd] != 10000000)
                                        {       
                                            if (nodeId[candidates[vqo][VN[vqo][dd]]]==0)                                            
                                                temp_sol[ao].push_back(candidates[vqo][VN[vqo][dd]]);                                        
                                        }
                                    }
                                }

                                ao--;
                                vqo = order[ao];
                                //UNPM += UNPMtemp;
                            }
         
                            
                            for (int aa=0;aa<qsiz;aa++){
                                if(nodeId[temp_sol[aa][0]]==0){                                    
                                    nodeId[temp_sol[aa][0]]=1;
                                    UNPM ++;
                                    tempPos[aa]=temp_sol[aa][0];
                                }
                                        
                                else if(temp_sol[aa].size()==1) 
                                    continue;  
                                else{
                                   int ja=temp_sol[aa].size()-1;
                                   while(ja>0){
                                    
                                    if (nodeId[temp_sol[aa][ja]]==0){
                                        nodeId[temp_sol[aa][ja]]=1;
                                        temp_sol[aa][0]=temp_sol[aa][ja];
                                        temp_sol[aa].pop_back();
                                        UNPM ++;
                                        break;
                                    }else{
                                        temp_sol[aa].pop_back();
                                        ja--;
                                    }

                                   }
                                }
                            }

                            if (UNPM > pq.top().first&&UNPM!=0)
                            {                                        
                                
                                if (heapcount >= ksize)
                                {
                                   auto [value, index] = pq.top();      
                                    pq.pop();
                                    pq.push({UNPM,index});
                                    addPos=index;
                                    
                                }else{
                                     pq.push({UNPM,heapcount});
                                     addPos=heapcount;            
                                }
                            for (int aa=0;aa<qsiz;aa++){
                                SolPos[addPos][aa]=temp_sol[aa][0];
                                //tempPos[aa]=temp_sol[aa][0];
                                    }
                                heapcount++;
                            }
                            //embedding_cnt += UNPM2 - 1;
                            int aa=0;
                            bool go=true;
                            int ii=1;
                            //reverse greedy
                            bool OG=true;
                            if(idx[cur_depth] == idx_count[cur_depth])
                            while (OG)
                            {   OG=false;
                                UNPM=0;
                                for (int aa=0;aa<qsiz;aa++){
                                if (temp_sol[aa].size()==1){
                                continue;
                                    }
                                int ja=temp_sol[aa].size()-1;
                                while(ja>0){
                                    
                                    if (nodeId[temp_sol[aa][ja]]==0){
                                        nodeId[temp_sol[aa][ja]]=1;
                                        temp_sol[aa][0]=temp_sol[aa][ja];
                                        temp_sol[aa].pop_back();
                                        UNPM++;
                                        break;
                                    }else{
                                        temp_sol[aa].pop_back();
                                        ja--;
                                    }
                                   }
                            }
                            if (UNPM > pq.top().first&&UNPM!=0)
                            {                                        
                                OG=true;
                                if (heapcount >= ksize)
                                {
                                   auto [value, index] = pq.top();      
                                    pq.pop();
                                    pq.push({UNPM,index});
                                    addPos=index;
                                    
                                }else{
                                     pq.push({UNPM,heapcount});
                                     addPos=heapcount;            
                                }
                            for (int aa=0;aa<qsiz;aa++){
                                SolPos[addPos][aa]=temp_sol[aa][0];
                                //tempPos[aa]=temp_sol[aa][0];
                                    }
                                heapcount++;
                            }

                            }
                            for (auto& inner_vec : temp_sol) {
                                inner_vec.clear();  // Clears each inner vector
                            }

                reverse_embedding.erase(embedding[u]);
                vec_failing_set[cur_depth].set();
                vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];


                if (embedding_cnt >= output_limit_num)
                {
                    goto EXIT;
                }
            }
            else
            {
                call_count += 1;
                cur_depth += 1;     // next depth
                idx[cur_depth] = 0; // set the element to 0
                // generateValidCandidateIndex(cur_depth, idx_embedding, idx_count, valid_candidate_idx, edge_matrix, bn,
                //                             bn_count, order, temp_buffer, candidates, candidates_count, query_graph);
                generateValidCandidateIndex(cur_depth, idx_embedding, idx_count, valid_candidate_idx, edge_matrix, bn,
                                            bn_count, order, temp_buffer);
                pruneCandidatesIndexBySymmetryBreakingC(cur_depth, embedding, order, candidates_count, candidates,
                                                       idx_count, valid_candidate_idx, ordered_constraints);

                // get the candidates for next depth

                if (idx_count[cur_depth] == 0)
                {
                    vec_failing_set[cur_depth - 1] = ancestors[order[cur_depth]];
                }
                else
                {
                    vec_failing_set[cur_depth - 1].reset();
                }

            }
        }

        cur_depth -= 1;
        LastMatc=-1;
        if (cur_depth == max_depth - 1)
        {
            continue;
        }

        if (cur_depth < 0)
            break;
        else
        {

            ui vq = order[cur_depth];
            if (VN[vq].size() > 1)
            { 
                int add = 0;
                for (int dd = 0; dd < VN[vq].size(); dd++)
                {

                    for (int kk = (idx[cur_depth]); kk < idx_count[cur_depth]; kk++)
                    {
                        ui valid_idxn = valid_candidate_idx[cur_depth][kk];
                        if (valid_idxn == 10000000)
                        {

                            continue;
                        }

                        if (valid_idxn == VN[vq][dd])
                        {
                            valid_candidate_idx[cur_depth][kk] = 10000000;
                            if(Match_BA[cur_depth]==true &&nodeId[candidates[order[cur_depth]][valid_idxn]]==0){
                                nodeId[candidates[order[cur_depth]][valid_idxn]]=1;                                for (int mm=0;mm<max_depth;mm++){
                                    if (mm<cur_depth)
                                    temp_sol[mm][0]=embedding[order[mm]];
                                    else
                                    temp_sol[mm][0]=candidates[order[cur_depth]][valid_idxn];
                                }
                             if (1 > pq.top().first)
                            {                                        
                                
                                if (heapcount >= ksize)
                                {
                                   auto [value, index] = pq.top();      
                                    pq.pop();
                                    pq.push({1,index});
                                    addPos=index;                                    
                                }else{
                                     pq.push({1,heapcount});
                                     addPos=heapcount;            
                                }
                            for (int aa=0;aa<qsiz;aa++){
                                SolPos[addPos][aa]=temp_sol[aa][0];
                                //tempPos[aa]=temp_sol[aa][0];
                                    }
                                heapcount++;
                            }       
                            }
                            break;
                        }
                    }
                }
            }
            VN[vq].clear();

            VertexID u = order[cur_depth];

            reverse_embedding.erase(embedding[u]);
            if (cur_depth != 0)
            {
                if (!vec_failing_set[cur_depth].test(u))
                {
                    vec_failing_set[cur_depth - 1] = vec_failing_set[cur_depth];
                    idx[cur_depth] = idx_count[cur_depth];
                }
                else
                {
                    vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
                }
            }

            visited_vertices[embedding[u]] = false;
            embd_vertices[embedding[u]] = false;
            counter--;
            // if (idx[cur_depth] == idx_count[cur_depth])
            // continue;

            // if (Match_BA[vq] == false && idx[cur_depth] < idx_count[cur_depth])
            //{ positive and negative removal now. dont need to count positive symmetric
        }
    }

    EXIT:
    releaseBuffer(max_depth, idx, idx_count, embedding, idx_embedding, temp_buffer, valid_candidate_idx,
                  visited_vertices,
                  bn, bn_count);


    std::priority_queue<std::pair<int, int>> maxHeap;
    int topkcounter = 0;
    int greedysum = 0;
    unordered_set<ui> EmbGreedy;
    while (!pq.empty())
    {   if(pq.top().first!=0)
        maxHeap.push(pq.top());
        pq.pop();
    }


    while (!maxHeap.empty())
    {   auto [value, index]= maxHeap.top();
        for (int aa=0;aa<qsiz;aa++){
            EmbGreedy.insert(SolPos[index][aa]);
        }
        maxHeap.pop();
        topkcounter++;
        switch (topkcounter)
        {

        case VALUE_10:
            s.topk[0] = EmbGreedy.size();
            break;
        case VALUE_50:
            s.topk[1] = EmbGreedy.size();
            break;
        case VALUE_100:
            s.topk[2] = EmbGreedy.size();
            break;
        case VALUE_500:
            s.topk[3] = EmbGreedy.size();
            break;
        case VALUE_1000:
            s.topk[4] = EmbGreedy.size();
            break;
        case VALUE_2500:
            s.topk[5] = EmbGreedy.size();
            break;    
        case VALUE_5000:
            s.topk[6] = EmbGreedy.size();
            break;
        case VALUE_10000:
            s.topk[7] = EmbGreedy.size();
            break;    
        default:
            break;
        }
    }
    for (int i = 0; i < 8; i++)
    {
        if (s.topk[i] == 0)
            s.topk[i] = EmbGreedy.size();
    }
    int true_cand_sum = 0;
    int countSol=0;
    for (int i=0;i<data_graph->getVerticesCount();i++){
        countSol+=nodeId[i];
    }
    s.Can_embed = countSol;
    //EmbSum.clear();
    s.embedding_cnt = embedding_cnt;
    s.candidate_true_count_sum = true_cand_sum;
    return s;
}
void EvaluateQuery::RestoreValidCans(const Graph *query_graph, const Graph *data_graph, bool *visited_u,
                                     VertexID last_u, VertexID last_v,
                                     std::vector<std::unordered_map<VertexID, ui>> &valid_cans)
{
    ui last_unbrs_count;
    const ui *last_unbrs = query_graph->getVertexNeighbors(last_u, last_unbrs_count);
    ui last_vnbrs_count;
    const ui *last_vnbrs = data_graph->getVertexNeighbors(last_v, last_vnbrs_count);

    for (ui i = 0; i < last_unbrs_count; i++)
    {
        ui last_unbr = last_unbrs[i];

        if (visited_u[last_unbr] == true)
            continue;
        for (ui j = 0; j < last_vnbrs_count; j++)
        {
            auto vertex = valid_cans[last_unbr].find(last_vnbrs[j]);
            if (vertex != valid_cans[last_unbr].end())
            {
                if (vertex->second == 1)
                {
                    valid_cans[last_unbr].erase(vertex);
                }
                else
                {
                    vertex->second--;
                }
            }
        }
    }
}

ui EvaluateQuery::generateNextU(const Graph *data_graph, const Graph *query_graph, ui **candidates, ui *candidates_count,
                                ui**valid_cans, ui*valid_cans_count, ui* extendable,  ui** nec,
                                ui depth, ui* embedding, Edges ***edge_matrix, bool *visited_vertices,
                                bool *visited_u, ui *order, TreeNode* tree) {
    
    ui query_vertices_num = query_graph->getVerticesCount();
    ui cur_vertex = -1;
    TreeNode &node = tree[order[depth - 1]];
    for (ui i = 0; i < node.fn_count_; ++i) {
        VertexID u = node.fn_[i];
        extendable[u] -= 1;
        if (extendable[u] == 0) { 
            ComputeValidCans(data_graph, query_graph, candidates, candidates_count, valid_cans,
                             valid_cans_count, embedding, u, visited_u);
        }
    }

    bool f_only_1 = true;
    
    for (ui i = 0; i < query_vertices_num; i++) {
        if (visited_u[i] == true || extendable[i] != 0) continue;
        ui u_count = 0;
        for (ui j = 0; j < valid_cans_count[i]; j++) {
            if (visited_vertices[valid_cans[i][j]] == false) u_count++;
        }
                
        if (nec[i] != NULL && *(nec[i]) == u_count) {
                (*(nec[i]))--;
                return i;

        }
        /*
        if (nec[i] != NULL && *(nec[i]) >= u_count) {
            if (*(nec[i]) > u_count) {
                
                return (ui)-1;
            } else {
                (*(nec[i]))--;
                return i;
            }
        }*/
        if (nec[i] != NULL && *(nec[i]) == u_count){
            (*(nec[i]))--;
                return i;
        }
        if (nec[i] != NULL) {
            if(cur_vertex == (ui)-1) cur_vertex = i;
        } else if (f_only_1 == true) {
            f_only_1 = false;
        }
    }

    if (f_only_1 == false) { 
        cur_vertex = (ui)-1;
        ui min_Um_u = (ui)-1;
        for (ui i = 0; i < query_vertices_num; i++) {
            if (visited_u[i] == true || nec[i] != 0 || extendable[i] != 0) continue;
            if ( min_Um_u > valid_cans_count[i]) {
            //if (valid_cans_count[i] != 0 && min_Um_u > valid_cans_count[i]) {
                cur_vertex = i;
                min_Um_u = valid_cans_count[i];
            }
            
        }
    }

  
    if (cur_vertex != (ui)-1 && nec[cur_vertex] != NULL) {
        (*nec[cur_vertex])--;
    }

    return cur_vertex;
}

void EvaluateQuery::ComputeValidCans(const Graph *data_graph, const Graph *query_graph, ui **candidates, ui *candidates_count,
                      ui**valid_cans, ui*valid_cans_count, ui* embedding, VertexID u, bool* visited_u) {
    ui unbrs_count;
    const ui* unbrs = query_graph->getVertexNeighbors(u, unbrs_count);
    valid_cans_count[u] = 0;
    
    
    
    
    
    for (ui i = 0; i < candidates_count[u]; i++) {
        VertexID v = candidates[u][i];
        bool flag = true;
        for (ui j = 0; j < unbrs_count; j++) {
            if (visited_u[unbrs[j]] == true && !data_graph->checkEdgeExistence(v, embedding[unbrs[j]])) {
                flag = false;
                break;
            }
        }
        if (flag == true) {
            valid_cans[u][valid_cans_count[u]++] = v;
        }
    }
    
    
    
    
    
}
void EvaluateQuery::computeNEC_Adaptive(const Graph *query_graph, Edges ***edge_matrix, ui *candidates_count,
                                        ui **candidates, std::vector<std::vector<ui>> &vec_index,
                                        std::vector<std::vector<ui>> &vec_set, ui i, ui j, ui vec_count)
{ // just give i and j inside and remove the dual loop.
    std::vector<ui> tmp_vec;

    tmp_vec.reserve(candidates_count[i]);
    ui unbrs_count;
    const ui *unbrs = query_graph->getVertexNeighbors(i, unbrs_count);
    vec_index[i][j] = vec_count;
    tmp_vec.push_back(candidates[i][j]);
    for (ui k = j + 1; k < candidates_count[i]; k++)
    {
        if (vec_index[i][k] != (ui)-1)
            continue;

        bool equ = true;
        for (ui u1 = 0; u1 < unbrs_count; u1++)
        {
            ui unbr = unbrs[u1];
            const Edges *edges = edge_matrix[i][unbr];
            if (edges->offset_[j + 1] - edges->offset_[j] != edges->offset_[k + 1] - edges->offset_[k])
            {
                equ = false;
                break;
            }

            for (ui u2 = 0; u2 < edges->offset_[j + 1] - edges->offset_[j]; u2++)
            {
                if (edges->edge_[u2 + edges->offset_[j]] != edges->edge_[u2 + edges->offset_[k]])
                {
                    equ = false;
                    break;
                }
            }
            if (equ == false)
                break;
        }
        if (equ == true)
        {
            tmp_vec.push_back(candidates[i][k]);
            vec_index[i][k] = vec_index[i][j];
        }
    }

    vec_set.push_back(tmp_vec);
    tmp_vec.clear();
}

void EvaluateQuery::computeNEC(const Graph *query_graph, Edges ***edge_matrix, ui *candidates_count,
                               ui**candidates, std::vector<std::vector<ui>>& vec_index,
                               std::vector<std::vector<ui>>& vec_set) {
    std::vector<ui> tmp_vec;
    ui vec_count = 0;
    for (ui i = 0; i < vec_index.size(); i++) {
        tmp_vec.reserve(candidates_count[i]);
        ui unbrs_count;
        const ui *unbrs = query_graph->getVertexNeighbors(i, unbrs_count);
        for (ui j = 0; j < candidates_count[i]; j++) {
            if (vec_index[i][j] != (ui)-1)
                continue;
            vec_index[i][j] = vec_count++;
            tmp_vec.push_back(candidates[i][j]);
            for (ui k = j + 1; k < candidates_count[i]; k++) {
                if (vec_index[i][k] != (ui)-1)
                    continue;
                
                bool equ = true;
                for (ui u1 = 0; u1 < unbrs_count; u1++) {
                    ui unbr = unbrs[u1];
                    const Edges* edges = edge_matrix[i][unbr];
                    if (edges->offset_[j+1]-edges->offset_[j] != edges->offset_[k+1] - edges->offset_[k]) {
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        equ = false;
                        break;
                    }
                    
                    
                    for (ui u2 = 0; u2 < edges->offset_[j+1] - edges->offset_[j]; u2++) {
                        if (edges->edge_[u2+edges->offset_[j]] != edges->edge_[u2+edges->offset_[k]]) {
                            
                            
                            
                            equ = false;
                            break;
                        }
                    }
                    if (equ == false) break;
                }
                if (equ == true) {
                    tmp_vec.push_back(candidates[i][k]);
                    vec_index[i][k] = vec_index[i][j];
                }
            }
            
            vec_set.push_back(tmp_vec);
            tmp_vec.clear();
        }
    }          
}

void EvaluateQuery::computeNEC(const Graph *query_graph, ui** nec) {
    ui query_vertices_num = query_graph->getVerticesCount();
    ui* nec_tmp = new ui[query_vertices_num];
    ui flag = true; 
    for (ui i = 0; i < query_vertices_num; i++) {
        
        if (query_graph->getVertexDegree(i) != 1 || nec[i] != NULL) {
            continue;
        }
        ui unbrs_num;
        const ui* unbrs = query_graph->getVertexNeighbors(i, unbrs_num);
        ui u_label = query_graph->getVertexLabel(i);
        ui* nec_count = new ui(0);
        nec_tmp[(*nec_count)++] = i;
        for (ui j = i + 1; j < query_vertices_num; j++) {
            if (query_graph->getVertexDegree(i) != 1 || nec[i] != NULL) {
                continue;
            }
            
            if (u_label != query_graph->getVertexLabel(j)) {
                continue;
            }
            ui u2_nbrs_num;
            const ui* u2_nbrs = query_graph->getVertexNeighbors(j, u2_nbrs_num);
            if (u2_nbrs_num != unbrs_num) {
                continue;
            }
            ui k1 = 0;
            for (k1 = 0; k1 < unbrs_num; k1++) {
                flag = false;
                for (ui k2 = 0; k2 < u2_nbrs_num; k2++) {
                    if (unbrs[k1] == u2_nbrs[k2]) {
                        flag = true;
                        break;
                    }
                }
                if (flag == false) break; 
            }
            if (k1 != unbrs_num) continue;
            nec_tmp[(*nec_count)++] = j;
        }
        
        
        for (ui j = 0; j < *nec_count; j++) {
            nec[nec_tmp[j]] = nec_count;
            
        }
        
    }
    delete []nec_tmp;
}

enumResult
EvaluateQuery::exploreRMStyle(const Graph *query_graph, const Graph *data_graph, catalog *&storage, Edges ***edge_matrix,
                              ui **candidates, ui *candidates_count, ui *order, size_t output_limit_num, size_t &call_count, int TimeL,
                              const std::unordered_map<VertexID, std::pair<std::set<VertexID>, std::set<VertexID>>> &order_constraints)
{
    enumResult s;
    auto start = std::chrono::high_resolution_clock::now();
    TimeL = TimeL * 1000;
    if (storage == NULL)
    {
        storage = new catalog(query_graph, data_graph);
        convertCans2Catalog(query_graph, candidates, edge_matrix, storage);
        storage->query_graph_->get2CoreSize();
        storage->data_graph_->getVerticesCount();
    }

    convert_to_encoded_relation(storage, order);
#ifdef SPARSE_BITMAP
    convert_encoded_relation_to_sparse_bitmap(storage, order);
#endif
    //if (getValue3() > MemSize)
   //     MemSize = getValue3();
    std::vector<ui> vorder;
    vorder.reserve(query_graph->getVerticesCount());
    vorder.insert(vorder.end(), order, order + query_graph->getVerticesCount());
    auto tree = execution_tree_generator::generate_single_node_execution_tree(vorder);
    auto end = std::chrono::high_resolution_clock::now();
    auto execution_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    // s.embedding_cnt = tree->execute(*storage, output_limit_num, call_count,TimeL-execution_time);
    s = tree->execute(*storage, output_limit_num, call_count, TimeL - execution_time);
       // if (getValue3() > MemSize)
       // MemSize = getValue3();
    //cout << "MemSize" << MemSize / 1000 << endl;
    return s;
}

void EvaluateQuery::convertCans2Catalog(const Graph *query_graph, ui **candidates, Edges ***edge_matrix, catalog *storage)
{

    for (ui u = 0; u < storage->num_sets_; u++)
    {
        ui unbrs_cnt;
        const ui *unbrs = query_graph->getVertexNeighbors(u, unbrs_cnt);
        for (ui i = 0; i < unbrs_cnt; i++)
        {
            ui v = unbrs[i];

            if (u > v)
                continue;
            std::vector<edge> tmp_edges;
            auto &edges = *edge_matrix[u][v];
            auto &relation = storage->edge_relations_[u][v];
            for (ui j = 0; j < edges.vertex_count_; j++)
            {
                ui src = candidates[u][j];
                for (ui k = edges.offset_[j]; k < edges.offset_[j + 1]; k++)
                {
                    ui dst = candidates[v][edges.edge_[k]];
                    tmp_edges.push_back(std::move(edge(src, dst)));
                }
            }
            relation.size_ = tmp_edges.size();
            relation.edges_ = new edge[relation.size_];
            memcpy(relation.edges_, tmp_edges.data(), sizeof(edge) * relation.size_);
            // std::cout << u << '|' << v << "(edge):" << tmp_edges.size() << std::endl;
        }
    }
}

void EvaluateQuery::convert_to_encoded_relation(catalog *storage, ui *order)
{
    auto &query_graph = storage->query_graph_;
    uint32_t core_vertices_cnt = query_graph->get2CoreSize();
    auto max_vertex_id = storage->data_graph_->getVerticesCount();

    auto projection_operator = new projection(max_vertex_id);
    for (uint32_t i = 0; i < core_vertices_cnt || i == 0; ++i)
    {
        uint32_t u = order[i];
        uint32_t nbr_cnt;
        const uint32_t *nbrs = query_graph->getVertexNeighbors(u, nbr_cnt);
        for (uint32_t j = 0; j < nbr_cnt; ++j)
        {
            uint32_t v = nbrs[j];
            uint32_t src = u;
            uint32_t dst = v;
            uint32_t kp = 0;
            if (src > dst)
            {
                std::swap(src, dst);
                kp = 1;
            }

            projection_operator->execute(&storage->edge_relations_[src][dst], kp, storage->candidate_sets_[u], storage->num_candidates_[u]);

            break;
        }
    }
    /*
    std::cout << "num_cans after encoded:";
    for (ui i = 0; i < query_graph->getVerticesCount(); i++)
    {
        std::cout << storage->num_candidates_[i] << ",";
    }
    std::cout << std::endl;
*/
    delete projection_operator;

    uint32_t n = query_graph->getVerticesCount();
    for (uint32_t i = 1; i < n; ++i)
    {
        uint32_t u = order[i];
        for (uint32_t j = 0; j < i; ++j)
        {
            uint32_t bn = order[j];
            if (query_graph->checkEdgeExistence(bn, u))
            {
                if (i < core_vertices_cnt)
                {
                    convert_to_encoded_relation(storage, bn, u);
                }
                else
                {
                    convert_to_hash_relation(storage, bn, u);
                }
            }
        }
    }
}

void EvaluateQuery::convert_to_encoded_relation(catalog *storage, uint32_t u, uint32_t v)
{
    uint32_t src = std::min(u, v);
    uint32_t dst = std::max(u, v);
    edge_relation &target_edge_relation = storage->edge_relations_[src][dst];
    edge *edges = target_edge_relation.edges_;
    uint32_t edge_size = target_edge_relation.size_;
    auto max_vertex_id = storage->data_graph_->getVerticesCount();
    assert(edge_size > 0);

    auto buffer = new uint32_t[max_vertex_id];
    memset(buffer, 0, sizeof(uint32_t) * max_vertex_id);

    uint32_t v_candidates_cnt = storage->get_num_candidates(v);
    uint32_t *v_candidates = storage->get_candidates(v);

    for (uint32_t i = 0; i < v_candidates_cnt; ++i)
    {
        uint32_t v_candidate = v_candidates[i];
        buffer[v_candidate] = i + 1;
    }

    uint32_t u_p = 0;
    uint32_t v_p = 1;
    if (u > v)
    {

        std::sort(edges, edges + edge_size, [](edge &l, edge &r) -> bool
                  {
            if (l.vertices_[1] == r.vertices_[1])
                return l.vertices_[0] < r.vertices_[0];
            return l.vertices_[1] < r.vertices_[1]; });
        u_p = 1;
        v_p = 0;
    }

    encoded_trie_relation &target_encoded_trie_relation = storage->encoded_trie_relations_[u][v];
    uint32_t edge_cnt = edge_size;
    uint32_t u_candidates_cnt = storage->get_num_candidates(u);
    uint32_t *u_candidates = storage->get_candidates(u);
    target_encoded_trie_relation.size_ = u_candidates_cnt;
    target_encoded_trie_relation.offset_ = new uint32_t[u_candidates_cnt + 1];
    target_encoded_trie_relation.children_ = new uint32_t[edge_size];

    uint32_t offset = 0;
    uint32_t edge_index = 0;

    for (uint32_t i = 0; i < u_candidates_cnt; ++i)
    {
        uint32_t u_candidate = u_candidates[i];
        target_encoded_trie_relation.offset_[i] = offset;
        uint32_t local_degree = 0;
        while (edge_index < edge_cnt)
        {
            uint32_t u0 = edges[edge_index].vertices_[u_p];
            uint32_t v0 = edges[edge_index].vertices_[v_p];
            if (u0 == u_candidate)
            {
                if (buffer[v0] > 0)
                {
                    target_encoded_trie_relation.children_[offset + local_degree] = buffer[v0] - 1;
                    local_degree += 1;
                }
            }
            else if (u0 > u_candidate)
            {
                break;
            }

            edge_index += 1;
        }

        offset += local_degree;

        if (local_degree > target_encoded_trie_relation.max_degree_)
        {
            target_encoded_trie_relation.max_degree_ = local_degree;
        }
    }

    target_encoded_trie_relation.offset_[u_candidates_cnt] = offset;

    for (uint32_t i = 0; i < v_candidates_cnt; ++i)
    {
        uint32_t v_candidate = v_candidates[i];
        buffer[v_candidate] = 0;
    }
}

void EvaluateQuery::convert_to_hash_relation(catalog *storage, uint32_t u, uint32_t v)
{

    uint32_t src = std::min(u, v);
    uint32_t dst = std::max(u, v);
    edge_relation &target_edge_relation = storage->edge_relations_[src][dst];
    hash_relation &target_hash_relation1 = storage->hash_relations_[u][v];
    auto max_vertex_id = storage->data_graph_->getVerticesCount();

    edge *edges = target_edge_relation.edges_;
    uint32_t edge_size = target_edge_relation.size_;

    assert(edge_size > 0);

    uint32_t u_key = 0;
    uint32_t v_key = 1;

    if (src != u)
    {
        std::swap(u_key, v_key);

        std::sort(edges, edges + edge_size, [](edge &l, edge &r) -> bool
                  {
            if (l.vertices_[1] == r.vertices_[1])
                return l.vertices_[0] < r.vertices_[0];
            return l.vertices_[1] < r.vertices_[1]; });
    }

    target_hash_relation1.children_ = new uint32_t[edge_size];

    uint32_t offset = 0;
    uint32_t local_degree = 0;
    uint32_t prev_u0 = max_vertex_id + 1;

    for (uint32_t i = 0; i < edge_size; ++i)
    {
        uint32_t u0 = edges[i].vertices_[u_key];
        uint32_t u1 = edges[i].vertices_[v_key];
        if (u0 != prev_u0)
        {
            if (prev_u0 != max_vertex_id + 1)
                target_hash_relation1.trie_->emplace(prev_u0, std::make_pair(local_degree, offset));

            offset += local_degree;

            if (local_degree > target_hash_relation1.max_degree_)
            {
                target_hash_relation1.max_degree_ = local_degree;
            }

            local_degree = 0;
            prev_u0 = u0;
        }

        target_hash_relation1.children_[offset + local_degree] = u1;
        local_degree += 1;
    }

    target_hash_relation1.cardinality_ = edge_size;
    target_hash_relation1.trie_->emplace(prev_u0, std::make_pair(local_degree, offset));
    if (local_degree > target_hash_relation1.max_degree_)
    {
        target_hash_relation1.max_degree_ = local_degree;
    }
}

void EvaluateQuery::convert_encoded_relation_to_sparse_bitmap(catalog *storage, ui *order)
{
    uint32_t core_vertices_cnt = storage->query_graph_->get2CoreSize();

    for (uint32_t i = 1; i < core_vertices_cnt; ++i)
    {
        uint32_t u = order[i];

        for (uint32_t j = 0; j < i; ++j)
        {
            uint32_t bn = order[j];
            if (storage->query_graph_->checkEdgeExistence(u, bn))
            {
                storage->bsr_relations_[bn][u].load(storage->encoded_trie_relations_[bn][u].get_size(),
                                                    storage->encoded_trie_relations_[bn][u].offset_,
                                                    storage->encoded_trie_relations_[bn][u].offset_,
                                                    storage->encoded_trie_relations_[bn][u].children_,
                                                    storage->max_num_candidates_per_vertex_, true);
            }
        }
    }
}

void EvaluateQuery::ComputeValidCans(const Graph *data_graph, const Graph *query_graph, ui **candidates, ui *candidates_count,
                                     ui **valid_cans, ui *valid_cans_count, ui *embedding, VertexID u, bool *visited_u, bool *visited_v)
{
    ui unbrs_count;
    const ui *unbrs = query_graph->getVertexNeighbors(u, unbrs_count);
    valid_cans_count[u] = 0;

    for (ui i = 0; i < candidates_count[u]; i++)
    {
        VertexID v = candidates[u][i];
        if (visited_v[v])
            continue;
        bool flag = true;
        for (ui j = 0; j < unbrs_count; j++)
        {
            if (visited_u[unbrs[j]] == true && !data_graph->checkEdgeExistence(v, embedding[unbrs[j]]))
            {
                flag = false;
                break;
            }
        }
        if (flag == true)
        {
            valid_cans[u][valid_cans_count[u]++] = v;
        }
    }
}

enumResult
EvaluateQuery::exploreKSSStyle(const Graph *query_graph, const Graph *data_graph, Edges ***edge_matrix,
                               ui **candidates, ui *candidates_count, ui *order, size_t output_limit_num,
                               size_t &call_count, int TimeL,
                               const std::unordered_map<VertexID, std::pair<std::set<VertexID>, std::set<VertexID>>> &order_constraints)
{
    auto start = std::chrono::high_resolution_clock::now();
    enumResult s;
    size_t embedding_cnt = 0;
    int cur_depth = 0;
    ui query_vertices_num = query_graph->getVerticesCount();
    ui data_vertices_num = data_graph->getVerticesCount();
    ui kernel_num = 0, shell_num = 0;
    ui *kernel = new ui[query_vertices_num];
    ui *shell = new ui[query_vertices_num];
    ui max_depth = query_vertices_num;
    unordered_set<ui> EmbSum;
    bool Match_BA[max_depth] = {false};
    TimeL = TimeL * 1000;
    double ens = 0;
    bool *kos = new bool[query_vertices_num];
    memset(kos, 0, sizeof(bool) * query_vertices_num);

    ui *degree = new ui[query_vertices_num];
    for (ui i = 0; i < query_vertices_num; i++)
    {
        degree[i] = query_graph->getVertexDegree(i);
    }
    for (ui i = 0; i < query_vertices_num; i++)
    {
        VertexID u = order[i];
        if (degree[u] == 0)
        {
            shell[shell_num++] = u;
            continue;
        }
        kernel[kernel_num++] = u;
        kos[u] = true;
        ui nbr_num = 0;
        const ui *nbrs = query_graph->getVertexNeighbors(u, nbr_num);
        for (ui j = 0; j < nbr_num; j++)
        {
            degree[nbrs[j]]--;
        }
    }

    ui *shell2kernel = new ui[query_vertices_num];
    memset(shell2kernel, 0, sizeof(ui) * query_vertices_num);
    for (ui i = 0; i < shell_num; i++)
    {
        VertexID v = shell[i];
        ui nbr_num = 0;
        const ui *nbrs = query_graph->getVertexNeighbors(v, nbr_num);
        for (ui j = 0; j < nbr_num; j++)
        {
            if (kos[nbrs[j]] == true)
            {
                shell2kernel[v]++;
            }
        }
    }

    ui *idx = new ui[kernel_num];
    ui *idx_count = new ui[kernel_num];
    ui *embedding = new ui[query_vertices_num];
    ui **valid_cans = new ui *[query_vertices_num];
    for (ui i = 0; i < query_vertices_num; i++)
    {
        valid_cans[i] = new ui[candidates_count[i]];
    }
    ui *valid_cans_cnt = new ui[query_vertices_num];
    bool *visited_v = new bool[data_vertices_num];
    memset(visited_v, 0, sizeof(bool) * data_vertices_num);
    bool *visited_u = new bool[query_vertices_num];
    memset(visited_u, 0, sizeof(bool) * query_vertices_num);

    VertexID start_vertex = kernel[0];
    visited_u[start_vertex] = true;
    idx[cur_depth] = 0;
    idx_count[cur_depth] = candidates_count[start_vertex];

    for (ui i = 0; i < idx_count[cur_depth]; ++i)
    {
        valid_cans[start_vertex][i] = candidates[start_vertex][i];
    }

    std::vector<VertexID> update;

    // std::cout << "query_vertices_num:" << query_vertices_num << std::endl;
    // std::cout << "kernel_num:" << kernel_num << std::endl;
    // std::cout << "shell_num:"  << shell_num << std::endl;

    while (true)
    {
        while (idx[cur_depth] < idx_count[cur_depth])
        {
            Match_BA[cur_depth] = false;
            auto end = std::chrono::high_resolution_clock::now();
            ens = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
            if (ens > TimeL)
            { // 1000 1 sec
                goto EXIT;
            }
            VertexID u = kernel[cur_depth];
            VertexID v = valid_cans[u][idx[cur_depth]];

            embedding[u] = v;
            visited_v[v] = true;
            idx[cur_depth] += 1;

            bool f_ok = true;
            updateShell2Kernel(query_graph, u, shell2kernel, kos, update);
            for (auto ushell : update)
            {
                ComputeValidCans(data_graph, query_graph, candidates, candidates_count, valid_cans,
                                 valid_cans_cnt, embedding, ushell, visited_u, visited_v);
                if (valid_cans_cnt[ushell] == 0)
                {
                    f_ok = false;
                    // std::cout << "break" << std::endl;
                    break;
                }
            }
            if (!f_ok)
            {
                visited_v[v] = false;
                restoreShell2Kernel(query_graph, u, shell2kernel, kos);
                continue;
            }

            if (cur_depth == kernel_num - 1)
            {
                for (int i = 0; i < shell_num; i++)
                {
                    VertexID u_shell = shell[i];

                    /*
                    std::cout << "candidate of u_shell " << u_shell << ": ";
                    for(int j = 0; j < valid_cans_cnt[u_shell]; j++) {
                        std::cout << valid_cans[u_shell][j] << ",";
                    }
                    std::cout << std::endl;*/
                }

                embedding_cnt += computeKSSEmbeddingNaive(shell_num, shell, valid_cans, valid_cans_cnt, visited_v);

                // std::cout << embedding_cnt << std::endl;
                visited_v[v] = false;
                if (embedding_cnt >= output_limit_num)
                {
                    goto EXIT;
                }

                restoreShell2Kernel(query_graph, u, shell2kernel, kos);
            }
            else
            {
                cur_depth += 1;

                VertexID next_u = kernel[cur_depth];
                idx[cur_depth] = 0;
                call_count += 1;
                ComputeValidCans(data_graph, query_graph, candidates, candidates_count, valid_cans,
                                 valid_cans_cnt, embedding, next_u, visited_u, visited_v);

                visited_u[next_u] = true;
                idx_count[cur_depth] = valid_cans_cnt[next_u];
            }
        }

        VertexID last_u = kernel[cur_depth];

        cur_depth -= 1;
        if (cur_depth < 0)
            break;
        visited_v[embedding[kernel[cur_depth]]] = false;
        visited_u[last_u] = false;

        restoreShell2Kernel(query_graph, kernel[cur_depth], shell2kernel, kos);
    }

EXIT:
    delete[] kernel;
    delete[] shell;
    delete[] kos;
    delete[] degree;
    delete[] shell2kernel;
    delete[] idx;
    delete[] idx_count;
    delete[] embedding;
    for (ui i = 0; i < query_vertices_num; i++)
    {
        delete[] valid_cans[i];
    }
    delete[] valid_cans;
    delete[] valid_cans_cnt;
    delete[] visited_u;
    delete[] visited_v;
    s.embedding_cnt = embedding_cnt;
    return s;
}

void EvaluateQuery::updateShell2Kernel(const Graph *query_graph, VertexID u, ui *shell2kernel, bool *kos, std::vector<VertexID> &update)
{
    ui nbr_num = 0;
    const ui *nbrs = query_graph->getVertexNeighbors(u, nbr_num);

    update.clear();
    for (ui i = 0; i < nbr_num; i++)
    {
        VertexID nbr = nbrs[i];
        if (kos[nbr] == false)
        {
            shell2kernel[nbr]--;

            if (shell2kernel[nbr] == 0)
                update.emplace_back(nbr);
        }
    }
}

void EvaluateQuery::restoreShell2Kernel(const Graph *query_graph, VertexID u, ui *shell2kernel, bool *kos)
{
    ui nbr_num = 0;

    const ui *nbrs = query_graph->getVertexNeighbors(u, nbr_num);
    for (ui i = 0; i < nbr_num; i++)
    {
        VertexID nbr = nbrs[i];
        if (kos[nbr] == false)
        {
            shell2kernel[nbr]++;
        }
    }
}

size_t
EvaluateQuery::computeKSSEmbeddingNaive(ui shell_num, ui *shell, ui **valid_cans, ui *valid_cans_count, bool *visited_v)
{

    return computeKSSEmbeddingNaiveImpl(0, shell_num, shell, valid_cans, valid_cans_count, visited_v);
}

size_t
EvaluateQuery::computeKSSEmbeddingNaiveImpl(ui depth, ui shell_num, ui *shell, ui **valid_cans, ui *valid_cans_count, bool *visited_v)
{
    VertexID u_shell = shell[depth];
    size_t embedding_cnt = 0;

    if (depth == shell_num - 1)
    {

        for (int i = 0; i < valid_cans_count[u_shell]; i++)
        {
            VertexID v_id = valid_cans[u_shell][i];
            if (!visited_v[v_id])
                embedding_cnt += 1;
        }
        return embedding_cnt;
    }

    else
    {
        for (int i = 0; i < valid_cans_count[u_shell]; i++)
        {
            VertexID v_id = valid_cans[u_shell][i];
            if (!visited_v[v_id])
            {
                visited_v[v_id] = true;
                embedding_cnt += computeKSSEmbeddingNaiveImpl(depth + 1, shell_num, shell, valid_cans, valid_cans_count, visited_v);
                visited_v[v_id] = false;
            }
        }
        return embedding_cnt;
    }
}

size_t
EvaluateQuery::computeKSSEmbeddingOpt1(ui shell_num, ui *shell, ui **valid_cans, ui *valid_cans_count, bool *visited_v)
{

    std::unordered_map<VertexID, ui> counter;

    for (int i = 0; i < shell_num; i++)
    {
        VertexID u_shell = shell[i];
        for (int j = 0; j < valid_cans_count[u_shell]; j++)
        {
            VertexID v_id = valid_cans[u_shell][j];
            if (counter.count(v_id) == 0)
            {
                counter[v_id] = 0;
            }
            else
            {
                counter[v_id] = i;
            }
        }
    }

    return computeKSSEmbeddingOpt1Impl(0, shell_num, shell, valid_cans, valid_cans_count, counter, visited_v);
}

size_t
EvaluateQuery::computeKSSEmbeddingOpt1Impl(ui depth, ui shell_num, ui *shell, ui **valid_cans, ui *valid_cans_count,
                                           std::unordered_map<VertexID, ui> &counter,
                                           bool *visited_v)
{
    size_t embedding_cnt = 0;
    VertexID u_shell = shell[depth];

    if (depth == shell_num - 1)
    {

        for (int i = 0; i < valid_cans_count[u_shell]; i++)
        {
            VertexID v_id = valid_cans[u_shell][i];
            if (!visited_v[v_id])
                embedding_cnt += 1;
        }
        return embedding_cnt;
    }

    else
    {

        ui no_conflict_cnt = 0;
        for (int i = 0; i < valid_cans_count[u_shell]; i++)
        {
            VertexID v_id = valid_cans[u_shell][i];

            if (!visited_v[v_id])
            {

                if (counter[v_id] == 0 || counter[v_id] < depth)
                {
                    no_conflict_cnt += 1;
                }
                else
                {
                    if (counter[v_id] < depth + 1)
                    {
                        no_conflict_cnt += 1;
                    }
                    else
                    {
                        visited_v[v_id] = true;
                        embedding_cnt += computeKSSEmbeddingOpt1Impl(depth + 1, shell_num, shell, valid_cans, valid_cans_count, counter, visited_v);
                        visited_v[v_id] = false;
                    }
                }
            }
        }

        embedding_cnt += no_conflict_cnt * computeKSSEmbeddingOpt1Impl(depth + 1, shell_num, shell, valid_cans, valid_cans_count, counter, visited_v);

        return embedding_cnt;
    }
}

size_t
EvaluateQuery::computeKSSEmbeddingOpt2(ui shell_num, ui *shell, ui **valid_cans, ui *valid_cans_count)
{

    std::unordered_map<VertexID, ui> counter;
    std::set<VertexID> *used = new std::set<VertexID>[shell_num];
    std::map<std::set<VertexID>, size_t> *memorized_table = new std::map<std::set<VertexID>, size_t>[shell_num];

    for (int i = 0; i < shell_num; i++)
    {
        VertexID u_shell = shell[i];
        for (int j = 0; j < valid_cans_count[u_shell]; j++)
        {
            VertexID v_id = valid_cans[u_shell][j];
            if (counter.count(v_id) == 0)
            {
                counter[v_id] = 0;
            }
            else
            {
                counter[v_id] = i;
            }
        }
    }

    return computeKSSEmbeddingOpt2Impl(0, shell_num, shell, valid_cans, valid_cans_count, counter, used, memorized_table);
}

size_t
EvaluateQuery::computeKSSEmbeddingOpt2Impl(ui depth, ui shell_num, ui *shell, ui **valid_cans, ui *valid_cans_count,
                                           std::unordered_map<VertexID, ui> &counter,
                                           std::set<VertexID> *used,
                                           std::map<std::set<VertexID>, size_t> *memorized_table)
{

    std::map<std::set<VertexID>, size_t> cur_table = memorized_table[depth];

    std::set<VertexID> &cur_used = used[depth];
    std::set<VertexID> &nxt_used = used[depth + 1];

    if (cur_table.count(cur_used))
    {
        return cur_table[cur_used];
    }

    size_t embedding_cnt = 0;
    VertexID u_shell = shell[depth];

    if (depth == shell_num - 1)
    {
        size_t dup = 0;
        return valid_cans_count[u_shell] - cur_used.size();
    }

    else
    {

        nxt_used.clear();
        for (auto &v_used : cur_used)
        {
            if (counter[v_used] > depth)
            {
                nxt_used.insert(v_used);
            }
        }

        ui no_conflict_cnt = 0;
        for (int i = 0; i < valid_cans_count[u_shell]; i++)
        {
            VertexID v_id = valid_cans[u_shell][i];

            if (counter[v_id] == 0 || counter[v_id] < depth)
            {
                no_conflict_cnt += 1;
            }
            else
            {
                if (cur_used.find(v_id) == cur_used.end())
                {
                    if (counter[v_id] == depth)
                    {
                        no_conflict_cnt += 1;
                    }
                    else
                    {
                        nxt_used.insert(v_id);
                        embedding_cnt += computeKSSEmbeddingOpt2Impl(depth + 1, shell_num, shell, valid_cans, valid_cans_count, counter, used, memorized_table);
                        nxt_used.erase(v_id);
                    }
                }
            }
        }

        embedding_cnt += computeKSSEmbeddingOpt2Impl(depth + 1, shell_num, shell, valid_cans, valid_cans_count, counter, used, memorized_table);

        cur_table[cur_used] = embedding_cnt;
        return embedding_cnt;
    }
}

void EvaluateQuery::pruneCandidatesIndexBySymmetryBreaking(ui depth, ui *embedding, ui *order,
                                                           ui *candidates_count, ui **candidates,
                                                           ui *idx_count, ui **valid_candidates_idx,
                                                           const std::unordered_map<VertexID, std::pair<std::set<VertexID>, std::set<VertexID>>> &ordered_constraints)
{
    VertexID u = order[depth];

    auto it = ordered_constraints.find(u);
    if (it != ordered_constraints.end())
    {

        VertexID min_id = 0;
        VertexID max_id = -1;

        const auto &u_cons_prior = it->second.first;
        const auto &u_cons_inferior = it->second.second;

        for (auto &u_p : u_cons_prior)
        {
            if (embedding[u_p] > min_id)
                min_id = embedding[u_p];
        }

        for (auto &u_i : u_cons_inferior)
        {
            if (embedding[u_i] < max_id)
                max_id = embedding[u_i];
        }

        int valid_cnt = 0;

        for (int i = 0; i < idx_count[depth]; i++)
        {
            ui valid_idx = valid_candidates_idx[depth][i];
            VertexID v = candidates[u][valid_idx];
            if (min_id <= v && v < max_id)
            {
                valid_candidates_idx[depth][valid_cnt++] = valid_idx;
            }
        }
        idx_count[depth] = valid_cnt;
    }
}
void EvaluateQuery::pruneCandidatesIndexBySymmetryBreakingC(ui depth, ui *embedding, ui *order,
                                                           ui *candidates_count, ui **candidates,
                                                           ui *idx_count, ui **valid_candidates_idx,
                                                           const std::unordered_map<VertexID, std::pair<std::set<VertexID>, std::set<VertexID>>> &ordered_constraints)
{
    VertexID u = order[depth];

    auto it = ordered_constraints.find(u);
    if (it != ordered_constraints.end())
    {

        VertexID min_id = 0;
        VertexID max_id = -1;

        const auto &u_cons_prior = it->second.first;
        const auto &u_cons_inferior = it->second.second;

        for (auto &u_p : u_cons_prior)
        {
            if (embedding[u_p] > min_id)
                min_id = embedding[u_p];
        }

        for (auto &u_i : u_cons_inferior)
        {
            if (embedding[u_i] < max_id)
                max_id = embedding[u_i];
        }

        int valid_cnt = 0;

        for (int i = 0; i < idx_count[depth]; i++)
        {
            ui valid_idx = valid_candidates_idx[depth][i];
            VertexID v = candidates[u][valid_idx];
            if (min_id <= v && v <= max_id)
            {
                valid_candidates_idx[depth][valid_cnt++] = valid_idx;
            }
        }
        idx_count[depth] = valid_cnt;
    }
}

void EvaluateQuery::pruneCandidatesBySymmetryBreaking(ui depth, ui *embedding, ui *order,
                                                      ui *idx_count, ui **valid_candidates,
                                                      const std::unordered_map<VertexID, std::pair<std::set<VertexID>, std::set<VertexID>>> &ordered_constraints)
{
    VertexID u = order[depth];

    auto it = ordered_constraints.find(u);
    if (it != ordered_constraints.end())
    {

        VertexID min_id = 0;
        VertexID max_id = -1;

        const auto &u_cons_prior = it->second.first;
        const auto &u_cons_inferior = it->second.second;

        for (auto &u_p : u_cons_prior)
        {
            if (embedding[u_p] > min_id)
                min_id = embedding[u_p];
        }

        for (auto &u_i : u_cons_inferior)
        {
            if (embedding[u_i] < max_id)
                max_id = embedding[u_i];
        }

        int valid_cnt = 0;
        for (int i = 0; i < idx_count[depth]; i++)
        {
            VertexID v = valid_candidates[depth][i];
            if (min_id <= v && v < max_id)
            {
                valid_candidates[depth][valid_cnt++] = v;
            }
        }
        idx_count[depth] = valid_cnt;
    }
}

bool EvaluateQuery::checkSingeVertexBySymmetryBreaking(ui depth, ui *embedding, ui *order, VertexID v,
                                                       const std::unordered_map<VertexID, std::pair<std::set<VertexID>, std::set<VertexID>>> &ordered_constraints)
{
    VertexID u = order[depth];

    auto it = ordered_constraints.find(u);
    if (it != ordered_constraints.end())
    {

        const auto &u_cons_prior = it->second.first;
        const auto &u_cons_inferior = it->second.second;

        for (auto &u_p : u_cons_prior)
        {
            if (embedding[u_p] > v)
                return false;
        }

        for (auto &u_i : u_cons_inferior)
        {
            if (embedding[u_i] < v)
                return false;
        }
    }
    return true;
}

bool EvaluateQuery::checkSingeVertexBySymmetryBreaking(ui depth, ui *embedding, ui *order, bool *visited_u, VertexID v,
                                                       const std::unordered_map<VertexID, std::pair<std::set<VertexID>, std::set<VertexID>>> &full_constraints,
                                                       bool *violated_symmetry_u,
                                                       ui max_depth)
{
    VertexID u = order[depth];
    bool is_violated = false;

    auto it = full_constraints.find(u);
    if (it != full_constraints.end())
    {

        std::fill(violated_symmetry_u, violated_symmetry_u + max_depth, false);

        const auto &u_cons_prior = it->second.first;
        const auto &u_cons_inferior = it->second.second;

        for (auto &u_p : u_cons_prior)
        {
            if (visited_u[u_p] && embedding[u_p] > v)
            {
                is_violated = true;
                violated_symmetry_u[u_p] = true;
            }
        }

        for (auto &u_i : u_cons_inferior)
        {
            if (visited_u[u_i] && embedding[u_i] < v)
            {
                is_violated = true;
                violated_symmetry_u[u_i] = true;
            }
        }
    }

    return is_violated;
}

void EvaluateQuery::AncestorWithSymmetryBreaking(const Graph *query_graph,
                                                 std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> &ancestors,
                                                 const std::unordered_map<VertexID, std::pair<std::set<VertexID>, std::set<VertexID>>> &ordered_constraints)
{
    ui query_vertices_num = query_graph->getVerticesCount();

    for (auto it = ordered_constraints.begin(); it != ordered_constraints.end(); it++)
    {
        VertexID u = it->first;
        const auto &u_cons_prior = it->second.first;
        const auto &u_cons_inferior = it->second.second;

        for (auto &u_p : u_cons_prior)
        {
            ancestors[u] |= ancestors[u_p];
        }

        for (auto &u_i : u_cons_inferior)
        {
            ancestors[u] |= ancestors[u_i];
        }
    }
}

void EvaluateQuery::generateValidCandidateIndex(ui depth, ui *idx_embedding, ui *idx_count, ui **valid_candidate_index,
                                                Edges ***edge_matrix, ui **bn, ui *bn_cnt, ui *order,
                                                ui *&temp_buffer, ui **candidates, ui *candidates_count, const Graph *query_graph)
{
    VertexID u = order[depth];
    VertexID previous_bn = bn[depth][0];
    ui previous_index_id = idx_embedding[previous_bn];
    ui valid_candidates_count = 0;

#if ENABLE_QFLITER == 1
    BSRGraph &bsr_graph = *qfliter_bsr_graph_[previous_bn][u];
    BSRSet &bsr_set = bsr_graph.bsrs[previous_index_id];

    if (bsr_set.size_ != 0)
    {
        offline_bsr_trans_uint(bsr_set.base_, bsr_set.states_, bsr_set.size_,
                               (int *)valid_candidate_index[depth]);
        // represent bsr size
        valid_candidates_count = bsr_set.size_;
    }

    if (bn_cnt[depth] > 0)
    {
        if (temp_bsr_base1_ == nullptr)
        {
            temp_bsr_base1_ = new int[1024 * 1024];
        }
        if (temp_bsr_state1_ == nullptr)
        {
            temp_bsr_state1_ = new int[1024 * 1024];
        }
        if (temp_bsr_base2_ == nullptr)
        {
            temp_bsr_base2_ = new int[1024 * 1024];
        }
        if (temp_bsr_state2_ == nullptr)
        {
            temp_bsr_state2_ = new int[1024 * 1024];
        }
        int *res_base_ = temp_bsr_base1_;
        int *res_state_ = temp_bsr_state1_;
        int *input_base_ = temp_bsr_base2_;
        int *input_state_ = temp_bsr_state2_;

        memcpy(input_base_, bsr_set.base_, sizeof(int) * bsr_set.size_);
        memcpy(input_state_, bsr_set.states_, sizeof(int) * bsr_set.size_);

        for (ui i = 1; i < bn_cnt[depth]; ++i)
        {
            VertexID current_bn = bn[depth][i];
            ui current_index_id = idx_embedding[current_bn];
            BSRGraph &cur_bsr_graph = *qfliter_bsr_graph_[current_bn][u];
            BSRSet &cur_bsr_set = cur_bsr_graph.bsrs[current_index_id];

            if (valid_candidates_count == 0 || cur_bsr_set.size_ == 0)
            {
                valid_candidates_count = 0;
                break;
            }

            valid_candidates_count = intersect_qfilter_bsr_b4_v2(cur_bsr_set.base_, cur_bsr_set.states_,
                                                                 cur_bsr_set.size_,
                                                                 input_base_, input_state_, valid_candidates_count,
                                                                 res_base_, res_state_);

            swap(res_base_, input_base_);
            swap(res_state_, input_state_);
        }

        if (valid_candidates_count != 0)
        {
            valid_candidates_count = offline_bsr_trans_uint(input_base_, input_state_, valid_candidates_count,
                                                            (int *)valid_candidate_index[depth]);
        }
    }

    idx_count[depth] = valid_candidates_count;

    // Debugging.
#ifdef YCHE_DEBUG
    Edges &previous_edge = *edge_matrix[previous_bn][u];

    auto gt_valid_candidates_count =
        previous_edge.offset_[previous_index_id + 1] - previous_edge.offset_[previous_index_id];
    ui *previous_candidates = previous_edge.edge_ + previous_edge.offset_[previous_index_id];
    ui *gt_valid_candidate_index = new ui[1024 * 1024];
    memcpy(gt_valid_candidate_index, previous_candidates, gt_valid_candidates_count * sizeof(ui));

    ui temp_count;
    for (ui i = 1; i < bn_cnt[depth]; ++i)
    {
        VertexID current_bn = bn[depth][i];
        Edges &current_edge = *edge_matrix[current_bn][u];
        ui current_index_id = idx_embedding[current_bn];

        ui current_candidates_count =
            current_edge.offset_[current_index_id + 1] - current_edge.offset_[current_index_id];
        ui *current_candidates = current_edge.edge_ + current_edge.offset_[current_index_id];

        ComputeSetIntersection::ComputeCandidates(current_candidates, current_candidates_count,
                                                  gt_valid_candidate_index, gt_valid_candidates_count, temp_buffer,
                                                  temp_count);

        std::swap(temp_buffer, gt_valid_candidate_index);
        gt_valid_candidates_count = temp_count;
    }
    assert(valid_candidates_count == gt_valid_candidates_count);

    cout << "Ret, Level:" << bn_cnt[depth] << ", BSR:"
         << pretty_print_array(valid_candidate_index[depth], valid_candidates_count)
         << "; GT: " << pretty_print_array(gt_valid_candidate_index, gt_valid_candidates_count) << "\n";

    for (auto i = 0; i < valid_candidates_count; i++)
    {
        assert(gt_valid_candidate_index[i] == valid_candidate_index[depth][i]);
    }
    delete[] gt_valid_candidate_index;
#endif
#else
    if (previous_bn > query_graph->getVerticesCount())
    {

        valid_candidate_index[depth] = new ui[candidates_count[u]];
        idx_count[depth] = candidates_count[u];
        for (ui i = 0; i < idx_count[depth]; ++i)
        {
            valid_candidate_index[depth][i] = i;
        }
    }
    else if (query_graph->checkEdgeExistence(previous_bn, u) == false)
    {

        valid_candidate_index[depth] = new ui[candidates_count[u]];
        idx_count[depth] = candidates_count[u];
        for (ui i = 0; i < idx_count[depth]; ++i)
        {
            valid_candidate_index[depth][i] = i;
        }
    }
    else
    {
        Edges &previous_edge = *edge_matrix[previous_bn][u];
        valid_candidates_count = previous_edge.offset_[previous_index_id + 1] - previous_edge.offset_[previous_index_id];
        ui *previous_candidates = previous_edge.edge_ + previous_edge.offset_[previous_index_id];

        memcpy(valid_candidate_index[depth], previous_candidates, valid_candidates_count * sizeof(ui));

        ui temp_count;
        for (ui i = 1; i < bn_cnt[depth]; ++i)
        {
            VertexID current_bn = bn[depth][i];
            Edges &current_edge = *edge_matrix[current_bn][u];
            ui current_index_id = idx_embedding[current_bn];

            ui current_candidates_count = current_edge.offset_[current_index_id + 1] - current_edge.offset_[current_index_id];
            ui *current_candidates = current_edge.edge_ + current_edge.offset_[current_index_id];

            ComputeSetIntersection::ComputeCandidates(current_candidates, current_candidates_count, valid_candidate_index[depth], valid_candidates_count,
                                                      temp_buffer, temp_count);
            std::swap(temp_buffer, valid_candidate_index[depth]);
            valid_candidates_count = temp_count;
        }

        idx_count[depth] = valid_candidates_count;
    }
#endif
}

enumResult EvaluateQuery::exploreVEQStyleOO(const Graph *data_graph, const Graph *query_graph,ui *&nodeId, TreeNode *tree,ui *order1,
                                            Edges ***edge_matrix, ui **candidates, ui *candidates_count,
                                            size_t output_limit_num, size_t &call_count,int TimeL,
                                            const std::unordered_map<VertexID, std::pair<std::set<VertexID>, std::set<VertexID>>> &full_constraints)
{   
    auto start = std::chrono::high_resolution_clock::now();
    enumResult s;
    ui query_vertices_num = query_graph->getVerticesCount();
    ui data_vertices_num = data_graph->getVerticesCount();
    ui** nec = new ui*[query_vertices_num];
    memset(nec, 0, sizeof(ui*)*query_vertices_num);
    computeNEC(query_graph, nec);
    int depth = 0;
    ui max_depth = query_graph->getVerticesCount();
    ui *embedding = new ui[query_vertices_num]; 
         
    ui embedding_cnt = 0;               
    ui *extendable = new ui[query_vertices_num]; 
    for (ui i = 0; i < query_vertices_num; ++i) {
        extendable[i] = tree[i].bn_count_;
    }
    ui *idx = new ui[max_depth];                   
    ui *idx_count = new ui[max_depth];             
    ui** valid_cans = new ui*[query_vertices_num]; 
    for (ui i = 0; i < query_vertices_num; i++) {
        valid_cans[i] = new ui[candidates_count[i]];
    }
    ui* valid_cans_count = new ui[query_vertices_num]; 
    bool* visited_u = new bool[query_vertices_num]; 
    memset(visited_u, 0, query_vertices_num*sizeof(bool));
    bool *visited_vertices = new bool[data_vertices_num]; 
    memset(visited_vertices, 0, data_vertices_num*sizeof(bool));
    bool *violated_symmetry_u = new bool[max_depth];
    ui **bn;
    ui *bn_count;
    generateBN(query_graph, order1, bn, bn_count);
    ui** TM = new ui*[query_vertices_num];   
    for (ui i = 0; i < query_vertices_num; i++) {
        TM[i] = new ui[candidates_count[i] + 1]; 
    }
    memset(TM[0], 0, sizeof(ui)*(candidates_count[0] + 1));
     
    std::vector<std::vector<ui>> vec_index(query_vertices_num);
    for (ui i = 0; i < query_vertices_num; i++) {
        vec_index[i].resize(candidates_count[i]);
        std::fill(vec_index[i].begin(), vec_index[i].end(), (ui)-1);
    }
    std::vector<std::vector<ui>> vec_set;  
    
    //computeNEC(query_graph, edge_matrix, candidates_count, candidates, vec_index, vec_set);
    
    std::vector<std::vector<ui>> pi_m_index(query_vertices_num);
    for (ui i = 0; i < query_vertices_num; i++) {
        pi_m_index[i].resize(candidates_count[i], (ui)-1);
    }
    std::vector<std::vector<ui>> pi_m;  
    ui* pi_m_count = new ui[max_depth]; 
    pi_m_count[0] = 0;
    ui vec_count = 0;
    std::vector<std::vector<ui>> dm(query_vertices_num);  
    //std::unordered_map<VertexID, VertexID> reverse_embedding;  
    #ifdef ENABLE_FAILING_SET
    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> ancestors;
    computeAncestor(query_graph, tree, order1, ancestors);
    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> vec_failing_set(max_depth);
    std::unordered_map<VertexID, VertexID> reverse_embedding;
    reverse_embedding.reserve(MAXIMUM_QUERY_GRAPH_SIZE * 2);
    #endif
    ui* order = new ui[max_depth]; 
    //reverse_embedding.reserve(query_vertices_num); 
    TimeL = TimeL * 1000;
    ui start_vertex = (ui)-1;
    for (ui i = 0; i < query_vertices_num; i++) {
        if (extendable[i] == 0) {
            start_vertex = i;
            break;
        }
    }
    assert(start_vertex != (ui)-1);
    order[depth] = start_vertex;
    visited_u[start_vertex] = true;
    idx[depth] = 0;
    valid_cans_count[start_vertex] = candidates_count[start_vertex];
    idx_count[depth] = valid_cans_count[start_vertex];
    for (ui i = 0; i < candidates_count[start_vertex]; i++) {
        valid_cans[start_vertex][i] = candidates[start_vertex][i];
    }
    std::fill(pi_m_index[start_vertex].begin(), pi_m_index[start_vertex].end(), (ui)-1);
    while (true) {
        while (idx[depth] < idx_count[depth]) {
            auto end = std::chrono::high_resolution_clock::now();
            auto ens = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
            if (ens > TimeL)
            { // 1000 1 sec
                goto EXIT;
            }
            VertexID u = order[depth];
            VertexID v = valid_cans[u][idx[depth]];
            ui v_index1 = 0;
            for (; v_index1 < candidates_count[u]; v_index1++)
            {
                if (candidates[u][v_index1] == v)
                    break;
            }
            if (vec_index[u][v_index1] == (ui)-1)
            {
                computeNEC_Adaptive(query_graph, edge_matrix, candidates_count, candidates, vec_index, vec_set, u, v_index1, vec_count);
                vec_count++;
            }
        #ifdef ENABLE_FAILING_SET
            if(depth==0){
            if (idx_count[depth] == 0)
            {
            vec_failing_set[depth] = ancestors[u];
            }
            else
            {
            vec_failing_set[depth].reset();
            }
            }
        #endif

            bool is_voilate_symmetry = false;
            if (!full_constraints.empty()) {
                is_voilate_symmetry = checkSingeVertexBySymmetryBreaking(depth, embedding, order, visited_u, v,
                                                                         full_constraints,
                                                                         violated_symmetry_u,
                                                                         max_depth);
            }

            if (visited_vertices[v] || is_voilate_symmetry) {
                TM[u][idx[depth]] = 0;

                if (visited_vertices[v]) {
                    VertexID con_u = reverse_embedding[v];

                    ui con_v_index = 0, v_index = 0;
                    for (; con_v_index < valid_cans_count[con_u]; con_v_index++) {
                        if (valid_cans[con_u][con_v_index] == v) break;
                    }
                    for (; v_index < candidates_count[u]; v_index++) {
                        if (candidates[u][v_index] == v) break;
                    } 
                    assert(con_v_index < valid_cans_count[con_u]);
                    assert(v_index < candidates_count[u]);                    
                    auto& con_uv_idx = pi_m_index[con_u][con_v_index];
                    for (ui i = 0; i < pi_m[con_uv_idx].size(); i++) {
                        bool f_in = false;
                        for (ui j = 0; j < vec_set[vec_index[u][v_index]].size(); j++) {
                            if (vec_set[vec_index[u][v_index]][j] == pi_m[con_uv_idx][i]) {
                                f_in = true;
                                break;
                            }
                        }
                        if (f_in == false) {
                            pi_m[con_uv_idx][i] = pi_m[con_uv_idx][pi_m[con_uv_idx].size() - 1];
                            pi_m[con_uv_idx].pop_back();
                            i--;
                        }
                    }
                }
                if (is_voilate_symmetry) {
                    for (ui i = 0; i < depth; ++i) {               
                        if (violated_symmetry_u[order[i]]) {
                            VertexID sym_u = order[i];
                            VertexID sym_v = embedding[sym_u];
                            ui sym_v_index = 0, v_index = 0;
                            for (; sym_v_index < valid_cans_count[sym_u]; sym_v_index++) {
                                if (valid_cans[sym_u][sym_v_index] == sym_v) break;
                            }
                            for (; v_index < candidates_count[u]; v_index++) {
                                if (candidates[u][v_index] == sym_v) break;
                            }

                            auto& sym_uv_idx = pi_m_index[sym_u][sym_v_index];

                            for (ui i = 0; i < pi_m[sym_uv_idx].size(); i++) {
                                bool f_in = false;
                                if ((vec_index[u][v_index]) == (ui)-1)
                                {
                                    computeNEC_Adaptive(query_graph, edge_matrix, candidates_count, candidates, vec_index, vec_set, u, v_index, vec_count);
                                    vec_count++;
                                }
                                for (ui j = 0; j < vec_set[vec_index[u][v_index]].size(); j++) {
                                    if (vec_set[vec_index[u][v_index]][j] == pi_m[sym_uv_idx][i]) {
                                        f_in = true;
                                        break;
                                    }
                                }
                                if (f_in == false) {
                                    pi_m[sym_uv_idx][i] = pi_m[sym_uv_idx][pi_m[sym_uv_idx].size() - 1];
                                    pi_m[sym_uv_idx].pop_back();
                                    i--;
                                }
                            }
                        }
                    }
                }
                idx[depth]++;
                #ifdef ENABLE_FAILING_SET
                vec_failing_set[depth] = ancestors[u];
                vec_failing_set[depth] |= ancestors[reverse_embedding[v]];
                vec_failing_set[depth - 1] |= vec_failing_set[depth];
                #endif
                continue;
            }
            if (pi_m_index[u][idx[depth]] != (ui)-1) {

                VertexID equ_v = pi_m[pi_m_index[u][idx[depth]]][0];
                ui equ_v_index = 0;
                for (; equ_v_index < idx_count[depth]; equ_v_index++) {
                    if (equ_v == valid_cans[u][equ_v_index]) break;
                }
                assert(equ_v_index < idx_count[depth]);                
                embedding_cnt += TM[u][equ_v_index];
                TM[u][candidates_count[u]] += TM[u][equ_v_index];
                #ifdef ENABLE_FAILING_SET
                if(depth - 1>=0){
                vec_failing_set[depth].set();
                vec_failing_set[depth - 1] |= vec_failing_set[depth];
                }else if (depth==0){
                    vec_failing_set[depth].set();
                }
                #endif
                idx[depth]++;
                continue;
            }
            reverse_embedding[v] = u;
            embedding[u] = v;
            visited_vertices[v] = true;
            ui cur_idx = idx[depth]++;
            pi_m_index[u][cur_idx] = pi_m_count[depth]++;
            ui v_index = 0;
            for (; v_index < candidates_count[u]; v_index++) {
                if (candidates[u][v_index] == v) break;
            }
            if ((vec_index[u][v_index]) == (ui)-1)
            {
                computeNEC_Adaptive(query_graph, edge_matrix, candidates_count, candidates, vec_index, vec_set, u, v_index, vec_count);
                vec_count++;
            } 
            auto& pi = vec_set[vec_index[u][v_index]];
            pi_m.push_back(pi); 
            assert(pi_m.size() == pi_m_count[depth]);
            dm[u].clear();                                  
            for (ui i = 0; i < depth; i++) {                
                bool va_in_pi = false, sec_empty = true;
                ui va_index = 0;
                VertexID ua = order[i];
                VertexID va = embedding[ua];
                for (; va_index < candidates_count[ua]; va_index++) {
                    if (candidates[ua][va_index] == va) break;
                }
                assert(va_index < candidates_count[ua]);
                if ((vec_index[ua][va_index]) == (ui)-1)
                {
                    computeNEC_Adaptive(query_graph, edge_matrix, candidates_count, candidates, vec_index, vec_set, ua, v_index, vec_count);
                    vec_count++;
                }
               
                auto& pi_a = vec_set[vec_index[ua][va_index]];
                for (ui j = 0; j < pi.size(); j++) { 
                    if (pi[j] == embedding[order[i]]) { 
                        va_in_pi = true;
                        break;
                    }
                    for (ui k = 0; k < pi_a.size(); k++) { 
                        if (pi_a[k] == pi[j]) {
                            sec_empty = false;
                            break;
                        }
                    }
                    if (sec_empty == false) break;
                }
                if (va_in_pi == false && sec_empty == false) { 
                    ui dm_size = dm[ua].size();
                    for (ui j = 0; j < pi.size(); j++) {
                        bool f_add = true;
                        for (ui k = 0; k < dm_size; k++) {
                            if (dm[ua][k] == pi[j]) {
                                f_add = false;
                                break;
                            }
                        }
                        if (f_add == true) {
                            dm[ua].push_back(pi[j]);
                        }
                    }
                }
            }
            if (depth == max_depth - 1) { 
                embedding_cnt += 1;
                visited_vertices[v] = false;
                visited_u[u] = false;
                                int ao = depth; //-1;
                ui vqo = order[ao];
                int UNPMtemp = 0;
                int UNPM2 = 0;
                while (ao >= 0)
                {
                    nodeId[embedding[vqo]]=1;

                        if (pi_m_index[vqo][idx[ao] - 1] != (ui)-1)
                        {
                            for (auto eq_node : pi_m[pi_m_index[vqo][idx[ao] - 1]])
                            {   
                                nodeId[eq_node]=1;
                            }
                        }
                    
                    ao--;
                    vqo = order[ao];

                }

                if (embedding_cnt >= output_limit_num)
                {
                    goto EXIT;
                }
                
                reverse_embedding.erase(embedding[u]);
                #ifdef ENABLE_FAILING_SET
                vec_failing_set[depth].set();
                vec_failing_set[depth - 1] |= vec_failing_set[depth];
                #endif
                TM[u][cur_idx] = 1;
                TM[u][candidates_count[u]]++;
                auto& uv_idx = pi_m_index[u][cur_idx];
                for (ui i = 0; i < pi_m[uv_idx].size(); i++) { 
                    bool f_del = false;
                    for (ui j = 0; j < dm[u].size(); j++) {
                        if (pi_m[uv_idx][i] == dm[u][j]) {
                            f_del = true;
                            break;
                        }
                    }
                    if (f_del == true) {
                        pi_m[uv_idx][i] = pi_m[uv_idx][pi_m[uv_idx].size() - 1];
                        pi_m[uv_idx].pop_back();
                        i--;
                    }
                }
                for (ui i = 1; i < pi_m[uv_idx].size(); i++) { 
                    ui v_equ_index = 0;
                    for (; v_equ_index < valid_cans_count[u]; v_equ_index++) {
                        if (valid_cans[u][v_equ_index] == pi_m[uv_idx][i])
                            break;
                    }
                    assert(v_equ_index < valid_cans_count[u]);
                    pi_m_index[u][v_equ_index] = uv_idx;
                }
                
                if (pi_m[uv_idx].size() > 1) {
                    auto pi_v_idx = std::find(pi_m[uv_idx].begin(), pi_m[uv_idx].end(), v);
                    assert(pi_v_idx != pi_m[uv_idx].end());
                    ui tmp = *pi_v_idx;
                    *pi_v_idx = pi_m[uv_idx][0];
                    pi_m[uv_idx][0] = tmp;
                }
                
            } else { 
                call_count += 1;
                depth += 1;
                order[depth] = generateNextU(data_graph, query_graph, candidates, candidates_count, valid_cans,
                                             valid_cans_count, extendable, nec, depth, embedding,
                                             edge_matrix, visited_vertices, visited_u, order, tree);
                
                if (order[depth] == (ui)-1) {
                     //vec_failing_set[depth - 1].reset();
                    break; 
                } else {
                    visited_u[order[depth]] = true;
                    idx[depth] =0 ;
                    idx_count[depth] = valid_cans_count[order[depth]];
                #ifdef ENABLE_FAILING_SET
                if (idx_count[depth] == 0)
                    {
                    vec_failing_set[depth - 1] = ancestors[order[depth]];
                    }
                    else
                    {
                    vec_failing_set[depth - 1].reset();
                    }
                    #endif
                    bool is_voilate_symmetry1=true;

                    while(is_voilate_symmetry1 && idx[depth]<idx_count[depth]&&false){
                    VertexID u = order[depth];
                    VertexID v = valid_cans[u][idx[depth]];
                    is_voilate_symmetry1 = checkSingeVertexBySymmetryBreaking(depth, embedding, order, visited_u, v,
                                                                         full_constraints,
                                                                         violated_symmetry_u,
                                                                         max_depth);
                    if(is_voilate_symmetry1)
                        idx[depth]++;
                    }
                }
                memset(TM[order[depth]], 0, sizeof(ui)*(candidates_count[order[depth]] + 1));
                std::fill(pi_m_index[order[depth]].begin(), pi_m_index[order[depth]].end(), (ui)-1);
                pi_m_count[depth] = pi_m_count[depth - 1];
            }
        }
        
        depth -= 1;
        
        if (depth < 0)
            break;
        VertexID u = order[depth];
        reverse_embedding.erase(embedding[u]);
        visited_vertices[embedding[u]] = false;
        ui cur_idx = idx[depth] - 1;
        restoreExtendableVertex(tree, u, extendable);
        

        if (order[depth + 1] != (ui)-1) { 
            VertexID last_u = order[depth + 1];          
            visited_u[last_u] = false;
            if (nec[last_u] != NULL) (*(nec[last_u]))++;
        }
                    #ifdef ENABLE_FAILING_SET
            if (depth != 0)
            {
                if (!vec_failing_set[depth].test(u))
                {   
                    vec_failing_set[depth - 1] = vec_failing_set[depth];
                    idx[depth] = idx_count[depth];
                    cur_idx=idx[depth];
                    continue;
                    //TM[u][cur_idx] = 0;
                }
                else
                {
                    vec_failing_set[depth - 1] |= vec_failing_set[depth];
                }
            }
        #endif
        if (order[depth + 1] != (ui)-1) {
            TM[u][cur_idx] = TM[order[depth + 1]][candidates_count[order[depth + 1]]];
        } else {
            TM[u][cur_idx] = 0;
        }
        TM[u][candidates_count[u]] += TM[u][cur_idx];
        auto& uv_idx = pi_m_index[u][cur_idx];
        if (TM[u][cur_idx] != 0) {
            
            for (ui i = 0; i < pi_m[uv_idx].size(); i++) {
                
            }
            
            for (ui i = 0; i < pi_m[uv_idx].size(); i++) {
                bool f_del = false;
                for (ui j = 0; j < dm[u].size(); j++) {
                    if (pi_m[uv_idx][i] == dm[u][j]) {
                        f_del = true;
                        break;
                    }
                }
                if (f_del == true) {
                    pi_m[uv_idx][i] = pi_m[uv_idx][pi_m[uv_idx].size() - 1];
                    pi_m[uv_idx].pop_back();
                    i--;
                }
            }
            
            for (ui i = 0; i < pi_m[uv_idx].size(); i++) {        
            }
            if (pi_m[uv_idx].size() > 1) {
                auto pi_v_idx = std::find(pi_m[uv_idx].begin(), pi_m[uv_idx].end(), embedding[u]);
                assert(pi_v_idx != pi_m[uv_idx].end());
                ui tmp = *pi_v_idx;
                *pi_v_idx = pi_m[uv_idx][0];
                pi_m[uv_idx][0] = tmp;
            }      
        }
        for (ui i = 1; i < pi_m[uv_idx].size(); i++) {
            ui v_equ_index = 0;
            for (; v_equ_index < valid_cans_count[u]; v_equ_index++) {
                if (valid_cans[u][v_equ_index] == pi_m[uv_idx][i])
                    break;
            }
            assert(v_equ_index < valid_cans_count[u]);
            
            pi_m_index[u][v_equ_index] = uv_idx;
        }
        pi_m.resize(pi_m_count[depth]);

       
    }

    EXIT:
    
    for (ui i = 0; i < query_vertices_num; i++) {
        
    }
    delete []nec;
    delete []embedding;
    delete []visited_u;
    delete []visited_vertices;
    delete []order1;
    delete []extendable;
    delete []idx;
    delete []idx_count;
    for (ui i = 0; i < query_vertices_num; i++) {
        delete []valid_cans[i];
    }
    delete []valid_cans;
    delete []valid_cans_count;
    delete[] TM;
    int countSol=0;
    for (int i=0;i<data_graph->getVerticesCount();i++){
        countSol+=nodeId[i];
    }
    s.embedding_cnt = embedding_cnt;
    s.Can_embed = countSol;
    return s;
}

enumResult EvaluateQuery::exploreVEQStyleOV(const Graph *data_graph, const Graph *query_graph, TreeNode *tree,
                                            Edges ***edge_matrix, ui **candidates, ui *candidates_count,
                                            size_t output_limit_num, size_t &call_count, int TimeL, ui *order,
                                             const std::unordered_map<VertexID, std::pair<std::set<VertexID>, std::set<VertexID>>> &ordered_constraints)
{

    ui query_vertices_num = query_graph->getVerticesCount();
    ui data_vertices_num = data_graph->getVerticesCount();
    ui **nec = new ui *[query_vertices_num];
    memset(nec, 0, sizeof(ui *) * query_vertices_num);
    computeNEC(query_graph, nec);
    enumResult s;

    int depth = 0;
    ui max_depth = query_graph->getVerticesCount();
    ui *embedding = new ui[query_vertices_num];

    ui embedding_cnt = 0;
    ui *extendable = new ui[query_vertices_num];
    for (ui i = 0; i < query_vertices_num; ++i)
    {
        extendable[i] = tree[i].bn_count_;
    }
    ui *idx = new ui[max_depth];
    ui *idx_count = new ui[max_depth];
    ui **valid_cans = new ui *[query_vertices_num];
    for (ui i = 0; i < query_vertices_num; i++)
    {
        valid_cans[i] = new ui[candidates_count[i]];
    }
    ui *valid_cans_count = new ui[query_vertices_num];
    bool *visited_u = new bool[query_vertices_num];
    memset(visited_u, 0, query_vertices_num * sizeof(bool));
    bool *visited_vertices = new bool[data_vertices_num];
    memset(visited_vertices, 0, data_vertices_num * sizeof(bool));
    bool *violated_symmetry_u = new bool[max_depth];

#ifdef ENABLE_EQUIVALENT_SET
    ui **TM = new ui *[query_vertices_num];
    for (ui i = 0; i < query_vertices_num; i++)
    {
        TM[i] = new ui[candidates_count[i] + 1];
    }
    memset(TM[0], 0, sizeof(ui) * (candidates_count[0] + 1));

    std::vector<std::vector<ui>> vec_index(query_vertices_num);
    for (ui i = 0; i < query_vertices_num; i++)
    {
        vec_index[i].resize(candidates_count[i]);
        std::fill(vec_index[i].begin(), vec_index[i].end(), (ui)-1);
    }
    std::vector<std::vector<ui>> vec_set;
    computeNEC(query_graph, edge_matrix, candidates_count, candidates, vec_index, vec_set);

    std::vector<std::vector<ui>> pi_m_index(query_vertices_num);
    for (ui i = 0; i < query_vertices_num; i++)
    {
        pi_m_index[i].resize(candidates_count[i], (ui)-1);
    }
    std::vector<std::vector<ui>> pi_m;
    ui *pi_m_count = new ui[max_depth];
    pi_m_count[0] = 0;

    std::vector<std::vector<ui>> dm(query_vertices_num);
    std::unordered_map<VertexID, VertexID> reverse_embedding;
    reverse_embedding.reserve(query_vertices_num);
#endif

    ui start_vertex = order[0];
    visited_u[start_vertex] = true;

    idx[depth] = 0;
    valid_cans_count[start_vertex] = candidates_count[start_vertex];
    idx_count[depth] = valid_cans_count[start_vertex];
    for (ui i = 0; i < candidates_count[start_vertex]; i++)
    {
        valid_cans[start_vertex][i] = candidates[start_vertex][i];
    }
#ifdef ENABLE_EQUIVALENT_SET
    std::fill(pi_m_index[start_vertex].begin(), pi_m_index[start_vertex].end(), (ui)-1);

#endif

    while (true)
    {
        while (idx[depth] < idx_count[depth])
        {

#ifdef ENABLE_DYNAMIC_CANS

            if (idx[order[depth]] != valid_cans[order[depth]].begin())
            {
                RestoreValidCans(query_graph, data_graph, visited_u, order[depth], embedding[order[depth]], valid_cans);
            }
#endif
            VertexID u = order[depth];
            VertexID v = valid_cans[u][idx[depth]];

            bool is_voilate_symmetry = false;


            if (visited_vertices[v])
            {
#ifdef ENABLE_EQUIVALENT_SET
                TM[u][idx[depth]] = 0;

                if (visited_vertices[v])
                {
                    VertexID con_u = reverse_embedding[v];

                    ui con_v_index = 0, v_index = 0;
                    for (; con_v_index < valid_cans_count[con_u]; con_v_index++)
                    {
                        if (valid_cans[con_u][con_v_index] == v)
                            break;
                    }
                    for (; v_index < candidates_count[u]; v_index++)
                    {
                        if (candidates[u][v_index] == v)
                            break;
                    }

                    assert(con_v_index < valid_cans_count[con_u]);
                    assert(v_index < candidates_count[u]);

                    auto &con_uv_idx = pi_m_index[con_u][con_v_index];

                    for (ui i = 0; i < pi_m[con_uv_idx].size(); i++)
                    {
                        bool f_in = false;
                        for (ui j = 0; j < vec_set[vec_index[u][v_index]].size(); j++)
                        {
                            if (vec_set[vec_index[u][v_index]][j] == pi_m[con_uv_idx][i])
                            {
                                f_in = true;
                                break;
                            }
                        }
                        if (f_in == false)
                        {
                            pi_m[con_uv_idx][i] = pi_m[con_uv_idx][pi_m[con_uv_idx].size() - 1];
                            pi_m[con_uv_idx].pop_back();
                            i--;
                        }
                    }
                }

#endif
                idx[depth]++;
                continue;
            }
#ifdef ENABLE_EQUIVALENT_SET
            if (pi_m_index[u][idx[depth]] != (ui)-1)
            {

                VertexID equ_v = pi_m[pi_m_index[u][idx[depth]]][0];
                ui equ_v_index = 0;
                for (; equ_v_index < idx_count[depth]; equ_v_index++)
                {
                    if (equ_v == valid_cans[u][equ_v_index])
                        break;
                }
                assert(equ_v_index < idx_count[depth]);

                embedding_cnt += TM[u][equ_v_index];
                TM[u][candidates_count[u]] += TM[u][equ_v_index];
                idx[depth]++;
                continue;
            }
            reverse_embedding[v] = u;
#endif
            embedding[u] = v;
            visited_vertices[v] = true;
            ui cur_idx = idx[depth]++;
#ifdef ENABLE_EQUIVALENT_SET
            pi_m_index[u][cur_idx] = pi_m_count[depth]++;
            ui v_index = 0;
            for (; v_index < candidates_count[u]; v_index++)
            {
                if (candidates[u][v_index] == v)
                    break;
            }

            auto &pi = vec_set[vec_index[u][v_index]];
            pi_m.push_back(pi);
            assert(pi_m.size() == pi_m_count[depth]);
            dm[u].clear();
            for (ui i = 0; i < depth; i++)
            {
                bool va_in_pi = false, sec_empty = true;
                ui va_index = 0;
                VertexID ua = order[i];
                VertexID va = embedding[ua];
                for (; va_index < candidates_count[ua]; va_index++)
                {
                    if (candidates[ua][va_index] == va)
                        break;
                }
                assert(va_index < candidates_count[ua]);
                auto &pi_a = vec_set[vec_index[ua][va_index]];
                for (ui j = 0; j < pi.size(); j++)
                {
                    if (pi[j] == embedding[order[i]])
                    {
                        va_in_pi = true;
                        break;
                    }
                    for (ui k = 0; k < pi_a.size(); k++)
                    {
                        if (pi_a[k] == pi[j])
                        {
                            sec_empty = false;
                            break;
                        }
                    }
                    if (sec_empty == false)
                        break;
                }
                if (va_in_pi == false && sec_empty == false)
                {
                    ui dm_size = dm[ua].size();
                    for (ui j = 0; j < pi.size(); j++)
                    {
                        bool f_add = true;
                        for (ui k = 0; k < dm_size; k++)
                        {
                            if (dm[ua][k] == pi[j])
                            {
                                f_add = false;
                                break;
                            }
                        }
                        if (f_add == true)
                        {
                            dm[ua].push_back(pi[j]);
                        }
                    }
                }
            }
#endif
            if (depth == max_depth - 1)
            {
                embedding_cnt += 1;
                visited_vertices[v] = false;
                visited_u[u] = false;

                if (embedding_cnt >= output_limit_num)
                {
                    goto EXIT;
                }
#ifdef ENABLE_EQUIVALENT_SET
                reverse_embedding.erase(v);
                TM[u][cur_idx] = 1;
                TM[u][candidates_count[u]]++;
                auto &uv_idx = pi_m_index[u][cur_idx];
                for (ui i = 0; i < pi_m[uv_idx].size(); i++)
                {
                    bool f_del = false;
                    for (ui j = 0; j < dm[u].size(); j++)
                    {
                        if (pi_m[uv_idx][i] == dm[u][j])
                        {
                            f_del = true;
                            break;
                        }
                    }
                    if (f_del == true)
                    {
                        pi_m[uv_idx][i] = pi_m[uv_idx][pi_m[uv_idx].size() - 1];
                        pi_m[uv_idx].pop_back();
                        i--;
                    }
                }
                for (ui i = 1; i < pi_m[uv_idx].size(); i++)
                {
                    ui v_equ_index = 0;
                    for (; v_equ_index < valid_cans_count[u]; v_equ_index++)
                    {
                        if (valid_cans[u][v_equ_index] == pi_m[uv_idx][i])
                            break;
                    }

                    //assert(v_equ_index < valid_cans_count[u]);
                    pi_m_index[u][v_equ_index] = uv_idx;
                }

                if (pi_m[uv_idx].size() > 1)
                {
                    auto pi_v_idx = std::find(pi_m[uv_idx].begin(), pi_m[uv_idx].end(), v);
                    assert(pi_v_idx != pi_m[uv_idx].end());
                    ui tmp = *pi_v_idx;
                    *pi_v_idx = pi_m[uv_idx][0];
                    pi_m[uv_idx][0] = tmp;
                }
#endif
            }
            else
            {
                call_count += 1;
                depth += 1;
                //
                ComputeValidCans(data_graph, query_graph, candidates, candidates_count, valid_cans,
                                 valid_cans_count, embedding, order[depth], visited_u);
                
                if (valid_cans_count==0)
                break;
                pruneCandidatesIndexBySymmetryBreaking(depth, embedding, order, candidates_count, candidates,
                                                       valid_cans_count, valid_cans, ordered_constraints);
                if (valid_cans_count==0)
                break;
                visited_u[order[depth]] = true;
                idx[depth] = 0;
                idx_count[depth] = valid_cans_count[order[depth]];

#ifdef ENABLE_EQUIVALENT_SET
                memset(TM[order[depth]], 0, sizeof(ui) * (candidates_count[order[depth]] + 1));
                std::fill(pi_m_index[order[depth]].begin(), pi_m_index[order[depth]].end(), (ui)-1);
                pi_m_count[depth] = pi_m_count[depth - 1];
#endif
            }
        }

        depth -= 1;

        if (depth < 0)
            break;
        VertexID u = order[depth];
        ui cur_idx = idx[depth] - 1;
        visited_vertices[embedding[u]] = false;

        if (order[depth + 1] != (ui)-1)
        {
            VertexID last_u = order[depth + 1];

            visited_u[last_u] = false;
            if (nec[last_u] != NULL)
                (*(nec[last_u]))++;
#ifdef ENABLE_DYNAMIC_CANS
            RestoreValidCans(query_graph, data_graph, visited_u, last_u, last_v, valid_cans);
#endif
        }
#ifdef ENABLE_EQUIVALENT_SET
        if (order[depth + 1] != (ui)-1)
        {
            TM[u][cur_idx] = TM[order[depth + 1]][candidates_count[order[depth + 1]]];
        }
        else
        {
            TM[u][cur_idx] = 0;
        }
        TM[u][candidates_count[u]] += TM[u][cur_idx];
        auto &uv_idx = pi_m_index[u][cur_idx];
        if (TM[u][cur_idx] != 0)
        {

            for (ui i = 0; i < pi_m[uv_idx].size(); i++)
            {
            }

            for (ui i = 0; i < pi_m[uv_idx].size(); i++)
            {
                bool f_del = false;
                for (ui j = 0; j < dm[u].size(); j++)
                {
                    if (pi_m[uv_idx][i] == dm[u][j])
                    {
                        f_del = true;
                        break;
                    }
                }
                if (f_del == true)
                {
                    pi_m[uv_idx][i] = pi_m[uv_idx][pi_m[uv_idx].size() - 1];
                    pi_m[uv_idx].pop_back();
                    i--;
                }
            }

            for (ui i = 0; i < pi_m[uv_idx].size(); i++)
            {
            }

            if (pi_m[uv_idx].size() > 1)
            {
                auto pi_v_idx = std::find(pi_m[uv_idx].begin(), pi_m[uv_idx].end(), embedding[u]);
                assert(pi_v_idx != pi_m[uv_idx].end());
                ui tmp = *pi_v_idx;
                *pi_v_idx = pi_m[uv_idx][0];
                pi_m[uv_idx][0] = tmp;
            }
        }
        for (ui i = 1; i < pi_m[uv_idx].size(); i++)
        {
            ui v_equ_index = 0;
            for (; v_equ_index < valid_cans_count[u]; v_equ_index++)
            {
                if (valid_cans[u][v_equ_index] == pi_m[uv_idx][i])
                    break;
            }
            assert(v_equ_index < valid_cans_count[u]);

            pi_m_index[u][v_equ_index] = uv_idx;
        }
        pi_m.resize(pi_m_count[depth]);
#endif
    }

EXIT:

    for (ui i = 0; i < query_vertices_num; i++)
    {
    }
    delete[] nec;
    delete[] embedding;
    delete[] visited_u;
    delete[] visited_vertices;
    // delete []order;
    delete[] extendable;
    delete[] idx;
    delete[] idx_count;
    for (ui i = 0; i < query_vertices_num; i++)
    {
        delete[] valid_cans[i];
    }
    delete[] valid_cans;
    delete[] valid_cans_count;
#ifdef ENABLE_EQUIVALENT_SET
    delete[] TM;
#endif
    s.embedding_cnt = embedding_cnt;
    s.Can_embed = 0;
    return s;
}


enumResult
EvaluateQuery::DSQLTOPK(const Graph *data_graph, const Graph *query_graph,ui *&nodeId, Edges ***edge_matrix, ui **candidates,
                    ui *candidates_count, ui *order, size_t output_limit_num, size_t &call_count, size_t &valid_vtx_cnt, int TimeL, int FairT,
                    const std::unordered_map<VertexID, std::pair<std::set<VertexID>, std::set<VertexID>>> &ordered_constraints)
    {

    auto start = std::chrono::high_resolution_clock::now();
    enumResult s;
    ui **bn;
    ui *bn_count;
    generateBN(query_graph, order, bn, bn_count);
    int UNPM = 0;
    ui *idx;
    ui *idx_count;
    ui *embedding;
    ui *idx_embedding;
    ui *temp_buffer;
    ui **valid_candidate_idx;
    int i=0;
    bool *visited_vertices;
        int count1=0;
    int count2=0;
    int K_size2=0;
        ui tempstore=0;
    ui nbr_cnt;
    ui nbr_cnt1;
    vector<int> indices(1);        
    //vector<int> indices1(max_depth); 
    allocateBuffer(data_graph, query_graph, candidates_count, idx, idx_count, embedding, idx_embedding,
                   temp_buffer, valid_candidate_idx, visited_vertices);
        int cc1=0;
    int cc2=0;
    size_t embedding_cnt = 0;
    int cur_depth = 0;
    ui max_depth = query_graph->getVerticesCount();
    VertexID start_vertex = order[0];
    ui order_matrix[max_depth];
    int index_comb[max_depth];
    bool is_used[max_depth];
    bool is_used_copy[max_depth];
    std::vector<std::vector<ui>> Cand_vec;
    Cand_vec.resize(max_depth);
    //unordered_set<ui> EmbSum;
    bool Match_BA[max_depth] = {false};
    bool conflicts[max_depth][max_depth] = {false};
    bool conflictT[max_depth][max_depth] = {false};
    std::vector<std::unordered_set<int>> nodeIDbad(max_depth);
        for (int i=0;i<max_depth;i++){
        for (int j=0;j<max_depth;j++)
        conflicts[i][j]=false;
    }
    for (int i=0;i<max_depth;i++){
    const VertexID *nbrs = query_graph->getVertexNeighbors(i, nbr_cnt);
            for (int j=0;j<nbr_cnt;j++){
                conflicts[i][nbrs[j]]=true;
                conflicts[nbrs[j]][i]=true;//not needed
            }
    }

    idx[cur_depth] = 0;
    idx_count[cur_depth] = candidates_count[start_vertex];

    for (ui i = 0; i < idx_count[cur_depth]; ++i)
    {
        valid_candidate_idx[cur_depth][i] = i;
    }    
    std::unordered_map<VertexID, VertexID> reverse_embedding;
    reverse_embedding.reserve(MAXIMUM_QUERY_GRAPH_SIZE * 2);
    bool hasNext = true;
    bool hasNext1 = true;
    ui K_size=1;
    TimeL = TimeL * 1000;
    vector<int> combination;
    double ens = 0;
    std::vector<std::vector<int>> vec(max_depth);
    ui rev_order[max_depth] = {0};


    int qsiz = query_graph->getVerticesCount();
    int ksize = 10000;
    int SolPos[ksize+1][qsiz];
    std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, std::greater<>> pq;
    int tempPos[qsiz];
    pq.push({0,0});
    int heapcount = 0;
                for (int i=0;i<max_depth;i++){
            is_used[i]=false;
            is_used_copy[i]=false;
            index_comb[i]=0;
            }
            is_used[0]=order[0];
            order_matrix[0] = order[0];
            is_used_copy[order[0]]=true;
            count1=1;
            count2=0;
            //count1++;
            int countend=max_depth-1;
            if (query_graph->getVertexDegree(order_matrix[0])==1){
                const VertexID *nbrs = query_graph->getVertexNeighbors(order_matrix[count2], nbr_cnt);
                ui vertexN=nbrs[0];
                is_used_copy[vertexN]=true;
                order_matrix[count1]=vertexN;
                count1++;
                count2++;
                while (query_graph->getVertexDegree(vertexN)==1){
                    const VertexID *nbrs = query_graph->getVertexNeighbors(order_matrix[count2], nbr_cnt);    
                    vertexN=nbrs[0];
                    is_used_copy[vertexN]=true;
                    order_matrix[count1]=vertexN;
                    count1++;
                    count2++;
                }
            }
        while(count1<=countend){
            if (is_used_copy[order_matrix[count2]]){
            const VertexID *nbrs = query_graph->getVertexNeighbors(order_matrix[count2], nbr_cnt);
             
            for (int i=0;i<nbr_cnt;i++){
                if(is_used_copy[nbrs[i]]==false){
                    is_used_copy[nbrs[i]]=true;
                    nbr_cnt1 = query_graph->getVertexDegree(order_matrix[count2]);
                    if(nbr_cnt1>1){                    
                        order_matrix[count1]=nbrs[i];
                        count1++;
                    }
                    
                    else{
                        //order_matrix[count1]=nbrs[i];
                        //count1++;
                        order_matrix[countend]=nbrs[i];
                        countend--;
                    }
                    
                    
                }
            }}count2++;
        }

    generateBNLM(query_graph, order_matrix, bn, bn_count);
        for (int i=0;i<max_depth;i++){
        rev_order[order_matrix[i]]=i;
    }
        for (int i=max_depth-2;i>0;i--){
        int pp=-1;
        const VertexID *nbrs = query_graph->getVertexNeighbors(order_matrix[i], nbr_cnt);
        for (int j=0;j<nbr_cnt;j++){
            if (rev_order[order_matrix[nbrs[j]]]>pp&&rev_order[order_matrix[nbrs[j]]]<i)
            pp=rev_order[order_matrix[nbrs[j]]];
        }
        if (pp!=-1)
        vec[pp].push_back(i-1);
    }
    while (true)
    {
        while (idx[cur_depth] < idx_count[cur_depth])
        {   
            ui valid_idx = valid_candidate_idx[cur_depth][idx[cur_depth]];
            VertexID u = order_matrix[cur_depth];
            VertexID v = candidates[u][valid_idx];
            auto end = std::chrono::high_resolution_clock::now();
            ens = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
            if (ens > TimeL)
            { // 1000 1 sec
                goto EXIT;
            }

            if (visited_vertices[v])
            {
                idx[cur_depth] += 1;
    
                conflictT[u][reverse_embedding[v]]=true;
                conflictT[reverse_embedding[v]][u]=true;
                //cout<<"conf"<<u<<","<<reverse_embedding[v]<<endl;
                continue;
            }
            if(nodeId[v]==1){
                idx[cur_depth] += 1;
                //Match_BA[cur_depth]=true;
                for (int b=max_depth-1;b>=0;b--){
                Match_BA[b]=true;
                }
                //conflictT[u][reverse_embedding[v]]=true;
                continue;
            }
            if((nodeIDbad[cur_depth].find(v) != nodeIDbad[cur_depth].end())){
                idx[cur_depth] += 1;
                continue;
            }

            Match_BA[cur_depth]=false;
            embedding[u] = v;
            idx_embedding[u] = valid_idx;
            visited_vertices[v] = true;
            idx[cur_depth] += 1;
            //if (cur_depth!=0&&cur_depth!=max_depth-1)

            reverse_embedding[v] = u;
            if (cur_depth == max_depth - 1)
            {
                UNPM = 0;
                if(nodeId[embedding[order_matrix[cur_depth]]]==0){
                    UNPM++;
                    nodeId[embedding[order_matrix[cur_depth]]]=1;
                }
                embedding_cnt += 1;
                visited_vertices[v] = false;
                int ao = cur_depth; //-1;
                Match_BA[cur_depth]=true;
                tempPos[ao]=embedding[order_matrix[cur_depth]];
                ui vqo = order_matrix[ao];
                nodeId[v]=1;
                while (ao > 0)
                {   
                    Match_BA[ao]=true;
                    if(nodeId[embedding[vqo]]==0){
                    UNPM++;
                    nodeId[embedding[vqo]]=1;
                    //nodeIDbad[ao].clear();
                }
                    tempPos[ao]=embedding[vqo];
                    idx[ao] = idx_count[ao];
                    ao--;
                    reverse_embedding.erase(embedding[vqo]);
                    visited_vertices[embedding[vqo]] = false;
                    vqo = order_matrix[ao];

                }
                Match_BA[ao]=true;
                tempPos[ao]=embedding[vqo];
                if(nodeId[embedding[vqo]]==0){
                    UNPM++;
                    nodeId[embedding[vqo]]=1;
                    }
                if (UNPM > pq.top().first)
                            {                                       
                                int addPos=0;
                                if (heapcount >= ksize)
                                {
                                   auto [value, index] = pq.top();      
                                    pq.pop();
                                    pq.push({UNPM,index});
                                    addPos=index;
                                    
                                }else{
                                     pq.push({UNPM,heapcount});
                                     addPos=heapcount;            
                                }
                            for (int aa=0;aa<qsiz;aa++){
                                SolPos[addPos][aa]=tempPos[aa];
                                    }
                                heapcount++;
                            }
                //reverse_embedding.erase(embedding[u]);
                //visited_vertices[embedding[vqo]] = false;
            }
            else
            {
                call_count += 1;
                cur_depth += 1;
                for (int i=0;i<max_depth;i++){
                conflictT[order_matrix[cur_depth]][order_matrix[i]]=conflicts[order_matrix[cur_depth]][order_matrix[i]];
                //conflictT[order[cur_depth]][order[i]]=false;
                }
                idx[cur_depth] = 0;
                
                generateValidCandidateIndex(cur_depth, idx_embedding, idx_count, valid_candidate_idx, edge_matrix, bn,
                                            bn_count, order_matrix, temp_buffer);

            //if(idx[cur_depth]!=idx_count[cur_depth]){
            //    Match_BA[cur_depth]=false;
            //}else{Match_BA[cur_depth]=true;
           // }


            //cout<<"hmm"<<endl;
           // ui valid_idx1 = valid_candidate_idx[cur_depth][idx[cur_depth]];
           // while(nodeId[candidates[order[cur_depth]][valid_idx1]]==1&&idx[cur_depth]<idx_count[cur_depth]){
          ///  idx[cur_depth]++;
           // valid_idx1 = valid_candidate_idx[cur_depth][idx[cur_depth]];    }
            }
        }


        cur_depth -= 1;

        if (cur_depth < 0)
            break;
        else
        {
            //get conflict level
            int testdepth=cur_depth-1;
            while(testdepth>=0){
            if(conflictT[order_matrix[cur_depth]][order_matrix[testdepth]]==false&&conflicts[order_matrix[cur_depth]][order_matrix[testdepth]]==false)
                break;
            testdepth--;
            }
            //testdepth=depth of conflict
        if (testdepth==0)    
        if(Match_BA[cur_depth]==false&&cur_depth!=max_depth - 1){
            nodeIDbad[cur_depth].insert(embedding[order_matrix[cur_depth]]);
        }
        else if(testdepth>0&&conflictT[order_matrix[testdepth]][order_matrix[testdepth-1]]==false&&conflicts[order_matrix[testdepth]][order_matrix[testdepth-1]]==false){
         if(Match_BA[cur_depth]==false&&cur_depth!=max_depth - 1){
            nodeIDbad[cur_depth].insert(embedding[order_matrix[cur_depth]]);
        }   
        }
        for (int mm=0;mm<vec[cur_depth].size();mm++){
            nodeIDbad[vec[cur_depth][mm]].clear();
        }
            VertexID u = order_matrix[cur_depth];
            reverse_embedding.erase(embedding[u]);

            if (cur_depth != 0&&Match_BA[cur_depth]==false&&cur_depth < max_depth-2)
            //if (cur_depth != 0&&Match_BA[cur_depth]==false)
            {
                //for (int k=cur_depth+1;k<max_depth-1;k++)
                //if(conflictT[order[cur_depth+1]][u]==false){
                bool Prun=true;
                ui test_depth=max_depth-1;
                while(test_depth>cur_depth){
                  if(conflictT[u][order_matrix[test_depth]]==false&&conflicts[u][order_matrix[test_depth]]==false)  
                        test_depth--;
                    else{
                        Prun=false;
                        break;
                    }
                }
                if (Prun==true)
                idx[cur_depth] = idx_count[cur_depth];
                //if(conflictT[u][order_matrix[cur_depth+1]]==false&&conflicts[u][order_matrix[cur_depth+1]]==false){
               //     idx[cur_depth] = idx_count[cur_depth];
               // }
            }
            visited_vertices[embedding[u]] = false;

        }
    }
   //goto EXIT;
    hasNext = true;
     K_size=1;
        //for (int i=0;i<max_depth;i++)
        //cout<<order[i];
    
    //std::vector<std::vector<ui>> Cand_vec;
       // 
    while (K_size<max_depth){
        indices.resize(K_size);
        for (ui i = 0; i < K_size; ++i) {
        indices[i] = i;
        }
        hasNext = true;
        
    for (int ik = 0; ik < max_depth; ik++) {
        
        Cand_vec[ik].clear();
        for (int i=0;i<candidates_count[ik];i++)
          if(nodeId[candidates[ik][i]]==1){
            //Cand_vec[ik].push_back(candidates[ik][i]);
            Cand_vec[ik].push_back(i);
                }
        }
        while (hasNext==true) {
            cc1++;
            
            count1=0;
            count2=0;
            K_size2=0;
            tempstore=0;
            nbr_cnt=0;
            for (int i=0;i<max_depth;i++){
            is_used[i]=false;
            is_used_copy[i]=false;
            index_comb[i]=0;
            }
            combination.clear();
            
            for (ui i = 0; i < K_size; ++i) {
                combination.push_back(order[indices[i]]);
                is_used[order[indices[i]]] = true;
            }

            order_matrix[0] = combination[0];
            is_used_copy[combination[0]]=true;
            count1=1;
            count2=0;
            //count1++;
            countend=max_depth-1;
            if (query_graph->getVertexDegree(order_matrix[0])==1){
                
                const VertexID *nbrs = query_graph->getVertexNeighbors(order_matrix[count2], nbr_cnt);
                ui vertexN=nbrs[0];
                is_used_copy[vertexN]=true;
                order_matrix[count1]=vertexN;
                count1++;
                count2++;
                while (query_graph->getVertexDegree(vertexN)==1){
                    const VertexID *nbrs = query_graph->getVertexNeighbors(order_matrix[count2], nbr_cnt);    
                    vertexN=nbrs[0];
                    is_used_copy[vertexN]=true;
                    order_matrix[count1]=vertexN;
                    count1++;
                    count2++;
                }
            }

        while(count1<max_depth){
            if (is_used_copy[order_matrix[count2]]){
            const VertexID *nbrs = query_graph->getVertexNeighbors(order_matrix[count2], nbr_cnt);
            for (int i=0;i<nbr_cnt;i++){
                if(is_used_copy[nbrs[i]]==false){
                    is_used_copy[nbrs[i]]=true;
                    nbr_cnt1 = query_graph->getVertexDegree(order_matrix[count2]);
                    if(nbr_cnt1>1){                        
                        order_matrix[count1]=nbrs[i];
                        count1++;
                    }
                    
                    else{
                        order_matrix[countend]=nbrs[i];
                        countend--;
                    }
                    
                    
                }
            }}count2++;
        }
        for (int i=0;i<max_depth;i++){
        rev_order[order_matrix[i]]=i;
        vec[i].clear();
    }
    for (int i=max_depth-2;i>0;i--){
        int pp=-1;
        const VertexID *nbrs = query_graph->getVertexNeighbors(order_matrix[i], nbr_cnt);
        for (int j=0;j<nbr_cnt;j++){
            if (rev_order[order_matrix[nbrs[j]]]>pp&&rev_order[order_matrix[nbrs[j]]]<i)
            pp=rev_order[order_matrix[nbrs[j]]];
        }
        if (pp!=-1)
        vec[pp].push_back(i-1);
    }
        count1=max_depth;
        count2=0;
        i=1;


        K_size2=0;
        i=0;
        while(is_used[order_matrix[i]]==true){
            K_size2++;
            i++;
        }
        
        i=0;
        
       // generateBN1(query_graph, order_matrix, bn, bn_count);
       generateBNLM(query_graph, order_matrix, bn, bn_count);
        allocateBufferLM(data_graph, query_graph, candidates_count, idx, idx_count, embedding, idx_embedding,
                   temp_buffer, valid_candidate_idx, visited_vertices);

        //std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> vec_failing_set(max_depth);
        while (hasNext1==true){
            cc2++;
            for (int aa1=0;aa1<max_depth;aa1++)
            nodeIDbad[aa1].clear();
            start_vertex = order_matrix[0];
            cur_depth=0;
            idx[cur_depth] = 0;
            idx_count[cur_depth] = candidates_count[start_vertex];
            Match_BA[0] = false;
        
        /*
        for (ui i = 0; i < idx_count[cur_depth]; ++i)
        { 
        if(candidates[start_vertex][i]==Cand_vec[start_vertex][index_comb[start_vertex]]) {
            valid_candidate_idx[cur_depth][0] = i;            
            break;
        }           
        }*/
            valid_candidate_idx[cur_depth][0]=Cand_vec[start_vertex][index_comb[start_vertex]];
            idx_count[cur_depth]=1;
            while (true)
                {
                 while (idx[cur_depth] < idx_count[cur_depth])
                {   

            
                    ui valid_idx = valid_candidate_idx[cur_depth][idx[cur_depth]];
                    VertexID u = order_matrix[cur_depth];
                    VertexID v = candidates[u][valid_idx];
                    auto end = std::chrono::high_resolution_clock::now();
                    ens = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
                    if (ens > TimeL)
                    { // 1000 1 sec
                        goto EXIT;
                    }
                    if (visited_vertices[v])
                    {
                     idx[cur_depth] += 1;
                    conflictT[u][reverse_embedding[v]]=1;
                    continue;
                    }
           /*if(nodeId[v]==0&&is_used[u]&&false){
                idx[cur_depth] += 1;
                conflictT[u][reverse_embedding[v]]=1;
                for (int b=max_depth-1;b>=0;b--)
                Match_BA[b]=true;
                continue;
            }*/ 
                if(nodeId[v]==1&& (!is_used[u])){
                    idx[cur_depth] += 1;
                    //conflictT[u][reverse_embedding[v]]=1;
                    //for (int b=max_depth-1;b>=0;b--)
                    //    Match_BA[b]=true;
                        continue;
                }
                if((nodeIDbad[cur_depth].find(v) != nodeIDbad[cur_depth].end())){
                idx[cur_depth] += 1;
                continue;
            }   
                embedding[u] = v;
                idx_embedding[u] = valid_idx;
                visited_vertices[v] = true;
                idx[cur_depth] += 1;
                reverse_embedding[v] = u;
            
                if (cur_depth == max_depth - 1)
                {   
                    UNPM = 0;
                    if(nodeId[embedding[order_matrix[cur_depth]]]==0){
                        UNPM++;
                        nodeId[embedding[order_matrix[cur_depth]]]=1;
                    }
                    embedding_cnt += 1;
                    visited_vertices[v] = false;
                    int ao = cur_depth;
                    ui vqo = order_matrix[ao];
                    tempPos[ao]=embedding[order_matrix[cur_depth]];
                    while (ao > 0)
                    {                
                        if(nodeId[embedding[vqo]]==0){
                        UNPM++;
                        nodeId[embedding[vqo]]=1;
                    }
                        tempPos[ao]=embedding[vqo];
                        if(ao>K_size2){
                            idx[ao] = idx_count[ao]; 
                            nodeIDbad[ao].clear();
                            reverse_embedding.erase(embedding[vqo]);           
                        }
                    
                    
                        Match_BA[ao]=true;
                        ao--;
                        vqo = order_matrix[ao];
                    }
                    Match_BA[ao]=true;
                    if(nodeId[embedding[vqo]]==0){
                    UNPM++;
                    nodeId[embedding[vqo]]=1;
                    //idx[ao] = idx_count[ao];
                }
                tempPos[ao]=embedding[vqo];
                if (UNPM > pq.top().first)
                            {                                        
                                int addPos=0;
                                if (heapcount >= ksize)
                                {
                                   auto [value, index] = pq.top();      
                                    pq.pop();
                                    pq.push({UNPM,index});
                                    addPos=index;
                                    
                                }else{
                                     pq.push({UNPM,heapcount});
                                     addPos=heapcount;            
                                }
                            for (int aa=0;aa<qsiz;aa++){
                                SolPos[addPos][aa]=tempPos[aa];
                                    }
                                heapcount++;
                            }
                
                    reverse_embedding.erase(embedding[u]);
                
                    if (embedding_cnt >= output_limit_num)
                    {
                        goto EXIT;
                    }
                }
                else
                {
                    call_count +=1;
                    cur_depth += 1;   
                    //Match_BA[cur_depth]=false;          
                    idx[cur_depth] = 0;
                    ui ut=order_matrix[cur_depth];
                    int Mem1=0;
                    if(is_used[ut]){                
                        generateValidCandidateIndexSS(cur_depth, idx_embedding, idx_count, valid_candidate_idx, edge_matrix, bn,
                                            bn_count, order_matrix, temp_buffer,Cand_vec[ut][index_comb[ut]],data_graph,candidates[ut][Cand_vec[ut][index_comb[ut]]],embedding);     

                //if(idx_count[cur_depth]==0)
                    }else{
                        generateValidCandidateIndex(cur_depth, idx_embedding, idx_count, valid_candidate_idx, edge_matrix, bn,
                                            bn_count, order_matrix, temp_buffer);

                }
                if (cur_depth!=max_depth-1)
                    for (int i=0;i<max_depth;i++){
                        conflictT[order_matrix[cur_depth]][i]=conflicts[order_matrix[cur_depth]][i];
                }
            
                }
            }

            cur_depth -= 1;

            if (cur_depth < 0)
                break;
            else
            {   
                int testdepth=cur_depth-1;
            while(testdepth>=0){
            if(conflictT[order[cur_depth]][order[testdepth]]==false)
                break;
            testdepth--;
            }
            //testdepth=depth of conflict
        if (testdepth==0)    
        if(Match_BA[cur_depth]==false&&cur_depth!=max_depth - 1){
            nodeIDbad[cur_depth].insert(embedding[order[cur_depth]]);
        }
        else if(testdepth>0&&conflictT[order[testdepth]][order[testdepth-1]]==false){
         if(Match_BA[cur_depth]==false&&cur_depth!=max_depth - 1){
            nodeIDbad[cur_depth].insert(embedding[order[cur_depth]]);
        }   
        }
        for (int mm=0;mm<vec[cur_depth].size();mm++){
            nodeIDbad[vec[cur_depth][mm]].clear();
        }
                VertexID u = order_matrix[cur_depth];
                reverse_embedding.erase(embedding[u]);

            if (cur_depth != 0&&Match_BA[cur_depth]==false&&cur_depth < max_depth-2)
            //if (cur_depth != 0&&Match_BA[cur_depth]==false)
            {
                //for (int k=cur_depth+1;k<max_depth-1;k++)
                //if(conflictT[order[cur_depth+1]][u]==false){
                bool Prun=true;
                ui test_depth=max_depth-1;
                while(test_depth>cur_depth){
                  if(conflictT[u][order_matrix[test_depth]]==false&&conflicts[u][order_matrix[test_depth]]==false)  
                        test_depth--;
                    else{
                        Prun=false;
                        break;
                    }
                }
                if (Prun==true)
                idx[cur_depth] = idx_count[cur_depth];
                //if(conflictT[u][order_matrix[cur_depth+1]]==false&&conflicts[u][order_matrix[cur_depth+1]]==false){
               //     idx[cur_depth] = idx_count[cur_depth];
               // }
            }
            visited_vertices[embedding[u]] = false;
           
        }
    }
    
    int mpe=K_size-1;
    while(mpe>=0){
        hasNext1=false;
        ui temp_u=combination[mpe];
        index_comb[temp_u]++;
        if(index_comb[temp_u]>=Cand_vec[temp_u].size()){
            index_comb[temp_u]=0;
            mpe--;
        }else{
            hasNext1=true;
            break;
        }

    }
    }   
        for (int i=0;i<max_depth;i++){
            index_comb[i]=0;
            }
            hasNext1=true;
        hasNext = nextCombination(indices, max_depth, K_size);
    }   
       hasNext1=true;
        K_size++;

    }


EXIT:
    releaseBuffer(max_depth, idx, idx_count, embedding, idx_embedding, temp_buffer, valid_candidate_idx,
                  visited_vertices,
                  bn, bn_count);

            std::priority_queue<std::pair<int, int>> maxHeap;

    int topkcounter = 0;
    int greedysum = 0;
    unordered_set<ui> EmbGreedy;
    while (!pq.empty())
    {   if(pq.top().first!=0)
        maxHeap.push(pq.top());
        pq.pop();
    }
    while (!maxHeap.empty())
    {   auto [value, index]= maxHeap.top();
        for (int aa=0;aa<qsiz;aa++){
            EmbGreedy.insert(SolPos[index][aa]);
        }
        maxHeap.pop();
        topkcounter++;
        switch (topkcounter)
        {

        case VALUE_10:
            s.topk[0] = EmbGreedy.size();
            break;
        case VALUE_50:
            s.topk[1] = EmbGreedy.size();
            break;
        case VALUE_100:
            s.topk[2] = EmbGreedy.size();
            break;
        case VALUE_500:
            s.topk[3] = EmbGreedy.size();
            break;
        case VALUE_1000:
            s.topk[4] = EmbGreedy.size();
            break;
        case VALUE_2500:
            s.topk[5] = EmbGreedy.size();
            break;    
        case VALUE_5000:
            s.topk[6] = EmbGreedy.size();
            break;
        case VALUE_10000:
        s.topk[7] = EmbGreedy.size();
            break;    
        default:
            break;
        }
    }
    for (int i = 0; i < 8; i++)
    {
        if (s.topk[i] == 0)
            s.topk[i] = EmbGreedy.size();
    }

        s.embedding_cnt = embedding_cnt;
        ui countSMU=0;
        for (int i=0;i<data_graph->getVerticesCount();i++){
            countSMU+=nodeId[i];
        }
        s.Can_embed = countSMU;
    // s.topk=greedysum;
        return s;
}

enumResult
EvaluateQuery::DSQL(const Graph *data_graph, const Graph *query_graph,ui *&nodeId, Edges ***edge_matrix, ui **candidates,
                    ui *candidates_count, ui *order, size_t output_limit_num, size_t &call_count, size_t &valid_vtx_cnt, int TimeL, int FairT,
                    const std::unordered_map<VertexID, std::pair<std::set<VertexID>, std::set<VertexID>>> &ordered_constraints)
    {

    auto start = std::chrono::high_resolution_clock::now();
    enumResult s;
    ui **bn;
    ui *bn_count;
    generateBN(query_graph, order, bn, bn_count);
    int UNPM = 0;
    ui *idx;
    ui *idx_count;
    ui *embedding;
    ui *idx_embedding;
    ui *temp_buffer;
    ui **valid_candidate_idx;
    int i=0;
    bool *visited_vertices;
        int count1=0;
    int count2=0;
    int K_size2=0;
        ui tempstore=0;
    ui nbr_cnt;
    ui nbr_cnt1;
    vector<int> indices(1);        
    //vector<int> indices1(max_depth); 
    allocateBuffer(data_graph, query_graph, candidates_count, idx, idx_count, embedding, idx_embedding,
                   temp_buffer, valid_candidate_idx, visited_vertices);
        int cc1=0;
    int cc2=0;
    size_t embedding_cnt = 0;
    int cur_depth = 0;
    ui max_depth = query_graph->getVerticesCount();
    VertexID start_vertex = order[0];
    ui order_matrix[max_depth];
    int index_comb[max_depth];
    bool is_used[max_depth];
    bool is_used_copy[max_depth];
    std::vector<std::vector<ui>> Cand_vec;
    Cand_vec.resize(max_depth);
    //unordered_set<ui> EmbSum;
    bool Match_BA[max_depth] = {false};
    bool conflicts[max_depth][max_depth] = {false};
    bool conflictT[max_depth][max_depth] = {false};
    std::vector<std::unordered_set<int>> nodeIDbad(max_depth);
        for (int i=0;i<max_depth;i++){
        for (int j=0;j<max_depth;j++)
        conflicts[i][j]=false;
    }
    for (int i=0;i<max_depth;i++){
    const VertexID *nbrs = query_graph->getVertexNeighbors(i, nbr_cnt);
            for (int j=0;j<nbr_cnt;j++){
                conflicts[i][nbrs[j]]=true;
                conflicts[nbrs[j]][i]=true;
            }
    }

    idx[cur_depth] = 0;
    idx_count[cur_depth] = candidates_count[start_vertex];

    for (ui i = 0; i < idx_count[cur_depth]; ++i)
    {
        valid_candidate_idx[cur_depth][i] = i;
    }    
    std::unordered_map<VertexID, VertexID> reverse_embedding;
    reverse_embedding.reserve(MAXIMUM_QUERY_GRAPH_SIZE * 2);
    bool hasNext = true;
    bool hasNext1 = true;
    ui K_size=1;
    TimeL = TimeL * 1000;
    vector<int> combination;
    std::vector<std::vector<int>> vec(max_depth);

    double ens = 0;
    ui rev_order[max_depth] = {0};
        for (int i=0;i<max_depth;i++){
            is_used[i]=false;
            is_used_copy[i]=false;
            index_comb[i]=0;
            }
            is_used[0]=order[0];
            order_matrix[0] = order[0];
            is_used_copy[order[0]]=true;
            count1=1;
            count2=0;
            //count1++;
            int countend=max_depth-1;
            if (query_graph->getVertexDegree(order_matrix[0])==1){
                const VertexID *nbrs = query_graph->getVertexNeighbors(order_matrix[count2], nbr_cnt);
                ui vertexN=nbrs[0];
                is_used_copy[vertexN]=true;
                order_matrix[count1]=vertexN;
                count1++;
                count2++;
                while (query_graph->getVertexDegree(vertexN)==1){
                    const VertexID *nbrs = query_graph->getVertexNeighbors(order_matrix[count2], nbr_cnt);    
                    vertexN=nbrs[0];
                    is_used_copy[vertexN]=true;
                    order_matrix[count1]=vertexN;
                    count1++;
                    count2++;
                }
            }
        while(count1<=countend){
            if (is_used_copy[order_matrix[count2]]){
            const VertexID *nbrs = query_graph->getVertexNeighbors(order_matrix[count2], nbr_cnt);
             
            for (int i=0;i<nbr_cnt;i++){
                if(is_used_copy[nbrs[i]]==false){
                    is_used_copy[nbrs[i]]=true;
                    nbr_cnt1 = query_graph->getVertexDegree(order_matrix[count2]);
                    if(nbr_cnt1>1){                    
                        order_matrix[count1]=nbrs[i];
                        count1++;
                    }
                    
                    else{
                        //order_matrix[count1]=nbrs[i];
                        //count1++;
                        order_matrix[countend]=nbrs[i];
                        countend--;
                    }
                    
                    
                }
            }}count2++;
        }

    generateBNLM(query_graph, order_matrix, bn, bn_count);
        for (int i=0;i<max_depth;i++){
        rev_order[order_matrix[i]]=i;
    }
        for (int i=max_depth-2;i>0;i--){
        int pp=-1;
        const VertexID *nbrs = query_graph->getVertexNeighbors(order_matrix[i], nbr_cnt);
        for (int j=0;j<nbr_cnt;j++){
            if (rev_order[order_matrix[nbrs[j]]]>pp&&rev_order[order_matrix[nbrs[j]]]<i)
            pp=rev_order[order_matrix[nbrs[j]]];
        }
        if (pp!=-1)
        vec[pp].push_back(i-1);
    }
    while (true)
    {
        while (idx[cur_depth] < idx_count[cur_depth])
        {   
            ui valid_idx = valid_candidate_idx[cur_depth][idx[cur_depth]];
            VertexID u = order_matrix[cur_depth];
            VertexID v = candidates[u][valid_idx];
            auto end = std::chrono::high_resolution_clock::now();
            ens = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
            if (ens > TimeL)
            { // 1000 1 sec
                goto EXIT;
            }

            if (visited_vertices[v])
            {
                idx[cur_depth] += 1;
    
                conflictT[u][reverse_embedding[v]]=true;
                conflictT[reverse_embedding[v]][u]=true;
                continue;
            }
            if(nodeId[v]==1){
                idx[cur_depth] += 1;
                
                //Match_BA[cur_depth]=true;
                for (int b=max_depth-1;b>=0;b--){
                Match_BA[b]=true;
                }
                //conflictT[u][reverse_embedding[v]]=true;
                continue;
            }
            if((nodeIDbad[cur_depth].find(v) != nodeIDbad[cur_depth].end())){
                idx[cur_depth] += 1;
                continue;
            }

            Match_BA[cur_depth]=false;
            embedding[u] = v;
            idx_embedding[u] = valid_idx;
            visited_vertices[v] = true;
            idx[cur_depth] += 1;

            reverse_embedding[v] = u;
            if (cur_depth == max_depth - 1)
            {
                UNPM = 0;
                if(nodeId[embedding[order_matrix[cur_depth]]]==0){
                    UNPM++;
                    nodeId[embedding[order_matrix[cur_depth]]]=1;
                }
                embedding_cnt += 1;
                visited_vertices[v] = false;
                int ao = cur_depth; //-1;
                ui vqo = order_matrix[ao];
                nodeId[v]=1;
                while (ao > 0)
                {   
                    Match_BA[ao]=true;
                    if(nodeId[embedding[vqo]]==0){
                    UNPM++;
                    nodeId[embedding[vqo]]=1;
                    nodeIDbad[ao].clear();
                }
                    idx[ao] = idx_count[ao];
                    ao--;
                    reverse_embedding.erase(embedding[vqo]);
                    visited_vertices[embedding[vqo]] = false;
                    vqo = order_matrix[ao];

                }
                Match_BA[ao]=true;
                nodeId[embedding[vqo]]=1;
                //reverse_embedding.erase(embedding[u]);
                //visited_vertices[embedding[vqo]] = false;
            }
            else
            {
                call_count += 1;
                cur_depth += 1;
                
                idx[cur_depth] = 0;
                
                generateValidCandidateIndex(cur_depth, idx_embedding, idx_count, valid_candidate_idx, edge_matrix, bn,
                                            bn_count, order_matrix, temp_buffer);
            //if (cur_depth!=0&&cur_depth!=max_depth-1)
            for (int i=0;i<max_depth;i++){
                conflictT[order_matrix[cur_depth]][order_matrix[i]]=conflicts[order_matrix[cur_depth]][order_matrix[i]];
            }
            //if(idx[cur_depth]!=idx_count[cur_depth]){
            //    Match_BA[cur_depth]=false;
            //}else{Match_BA[cur_depth]=true;
           // }


            //cout<<"hmm"<<endl;
           // ui valid_idx1 = valid_candidate_idx[cur_depth][idx[cur_depth]];
           // while(nodeId[candidates[order[cur_depth]][valid_idx1]]==1&&idx[cur_depth]<idx_count[cur_depth]){
          ///  idx[cur_depth]++;
           // valid_idx1 = valid_candidate_idx[cur_depth][idx[cur_depth]];    }
            }
        }


        cur_depth -= 1;

        if (cur_depth < 0)
            break;
        else
        {
            //get conflict level
            int testdepth=cur_depth-1;
            while(testdepth>=0){
            if(conflictT[order_matrix[cur_depth]][order_matrix[testdepth]]==false)
                break;
            testdepth--;
            }
            //testdepth=depth of conflict
        if (testdepth==0)    
        if(Match_BA[cur_depth]==false&&cur_depth!=max_depth - 1){
            nodeIDbad[cur_depth].insert(embedding[order_matrix[cur_depth]]);
        }
        else if(testdepth>0&&conflictT[order_matrix[testdepth]][order_matrix[testdepth-1]]==false){
         if(Match_BA[cur_depth]==false&&cur_depth!=max_depth - 1){
            nodeIDbad[cur_depth].insert(embedding[order_matrix[cur_depth]]);
        }   
        }
        for (int mm=0;mm<vec[cur_depth].size();mm++){
            nodeIDbad[vec[cur_depth][mm]].clear();
        }
            VertexID u = order_matrix[cur_depth];
            reverse_embedding.erase(embedding[u]);
            if (cur_depth != 0&&Match_BA[cur_depth]==false&&cur_depth < max_depth-2)
            //if (cur_depth != 0&&Match_BA[cur_depth]==false)
            {
                //for (int k=cur_depth+1;k<max_depth-1;k++)
                //if(conflictT[order[cur_depth+1]][u]==false){
                bool Prun=true;
                ui test_depth=max_depth-1;
                while(test_depth>cur_depth){
                  if(conflictT[u][order_matrix[test_depth]]==false&&conflicts[u][order_matrix[test_depth]]==false)  
                        test_depth--;
                    else{
                        Prun=false;
                        break;
                    }
                }
                if (Prun==true)
                idx[cur_depth] = idx_count[cur_depth];
                //if(conflictT[u][order_matrix[cur_depth+1]]==false&&conflicts[u][order_matrix[cur_depth+1]]==false){
               //     idx[cur_depth] = idx_count[cur_depth];
               // }
            }
            visited_vertices[embedding[u]] = false;

        }
    }
   //goto EXIT;
    hasNext = true;
     K_size=1;
        //for (int i=0;i<max_depth;i++)
        //cout<<order[i];
    
    //std::vector<std::vector<ui>> Cand_vec;
       // 
    while (K_size<max_depth){
        indices.resize(K_size);
        for (ui i = 0; i < K_size; ++i) {
        indices[i] = i;
        }
        hasNext = true;
        
    for (int ik = 0; ik < max_depth; ik++) {
        
        Cand_vec[ik].clear();
        for (int i=0;i<candidates_count[ik];i++)
          if(nodeId[candidates[ik][i]]==1){
            //Cand_vec[ik].push_back(candidates[ik][i]);
            Cand_vec[ik].push_back(i);
                }
        }

        while (hasNext==true) {
            cc1++;
            
            count1=0;
            count2=0;
            K_size2=0;
            tempstore=0;
            nbr_cnt=0;
            for (int i=0;i<max_depth;i++){
            is_used[i]=false;
            is_used_copy[i]=false;
            index_comb[i]=0;
            }
            combination.clear();
            
            for (ui i = 0; i < K_size; ++i) {
                combination.push_back(order[indices[i]]);
                is_used[order[indices[i]]] = true;
            }

            order_matrix[0] = combination[0];
            is_used_copy[combination[0]]=true;
            count1=1;
            count2=0;
            //count1++;
            int countend=max_depth-1;
            if (query_graph->getVertexDegree(order_matrix[0])==1){
                
                const VertexID *nbrs = query_graph->getVertexNeighbors(order_matrix[count2], nbr_cnt);
                ui vertexN=nbrs[0];
                is_used_copy[vertexN]=true;
                order_matrix[count1]=vertexN;
                count1++;
                count2++;
                while (query_graph->getVertexDegree(vertexN)==1){
                    const VertexID *nbrs = query_graph->getVertexNeighbors(order_matrix[count2], nbr_cnt);    
                    vertexN=nbrs[0];
                    is_used_copy[vertexN]=true;
                    order_matrix[count1]=vertexN;
                    count1++;
                    count2++;
                }
            }
        while(count1<max_depth){
            if (is_used_copy[order_matrix[count2]]){
            const VertexID *nbrs = query_graph->getVertexNeighbors(order_matrix[count2], nbr_cnt);
            for (int i=0;i<nbr_cnt;i++){
                if(is_used_copy[nbrs[i]]==false){
                    is_used_copy[nbrs[i]]=true;
                    nbr_cnt1 = query_graph->getVertexDegree(order_matrix[count2]);
                    if(nbr_cnt1>1){                        
                        order_matrix[count1]=nbrs[i];
                        count1++;
                    }
                    
                    else{
                        order_matrix[countend]=nbrs[i];
                        countend--;
                    }
                    
                    
                }
            }}count2++;
        }
        for (int i=0;i<max_depth;i++){
        rev_order[order_matrix[i]]=i;
        vec[i].clear();
    }
    for (int i=max_depth-2;i>0;i--){
        int pp=-1;
        const VertexID *nbrs = query_graph->getVertexNeighbors(order_matrix[i], nbr_cnt);
        for (int j=0;j<nbr_cnt;j++){
            if (rev_order[order_matrix[nbrs[j]]]>pp&&rev_order[order_matrix[nbrs[j]]]<i)
            pp=rev_order[order_matrix[nbrs[j]]];
        }
        if (pp!=-1)
        vec[pp].push_back(i-1);
    }

        count1=max_depth;
        count2=0;
        i=1;


        K_size2=0;
        i=0;
        while(is_used[order_matrix[i]]==true){
            K_size2++;
            i++;
        }
        
        i=0;
        
       // generateBN1(query_graph, order_matrix, bn, bn_count);
       generateBNLM(query_graph, order_matrix, bn, bn_count);
        allocateBufferLM(data_graph, query_graph, candidates_count, idx, idx_count, embedding, idx_embedding,
                   temp_buffer, valid_candidate_idx, visited_vertices);

        //std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> vec_failing_set(max_depth);
        while (hasNext1==true){
            cc2++;
            for (int aa1=0;aa1<max_depth;aa1++)
            nodeIDbad[aa1].clear();
            start_vertex = order_matrix[0];
            cur_depth=0;
            idx[cur_depth] = 0;
            idx_count[cur_depth] = candidates_count[start_vertex];
            Match_BA[0] = false;
        
        /*
        for (ui i = 0; i < idx_count[cur_depth]; ++i)
        { 
        if(candidates[start_vertex][i]==Cand_vec[start_vertex][index_comb[start_vertex]]) {
            valid_candidate_idx[cur_depth][0] = i;            
            break;
        }           
        }*/
            
            valid_candidate_idx[cur_depth][0]=Cand_vec[start_vertex][index_comb[start_vertex]];
            idx_count[cur_depth]=1;
            while (true)
                {
                 while (idx[cur_depth] < idx_count[cur_depth])
                {   

            
                    ui valid_idx = valid_candidate_idx[cur_depth][idx[cur_depth]];
                    VertexID u = order_matrix[cur_depth];
                    VertexID v = candidates[u][valid_idx];
                    auto end = std::chrono::high_resolution_clock::now();
                    ens = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
                    if (ens > TimeL)
                    { // 1000 1 sec
                        goto EXIT;
                    }
                    if (visited_vertices[v])
                    {
                     idx[cur_depth] += 1;
                    conflictT[u][reverse_embedding[v]]=1;
                    continue;
                    }
           /*if(nodeId[v]==0&&is_used[u]&&false){
                idx[cur_depth] += 1;
                conflictT[u][reverse_embedding[v]]=1;
                for (int b=max_depth-1;b>=0;b--)
                Match_BA[b]=true;
                continue;
            }*/ 
                if(nodeId[v]==1&& (!is_used[u])){
                    idx[cur_depth] += 1;
                    //conflictT[u][reverse_embedding[v]]=1;
                    //for (int b=max_depth-1;b>=0;b--)
                    //    Match_BA[b]=true;
                        continue;
                }
                if((nodeIDbad[cur_depth].find(v) != nodeIDbad[cur_depth].end())){
                idx[cur_depth] += 1;
                continue;
            }
                embedding[u] = v;
                idx_embedding[u] = valid_idx;
                visited_vertices[v] = true;
                idx[cur_depth] += 1;
                reverse_embedding[v] = u;
            
                if (cur_depth == max_depth - 1)
                {   
                    UNPM = 0;
                    if(nodeId[embedding[order_matrix[cur_depth]]]==0){
                        UNPM++;
                        nodeId[embedding[order_matrix[cur_depth]]]=1;
                    }
                    embedding_cnt += 1;
                    visited_vertices[v] = false;
                    int ao = cur_depth;
                    ui vqo = order_matrix[ao];

                    while (ao > 0)
                    {                
                        if(nodeId[embedding[vqo]]==0){
                        UNPM++;
                        nodeId[embedding[vqo]]=1;
                    }
                        if(ao>K_size2){
                            idx[ao] = idx_count[ao]; 
                            nodeIDbad[ao].clear();
                            reverse_embedding.erase(embedding[vqo]);           
                        }
                    
                    
                        Match_BA[ao]=true;
                        ao--;
                        vqo = order_matrix[ao];
                    }
                    Match_BA[ao]=true;
                    nodeId[embedding[vqo]]=1;
                
                    reverse_embedding.erase(embedding[u]);
                
                    if (embedding_cnt >= output_limit_num)
                    {
                        goto EXIT;
                    }
                }
                else
                {
                    call_count +=1;
                    cur_depth += 1;   
                    //Match_BA[cur_depth]=false;          
                    idx[cur_depth] = 0;
                    ui ut=order_matrix[cur_depth];
                    int Mem1=0;
                    if(is_used[ut]){                
                        generateValidCandidateIndexSS(cur_depth, idx_embedding, idx_count, valid_candidate_idx, edge_matrix, bn,
                                            bn_count, order_matrix, temp_buffer,Cand_vec[ut][index_comb[ut]],data_graph,candidates[ut][Cand_vec[ut][index_comb[ut]]],embedding);     

                //if(idx_count[cur_depth]==0)
                    }else{
                        generateValidCandidateIndex(cur_depth, idx_embedding, idx_count, valid_candidate_idx, edge_matrix, bn,
                                            bn_count, order_matrix, temp_buffer);

                }
                if (cur_depth!=max_depth-1)
                    for (int i=0;i<max_depth;i++){
                        conflictT[order_matrix[cur_depth]][i]=conflicts[order_matrix[cur_depth]][i];
                }
            
                }
            }

            cur_depth -= 1;

            if (cur_depth < 0)
                break;
            else
            {   
                int testdepth=cur_depth-1;
            while(testdepth>=0){
            if(conflictT[order[cur_depth]][order[testdepth]]==false)
                break;
            testdepth--;
            }
            //testdepth=depth of conflict
        if (testdepth==0)    
        if(Match_BA[cur_depth]==false&&cur_depth!=max_depth - 1){
            nodeIDbad[cur_depth].insert(embedding[order[cur_depth]]);
        }
        else if(testdepth>0&&conflictT[order[testdepth]][order[testdepth-1]]==false){
         if(Match_BA[cur_depth]==false&&cur_depth!=max_depth - 1){
            nodeIDbad[cur_depth].insert(embedding[order[cur_depth]]);
        }   
        }
        for (int mm=0;mm<vec[cur_depth].size();mm++){
            nodeIDbad[vec[cur_depth][mm]].clear();
        }
                VertexID u = order_matrix[cur_depth];
                reverse_embedding.erase(embedding[u]);

if (cur_depth != 0&&Match_BA[cur_depth]==false&&cur_depth < max_depth-2)
            //if (cur_depth != 0&&Match_BA[cur_depth]==false)
            {
                //for (int k=cur_depth+1;k<max_depth-1;k++)
                //if(conflictT[order[cur_depth+1]][u]==false){
                bool Prun=true;
                ui test_depth=max_depth-1;
                while(test_depth>cur_depth){
                  if(conflictT[u][order_matrix[test_depth]]==false&&conflicts[u][order_matrix[test_depth]]==false)  
                        test_depth--;
                    else{
                        Prun=false;
                        break;
                    }
                }
                if (Prun==true)
                idx[cur_depth] = idx_count[cur_depth];
                //if(conflictT[u][order_matrix[cur_depth+1]]==false&&conflicts[u][order_matrix[cur_depth+1]]==false){
               //     idx[cur_depth] = idx_count[cur_depth];
               // }
            }
            visited_vertices[embedding[u]] = false;
           
        }
    }
    
    int mpe=K_size-1;
    while(mpe>=0){
        hasNext1=false;
        ui temp_u=combination[mpe];
        index_comb[temp_u]++;
        if(index_comb[temp_u]>=Cand_vec[temp_u].size()){
            index_comb[temp_u]=0;
            mpe--;
        }else{
            hasNext1=true;
            break;
        }

    }
    }   
        for (int i=0;i<max_depth;i++){
            index_comb[i]=0;
            }
            hasNext1=true;
        hasNext = nextCombination(indices, max_depth, K_size);
    }
       hasNext1=true;
        K_size++;

    }


EXIT:
    releaseBuffer(max_depth, idx, idx_count, embedding, idx_embedding, temp_buffer, valid_candidate_idx,
                  visited_vertices,
                  bn, bn_count);

        s.embedding_cnt = embedding_cnt;
        ui countSMU=0;
        for (int i=0;i<data_graph->getVerticesCount();i++){
            countSMU+=nodeId[i];
        }
        s.Can_embed = countSMU;
    // s.topk=greedysum;
        return s;
}

void EvaluateQuery::generateValidCandidateIndexSS(ui depth, ui *idx_embedding, ui *idx_count, ui **valid_candidate_index,
                                                Edges ***edge_matrix, ui **bn, ui *bn_cnt, ui *order,
                                                ui *&temp_buffer,ui CN,const Graph *data_graph,ui dataV,ui *embedding)
{   

    VertexID u = order[depth];
    VertexID previous_bn = bn[depth][0];
    ui previous_index_id = idx_embedding[previous_bn];
    ui valid_candidates_count = 0;
    Edges &previous_edge = *edge_matrix[previous_bn][u];
    valid_candidates_count=1;
    std::fill(valid_candidate_index[depth], valid_candidate_index[depth] + valid_candidates_count, CN);
    ui temp_count;
    bool okay=false;
    for (ui i = 0; i < bn_cnt[depth]; ++i)
    {     
        okay=data_graph->checkEdgeExistence(dataV,embedding[bn[depth][i]]);
        if(okay==false){
            valid_candidates_count=0;
            break;
        }
    }
    idx_count[depth] = valid_candidates_count;
}