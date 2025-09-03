#include "leapfrogtriejoin.h"
#include "../computesetintersection.h"
#include <unordered_set>
#include <cstring>
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

enumResult leapfrogtriejoin::execute(size_t& call_count) {
    int MemSize=0;
        //if (getValue1() > MemSize)
      //  MemSize = getValue1();

    //cout << "MemSize" << MemSize / 1000 << endl;
    
    //nodeId = new ui[catalog_->data_graph_->getVerticesCount()];
    
    auto start = std::chrono::high_resolution_clock::now();
    uint32_t start_depth = input_->get_tuple_length();
    uint32_t max_depth = num_vertex_;
    uint32_t core_vertex_depth = num_core_vertex_;
    double ens = 0;
    uint32_t last_vertex = vertex_ordering_[max_depth - 1];
    uint32_t* last_vertex_candidate_sets = catalog_->candidate_sets_[last_vertex];
    intermediate_result_count_[0] = input_->get_size();

    #ifdef TOPKGREEDY
    //pq.push(0);
    
    //int **SolPos = new int*[ksize+1];
    SolPos = new int*[ksize+1];
    for (int i = 0; i < (ksize+1); ++i) {
        SolPos[i] = new int[max_depth];
    }
    //int *tempPos= new int [max_depth];
    tempPos= new int [max_depth];
    pq.push({0,0});
    #endif
    for (uint64_t i = 0; i < input_->get_size(); ++i) {

        memcpy(embedding_, input_->get_tuple(i), sizeof(uint32_t) * start_depth);

        set_idx(start_depth);

        for (uint32_t j = 0; j < start_depth; ++j) {
            visited_[embedding_[j]] = true;
        }

#ifdef FAILING_SET_PRUNING
        for (uint32_t j = 0; j < start_depth; ++j) {
            reverse_embedding_[embedding_[j]] = vertex_ordering_[j];
        }
#endif
                
                
                auto end = std::chrono::high_resolution_clock::now();
                    ens = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
                    if (ens > timeL_)
                    { // 1000 1 sec
                  goto EXIT;
                        //return count_;
                    }
        if (start_depth >= num_core_vertex_) {
            uint32_t temp_u = vertex_ordering_[start_depth];
            uint32_t temp_bn = get_backward_neighbors(start_depth)[0];
            uint32_t temp_bns_depth = get_backward_neighbors_depth(start_depth)[0];
            uint32_t temp_key = embedding_[temp_bns_depth];
            local_candidates_[start_depth] = catalog_->get_non_core_relation_children(temp_bn, temp_u, temp_key, num_local_candidates_[start_depth]);
            for (uint32_t j = 0; j < num_local_candidates_[start_depth]; ++j) {
                si_buffer_[j] = local_candidates_[start_depth][j];
            }

#ifdef FAILING_SET_PRUNING
            if (num_local_candidates_[start_depth] == 0) {
                vec_failing_set_[start_depth - 1] = ancestors_[vertex_ordering_[start_depth]];
            } else {
                vec_failing_set_[start_depth - 1].reset();
            }
#endif
            
            //end = std::chrono::high_resolution_clock::now();
            //ens = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
            
            bool exit = enumerate_non_core_results(start_depth, max_depth,ens);
            if (exit) {
                  goto EXIT;
            }
               
        }
        else {
            uint32_t cur_depth = start_depth;
            idx_[cur_depth] = 0;
            num_local_candidates_[cur_depth] = compute_local_candidates(cur_depth);    
#ifdef FAILING_SET_PRUNING
            if (num_local_candidates_[cur_depth] == 0) {
                vec_failing_set_[cur_depth - 1] = ancestors_[vertex_ordering_[cur_depth]];
            } else {
                vec_failing_set_[cur_depth - 1].reset();
            }
#endif

            intermediate_result_count_[cur_depth] += num_local_candidates_[cur_depth];

            while (true) {
                while (idx_[cur_depth] < num_local_candidates_[cur_depth]) {
                    
                    //first here
                    auto end = std::chrono::high_resolution_clock::now();
                    ens = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
                    if (ens > timeL_)
                    { // 1000 1 sec
                     goto EXIT;
                    }

                    uint32_t u = vertex_ordering_[cur_depth];
                    uint32_t encoded_id = local_candidates_[cur_depth][idx_[cur_depth]++];
                    uint32_t v = catalog_->candidate_sets_[u][encoded_id];
                    idx_embedding_[cur_depth] = encoded_id;

                    if (visited_[v]) {
#ifdef FAILING_SET_PRUNING
                        vec_failing_set_[cur_depth] = ancestors_[u];
                        vec_failing_set_[cur_depth] |= ancestors_[reverse_embedding_[v]];
                        vec_failing_set_[cur_depth - 1] |= vec_failing_set_[cur_depth];
#endif
                        continue;
                    }
                    visited_[v] = true;

#ifdef FAILING_SET_PRUNING
                    reverse_embedding_[v] = u;
#endif

                    embedding_[cur_depth] = v;


                    num_local_candidates_[cur_depth + 1] = compute_local_candidates(cur_depth + 1);
#ifdef FAILING_SET_PRUNING
                    if (num_local_candidates_[cur_depth + 1] == 0) {
                        vec_failing_set_[cur_depth] = ancestors_[vertex_ordering_[cur_depth + 1]];
                    } else {
                        vec_failing_set_[cur_depth].reset();
                    }
#endif

                    intermediate_result_count_[cur_depth + 1] += num_local_candidates_[cur_depth + 1];

                    if (cur_depth == max_depth - 2) {
                        bool doit=true;
                        uint32_t *temp_buffer = local_candidates_[cur_depth + 1];
                        uint32_t loop_count = num_local_candidates_[cur_depth + 1];
                        for (uint32_t j = 0; j < loop_count; ++j) {
                            uint32_t temp_v = last_vertex_candidate_sets[temp_buffer[j]];
                    
                    end = std::chrono::high_resolution_clock::now();
                    ens = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
                    if (ens > timeL_)
                    { 
                   goto EXIT;
                        //return count_;
                    }
                            if (visited_[temp_v]) {
#ifdef FAILING_SET_PRUNING
                                vec_failing_set_[cur_depth + 1] = ancestors_[last_vertex];
                                vec_failing_set_[cur_depth + 1] |= ancestors_[reverse_embedding_[temp_v]];
                                vec_failing_set_[cur_depth] |= vec_failing_set_[cur_depth + 1];
#endif
                                continue;
                            }

#ifdef FAILING_SET_PRUNING
                            vec_failing_set_[cur_depth + 1].set();
                            vec_failing_set_[cur_depth] |= vec_failing_set_[cur_depth + 1];
#endif                      
                            
                           
                            UNPM=0; 
                            if(nodeId[temp_v]==0){
                                nodeId[temp_v]=1;
                                UNPM++;
                            } 
                    #ifdef TOPKGREEDY
                    tempPos[cur_depth+1]=temp_v;
                    #endif  
                            int ao=cur_depth;
                        while (ao >= 0)
                        {    if(nodeId[embedding_[ao]]==0){
                                nodeId[embedding_[ao]]=1;
                                UNPM++;
                            }
                    #ifdef TOPKGREEDY
                    tempPos[ao]=embedding_[ao];
                    #endif
                                ao--;
                        }//doit=false;
                    
                    last_level_candidates_[j] = temp_v;
                    count_ += 1;
                
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
                            for (int aa=0;aa<max_depth;aa++){
                                SolPos[addPos][aa]=tempPos[aa];
                                    }
                                heapcount++;
                            }
#endif
                
                            
                        }
#ifdef FAILING_SET_PRUNING
                        if (!vec_failing_set_[cur_depth].test(vertex_ordering_[cur_depth])) {
                            vec_failing_set_[cur_depth - 1] = vec_failing_set_[cur_depth];
                            idx_[cur_depth] = num_local_candidates_[cur_depth];
                        } else {
                            vec_failing_set_[cur_depth - 1] |= vec_failing_set_[cur_depth];
                        }

                        reverse_embedding_.erase(v);
#endif

                        visited_[v] = false;


                        if (count_ >= output_count_limit_) {
                             goto EXIT;                       
                    //return count_;
                        }
                    }
                    else if (cur_depth == core_vertex_depth - 2) {
                        call_count += 1;
                        uint32_t next_depth = cur_depth + 1;
                        uint32_t next_vertex = vertex_ordering_[next_depth];

                        for (uint32_t j = 0; j < num_local_candidates_[next_depth]; ++j) {
                            uint32_t temp_idx = local_candidates_[next_depth][j];
                            
                            si_buffer_[j] = catalog_->candidate_sets_[next_vertex][temp_idx];
                        }
                        auto end = std::chrono::high_resolution_clock::now();
                        ens = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
                        bool exit = enumerate_non_core_results(next_depth, max_depth,ens);
                        if (exit) {   
                            goto EXIT;
                       // return count_;
                        } 
#ifdef FAILING_SET_PRUNING
                        if (!vec_failing_set_[cur_depth].test(vertex_ordering_[cur_depth])) {
                            vec_failing_set_[cur_depth - 1] = vec_failing_set_[cur_depth];
                            idx_[cur_depth] = num_local_candidates_[cur_depth];
                        } else {
                            vec_failing_set_[cur_depth - 1] |= vec_failing_set_[cur_depth];
                        }

                        reverse_embedding_.erase(v);
#endif

                        visited_[v] = false;
                    }
                    else {
                        call_count += 1;
                        cur_depth += 1;
                        idx_[cur_depth] = 0;
                    }
                }

                cur_depth -= 1;


                if (cur_depth < start_depth) {
                    break;
                } else {

#ifdef FAILING_SET_PRUNING
                    if (!vec_failing_set_[cur_depth].test(vertex_ordering_[cur_depth])) {
                        vec_failing_set_[cur_depth - 1] = vec_failing_set_[cur_depth];
                        idx_[cur_depth] = num_local_candidates_[cur_depth];
                    } else {
                        vec_failing_set_[cur_depth - 1] |= vec_failing_set_[cur_depth];
                    }

                    reverse_embedding_.erase(embedding_[cur_depth]);
#endif

                    visited_[embedding_[cur_depth]] = false;
                }
            }
        }
        for (uint32_t j = 0; j < start_depth; ++j) {
            visited_[embedding_[j]] = false;
        }

#ifdef FAILING_SET_PRUNING
        reverse_embedding_.clear();
#endif
    }

    EXIT:
        //if (getValue1() > MemSize)
       // MemSize = getValue1();
    enumOutput.embedding_cnt=count_;
    enumOutput.Can_embed=0;
    ui countAA=0;
    for (int i=0;i<catalog_->data_graph_->getVerticesCount();i++){
        countAA=countAA+nodeId[i];
    }
    enumOutput.Can_embed=countAA;
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
        for (int aa=0;aa<max_depth;aa++){
            EmbGreedy.insert(SolPos[index][aa]);
        }
        maxHeap.pop();
        topkcounter++;
        switch (topkcounter)
        {
        case VALUE_5:
            enumOutput.topk[0] = EmbGreedy.size();
            break;
        case VALUE_10:
            enumOutput.topk[1] = EmbGreedy.size();
            break;
        case VALUE_50:
            enumOutput.topk[2] = EmbGreedy.size();
            break;
        case VALUE_100:
            enumOutput.topk[3] = EmbGreedy.size();
            break;
        case VALUE_250:
            enumOutput.topk[4] = EmbGreedy.size();
            break;
        case VALUE_500:
            enumOutput.topk[5] = EmbGreedy.size();
            break;
        case VALUE_750:
            enumOutput.topk[6] = EmbGreedy.size();
            break;
        case VALUE_1000:
            enumOutput.topk[7] = EmbGreedy.size();
            break;
        default:
            break;
        }
    }
    for (int i = 0; i < 8; i++)
    {
        if (enumOutput.topk[i] == 0)
            enumOutput.topk[i] = EmbGreedy.size();
    }

#endif
    
    return enumOutput;
   // return count_; //
    
}

bool leapfrogtriejoin::enumerate_non_core_results(uint32_t start_depth, uint32_t max_depth,double TimeL) {
    auto start = std::chrono::high_resolution_clock::now();
    double ens = 0;
    for (idx_[start_depth] = 0; idx_[start_depth] < num_local_candidates_[start_depth]; ++idx_[start_depth]) {
        uint32_t v = si_buffer_[idx_[start_depth]];
       
        auto end = std::chrono::high_resolution_clock::now();
                    ens = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
                    if ((ens+TimeL) > timeL_)
                    { // 1000 1 sec
                        return true;
                    }
        if (visited_[v]) {
#ifdef FAILING_SET_PRUNING
            vec_failing_set_[start_depth] = ancestors_[vertex_ordering_[start_depth]];
            vec_failing_set_[start_depth] |= ancestors_[reverse_embedding_[v]];
            vec_failing_set_[start_depth - 1] |= vec_failing_set_[start_depth];
#endif
            continue;
        }

#ifdef FAILING_SET_PRUNING
        reverse_embedding_[v] = vertex_ordering_[start_depth];
#endif
        visited_[v] = true;

        embedding_[start_depth] = v;

        uint32_t cur_depth = start_depth + 1;
        uint32_t u = vertex_ordering_[cur_depth];
        uint32_t bn = get_backward_neighbors(cur_depth)[0];
        uint32_t bns_depth = get_backward_neighbors_depth(cur_depth)[0];
        uint32_t key = embedding_[bns_depth];

        local_candidates_[cur_depth] = catalog_->get_non_core_relation_children(bn, u, key, num_local_candidates_[cur_depth]);
        intermediate_result_count_[cur_depth] += num_local_candidates_[cur_depth];
#ifdef FAILING_SET_PRUNING
        if (num_local_candidates_[cur_depth] == 0) {
            vec_failing_set_[cur_depth - 1] = ancestors_[vertex_ordering_[cur_depth]];
        } else {
            vec_failing_set_[cur_depth - 1].reset();
        }
#endif
        if (cur_depth == max_depth - 1) {
            bool doit=true;
            uint32_t temp_count = num_local_candidates_[cur_depth];
            uint32_t* temp_buffer = local_candidates_[cur_depth];
            for (uint32_t j = 0; j < temp_count; ++j) {
                if (visited_[temp_buffer[j]]) {
#ifdef FAILING_SET_PRUNING
                    vec_failing_set_[cur_depth] = ancestors_[vertex_ordering_[cur_depth]];
                    vec_failing_set_[cur_depth] |= ancestors_[reverse_embedding_[temp_buffer[j]]];
                    vec_failing_set_[cur_depth - 1] |= vec_failing_set_[cur_depth];
#endif
                    continue;
                }

#ifdef FAILING_SET_PRUNING
                vec_failing_set_[cur_depth].set();
                vec_failing_set_[cur_depth - 1] |= vec_failing_set_[cur_depth];
#endif              
                    UNPM=0; 
                    if(nodeId[temp_buffer[j]]==0){
                          UNPM++;  
                          nodeId[temp_buffer[j]]=1;
                        }
                        #ifdef TOPKGREEDY
                    tempPos[cur_depth]=temp_buffer[j];
                    #endif
                    //UNPM+=EmbSum.insert(temp_buffer[j]).second;
                            int ao=cur_depth-1;
                        while (ao >=0 && doit == true)
                        {  
                            if(nodeId[embedding_[ao]]==0){
                          UNPM++;  
                          nodeId[embedding_[ao]]=1;
                        }                    
                        #ifdef TOPKGREEDY
                        tempPos[ao]=embedding_[ao];
                    #endif
                            //UNPM+=EmbSum.insert(embedding_[ao]).second;
                                ao--;
                }//doit=false;
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
                            for (int aa=0;aa<max_depth;aa++){
                                SolPos[addPos][aa]=tempPos[aa];
                                    }
                                heapcount++;
                            }
                #endif
                last_level_candidates_[j] = temp_buffer[j];
                count_ += 1;


            }


#ifdef FAILING_SET_PRUNING
            if (!vec_failing_set_[cur_depth - 1].test(vertex_ordering_[cur_depth - 1])) {
                vec_failing_set_[cur_depth - 2] = vec_failing_set_[cur_depth - 1];
                idx_[cur_depth - 1] = num_local_candidates_[cur_depth - 1];
            } else {
                vec_failing_set_[cur_depth - 2] |= vec_failing_set_[cur_depth - 1];
            }

            reverse_embedding_.erase(v);
#endif

            visited_[v] = false;
            if (count_ >= output_count_limit_) {
                return true;
            }
        }
        else {
            idx_[cur_depth] = 0;
            while (true) {
                while (idx_[cur_depth] < num_local_candidates_[cur_depth]) {

                    auto end = std::chrono::high_resolution_clock::now();
                    ens = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
                    if ((ens+TimeL) > timeL_)
                    { // 1000 1 sec
                        return true;
                    }
                    v = local_candidates_[cur_depth][idx_[cur_depth]++];

                    if (visited_[v]) {
#ifdef FAILING_SET_PRUNING
                        vec_failing_set_[cur_depth] = ancestors_[vertex_ordering_[cur_depth]];
                        vec_failing_set_[cur_depth] |= ancestors_[reverse_embedding_[v]];
                        vec_failing_set_[cur_depth - 1] |= vec_failing_set_[cur_depth];
#endif
                        continue;
                    }

#ifdef FAILING_SET_PRUNING
                    reverse_embedding_[v] = vertex_ordering_[cur_depth];
#endif
                    visited_[v] = true;
                    embedding_[cur_depth] = v;
                    uint32_t next_depth = cur_depth + 1;
                    u = vertex_ordering_[next_depth];

                    bn = get_backward_neighbors(next_depth)[0];
                    bns_depth = get_backward_neighbors_depth(next_depth)[0];

                    key = embedding_[bns_depth];
                    local_candidates_[next_depth] = catalog_->get_non_core_relation_children(bn, u, key,
                                                                                             num_local_candidates_[next_depth]);
                    intermediate_result_count_[next_depth] += num_local_candidates_[next_depth];
#ifdef FAILING_SET_PRUNING
                    if (num_local_candidates_[next_depth] == 0) {
                        vec_failing_set_[next_depth - 1] = ancestors_[vertex_ordering_[next_depth]];
                    } else {
                        vec_failing_set_[next_depth - 1].reset();
                    }
#endif

                    if (cur_depth == max_depth - 2) {
                        bool doit=true;
                        uint32_t temp_count = num_local_candidates_[next_depth];
                        uint32_t *temp_buffer = local_candidates_[next_depth];
                        for (uint32_t j = 0; j < temp_count; ++j) {
                            if (visited_[temp_buffer[j]]) {
#ifdef FAILING_SET_PRUNING
                                vec_failing_set_[next_depth] = ancestors_[vertex_ordering_[next_depth]];
                                vec_failing_set_[next_depth] |= ancestors_[reverse_embedding_[temp_buffer[j]]];
                                vec_failing_set_[next_depth - 1] |= vec_failing_set_[next_depth];
#endif
                                continue;
                            }

#ifdef FAILING_SET_PRUNING
                            vec_failing_set_[next_depth].set();
                            vec_failing_set_[next_depth - 1] |= vec_failing_set_[next_depth];
#endif              
            
                        UNPM=0;  
                        if(nodeId[temp_buffer[j]]==0){
                          UNPM++;  
                          nodeId[temp_buffer[j]]=1;
                        }  
                        
                            #ifdef TOPKGREEDY
                            tempPos[next_depth]=temp_buffer[j];
                            #endif
                        
                            int ao=cur_depth;
                        while (ao >= 0 )
                        {     
                            if(nodeId[embedding_[ao]]==0){
                          UNPM++;  
                          nodeId[embedding_[ao]]=1;
                        }                     
                            #ifdef TOPKGREEDY
                            tempPos[ao]=embedding_[ao];
                            #endif
                                ao--;
                            }//doit=false;
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
                            for (int aa=0;aa<max_depth;aa++){
                                 
                                SolPos[addPos][aa]=tempPos[aa];
                                    }
                                heapcount++;
                            }
                #endif
                
                            last_level_candidates_[j] = temp_buffer[j];
                            count_ += 1;
                        }

#ifdef FAILING_SET_PRUNING
                        if (!vec_failing_set_[cur_depth].test(vertex_ordering_[cur_depth])) {
                            vec_failing_set_[cur_depth - 1] = vec_failing_set_[cur_depth];
                            idx_[cur_depth] = num_local_candidates_[cur_depth];
                        } else {
                            vec_failing_set_[cur_depth - 1] |= vec_failing_set_[cur_depth];
                        }

                        reverse_embedding_.erase(v);
#endif

                        visited_[v] = false;
                        if (count_ >= output_count_limit_) {
                            return true;
                        }
                    } else {
                        cur_depth += 1;
                        idx_[cur_depth] = 0;
                    }
                }

                cur_depth -= 1;

                if (cur_depth <= start_depth) {
                    break;
                } else {
                    visited_[embedding_[cur_depth]] = false;

#ifdef FAILING_SET_PRUNING
                    reverse_embedding_.erase(embedding_[cur_depth]);
                    if (!vec_failing_set_[cur_depth].test(vertex_ordering_[cur_depth])) {
                        vec_failing_set_[cur_depth - 1] = vec_failing_set_[cur_depth];
                        idx_[cur_depth] = num_local_candidates_[cur_depth];
                    } else {
                        vec_failing_set_[cur_depth - 1] |= vec_failing_set_[cur_depth];
                    }
#endif
                }

            }
            visited_[embedding_[start_depth]] = false;

#ifdef FAILING_SET_PRUNING
            reverse_embedding_.erase(embedding_[start_depth]);
            if (!vec_failing_set_[start_depth].test(vertex_ordering_[start_depth])) {
                vec_failing_set_[start_depth - 1] = vec_failing_set_[start_depth];
                idx_[start_depth] = num_local_candidates_[start_depth];
            } else {
                vec_failing_set_[start_depth - 1] |= vec_failing_set_[start_depth];
            }
#endif
        }
    }

    return false;
}

uint32_t leapfrogtriejoin::compute_local_candidates(uint32_t depth) {
    uint32_t u = vertex_ordering_[depth];
    uint32_t bn_count = get_backward_neighbors_count(depth);
    uint32_t* bns = get_backward_neighbors(depth);
    uint32_t* bns_depth = get_backward_neighbors_depth(depth);
    uint32_t lc_count = 0;
#ifdef INTERSECTION_CACHE
    bool recompute = false;
#endif
    if (bn_count == 1) {
        uint32_t bn = bns[0];
        uint32_t key = idx_embedding_[bns_depth[0]];
        local_candidates_[depth] = catalog_->get_core_relation_children(bn, u, key, lc_count);
    }
    else {
#ifndef INTERSECTION_CACHE
        uint32_t bn1 = bns[0];
        uint32_t bn2 = bns[1];
        uint32_t key1 = idx_embedding_[bns_depth[0]];
        uint32_t key2 = idx_embedding_[bns_depth[1]];

#ifndef SPARSE_BITMAP
        uint32_t lc_count1 = 0;
        uint32_t lc_count2 = 0;

        uint32_t *lc1 = catalog_->get_core_relation_children(bn1, u, key1, lc_count1);
        uint32_t *lc2 = catalog_->get_core_relation_children(bn2, u, key2, lc_count2);
        ComputeSetIntersection::ComputeCandidates(lc1, lc_count1, lc2, lc_count2, local_candidates_[depth],
                                                  lc_count);

#else
        BSRSet& bsr_set1 = catalog_->bsr_relations_[bn1][u].bsrs[key1];
        BSRSet& bsr_set2 = catalog_->bsr_relations_[bn2][u].bsrs[key2];

        cached_bsr_sets_[depth].size_ = intersect_qfilter_bsr_hybrid(bsr_set1.base_, bsr_set1.states_, bsr_set1.size_,
                                                                     bsr_set2.base_, bsr_set2.states_, bsr_set2.size_,
                                                                     cached_bsr_sets_[depth].base_, cached_bsr_sets_[depth].states_);
#endif

        for (uint32_t i = 2; i < bn_count; ++i) {
            uint32_t bn = bns[i];
            uint32_t key = idx_embedding_[bns_depth[i]];

#ifndef SPARSE_BITMAP
            lc_count1 = 0;
            lc1 = catalog_->get_core_relation_children(bn, u, key, lc_count1);

            uint32_t temp_lc_count = 0;

            ComputeSetIntersection::ComputeCandidates(lc1, lc_count1, local_candidates_[depth], lc_count,
                                                      si_buffer_, temp_lc_count);
            std::swap(local_candidates_[depth], si_buffer_);
            lc_count = temp_lc_count;

            if (lc_count == 0) {
                goto EXIT_EMPTY_SET;
            }
#else
            BSRSet& bsr_set = catalog_->bsr_relations_[bn][u].bsrs[key];



            buffer_bsr_set1_.size_ = intersect_qfilter_bsr_hybrid(bsr_set.base_, bsr_set.states_, bsr_set.size_,
                                                                  cached_bsr_sets_[depth].base_, cached_bsr_sets_[depth].states_,
                                                                  cached_bsr_sets_[depth].size_, buffer_bsr_set1_.base_,
                                                                  buffer_bsr_set1_.states_);

            std::swap(cached_bsr_sets_[depth], buffer_bsr_set1_);
#endif
        }

#ifdef SPARSE_BITMAP
        lc_count = (uint32_t)offline_bsr_trans_uint(cached_bsr_sets_[depth].base_, cached_bsr_sets_[depth].states_,
                                                    cached_bsr_sets_[depth].size_, (int*) local_candidates_[depth]);
#endif

#else

        if (intersection_cache_type_[depth] == INTERSECTION_CACHE_TYPE::none || !is_intersection_cache_valid(depth)) {
            uint32_t cached_bn_count = cached_bns_offset_[depth + 1] - cached_bns_offset_[depth];
            uint32_t *cached_bns = cached_bns_ + cached_bns_offset_[depth];
            uint32_t *cached_bns_depth = cached_bns_depth_ + cached_bns_offset_[depth];
            uint32_t cached_bn1 = cached_bns[0];
            uint32_t cached_bn2 = cached_bns[1];

            uint32_t key1 = idx_embedding_[cached_bns_depth[0]];
            uint32_t key2 = idx_embedding_[cached_bns_depth[1]];

#ifndef SPARSE_BITMAP
            uint32_t cached_candidates_count = 0;
            uint32_t lc_count1 = 0;

            uint32_t lc_count2 = 0;

            uint32_t *lc1 = catalog_->get_core_relation_children(cached_bn1, u, key1, lc_count1);
            uint32_t *lc2 = catalog_->get_core_relation_children(cached_bn2, u, key2, lc_count2);

            ComputeSetIntersection::ComputeCandidates(lc1, lc_count1, lc2, lc_count2, cached_candidates_[depth],
                                                      cached_candidates_count);

#else
            BSRSet& bsr_set1 = catalog_->bsr_relations_[cached_bn1][u].bsrs[key1];
            BSRSet& bsr_set2 = catalog_->bsr_relations_[cached_bn2][u].bsrs[key2];

            cached_bsr_sets_[depth].size_ = intersect_qfilter_bsr_hybrid(bsr_set1.base_, bsr_set1.states_, bsr_set1.size_,
                                                                        bsr_set2.base_, bsr_set2.states_, bsr_set2.size_,
                                                                        cached_bsr_sets_[depth].base_, cached_bsr_sets_[depth].states_);
#endif
            for (uint32_t i = 2; i < cached_bn_count; ++i) {
                uint32_t cached_bn = cached_bns[i];
                uint32_t key = idx_embedding_[cached_bns_depth[i]];

#ifndef SPARSE_BITMAP
                lc_count1 = 0;
                lc1 = catalog_->get_core_relation_children(cached_bn, u, key, lc_count1);

                ComputeSetIntersection::ComputeCandidates(lc1, lc_count1, cached_candidates_[depth], cached_candidates_count,
                                                          si_buffer_, cached_candidates_count);
                std::swap(cached_candidates_[depth], si_buffer_);
#else
                BSRSet& bsr_set = catalog_->bsr_relations_[cached_bn][u].bsrs[key];

                buffer_bsr_set1_.size_ = intersect_qfilter_bsr_hybrid(bsr_set.base_, bsr_set.states_, bsr_set.size_,
                                                                      cached_bsr_sets_[depth].base_, cached_bsr_sets_[depth].states_,
                                                                      cached_bsr_sets_[depth].size_, buffer_bsr_set1_.base_,
                                                                      buffer_bsr_set1_.states_);

                std::swap(cached_bsr_sets_[depth], buffer_bsr_set1_);
#endif
            }

#ifndef SPARSE_BITMAP
            cached_candidates_count_[depth] = cached_candidates_count;
#endif
            recompute = true;
        }
        if (intersection_cache_type_[depth] != INTERSECTION_CACHE_TYPE::partial) {
#ifndef SPARSE_BITMAP
            lc_count = cached_candidates_count_[depth];
#else
            lc_count = num_local_candidates_[depth];
#endif
            if (recompute) {
#ifndef SPARSE_BITMAP
                std::swap(cached_candidates_[depth], local_candidates_[depth]);
#else
                lc_count = (uint32_t)offline_bsr_trans_uint(cached_bsr_sets_[depth].base_, cached_bsr_sets_[depth].states_,
                                                                cached_bsr_sets_[depth].size_, (int*) local_candidates_[depth]);
#endif
            }
        }
        else {
            uint32_t no_cached_bn_count = no_cached_bns_offset_[depth + 1] - no_cached_bns_offset_[depth];
            uint32_t* no_cached_bns = no_cached_bns_ + no_cached_bns_offset_[depth];
            uint32_t* no_cached_bns_depth = no_cached_bns_depth_ + no_cached_bns_offset_[depth];

            uint32_t no_cached_bn = no_cached_bns[0];
            uint32_t key = idx_embedding_[no_cached_bns_depth[0]];


#ifndef SPARSE_BITMAP
            uint32_t lc_count1 = 0;
            uint32_t* lc1 = catalog_->get_core_relation_children(no_cached_bn, u, key, lc_count1);

            ComputeSetIntersection::ComputeCandidates(lc1, lc_count1, cached_candidates_[depth], cached_candidates_count_[depth],
                                                      local_candidates_[depth], lc_count);

#else
            BSRSet& bsr_set1 = catalog_->bsr_relations_[no_cached_bn][u].bsrs[key];

            buffer_bsr_set1_.size_ = intersect_qfilter_bsr_hybrid(bsr_set1.base_, bsr_set1.states_, bsr_set1.size_,
                                                                 cached_bsr_sets_[depth].base_, cached_bsr_sets_[depth].states_, cached_bsr_sets_[depth].size_,
                                                                 buffer_bsr_set1_.base_, buffer_bsr_set1_.states_);
#endif
            for (uint32_t i = 1; i < no_cached_bn_count; ++i) {
                no_cached_bn = no_cached_bns[i];

                key = idx_embedding_[no_cached_bns_depth[i]];


#ifndef SPARSE_BITMAP
                lc_count1 = 0;
                lc1 = catalog_->get_core_relation_children(no_cached_bn, u, key, lc_count1);

                ComputeSetIntersection::ComputeCandidates(lc1, lc_count1, local_candidates_[depth], lc_count,
                                                          si_buffer_, lc_count);
                std::swap(local_candidates_[depth], si_buffer_);
#else
                BSRSet& bsr_set = catalog_->bsr_relations_[no_cached_bn][u].bsrs[key];

                buffer_bsr_set2_.size_ = intersect_qfilter_bsr_hybrid(bsr_set.base_, bsr_set.states_, bsr_set.size_,
                                                                     buffer_bsr_set1_.base_, buffer_bsr_set1_.states_, buffer_bsr_set1_.size_,
                                                                     buffer_bsr_set2_.base_, buffer_bsr_set2_.states_);

                std::swap(buffer_bsr_set1_, buffer_bsr_set2_);
#endif
            }

#ifdef SPARSE_BITMAP
            lc_count = (uint32_t)offline_bsr_trans_uint(buffer_bsr_set1_.base_, buffer_bsr_set1_.states_,
                                                       buffer_bsr_set1_.size_, (int*) local_candidates_[depth]);
#endif
        }
#endif
    }


    if (lc_count == 0){
        EXIT_EMPTY_SET:
        return 0;
    }

    return lc_count;
}

void leapfrogtriejoin::initialize() {
    
    fail_si_count_ = 0;
    iso_conflict_count_ = 0;
    invalid_core_pr_count_ = 0;
    invalid_leaf_pr_count_ = 0;
    enter_result_count_ = new uint64_t[num_vertex_];
    std::fill(enter_result_count_, enter_result_count_ + num_vertex_, 0);

    idx_ = new uint32_t[num_vertex_];
    num_local_candidates_ = new uint32_t[num_vertex_];
    embedding_ = new uint32_t[num_vertex_];
    si_buffer_ = new uint32_t[catalog_->max_num_candidates_per_vertex_];
    last_level_candidates_ = new uint32_t[catalog_->max_num_candidates_per_vertex_];

    bns_offset_ = new uint32_t[num_vertex_ + 1];
    fns_offset_ = new uint32_t[num_vertex_ + 1];
    bns_ = new uint32_t[num_vertex_ * num_vertex_];
    bns_depth_ = new uint32_t[num_vertex_ * num_vertex_];
    fns_ = new uint32_t[num_vertex_ * num_vertex_];
    fns_depth_ = new uint32_t[num_vertex_ * num_vertex_];
    nodeId = new ui[catalog_->data_graph_->getVerticesCount()];   
     set_intersection_cost_ = new uint64_t[num_vertex_];
    set_intersection_count_ = new uint64_t[num_vertex_];
    intermediate_result_count_ = new uint64_t[num_vertex_];
    std::fill(set_intersection_cost_, set_intersection_cost_ + num_vertex_, 0);
    std::fill(set_intersection_count_, set_intersection_count_ + num_vertex_, 0);
    std::fill(intermediate_result_count_, intermediate_result_count_ + num_vertex_, 0);
    std::fill(nodeId, nodeId + catalog_->data_graph_->getVerticesCount(), 0);


#ifdef INTERSECTION_CACHE
    cached_bns_offset_ = new uint32_t[num_vertex_ + 1];
    cached_bns_ = new uint32_t[num_vertex_ * num_vertex_];
    cached_bns_depth_ = new uint32_t[num_vertex_ * num_vertex_];
    no_cached_bns_offset_ = new uint32_t[num_vertex_ + 1];
    no_cached_bns_ = new uint32_t[num_vertex_ * num_vertex_];
    no_cached_bns_depth_ = new uint32_t[num_vertex_ * num_vertex_];
    intersection_cache_type_ = new INTERSECTION_CACHE_TYPE[num_vertex_];
    mapping_to_cached_bns_ = new uint32_t[num_vertex_ * num_vertex_];
    std::fill(mapping_to_cached_bns_, mapping_to_cached_bns_ + num_vertex_ * num_vertex_, catalog_->max_data_vertex_id_ + 1);
#endif

    visited_ = new bool[catalog_->max_data_vertex_id_ + 1];
    memset(visited_, 0, sizeof(bool) * (catalog_->max_data_vertex_id_ + 1));

    idx_embedding_ = new uint32_t[num_vertex_];

#ifdef SPARSE_BITMAP
    cached_bsr_sets_ = new BSRSet[num_vertex_];
    for (uint32_t i = 0; i < num_vertex_; ++i) {
        cached_bsr_sets_[i].states_ = new int[catalog_->max_num_candidates_per_vertex_];
        cached_bsr_sets_[i].base_ = new int[catalog_->max_num_candidates_per_vertex_];
    }
    buffer_bsr_set1_.states_ = new int[catalog_->max_num_candidates_per_vertex_];
    buffer_bsr_set1_.base_ = new int[catalog_->max_num_candidates_per_vertex_];
    buffer_bsr_set2_.states_ = new int[catalog_->max_num_candidates_per_vertex_];
    buffer_bsr_set2_.base_ = new int[catalog_->max_num_candidates_per_vertex_];
#endif

    if (output_ != nullptr) {
        materialize_ = true;
    }
}

void leapfrogtriejoin::clear() {
    delete[] idx_;
    delete[] num_local_candidates_;
    delete[] embedding_;
    delete[] si_buffer_;
    delete[] last_level_candidates_;

    for (uint32_t i = 0; i < num_vertex_; ++i) {
        uint32_t bn_count = bns_offset_[i + 1] - bns_offset_[i];
        if (bn_count > 1) {
            delete[] local_candidates_[i];
        }
    }
    delete[] local_candidates_;

    delete[] bns_offset_;
    delete[] fns_offset_;
    delete[] bns_;
    delete[] fns_;
    delete[] bns_depth_;
    delete[] fns_depth_;

    delete[] set_intersection_cost_;
    delete[] intermediate_result_count_;
    delete[] set_intersection_count_;

#ifdef INTERSECTION_CACHE
    for (uint32_t i = 0; i < num_vertex_; ++i) {
        delete[] cached_candidates_[i];
    }

    delete[] cached_candidates_;
    delete[] cached_candidates_count_;
    delete[] cached_bns_offset_;
    delete[] cached_bns_;
    delete[] cached_bns_depth_;
    delete[] no_cached_bns_offset_;
    delete[] no_cached_bns_;
    delete[] no_cached_bns_depth_;
    delete[] intersection_cache_type_;
    delete[] mapping_to_cached_bns_;
#endif

    delete[] visited_;

#ifdef SPARSE_BITMAP
    delete[] cached_bsr_sets_;
#endif

    delete[] idx_embedding_;
}

void leapfrogtriejoin::set_idx(uint32_t depth) {
    for (uint32_t i = 0; i < depth; ++i) {
        uint32_t u = vertex_ordering_[i];
        if (get_forward_neighbors_count(i) != 0) {
            idx_embedding_[i] = catalog_->get_candidate_index(u, embedding_[i]);
        }
    }
}

void leapfrogtriejoin::initialize_bn_fn() {
    uint32_t bn_offset = 0;
    uint32_t fn_offset = 0;
    for (uint32_t i = 0; i < num_vertex_; ++i) {
        uint32_t u = vertex_ordering_[i];
        bns_offset_[i] = bn_offset;
        fns_offset_[i] = fn_offset;
        for (uint32_t j = 0; j < num_vertex_; ++j) {
            uint32_t v = vertex_ordering_[j];
            if (catalog_->is_adjacent(u, v)) {
                if (j < i) {
                    bns_[bn_offset] = v;
                    bns_depth_[bn_offset++] = j;
                }
                else if (j > i) {
                    fns_[fn_offset] = v;
                    fns_depth_[fn_offset++] = j;
                }
            }
        }
    }
    bns_offset_[num_vertex_] = bn_offset;
    fns_offset_[num_vertex_] = fn_offset;

    local_candidates_ = new uint32_t*[num_vertex_];
    for (uint32_t i = 0; i < num_vertex_; ++i) {
        uint32_t bn_count = bns_offset_[i + 1] - bns_offset_[i];

        if (i >= input_->get_tuple_length() && bn_count > 1) {
            local_candidates_[i] = new uint32_t[catalog_->max_num_candidates_per_vertex_];
        }
        else {
            local_candidates_[i] = nullptr;
        }
    }
}

#ifdef INTERSECTION_CACHE

void leapfrogtriejoin::initialize_intersection_cache() {
    uint32_t cached_bn_offset = 0;
    uint32_t no_cached_bn_offset = 0;

    for (uint32_t i = 0; i < num_vertex_; ++i) {
        cached_bns_offset_[i] = cached_bn_offset;
        no_cached_bns_offset_[i] = no_cached_bn_offset;
        intersection_cache_type_[i] = INTERSECTION_CACHE_TYPE::none;

        
        if (bns_offset_[i + 1] - bns_offset_[i] >= 2) {
            for (uint32_t j = bns_offset_[i]; j < bns_offset_[i + 1]; ++j) {
                if (bns_depth_[j] < i - 1) {
                    cached_bns_depth_[cached_bn_offset] = bns_depth_[j];
                    cached_bns_[cached_bn_offset++] = bns_[j];
                } else {
                    no_cached_bns_depth_[no_cached_bn_offset] = bns_depth_[j];
                    no_cached_bns_[no_cached_bn_offset++] = bns_[j];
                }
            }

            
            
            if (cached_bn_offset - cached_bns_offset_[i] == 1) {
                no_cached_bn_offset = no_cached_bns_offset_[i];
                cached_bns_depth_[cached_bn_offset] = no_cached_bns_depth_[no_cached_bn_offset];
                cached_bns_[cached_bn_offset++] = no_cached_bns_[no_cached_bn_offset];
            } else {
                if (no_cached_bn_offset - no_cached_bns_offset_[i] > 0) {
                    intersection_cache_type_[i] = INTERSECTION_CACHE_TYPE::partial;
                } else {
                    intersection_cache_type_[i] = INTERSECTION_CACHE_TYPE::full;
                }
            }
        }
    }
    cached_bns_offset_[num_vertex_] = cached_bn_offset;
    no_cached_bns_offset_[num_vertex_] = no_cached_bn_offset;

    cached_candidates_count_ = new uint32_t[num_vertex_];
    cached_candidates_ = new uint32_t*[num_vertex_];
    for (uint32_t i = 0; i < num_vertex_; ++i) {
        if (i >= input_->get_tuple_length() && bns_offset_[i + 1] - bns_offset_[i] >= 2) {
            cached_candidates_[i] = new uint32_t[catalog_->max_num_candidates_per_vertex_];
        }
        else {
            cached_candidates_[i] = nullptr;
        }
    }
}

bool leapfrogtriejoin::is_intersection_cache_valid(uint32_t depth) {
    bool valid = true;
    for (uint32_t i = cached_bns_offset_[depth]; i < cached_bns_offset_[depth + 1]; ++i) {
        if (mapping_to_cached_bns_[i] != embedding_[cached_bns_depth_[i]]) {
            valid = false;
            mapping_to_cached_bns_[i] = embedding_[cached_bns_depth_[i]];
        }
    }
    return valid;
}

#endif

#ifdef FAILING_SET_PRUNING

void leapfrogtriejoin::initialize_failing_set_pruning() {
    ancestors_.resize(num_vertex_);
    vec_failing_set_.resize(num_vertex_);
    reverse_embedding_.reserve(MAXIMUM_QUERY_GRAPH_SIZE * 2);
    compute_ancestors();
}

void leapfrogtriejoin::compute_ancestors() {
    
    for (uint32_t i = 0; i < num_vertex_; ++i) {
        uint32_t u = vertex_ordering_[i];
        ancestors_[u].set(u);

        for (uint32_t j = bns_offset_[i]; j < bns_offset_[i + 1]; ++j) {
            uint32_t u_bn = bns_[j];
            ancestors_[u] |= ancestors_[u_bn];
        }
    }
}

#endif
