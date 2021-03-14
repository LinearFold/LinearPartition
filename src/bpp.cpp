/*
 *bpp.cpp*
 The main code for base pair probability calculation.

 author: He Zhang
 created by: 04/2019
*/

#include <stdio.h> 
#include <set>
#include <algorithm>
#include "LinearPartition.h"

using namespace std;

void BeamCKYParser::output_to_file(string file_name, const char * type) {
    if(!file_name.empty()) {
        printf("Outputing base pairing probability matrix to %s...\n", file_name.c_str()); 
        FILE *fptr = fopen(file_name.c_str(), type); 
        if (fptr == NULL) { 
            printf("Could not open file!\n"); 
            return; 
        }

        int turn = no_sharp_turn?3:0;
        for (int i = 1; i <= seq_length; i++) {
            for (int j = i + turn + 1; j <= seq_length; j++) {
                pair<int, int> key = make_pair(i,j);
                auto got = Pij.find(key);
                if (got != Pij.end()){
                    fprintf(fptr, "%d %d %.4e\n", i, j, got->second);
                }
            }
        }
        fprintf(fptr, "\n");
        fclose(fptr); 
        printf("Done!\n"); 
    }

    return;
}

void BeamCKYParser::output_to_file_threshknot_bpseq(string file_name, const char * type, map<int,int>& pairs, string & seq) {

    int i,j;
    char nuc;
    if(!file_name.empty()) {
        printf("Outputing ThreshKnot base pairing matrix to %s...\n", file_name.c_str()); 
        FILE *fptr = fopen(file_name.c_str(), type); 
        if (fptr == NULL) { 
            printf("Could not open file!\n"); 
            return; 
        }


        for (int i = 1; i <= seq_length; i++) {
            if (pairs.find(i) != pairs.end()){
                j = pairs[i];
            }
            else{
                j = 0;
            }
            nuc = seq[i-1];
            fprintf(fptr, "%d %c %d\n", i, nuc, j);
        }

        fprintf(fptr, "\n");
        fclose(fptr); 
        printf("Done!\n"); 
    }
    else{
        for (int i = 1; i <= seq_length; i++) {
            if (pairs.find(i) != pairs.end()){
                j = pairs[i];
            }
            else{
                j = 0;
            }
            nuc = seq[i-1];
            printf("%d %c %d\n", i, nuc, j);
        }
        printf("\n");
    }

}


void BeamCKYParser::cal_PairProb(State& viterbi) {

    for(int j=0; j<seq_length; j++){
        for(auto &item : bestP[j]){
            int i = item.first;
            State state = item.second;
            
            float temp_prob_inside = state.alpha + state.beta - viterbi.alpha;
            if (temp_prob_inside > float(-9.91152)) {
                float prob = Fast_Exp(temp_prob_inside);
                if(prob > 1.0) prob = 1.0;
                if(prob < bpp_cutoff) continue;
                Pij[make_pair(i+1, j+1)] = prob;
            }
        }
    }

    // -o mode: output to a single file with user specified name;
    // bpp matrices for different sequences are separated with empty lines
    if (!bpp_file.empty()){
        output_to_file(bpp_file, "a");
    } 

    // -prefix mode: output to multiple files with user specified prefix;
    else if (!bpp_file_index.empty()) {
        output_to_file(bpp_file_index, "w");
    }

    return;
}




string BeamCKYParser::back_trace(const int i, const int j, const vector<vector<int> >& back_pointer){
    // cout<<"line 74 "<<i<<' '<<j<<endl;

    if (i>j){
        return "";
    }

    if (back_pointer[i][j] == -1){
        if (i == j){
            return ".";
        }

        else{

            return "." + back_trace(i+1,j, back_pointer);
        }
    }

    else if (back_pointer[i][j] != 0){
        int k = back_pointer[i][j];

        assert(k + 1 > 0 && k + 1 <= seq_length);

        string temp;

        if (k == j){
            temp = "";
        }
        else{
            temp = back_trace(k+1,j, back_pointer);
        }

        return "(" + back_trace(i+1,k-1, back_pointer) + ")" + temp;
    }


    // return "";
    assert(false);
    return "";
}




void BeamCKYParser::ThreshKnot(string & seq){
    
    map<int, float> rowprob;
    vector<tuple<int, int, float> > prob_list;

    map<int, int> pairs;
    set<int> visited;

    for(auto& pij : Pij){
        auto i = pij.first.first; //index starts from 1
        auto j = pij.first.second; 
        auto score = pij.second;

        if (score < threshknot_threshold){
            continue;
        }

        prob_list.push_back(make_tuple(i,j,score));

        rowprob[i] = max(rowprob[i], score);
        rowprob[j] = max(rowprob[j], score);

    }


    for(auto& elem : prob_list){

        auto i = std::get<0>(elem);
        auto j = std::get<1>(elem);
        auto score =  std::get<2>(elem);



        if (score == rowprob[i] && score == rowprob[j]){

            if ((visited.find(i) != visited.end()) || (visited.find(j) != visited.end())){
                continue;
            }
            visited.insert(i);
            visited.insert(j);

            pairs[i] = j;
            pairs[j] = i;


        }
    }


    output_to_file_threshknot_bpseq(threshknot_file_index, "w", pairs, seq);


}


pair<float, string> BeamCKYParser::PairProb_MEA() {
    
    vector<vector<float> > OPT;
    OPT.resize(seq_length);

    for (int i = 0; i < seq_length; ++i){
        OPT[i].resize(seq_length);
    }

    vector<vector<float> > P;
    P.resize(seq_length);

    for (int i = 0; i < seq_length; ++i){
        P[i].resize(seq_length);
    }


    vector<vector<int> > back_pointer;
    back_pointer.resize(seq_length);

    for (int i = 0; i < seq_length; ++i){
        back_pointer[i].resize(seq_length);
    }

    vector<vector<int>> paired;

    paired.resize(seq_length);


    vector<float> Q;
    for (int i = 0; i < seq_length; ++i){
        Q.push_back(float(1.0));
    }



    int count = 0;
    for(auto& pij : Pij){
        count+=1;
        auto i = pij.first.first-1;
        auto j = pij.first.second-1;
        auto score = pij.second;

        P[i][j] = score;

        // if (score < 0.000001)
        //     continue;

        paired[i].push_back(j);
        // assert(i<j);
        Q[i] -= score;
        Q[j] -= score;

    }


    for (int i = 0; i < seq_length; ++i){
        std::sort (paired[i].begin(), paired[i].end());
    }


    // clock_t tt = clock();


    for (int l = 0; l< seq_length; l++){

        for (int i = 0; i<seq_length - l; i++){
            int j = i + l;

            if (i == j){
                OPT[i][j] = Q[i];
                back_pointer[i][j] = -1;
                continue;
            }

            OPT[i][j] = OPT[i][i] + OPT[i+1][j];
            back_pointer[i][j] = -1;

            for (int k : paired[i]){

                if (k>j){
                    break;
                }

                float temp_OPT_k1_j;
                if (k<j){
                    temp_OPT_k1_j = OPT[k+1][j];
                }
                else{

                    temp_OPT_k1_j = float(0.);
                }
 
                auto temp_score = 2 * gamma * P[i][k] + OPT[i+1][k-1] + temp_OPT_k1_j;
                if (OPT[i][j] < temp_score){
                    OPT[i][j] = temp_score;

                    back_pointer[i][j] = k;

                }

            }

        }
    
    }



    auto structure = back_trace(0,seq_length-1, back_pointer);
    return make_pair(OPT[0][seq_length-1], structure);

}









void BeamCKYParser::outside(vector<int> next_pair[]){
      
    struct timeval parse_starttime, parse_endtime;

    gettimeofday(&parse_starttime, NULL);

    bestC[seq_length-1].beta = 0.0;

    // from right to left
    value_type newscore;
    for(int j = seq_length-1; j > 0; --j) {
        int nucj = nucs[j];
        int nucj1 = (j+1) < seq_length ? nucs[j+1] : -1;

        unordered_map<int, State>& beamstepH = bestH[j];
        unordered_map<int, State>& beamstepMulti = bestMulti[j];
        unordered_map<int, State>& beamstepP = bestP[j];
        unordered_map<int, State>& beamstepM2 = bestM2[j];
        unordered_map<int, State>& beamstepM = bestM[j];
        State& beamstepC = bestC[j];

        // beam of C
        {
            // C = C + U
            if (j < seq_length-1) {
#ifdef lpv
            Fast_LogPlusEquals(beamstepC.beta, (bestC[j+1].beta));
                    
#else
            newscore = score_external_unpaired(j+1, j+1);
            Fast_LogPlusEquals(beamstepC.beta, bestC[j+1].beta + newscore);
#endif
            }
        }
    
        // beam of M
        {
            for(auto& item : beamstepM) {
                int i = item.first;
                State& state = item.second;
                if (j < seq_length-1) {
#ifdef lpv
                    Fast_LogPlusEquals(state.beta, bestM[j+1][i].beta);
#else
                    newscore = score_multi_unpaired(j + 1, j + 1);
                    Fast_LogPlusEquals(state.beta, bestM[j+1][i].beta + newscore);
#endif
                }
            }
        }

        // beam of M2
        {
            for(auto& item : beamstepM2) {
                int i = item.first;
                State& state = item.second;

                // 1. multi-loop
                {
                    for (int p = i-1; p >= std::max(i - SINGLE_MAX_LEN, 0); --p) {
                        int nucp = nucs[p];
                        int q = next_pair[nucp][j];
                        if (q != -1 && ((i - p - 1) <= SINGLE_MAX_LEN)) {
#ifdef lpv
                            Fast_LogPlusEquals(state.beta, bestMulti[q][p].beta);
#else
                            newscore = score_multi_unpaired(p+1, i-1) +
                                    score_multi_unpaired(j+1, q-1);
                            Fast_LogPlusEquals(state.beta, bestMulti[q][p].beta + newscore);
#endif
                        }
                    }
                }

                // 2. M = M2
                Fast_LogPlusEquals(state.beta, beamstepM[i].beta);
            }
        }

        // beam of P
        {  
            for(auto& item : beamstepP) {
                int i = item.first;
                State& state = item.second;
                int nuci = nucs[i];
                int nuci_1 = (i-1>-1) ? nucs[i-1] : -1;

                if (i >0 && j<seq_length-1) {
#ifndef lpv
                    value_type precomputed = score_junction_B(j, i, nucj, nucj1, nuci_1, nuci);
#endif
                    for (int p = i - 1; p >= std::max(i - SINGLE_MAX_LEN, 0); --p) {
                        int nucp = nucs[p];
                        int nucp1 = nucs[p + 1]; 
                        int q = next_pair[nucp][j];
                        while (q != -1 && ((i - p) + (q - j) - 2 <= SINGLE_MAX_LEN)) {
                            int nucq = nucs[q];
                            int nucq_1 = nucs[q - 1];

                            if (p == i - 1 && q == j + 1) {
                                // helix
#ifdef lpv
                                int score_single = -v_score_single(p,q,i,j, nucp, nucp1, nucq_1, nucq,
                                                             nuci_1, nuci, nucj, nucj1);
                                Fast_LogPlusEquals(state.beta, (bestP[q][p].beta + score_single/kT));
#else
                                newscore = score_helix(nucp, nucp1, nucq_1, nucq);
                                Fast_LogPlusEquals(state.beta, bestP[q][p].beta + newscore);
#endif
                            } else {
                                // single branch
#ifdef lpv
                                int score_single = - v_score_single(p,q,i,j, nucp, nucp1, nucq_1, nucq,
                                                   nuci_1, nuci, nucj, nucj1);
                                Fast_LogPlusEquals(state.beta, (bestP[q][p].beta + score_single/kT));
#else
                                newscore = score_junction_B(p, q, nucp, nucp1, nucq_1, nucq) +
                                        precomputed + 
                                        score_single_without_junctionB(p, q, i, j, nuci_1, nuci, nucj, nucj1);
                                Fast_LogPlusEquals(state.beta, bestP[q][p].beta + newscore);
#endif
                            }
                            q = next_pair[nucp][q];
                        }
                    }
                }

                // 2. M = P
                if(i > 0 && j < seq_length-1){
#ifdef lpv
                        int score_M1 = - v_score_M1(i, j, j, nuci_1, nuci, nucj, nucj1, seq_length);
                        Fast_LogPlusEquals(state.beta, (beamstepM[i].beta + score_M1/kT));
#else
                        newscore = score_M1(i, j, j, nuci_1, nuci, nucj, nucj1, seq_length);
                        Fast_LogPlusEquals(state.beta, beamstepM[i].beta + newscore);
#endif
                }

                // 3. M2 = M + P
                int k = i - 1;
                if ( k > 0 && !bestM[k].empty()) {
#ifdef lpv
                    int M1_score = - v_score_M1(i, j, j, nuci_1, nuci, nucj, nucj1, seq_length);
                    float m1_alpha = M1_score/kT;
#else
                    newscore = score_M1(i, j, j, nuci_1, nuci, nucj, nucj1, seq_length);
                    float m1_alpha = newscore;
#endif
                    float m1_plus_P_alpha = state.alpha + m1_alpha;
                    for (auto &m : bestM[k]) {
                        int newi = m.first;
                        State& m_state = m.second;
                        Fast_LogPlusEquals(state.beta, (beamstepM2[newi].beta + m_state.alpha + m1_alpha));
                        Fast_LogPlusEquals(m_state.beta, (beamstepM2[newi].beta + m1_plus_P_alpha));
                    }
                }

                // 4. C = C + P
                {
                    int k = i - 1;
                    if (k >= 0) {
                        int nuck = nuci_1;
                        int nuck1 = nuci;
#ifdef lpv
                        int score_external_paired = - v_score_external_paired(k+1, j, nuck, nuck1,
                                                                 nucj, nucj1, seq_length);
                        float external_paired_alpha_plus_beamstepC_beta = beamstepC.beta + score_external_paired/kT;

#else
                        newscore = score_external_paired(k+1, j, nuck, nuck1, nucj, nucj1, seq_length);
                        float external_paired_alpha_plus_beamstepC_beta = beamstepC.beta + newscore;
#endif
                        Fast_LogPlusEquals(bestC[k].beta, state.alpha + external_paired_alpha_plus_beamstepC_beta);
                        Fast_LogPlusEquals(state.beta, bestC[k].alpha + external_paired_alpha_plus_beamstepC_beta);
                    } else {
                        // value_type newscore;
#ifdef lpv
                        int score_external_paired = - v_score_external_paired(0, j, -1, nucs[0],
                                                                 nucj, nucj1, seq_length);
                        Fast_LogPlusEquals(state.beta, (beamstepC.beta + score_external_paired/kT));
#else
                        newscore = score_external_paired(0, j, -1, nucs[0],
                                                             nucj, nucj1, seq_length);
                        Fast_LogPlusEquals(state.beta, beamstepC.beta + newscore);
#endif
                    }
                }
            }
        }

        // beam of Multi
        {
            for(auto& item : beamstepMulti) {
                int i = item.first;
                State& state = item.second;

                int nuci = nucs[i];
                int nuci1 = nucs[i+1];
                int jnext = next_pair[nuci][j];

                // 1. extend (i, j) to (i, jnext)
                {
                    if (jnext != -1) {
#ifdef lpv
                        Fast_LogPlusEquals(state.beta, (bestMulti[jnext][i].beta));
#else
                        newscore = score_multi_unpaired(j, jnext - 1);
                        Fast_LogPlusEquals(state.beta, bestMulti[jnext][i].beta + newscore);
#endif
                    }
                }

                // 2. generate P (i, j)
                {
#ifdef lpv
                    int score_multi = - v_score_multi(i, j, nuci, nuci1, nucs[j-1], nucj, seq_length);
                    Fast_LogPlusEquals(state.beta, (beamstepP[i].beta + score_multi/kT));
#else
                    newscore = score_multi(i, j, nuci, nuci1, nucs[j-1], nucj, seq_length);
                    Fast_LogPlusEquals(state.beta, beamstepP[i].beta + newscore);
#endif
                }
            }
        }
    }  // end of for-loo j

    gettimeofday(&parse_endtime, NULL);
    double parse_elapsed_time = parse_endtime.tv_sec - parse_starttime.tv_sec + (parse_endtime.tv_usec-parse_starttime.tv_usec)/1000000.0;

// #ifdef testtime
    if(is_verbose) printf("Base Pairing Probabilities Calculation Time: %.2f seconds.\n", parse_elapsed_time);
// #endif

    fflush(stdout);

    return;
}

