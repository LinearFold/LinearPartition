/*
 *LinearPartition.cpp*
 The main code for LinearPartition: Linear-Time Approximation of 
                                    RNA Folding Partition Function 
                                    and Base Pairing Probabilities

 author: He Zhang
 created by: 03/2019
*/

#include <fstream>
#include <iostream>
#include <sys/time.h>
#include <stack>
#include <tuple>
#include <cassert>
#include <unordered_map>
#include <algorithm>
#include <string>
#include <map>

#include "LinearPartition.h"
#include "Utils/utility.h"
#include "Utils/utility_v.h"
#include "bpp.cpp"

#define SPECIAL_HP


#define pf

using namespace std;

unsigned long quickselect_partition(vector<pair<float, int>>& scores, unsigned long lower, unsigned long upper) {
    float pivot = scores[upper].first;
    while (lower < upper) {
        while (scores[lower].first < pivot) ++lower;
        while (scores[upper].first > pivot) --upper;
        if (scores[lower].first == scores[upper].first) ++lower;
        else if (lower < upper) swap(scores[lower], scores[upper]);
    }
    return upper;
}

// in-place quick-select
float quickselect(vector<pair<float, int>>& scores, unsigned long lower, unsigned long upper, unsigned long k) {
    if ( lower == upper ) return scores[lower].first;
    unsigned long split = quickselect_partition(scores, lower, upper);
    unsigned long length = split - lower + 1;
    if (length == k) return scores[split].first;
    else if (k  < length) return quickselect(scores, lower, split-1, k);
    else return quickselect(scores, split+1, upper, k - length);
}


float BeamCKYParser::beam_prune(std::unordered_map<int, State> &beamstep) {
    scores.clear();
    for (auto &item : beamstep) {
        int i = item.first;
        State &cand = item.second;
        int k = i - 1;
        float newalpha = (k >= 0 ? bestC[k].alpha : 0.0) + cand.alpha;
        scores.push_back(make_pair(newalpha, i));
    }
    if (scores.size() <= beam) return VALUE_MIN;
    float threshold = quickselect(scores, 0, scores.size() - 1, scores.size() - beam);
    for (auto &p : scores) {
        if (p.first < threshold) beamstep.erase(p.second);
    }

    return threshold;
}


void BeamCKYParser::prepare(unsigned len) {
    seq_length = len;

    bestH.clear();
    bestH.resize(seq_length);
    bestP.clear();
    bestP.resize(seq_length);
    bestM2.clear();
    bestM2.resize(seq_length);
    bestM.clear();
    bestM.resize(seq_length);
    bestC.clear();
    bestC.resize(seq_length);
    bestMulti.clear();
    bestMulti.resize(seq_length);

    nucs.clear();
    nucs.resize(seq_length);
    scores.reserve(seq_length);
}


BeamCKYParser::DecoderResult BeamCKYParser::parse(string& seq) {
      
    struct timeval parse_starttime, parse_endtime;

    gettimeofday(&parse_starttime, NULL);

    prepare(static_cast<unsigned>(seq.length()));

    for (int i = 0; i < seq_length; ++i)
        nucs[i] = GET_ACGU_NUM(seq[i]);

    vector<int> next_pair[NOTON];
    {
        for (int nuci = 0; nuci < NOTON; ++nuci) {
            // next_pair
            next_pair[nuci].resize(seq_length, -1);
            int next = -1;
            for (int j = seq_length-1; j >=0; --j) {
                next_pair[nuci][j] = next;
                if (_allowed_pairs[nuci][nucs[j]]) next = j;
            }
        }
    }


#ifdef SPECIAL_HP
#ifdef lpv
    v_init_tetra_hex_tri(seq, seq_length, if_tetraloops, if_hexaloops, if_triloops);
#else // only for eval
    if (is_verbose)
        v_init_tetra_hex_tri(seq, seq_length, if_tetraloops, if_hexaloops, if_triloops);
#endif
#endif

#ifdef lpv
        if(seq_length > 0) bestC[0].alpha = 0.0;
        if(seq_length > 1) bestC[1].alpha = 0.0;
#else
        if(seq_length > 0) Fast_LogPlusEquals(bestC[0].alpha, score_external_unpaired(0, 0));
        if(seq_length > 1) Fast_LogPlusEquals(bestC[1].alpha, score_external_unpaired(0, 1));
#endif

    value_type newscore;
    for(int j = 0; j < seq_length; ++j) {
        int nucj = nucs[j];
        int nucj1 = (j+1) < seq_length ? nucs[j+1] : -1;

        unordered_map<int, State>& beamstepH = bestH[j];
        unordered_map<int, State>& beamstepMulti = bestMulti[j];
        unordered_map<int, State>& beamstepP = bestP[j];
        unordered_map<int, State>& beamstepM2 = bestM2[j];
        unordered_map<int, State>& beamstepM = bestM[j];
        State& beamstepC = bestC[j];

        // beam of H
        {
            if (beam > 0 && beamstepH.size() > beam) beam_prune(beamstepH);

            {
                // for nucj put H(j, j_next) into H[j_next]
                int jnext = next_pair[nucj][j];
                if (no_sharp_turn) while (jnext - j < 4 && jnext != -1) jnext = next_pair[nucj][jnext];
                if (jnext != -1) {
                    int nucjnext = nucs[jnext];
                    int nucjnext_1 = (jnext - 1) > -1 ? nucs[jnext - 1] : -1;
#ifdef lpv
                        int tetra_hex_tri = -1;
#ifdef SPECIAL_HP
                        if (jnext-j-1 == 4) // 6:tetra
                            tetra_hex_tri = if_tetraloops[j];
                        else if (jnext-j-1 == 6) // 8:hexa
                            tetra_hex_tri = if_hexaloops[j];
                        else if (jnext-j-1 == 3) // 5:tri
                            tetra_hex_tri = if_triloops[j];
#endif
                        newscore = - v_score_hairpin(j, jnext, nucj, nucj1, nucjnext_1, nucjnext, tetra_hex_tri);
                        Fast_LogPlusEquals(bestH[jnext][j].alpha, newscore/kT);
#else
                        newscore = score_hairpin(j, jnext, nucj, nucj1, nucjnext_1, nucjnext);
                        Fast_LogPlusEquals(bestH[jnext][j].alpha, newscore);
#endif
                }
            }

            {
                // for every state h in H[j]
                //   1. extend h(i, j) to h(i, jnext)
                //   2. generate p(i, j)
                for (auto &item : beamstepH) {
                    int i = item.first;
                    State &state = item.second;
                    int nuci = nucs[i];
                    int jnext = next_pair[nuci][j];

                    if (jnext != -1) {
                        int nuci1 = (i + 1) < seq_length ? nucs[i + 1] : -1;
                        int nucjnext = nucs[jnext];
                        int nucjnext_1 = (jnext - 1) > -1 ? nucs[jnext - 1] : -1;

                        // 1. extend h(i, j) to h(i, jnext)
                        // value_type newscore;

#ifdef lpv
                            int tetra_hex_tri = -1;
#ifdef SPECIAL_HP
                            if (jnext-i-1 == 4) // 6:tetra
                                tetra_hex_tri = if_tetraloops[i];
                            else if (jnext-i-1 == 6) // 8:hexa
                                tetra_hex_tri = if_hexaloops[i];
                            else if (jnext-i-1 == 3) // 5:tri
                                tetra_hex_tri = if_triloops[i];
#endif
                            newscore = - v_score_hairpin(i, jnext, nuci, nuci1, nucjnext_1, nucjnext, tetra_hex_tri);
                            Fast_LogPlusEquals(bestH[jnext][i].alpha, (newscore/kT));
#else
                            newscore = score_hairpin(i, jnext, nuci, nuci1, nucjnext_1, nucjnext);
                            Fast_LogPlusEquals(bestH[jnext][i].alpha, newscore);
#endif
                    }

                    // 2. generate p(i, j)
                    Fast_LogPlusEquals(beamstepP[i].alpha, state.alpha);
                }
            }
        }
        if (j == 0) continue;

        // beam of Multi
        {
            if (beam > 0 && beamstepMulti.size() > beam) beam_prune(beamstepMulti);

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
                            Fast_LogPlusEquals(bestMulti[jnext][i].alpha, (state.alpha));
#else
                            newscore = score_multi_unpaired(j, jnext - 1);
                            Fast_LogPlusEquals(bestMulti[jnext][i].alpha, state.alpha + newscore);
#endif
                    }
                }

                // 2. generate P (i, j)
                {
#ifdef lpv
                    value_type score_multi = - v_score_multi(i, j, nuci, nuci1, nucs[j-1], nucj, seq_length);
                    Fast_LogPlusEquals(beamstepP[i].alpha, (state.alpha + score_multi/kT));
#else
                    newscore = score_multi(i, j, nuci, nuci1, nucs[j-1], nucj, seq_length);
                    Fast_LogPlusEquals(beamstepP[i].alpha, state.alpha + newscore);
#endif
                }
            }
        }

        // beam of P
        {   
            if (beam > 0 && beamstepP.size() > beam) beam_prune(beamstepP);

            // for every state in P[j]
            //   1. generate new helix/bulge
            //   2. M = P
            //   3. M2 = M + P
            //   4. C = C + P
            for(auto& item : beamstepP) {
                int i = item.first;
                State& state = item.second;
                int nuci = nucs[i];
                int nuci_1 = (i-1>-1) ? nucs[i-1] : -1;

                // 1. generate new helix / single_branch
                // new state is of shape p..i..j..q
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
                                    Fast_LogPlusEquals(bestP[q][p].alpha, (state.alpha + score_single/kT));
#else
                                    newscore = score_helix(nucp, nucp1, nucq_1, nucq);
                                    Fast_LogPlusEquals(bestP[q][p].alpha, state.alpha + newscore);

#endif
                            } else {
                                // single branch
#ifdef lpv
                                    int score_single = - v_score_single(p,q,i,j, nucp, nucp1, nucq_1, nucq,
                                                   nuci_1, nuci, nucj, nucj1);
                                    Fast_LogPlusEquals(bestP[q][p].alpha, (state.alpha + score_single/kT));
#else
                                    newscore = score_junction_B(p, q, nucp, nucp1, nucq_1, nucq) +
                                        precomputed +
                                        score_single_without_junctionB(p, q, i, j,
                                                                       nuci_1, nuci, nucj, nucj1);
                                    Fast_LogPlusEquals(bestP[q][p].alpha, state.alpha + newscore);
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
                        Fast_LogPlusEquals(beamstepM[i].alpha, (state.alpha + score_M1/kT));
#else
                        newscore = score_M1(i, j, j, nuci_1, nuci, nucj, nucj1, seq_length);
                        Fast_LogPlusEquals(beamstepM[i].alpha, state.alpha + newscore);
#endif
                }

                // 3. M2 = M + P
                int k = i - 1;
                if ( k > 0 && !bestM[k].empty()) {
#ifdef lpv
                    int M1_score = - v_score_M1(i, j, j, nuci_1, nuci, nucj, nucj1, seq_length);
                    float m1_alpha = state.alpha + M1_score/kT;
#else
                    value_type M1_score = score_M1(i, j, j, nuci_1, nuci, nucj, nucj1, seq_length);
                    float m1_alpha = state.alpha + M1_score;
#endif
                    for (auto &m : bestM[k]) {
                        int newi = m.first;
                        State& m_state = m.second;
                        Fast_LogPlusEquals(beamstepM2[newi].alpha, m_state.alpha + m1_alpha);
                    }
                }

                // 4. C = C + P
                {
                    int k = i - 1;
                    if (k >= 0) {
                      State& prefix_C = bestC[k];
                        int nuck = nuci_1;
                        int nuck1 = nuci;
#ifdef lpv
                        int score_external_paired = - v_score_external_paired(k+1, j, nuck, nuck1,
                                                             nucj, nucj1, seq_length);
                        Fast_LogPlusEquals(beamstepC.alpha, prefix_C.alpha + state.alpha + score_external_paired/kT);      
#else
                        newscore = score_external_paired(k+1, j, nuck, nuck1,
                                                             nucj, nucj1, seq_length);
                        Fast_LogPlusEquals(beamstepC.alpha, prefix_C.alpha + state.alpha + newscore);
#endif
                    } else {
#ifdef lpv
                        int score_external_paired = - v_score_external_paired(0, j, -1, nucs[0],
                                                                 nucj, nucj1, seq_length);
                        Fast_LogPlusEquals(beamstepC.alpha, state.alpha + score_external_paired/kT);       
#else
                        newscore = score_external_paired(0, j, -1, nucs[0],
                                                             nucj, nucj1, seq_length);
                        Fast_LogPlusEquals(beamstepC.alpha, state.alpha + newscore);
#endif
                    }
                }
            }
        }


        // beam of M2
        {
            if (beam > 0 && beamstepM2.size() > beam) beam_prune(beamstepM2);

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
                        Fast_LogPlusEquals(bestMulti[q][p].alpha, state.alpha);      

#else
                        newscore = score_multi_unpaired(p+1, i-1) +
                                    score_multi_unpaired(j+1, q-1);
                        Fast_LogPlusEquals(bestMulti[q][p].alpha, state.alpha + newscore);      
#endif
                        }
                    }
                }

                // 2. M = M2
                Fast_LogPlusEquals(beamstepM[i].alpha, state.alpha);  
            }
        }

        // beam of M
        {
            // float threshold = VALUE_MIN;
            // if (beam > 0 && beamstepM.size() > beam) threshold = beam_prune(beamstepM);
            if (beam > 0 && beamstepM.size() > beam) beam_prune(beamstepM);

            for(auto& item : beamstepM) {
                int i = item.first;
                State& state = item.second;
                if (j < seq_length-1) {
#ifdef lpv
                    Fast_LogPlusEquals(bestM[j+1][i].alpha, state.alpha); 
#else
                    newscore = score_multi_unpaired(j + 1, j + 1);
                    Fast_LogPlusEquals(bestM[j+1][i].alpha, state.alpha + newscore); 
#endif
                }
            }
        }

        // beam of C
        {
            // C = C + U
            if (j < seq_length-1) {
#ifdef lpv
                Fast_LogPlusEquals(bestC[j+1].alpha, beamstepC.alpha); 
                    
#else
                newscore = score_external_unpaired(j+1, j+1);
                Fast_LogPlusEquals(bestC[j+1].alpha, beamstepC.alpha + newscore); 
#endif
            }
        }

    }  // end of for-loo j

    State& viterbi = bestC[seq_length-1];

    gettimeofday(&parse_endtime, NULL);
    double parse_elapsed_time = parse_endtime.tv_sec - parse_starttime.tv_sec + (parse_endtime.tv_usec-parse_starttime.tv_usec)/1000000.0;

    // unsigned long nos_tot = nos_H + nos_P + nos_M2 + nos_Multi + nos_M + nos_C;

#ifdef lpv
    printf("Free Energy of Ensumble: %.2f kcal/mol\n", -kT * viterbi.alpha / 100.0);
#else
    printf("Log Partition Coefficient: %.5f\n", viterbi.alpha);
#endif

    if(is_verbose) printf("Partition Function Calculation Time: %f seconds.\n", parse_elapsed_time);

    if(!pf_only){
        outside(next_pair);
        cal_PairProb(viterbi);
    }

    fflush(stdout);

    // return {viterbi.alpha, nos_tot, parse_elapsed_time};
    return {viterbi.alpha, parse_elapsed_time};
}

BeamCKYParser::BeamCKYParser(int beam_size,
                             bool nosharpturn,
                             bool verbose,
                             string bppfile,
                             string bppfileindex,
                             bool pfonly)
    : beam(beam_size), 
      no_sharp_turn(nosharpturn), 
      is_verbose(verbose),
      bpp_file(bppfile),
      bpp_file_index(bppfileindex),
      pf_only(pfonly) {
#ifdef lpv
        initialize();
#else
        initialize();
        initialize_cachesingle();
#endif
}

int main(int argc, char** argv){

    struct timeval total_starttime, total_endtime;

    gettimeofday(&total_starttime, NULL);


    int beamsize = 100;
    bool sharpturn = false;
    bool is_verbose = false;
    string bpp_file;
    string bpp_prefix;
    bool pf_only = false;

    if (argc > 1) {
        beamsize = atoi(argv[1]);
        sharpturn = atoi(argv[2]) == 1;
        is_verbose = atoi(argv[3]) == 1;
        bpp_file = argv[4];
        bpp_prefix = argv[5];
        pf_only = atoi(argv[6]) == 1;
    }

    if (is_verbose) printf("beam size: %d\n", beamsize);

    // variables for decoding
    int num=0, total_len = 0;
    unsigned long long total_states = 0;
    double total_score = .0;
    double total_time = .0;

    int seq_index = 0;
    string bpp_file_index = "";
    for (string seq; getline(cin, seq);) {
        if (seq.length() == 0)
            continue;

        if (seq[0] == ';' || seq[0] == '>') {
            printf("%s\n", seq.c_str());
            continue;
        }

        if (!isalpha(seq[0])){
            printf("Unrecognized sequence: %s\n", seq.c_str());
            continue;
        }

        if (!bpp_prefix.empty()) {
            seq_index ++;
            bpp_file_index = bpp_prefix + to_string(seq_index);
        }

        // if (seq_index > 1 && bpp_prefix.empty() && bpp_file.empty()){
        //     printf("\nWARNING: use --prefix mode for multiple sequences\n");
        //     return 0;
        // }
        printf("%s\n", seq.c_str());
        
        // convert to uppercase
        transform(seq.begin(), seq.end(), seq.begin(), ::toupper);

        // convert T to U
        replace(seq.begin(), seq.end(), 'T', 'U');

        // lhuang: moved inside loop, fixing an obscure but crucial bug in initialization
        BeamCKYParser parser(beamsize, !sharpturn, is_verbose, bpp_file, bpp_file_index, pf_only);

        BeamCKYParser::DecoderResult result = parser.parse(seq);
    }

    gettimeofday(&total_endtime, NULL);
    double total_elapsed_time = total_endtime.tv_sec - total_starttime.tv_sec + (total_endtime.tv_usec-total_starttime.tv_usec)/1000000.0;

    // if(is_verbose) printf("Total Time: %f\n", total_elapsed_time);

    return 0;
}
