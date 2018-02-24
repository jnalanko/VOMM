#ifndef LMA_SUPPORT_TESTS
#define LMA_SUPPORT_TESTS

#include "LMA_Support.hh"
#include <iostream>

using namespace std;

int64_t LMA_naive(int64_t open, sdsl::bp_support_g<>& bpr_bps, sdsl::bit_vector& marks){
    while(!marks[open]){
        if(open == 0) return -1; // Unmarked root
        open = bpr_bps.enclose(open);
    }
    return open;
}

void test_all_queries(sdsl::bit_vector& bpr, sdsl::bit_vector& marks){
    // Test all queries
    LMA_Support nmas(bpr,marks);
    sdsl::bp_support_g<> bpr_bps;
    sdsl::util::init_support(bpr_bps, &bpr);
    for(int64_t open = 0; open < bpr.size(); open++){
        if(bpr[open] == 1){
            int64_t correct = LMA_naive(open, bpr_bps, marks);
            int64_t test = nmas.LMA(open);
            assert(test == correct);
        }
    }
}

// Random double between fMin and fMax
double dRand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

// Test for a random input with N pairs of parenthesis such
// that each is marked with probility marked_prob

void test_random(int64_t N, double marked_prob){
    
    sdsl::bit_vector bpr(2*N,0);
    bpr[0] = 1; // Open root
    bpr[2*N-1] = 0; // Close root
    
    // Generate a sequnece of N-1 pairs inside the root i.e. in bpr[1..2N-2)
    int64_t balance = 0;
    int64_t n_open = 0;
    for(int64_t i = 1; i <= 2*N-2; i++){
        if(balance == 0){
            bpr[i] = 1;
        } else if(n_open == N-1){
            bpr[i] = 0;
        }
        else{
            if(dRand(0,1) > 0.5) bpr[i] = 1;
            else bpr[i] = 0;
        }
        
        if(bpr[i] == 1){
            n_open++;
            balance++;
        } else balance--;
    }
    
    sdsl::bp_support_g<> bpr_bps;
    sdsl::util::init_support(bpr_bps, &bpr);
    
    // Generate marks
    sdsl::bit_vector marks(2*N,0);
    for(int64_t i = 0; i < 2*N; i++){
        if(bpr[i] == 1 && dRand(0,1) < marked_prob){
            marks[i] = 1;
            marks[bpr_bps.find_close(i)] = 1;
        }
    }
    
    // Test
    test_all_queries(bpr,marks);
    
    cout << "Lowest marked ancestor " << N << " " << marked_prob << " test OK" << endl;
    
}

void LMA_testcase1(){
    //                        0 1 2 3 4 5 6 7 8 9 101112131415161718192021222324252627
    //                        ( ( ( ( ( ) ( ) ) ( ( ) ( ( ) ( ( ) ) ) ) ) ( ( ) ) ) )
    sdsl::bit_vector bpr =   {1,1,1,1,1,0,1,0,0,1,1,0,1,1,0,1,1,0,0,0,0,0,1,1,0,0,0,0};
    sdsl::bit_vector marks = {1,0,0,1,0,0,0,0,1,1,1,1,0,0,0,1,0,0,1,0,1,0,1,0,0,1,0,1};
    
    LMA_Support nmas(bpr,marks);
    assert(nmas.LMA(13) == 9); // Checked manually
    test_all_queries(bpr,marks);
}

void LMA_testcase2(){
    string bpr_string = "11011011011010001110111010010011010001110110100010100111011010001110110110100001000011101101000110101000011101011101101001001110110100010011010001110110100011010100001110111010010110100001010011101101101000011110100110110100001010000111101001110100111010011010000011110100111010011010000110100110100000";
    string marks_string = "10010010000000011100110000100110000111000000000000001110010000111100100000000110011100000000000000000000111000011001000010011100100001100100000010100100001110000001011100110000100000000110000111000000000000111100001000000000010000111111000011100001010000110000101111000000110000100000011100001100001111";
    
    sdsl::bit_vector bpr(bpr_string.size());
    sdsl::bit_vector marks(marks_string.size());
    
    for(int64_t i = 0; i < bpr.size(); i++){
        bpr[i] = bpr_string[i] == '1' ? 1 : 0;
        marks[i] = marks_string[i] == '1' ? 1 : 0;
    }
    
    LMA_Support nmas(bpr,marks);
    test_all_queries(bpr,marks);

}

void LMA_Support_Tests(){

    
    LMA_testcase1();
    LMA_testcase2();
    
    // Test random inputs
    srand(234789234);
    for(int64_t i = 0; i < 100; i++){
        test_random(10000, 0.01 * i);
    }
    cout << "Nearest marked ancestor tests OK" << endl;
}

#endif
