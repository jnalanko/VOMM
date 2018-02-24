#include "bwt.hh"
#include <cstring>
#include <iostream>
#include <fstream>
#include <string>

extern "C" { 
#include "dbwt.h"
}

using namespace std;

uint8_t* bwt_dbwt(uint8_t* text, int64_t length, uint8_t end_char){

    unsigned int last;
    int64_t n = length;
    uint8_t* d = dbwt_bwt(text, n, &last, 0);
    d = (uint8_t*)realloc(d, (n + 2) * sizeof(uint8_t));
    d[last] = end_char;
    d[n + 1] = 0;

    return d;
}
