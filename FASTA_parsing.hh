#ifndef FASTA_PARSING_HH
#define FASTA_PARSING_HH

#include <utility>
#include <string>
#include <vector>
#include <algorithm>
#include <fstream>

bool is_invalid_read(const std::pair<std::string, std::string>& read){
    // Check that the read has only valid characters
    const std::vector<char> alphabet = {'A','C','G','T'};
    for(char c : read.first){
        c = toupper(c);
        if(find(alphabet.begin(), alphabet.end(), c) == alphabet.end())
            return true;
    }
    return false;
}

// Vector of (read, header) pairs
std::vector<std::pair<std::string, std::string> > parse_FASTA(std::string filename){
    std::ifstream input(filename);
    if(!input.good()){
        cerr << "Error opening file " << filename << endl;
        exit(-1);
    }

    std::vector<std::pair<std::string,std::string> > reads;

    std::string line;
    while(std::getline(input,line)){
        while(line.size() > 0 && isspace(line.back()))
            line.pop_back(); // Trim trailing whitespace just in case

        if(line.size() == 0 && !input.eof()) continue; // Ignore empty lines, just in case

        if(line[0] == '>')
            reads.push_back({"", line}); // Start new read
        else 
            reads.back().first += line; // Append base pairs to read
    }

    //Erase invalid reads
    //reads.erase(std::remove_if(reads.begin(), reads.end(), is_invalid_read), reads.end());
    return reads;
}


#endif
