#ifndef INPUT_READING_HH
#define INPUT_READING_HH

#include <utility>
#include <string>
#include <vector>
#include <algorithm>
#include <fstream>

class Raw_file_stream{
private:
    
    std::ifstream file;
    
public:
    
    Raw_file_stream(string filename) : file(filename, ios::in | ios::binary) {
        if(file.bad()){
            cerr << "Error opening file " << filename << endl;
            exit(-1);
        }
    }
    
    bool getchar(char& c){
        if(file.eof()) return false; // End of stream
        file.read(&c,1); // Read 1 byte
        if(file.eof()) return false;
        if(c == '\n' || c == '\r')
            std::cerr << "Warning: file contains a newline character" << std::endl;
        return true;
    }
    
};

class Read_stream{
    
private:
    
    bool end_of_read;
    std::ifstream* file;
    
public:
    
    Read_stream(std::ifstream* file) : end_of_read(false), file(file) {
    
    }

    // Behaviour: Tries to read a byte to c. If eof or '>' was read,
    // return false and return false from here on. If c
    // is '\n' or '\r', read another byte recursively.
    bool getchar(char& c){
        if(end_of_read) return false;
        file->read(&c,1);
        if(file->eof()) end_of_read = true;
        else{
            // Read was successful
            if(c == '\n' || c == '\r') return getchar(c);
            if(c == '>') end_of_read = true;
        }
        return !end_of_read;
    }

};



class FASTA_reader{
    
private:
    
    std::ifstream file;
    
public:
    
    FASTA_reader(string filename) : file(filename, ios::in | ios::binary) {
        if(file.bad()){
            cerr << "Error opening file " << filename << endl;
            exit(-1);
        }
    }
    
    bool done(){
        return file.eof();
    }
    
    Read_stream get_next_query_stream(){
        string line;
        getline(file, line); // Discard the header
        Read_stream rs(&file);
        return rs;
    }
};

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
