#ifndef LOGGING_HH
#define LOGGING_HH

#include <iostream>

string getTimeString(){
    std::time_t result = std::time(NULL);
    string time = std::asctime(std::localtime(&result));
    return time.substr(0,time.size() - 1); // Trim the trailing newline
}

static bool logging_enabled = true;

void enable_logging(){
    logging_enabled = true;
}

void disable_logging(){
    logging_enabled = false;
}

void write_log(string message){
    if(logging_enabled){
        std::cerr << getTimeString() << " " << message << std::endl;
    }
}

#endif
