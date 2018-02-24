#ifndef LOGGING_HH
#define LOGGING_HH

string getTimeString(){
    std::time_t result = std::time(NULL);
    string time = std::asctime(std::localtime(&result));
    return time.substr(0,time.size() - 1); // Trim the trailing newline
}

void write_log(string message){
    cerr << getTimeString() << " " << message << endl;
}

#endif
