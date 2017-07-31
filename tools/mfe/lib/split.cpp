#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

vector<string> split(const string& str)
{
    string buf; // Have a buffer string
    stringstream ss(str); // Insert the string into a stream

    vector<string> tokens; // Create vector to hold the words

    while (ss >> buf) tokens.push_back(buf);
    
    return tokens;
}

vector<string> split(const string& str, const char* delimiters)
{
    string buf = str; /* buffer */
    vector<string> tokens; // Create vector to hold the words

    while (buf.size()) {
        string token = buf.substr(0, buf.find_first_of(delimiters));
        
        tokens.push_back(token);
        if (buf.size() > token.size()) {
                buf = buf.substr(token.length()+1);
                if (buf.size() ==0) {
                    token = "";
                    tokens.push_back(token);
                }
        }
        else {
            buf.clear();
        }
    }
    return tokens;
}

string rm_comment(const string& str)
{
    string ret_val = str;
    ret_val = ret_val.substr(0, ret_val.find("REMARK"));    /* return the string before "REMARK */
    ret_val = ret_val.substr(0, ret_val.find_first_of("#\n"));
    
    return ret_val;
}

vector<string> file2str(ifstream& infile)
{
    vector<string> str_vec;
    string line;
    
    while(getline(infile, line)) {
        str_vec.push_back(line);
    }
    return str_vec;
}

vector<string> file2str(const char* file_name)
{
    ifstream infile;
    vector<string> str_vec;
    string msg;

    infile.open(file_name, ios::in);
    if (!infile) {
        cerr << "Error: can't open file " << file_name << " for reading." << endl;
        throw msg = "read_file_error";
    }
    str_vec = file2str(infile);
    infile.close();
    
    return str_vec;
}

bool equ_str(const string& str1, const string& str2)
{
    if (str1.size() != str2.size()) return false;
    int i;
    for (i=0; i<str1.size(); ++i) {
        if (str1[i] == '.') continue;
        if (str2[i] == '.') continue;
        if (str1[i] == str2[i]) continue;
        
        break;
    }
    if (i<str1.size()) return false;
    else return true;
}

string strip(const string& str)
{
    int i_start, i_end;
    for (i_start=0; i_start<str.size(); i_start++) {
        if (str[i_start] != ' ' && str[i_start] != '\n') {
            break;
        }
    }
    for (i_end=str.size()-1; i_end>=0; i_end--) {
        if (str[i_end] != ' ' && str[i_end] != '\n') {
            break;
        }
    }
    return str.substr(i_start, i_end-i_start+1);
}

string strip_tail(const string& str)
{
    int i_end;
    i_end = str.find_last_not_of(" \n");
    return str.substr(0,i_end+1);
}

string add_under(string str_val)
{
    string return_val = strip(str_val);
    for (int i=0; i<return_val.size(); ++i) {
        if (return_val[i] == ' ') {
            return_val[i] = '_';
        }
    }
    return return_val;
}

string rm_under(string str_val)
{
    string return_val = strip(str_val);
    for (int i=0; i<return_val.size(); ++i) {
        if (return_val[i] == '_') {
            return_val[i] = ' ';
        }
    }
    return return_val;
}

string upper(const string str_in) {
    string str_out;
    str_out = str_in;
    for (int i_char = 0; i_char <str_in.length(); i_char++) {
        if (str_in[i_char] >= 'a' && str_in[i_char] <= 'z') {
            str_out[i_char] = str_in[i_char] - 'a' + 'A';
        }
    }
    return str_out;
}
