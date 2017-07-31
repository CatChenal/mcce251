#include <sys/types.h>
#include <ctype.h>
#include <dirent.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

using namespace std;

vector <string> get_files(const char *dir_name)
{
    vector <string> files;
    DIR     *dir(opendir(dir_name));;              /* DIR returned by opendir() */
    
    if (dir != 0) {
        struct dirent *dirent;
        while(( dirent = readdir(dir)) !=0) {
            string name(dirent->d_name);
            files.push_back(name);
        }
        closedir(dir);
    }
    else {
        cerr << "Error opening directory \"" << dir_name << "\"" << endl;
    }
    
    return files;
}

