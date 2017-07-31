#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include "mcce.hpp"

using namespace std;

int main(int argc, char *argv[]) {
    if (argc < 3) {
        cout << "conv2o in.pdb out.pdb\n" ;
    }
    
    vector <string> lines = file2str(argv[1]);
    for (int i_line = 0; i_line<lines.size(); i_line++) {
        if (lines[i_line].length() > 83) {
            lines[i_line][82] = 'O';
        }
    }
    
    ofstream output;
    output.open(argv[2], ios::out);
    for (int i_line = 0; i_line<lines.size(); i_line++) {
        output << lines[i_line] << endl;
    }
    
    return 0;
        
}
