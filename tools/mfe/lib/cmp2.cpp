int main(int argc, char** argv) {
    vector<string> lines_1, lines_2;
    
    lines_1 = file2str(argv[1]);
    lines_2 = file2str(argv[2]);
    
    if (lines_1[0] != lines_2[0]) {
        cout << "###Warning! the head line doesn't match\n";
    }
    
    j_line = 1;
    for (i_line = 1; i_line < lines_1.size(); ++i_line) {
        vector<string> fields = split(lines_1[i_line]);
        while (1) {
            ++j_line;
        }
        
    }
    return 0;
}
