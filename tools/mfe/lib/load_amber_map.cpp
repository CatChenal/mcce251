#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include "mcce.hpp"

using namespace std;

vector <Res> load_amber_map(const char* file_name);
vector <Res> load_amber_map(ifstream& dat_file);
map< string, vector<float> > load_amber_param(const char* file_name);
map< string, vector<float> > load_amber_param(ifstream& dat_file);
map < string, vector<string> > load_mcce2amber(const char* file_name);
map < string, vector<string> > load_mcce2amber(ifstream& dat_file);
vector <vector <string> > get_all_atoms(const char *dir_name);

vector<string> split_w_quote(const string& str);

int main(int argc, char **argv)
{
    if (argc < 4) {
        cout << "Usage: \n";
        cout << "   test_amber file.lib file.dat param_dir\n";
        return -1;
    }
    
    vector <Res> amber_map = load_amber_map(argv[1]);
    map< string, vector<float> > amber_param = load_amber_param(argv[2]);
    ofstream mcce2amber_file;
    ofstream tpl_file;

    for (int i_res =0; i_res<amber_map.size(); ++i_res) {
        for (int i_atom =0; i_atom<amber_map[i_res].conf[0].atom.size(); ++i_atom) {
            string atomType = amber_map[i_res].conf[0].atom[i_atom].atomType;
            if (amber_param.find(atomType) == amber_param.end()) {
                cout << "Warning: no vdw parameter for atom type " << atomType << endl;
                continue;
            }
            
            vector <float> vdw_param = amber_param[atomType];
            amber_map[i_res].conf[0].atom[i_atom].vdw_rad = vdw_param[0];
            amber_map[i_res].conf[0].atom[i_atom].vdw_eps = vdw_param[1];
            
            /*
            cout << "VDW_RAD   " << amber_map[i_res].resName << "  "
            << amber_map[i_res].conf[0].atom[i_atom].atomName << "  "
            << amber_map[i_res].conf[0].atom[i_atom].vdw_rad << endl;
            cout << "VDW_EPS   " << amber_map[i_res].resName << "  "
            << amber_map[i_res].conf[0].atom[i_atom].atomName << "  "
            << amber_map[i_res].conf[0].atom[i_atom].vdw_eps << endl;
            cout << "CHARGE    " << amber_map[i_res].resName << "  "
            << amber_map[i_res].conf[0].atom[i_atom].atomName << "  "
            << amber_map[i_res].conf[0].atom[i_atom].crg << endl;
            */
        }
    }
    
    map < string, vector<string> > mcce2amber;
    ifstream test_file("mcce2amber.txt");
    if (test_file) {
        test_file.close();
        mcce2amber = load_mcce2amber("mcce2amber.txt");
    }
    
    map < string, bool > exported_confs;
    
    vector <vector <string> > all_atoms = get_all_atoms(argv[3]);
    for (int i_a = 0; i_a < all_atoms.size(); ++i_a) {
        cout << all_atoms[i_a][0] << " " << all_atoms[i_a][1] << " " << endl;
    }
    for (int i_a = 0; i_a < all_atoms.size(); ++i_a) {
        string key = strip(all_atoms[i_a][0]) + "|" + strip(all_atoms[i_a][1]);
        
        if (mcce2amber.find(key) != mcce2amber.end()) {
            int k_res;
            for (k_res =0; k_res<amber_map.size(); ++k_res) {
                if (mcce2amber[key][0] == amber_map[k_res].resName) break;
            }
            if (k_res >= amber_map.size()) {
                mcce2amber.erase(key);
            }
            else {
                int k_atom;
                for (k_atom =0; k_atom<amber_map[k_res].conf[0].atom.size(); ++k_atom) {
                    if (mcce2amber[key][1] == amber_map[k_res].conf[0].atom[k_atom].atomName) break;
                }
                if (k_atom >= amber_map[k_res].conf[0].atom.size()) {
                    mcce2amber.erase(key);
                }
            }
        }
        while (mcce2amber.find(key) == mcce2amber.end()) {
            cout << "Can't find atom " << all_atoms[i_a][0] << " " << all_atoms[i_a][1] << " in the parameter list\n";
            cout << "Choose the closest residue from the list:\n";
            cout << "[-1] skip this\n";
            for (int j_res =0; j_res<amber_map.size(); ++j_res) {
                cout << "[" << j_res+1 << "] " << amber_map[j_res].resName << " ";
                if (fmod(j_res+1., 10.) < 1e-4) { cout << endl;}
            }
            cout << endl;
            
            int k_res;
            cin >> k_res;

            if (k_res == -1) {
                break;
            }
            else if (k_res >0 && k_res <= amber_map.size()) {
                k_res = k_res-1;
            }
            else {
                continue;
            }
            
            cout << "Choose the closest atom from the list:\n";
            cout << "[-1] skip this\n";
            for (int j_atom =0; j_atom<amber_map[k_res].conf[0].atom.size(); ++j_atom) {
                cout << "[" << j_atom+1 << "] " << amber_map[k_res].conf[0].atom[j_atom].atomName << " ";
                if (fmod(j_atom+1., 10.) < 1e-4) { cout << endl;}
            }
            cout << endl;
            
            int k_atom;
            cin >> k_atom;

            if (k_atom == -1) {
                break;
            }
            else if (k_atom >0 && k_atom <= amber_map[k_res].conf[0].atom.size()) {
                k_atom = k_atom-1;
                vector <string> amber_names(2);
                amber_names[0] = amber_map[k_res].resName;
                amber_names[1] = amber_map[k_res].conf[0].atom[k_atom].atomName;
                
                mcce2amber[key] = amber_names;
                
                mcce2amber_file.open("mcce2amber.txt", ios::out);
                for (int k_a = 0; k_a < all_atoms.size(); ++k_a) {
                    string key_print = strip(all_atoms[k_a][0]) + "|" + strip(all_atoms[k_a][1]);
                    if (mcce2amber.find(key_print) != mcce2amber.end()) {
                        mcce2amber_file << all_atoms[k_a][0] << " " << all_atoms[k_a][1] << " "
                        << mcce2amber[key_print][0] << " " << mcce2amber[key_print][1] << endl;
                    }
                }
                //mcce2amber_file << all_atoms[i_a][0] << " " << all_atoms[i_a][1] << " "
                //<< amber_names[0] << " " << amber_names[1] << endl;
                mcce2amber_file.close();

                tpl_file.open("amber.tpl", ios::out);
                for (int k_a = 0; k_a < all_atoms.size(); ++k_a) {
                    string key_print = strip(all_atoms[k_a][0]) + "|" + strip(all_atoms[k_a][1]);
                    if (mcce2amber.find(key_print) != mcce2amber.end()) {
                        int j_res;
                        for (j_res =0; j_res<amber_map.size(); ++j_res) {
                            cout << mcce2amber[key_print][0] << " " << amber_map[j_res].resName <<endl;
                            if (mcce2amber[key_print][0] == amber_map[j_res].resName) break;
                        }
                        if (j_res<amber_map.size()) {
                            int j_atom;
                            for (j_atom =0; j_atom<amber_map[j_res].conf[0].atom.size(); ++j_atom) {
                                if (mcce2amber[key_print][1] == amber_map[j_res].conf[0].atom[j_atom].atomName) break;
                            }
                            if (j_atom<amber_map[j_res].conf[0].atom.size()) {
                                tpl_file << "VDW_RAD   " << all_atoms[k_a][0] << " " << all_atoms[k_a][1] << " "
                                << amber_map[j_res].conf[0].atom[j_atom].vdw_rad << endl;
                                tpl_file << "VDW_EPS   " << all_atoms[k_a][0] << " " << all_atoms[k_a][1] << " "
                                << amber_map[j_res].conf[0].atom[j_atom].vdw_eps << endl;
                            }
                        }
                    }
                }
                //tpl_file << all_atoms[i_a][0] << " " << all_atoms[i_a][1] << " "
                //<< amber_names[0] << " " << amber_names[1] << endl;
                tpl_file.close();
            }
            else {
                continue;
            }
        }
    }
    /*
    if (exported_confs.find(prot.res[i_res].conf[i_conf].confName) != exported_confs.end())
        continue;
    exported_confs[prot.res[i_res].conf[i_conf].confName] = true;
    
    while (mcce2amber.find(prot.res[i_res].conf[i_conf].confName) == mcce2amber.end()) {
        cout << "Can't find conformer " << prot.res[i_res].conf[i_conf].confName << " in the parameter list\n";
        cout << "Choose the closest one from the list:\n";
        cout << "[-1] skip this conformer\n";
        for (int j_res =0; j_res<amber_map.size(); ++j_res) {
            cout << "[" << j_res+1 << "] " << amber_map[j_res].resName << " ";
            if (fmod(j_res+1., 10.) < 1e-4) { cout << endl;}
        }
        cout << endl;
        
        int k_res;
        cin >> k_res;
        
        if (k_res == -1) {
            break;
        }
        else if (k_res >=0 && k_res < amber_map.size()) {
            k_res = k_res-1;
            vector <string> amber_names(1);
            amber_names[0] = amber_map[k_res].resName;
            mcce2amber[prot.res[i_res].conf[i_conf].confName] = amber_names;
            
            ofstream mcce2amber_file;
            mcce2amber_file.open("mcce2amber.txt", ios::app);
            mcce2amber_file << prot.res[i_res].conf[i_conf].confName << "    * " << amber_map[k_res].resName << "    *\n";
            mcce2amber_file.close();
        }
        else {
            continue;
        }
    }
    
    for (int i_atom = 0; i_atom<prot.res[i_res].conf[i_conf].atom.size(); ++i_atom) {
        string key = prot.res[i_res].conf[i_conf].confName + "|" + prot.res[i_res].conf[i_conf].atom[i_atom].atomName;
        if (mcce2amber.find(prot.res[i_res].conf[i_conf].confName) == mcce2amber.end()) {
            
        }
        vector <string> value(2);
        value[0] = tokens[2];
        value[1] = tokens[3];
        
    }
        }
    }
    */
    tpl_file.open("amber.tpl", ios::out);
    for (int k_a = 0; k_a < all_atoms.size(); ++k_a) {
        string key_print = strip(all_atoms[k_a][0]) + "|" + strip(all_atoms[k_a][1]);
        if (mcce2amber.find(key_print) != mcce2amber.end()) {
            int j_res;
            for (j_res =0; j_res<amber_map.size(); ++j_res) {
                cout << mcce2amber[key_print][0] << " " << amber_map[j_res].resName <<endl;
                if (mcce2amber[key_print][0] == amber_map[j_res].resName) break;
            }
            if (j_res<amber_map.size()) {
                int j_atom;
                for (j_atom =0; j_atom<amber_map[j_res].conf[0].atom.size(); ++j_atom) {
                    if (mcce2amber[key_print][1] == amber_map[j_res].conf[0].atom[j_atom].atomName) break;
                }
                if (j_atom<amber_map[j_res].conf[0].atom.size()) {
                    tpl_file << "VDW_RAD   " << all_atoms[k_a][0] << " " << all_atoms[k_a][1] << " "
                    << amber_map[j_res].conf[0].atom[j_atom].vdw_rad << endl;
                    tpl_file << "VDW_EPS   " << all_atoms[k_a][0] << " " << all_atoms[k_a][1] << " "
                    << amber_map[j_res].conf[0].atom[j_atom].vdw_eps << endl;
                }
            }
        }
    }
    //tpl_file << all_atoms[i_a][0] << " " << all_atoms[i_a][1] << " "
    //<< amber_names[0] << " " << amber_names[1] << endl;
    tpl_file.close();
    return 0;
}

vector <Res> load_amber_map(const char* file_name)
{
    vector <Res> res_list;
    ifstream dat_file;
    string msg;
    
    dat_file.open(file_name, ios::in);
    if (!dat_file) {
        cerr << "Error: can't open file " << file_name << " for reading." << endl;
        throw msg = "read_file_error";
    }
    res_list = load_amber_map(dat_file);
    dat_file.close();
    return res_list;
}

vector <Res> load_amber_map(ifstream& dat_file)
{
    string line;
    vector <Res> res_list;
    int k_res = -1; /* telling the program which residue to store the atom:atomType map information */
    
    while (getline(dat_file, line)) {
        if (line.substr(0,2) == "!!") {
            /* index array */
            k_res = -1;
        }
        else if (line.substr(0,1) == "!") {
            /* starting a new section,
            all entries in this section belongs to one residue */
            if (line.find("atoms ") == string::npos) {
                k_res = -1;
                continue;
            }
            
            string resName = line.substr(line.find_first_of("."));
            resName = resName.substr(1);
            resName = resName.substr(0, resName.find_first_of("."));
            
            /* find the matching residue */
            for (k_res = res_list.size() -1; k_res >= 0; --k_res) {
                if (resName == res_list[k_res].resName) break;
            }
            
            /* If this residue is new, add one */
            if (k_res < 0) {
                res_list.resize(res_list.size()+1);
                k_res = res_list.size()-1;
                res_list[k_res].conf.resize(1);
            }
            res_list[k_res].resName = resName;
        }
        else {
            if (k_res == -1) continue;
            string atomName;
            string atomType;
            double crg;
            int k_atom;
            
            string resName = res_list[k_res].resName;
            
            int word_start_pos, word_end_pos, word_size;
            
            word_start_pos = line.find("\"");
            word_end_pos = line.find("\"", word_start_pos+1);
            word_size = word_end_pos - word_start_pos - 1;
            atomName = line.substr(word_start_pos+1, word_size);
            
            word_start_pos = line.find("\"", word_end_pos+1);
            word_end_pos = line.find("\"", word_start_pos+1);
            word_size = word_end_pos - word_start_pos - 1;
            atomType = line.substr(word_start_pos+1, word_size);
            
            stringstream line_ss(strip(line.substr(word_end_pos+1)));
            for (int i=0;i<6;++i) line_ss >> crg;
            
            for (k_atom = res_list[k_res].conf[0].atom.size() -1; k_atom >= 0; --k_atom) {
                if (res_list[k_res].conf[0].atom[k_atom].atomName == atomName) break;
            }
            if (k_atom < 0) {
                int n_atom = res_list[k_res].conf[0].atom.size() + 1;
                res_list[k_res].conf[0].atom.resize(n_atom);
                k_atom = n_atom -1;
                res_list[k_res].conf[0].atom[k_atom].atomName = atomName;
                res_list[k_res].conf[0].atom[k_atom].atomType = atomType;
                res_list[k_res].conf[0].atom[k_atom].crg    = crg;
            }
            
            /* 
            cout << res_list[k_res].resName << " " << atomName << " " << k_atom << " "
            << res_list[k_res].conf[0].atom[k_atom].atomName << " "
            << res_list[k_res].conf[0].atom[k_atom].atomType << " "
            << res_list[k_res].conf[0].atom[k_atom].crg << " "
            << endl;
            */
        }
    }
    return res_list;
}

map< string, vector<float> > load_amber_param(const char* file_name)
{
    map< string, vector<float> > vdw_param_map;
    ifstream dat_file;
    string msg;
    
    dat_file.open(file_name, ios::in);
    if (!dat_file) {
        cerr << "Error: can't open file " << file_name << " for reading." << endl;
        throw msg = "read_file_error";
    }
    vdw_param_map = load_amber_param(dat_file);
    dat_file.close();
    return vdw_param_map;
}

map< string, vector<float> > load_amber_param(ifstream& dat_file)
{
    vector <string> lines;
    map< string, vector<float> > vdw_param_map;
    int k_line;
    lines = file2str(dat_file);
    
    /* Jump to VDW section (looking for a line says "MOD4  RE" */
    for (k_line=0; k_line<lines.size(); ++k_line) {
        string line = lines[k_line];
        if (line.find("MOD4 ") != string::npos) {
            if (line.find("RE") != string::npos) {
                break;
            }
        }
    }
    
    /* read this section until reaching a blank line */
    for (int i_line=k_line+1; i_line<lines.size(); ++i_line) {
        string line = lines[i_line];
        if (line.size() == 0) break;
        
        stringstream line_ss(line);
        string atomType;
        vector <float> vdw_param;
        float r, epsilon;
        
        line_ss >> atomType;
        line_ss >> r;
        line_ss >> epsilon;
        
        vdw_param.resize(2);
        vdw_param[0] = r;
        vdw_param[1] = epsilon;
        
        vdw_param_map[atomType] = vdw_param;
    }
    
    /* read previous section for type alias information */
    for (int i_line=k_line-2; i_line>=0; --i_line) {
        string line = lines[i_line];
        if (line.size() == 0) break;
        
        vector <string> alias;
        alias = split(line);
        for (int i_alias = 1; i_alias < alias.size(); ++i_alias) {
            vdw_param_map[alias[i_alias]] = vdw_param_map[alias[0]];
        }
    }
    return vdw_param_map;
}

map < string, vector<string> > load_mcce2amber(const char* file_name)
{
    map < string, vector<string> > mcce2amber;
    ifstream dat_file;
    string msg;
    
    dat_file.open(file_name, ios::in);
    if (!dat_file) {
        cerr << "Error: can't open file " << file_name << " for reading." << endl;
        throw msg = "read_file_error";
    }
    mcce2amber = load_mcce2amber(dat_file);
    dat_file.close();
    return mcce2amber;

}

map < string, vector<string> > load_mcce2amber(ifstream& dat_file)
{
    map < string, vector<string> > mcce2amber;
    string line;
    
    while (getline(dat_file, line)) {
        line = rm_comment(line);
        vector <string> tokens = split_w_quote(line);
        if (tokens.size() <4) continue;
        
        if (tokens[1] == "*") {
            string key = tokens[0];
            vector <string> value(1); /* make a string vector with 1 element */
            value[0] = tokens[2];
            
            mcce2amber[key] = value;
        }
        else {
            string key = tokens[0] + "|" + tokens[1];
            vector <string> value(2); /* make a string vector with 2 elements */
            value[0] = tokens[2];
            value[1] = tokens[3];
            
            mcce2amber[key] = value;
        }
    }
    return mcce2amber;
}

vector<string> split_w_quote(const string& str)
{
    /* break a line where quote marks can be used to enclose one token especially
    when spaces occur in tokens */
    vector<string> tokens; // Create vector to hold the words
    bool quote_on = false;
    int last_break = -1;
    
    /* when reading character by character, each can be either a space, or a quote mark
    or a regular character. spaces will break the tokens unless it is enclosed in quote marks */
    for (int i_char=0; i_char<str.size(); ++i_char) {
        if (str[i_char] == ' ') {
            if (quote_on) continue;
            
            if ((i_char-last_break) >1) {
                tokens.push_back(str.substr(last_break+1, i_char-last_break-1));
            }
            last_break = i_char;
        }
        else if (str[i_char] == '"') {
            quote_on = not quote_on;
            if ((i_char-last_break) >1 or not quote_on) {
                tokens.push_back(str.substr(last_break+1, i_char-last_break-1));
            }
            last_break = i_char;
        }
        else {
            continue;
        }
    }
    if (last_break < str.size()-1) {
        tokens.push_back(str.substr(last_break+1, str.size()-last_break-1));
    }
    
    return tokens;
}

vector <vector <string> > get_all_atoms(const char *dir_name)
{
    vector <vector <string> > atoms;
    vector <string> file_names = get_files(dir_name);
    for (int i_file = 0; i_file < file_names.size(); ++i_file) {
        if (file_names[i_file].size() < 5) continue;
        if (file_names[i_file].substr(file_names[i_file].size()-4,4) != ".tpl") continue;
        
        string file_path = dir_name;
        file_path += '/';
        file_path += file_names[i_file];
        
        vector <string> lines = file2str(file_path.c_str());
        for (int i_line=0; i_line<lines.size(); ++i_line) {
            if (lines[i_line].size() < 20) continue;
            if (lines[i_line].substr(0, 7) != "CONNECT") continue;
            
            vector <string> params(2);
            params[0] = lines[i_line].substr(9, 5);
            params[1] = lines[i_line].substr(15, 4);
            
            atoms.push_back(params);
        }
    }
    return atoms;
}

/* amber_map
string mcce_sch_str = prot.res[i_res].conf[i_conf].confName + "|" + prot.res[i_res].conf[i_conf].atom[i_atom].atomName;
while (mcce2amber.find(mcce_sch_str) == mcce2amber.end()) {
    if (mcce2amber.find(prot.res[i_res].conf[i_conf].confName) != mcce2amber.end()) {
        for (int j_atom =0; i_atom<amber_map[i_res].conf[0].atom.size(); ++i_atom) {
            if () {
                
            }
        }
        else {
            break;
        }
    }
    
    */
