#include <cstdlib>
#include <cstring>
#include <vector>
#include <cmath>
#include "mcce.hpp"
#include "math3d.hpp"

using namespace std;
using namespace Math3D;

vector<float> sum_2vec(vector<float> vec1, vector<float> vec2)
{
    /* return summation of 2 vectors.
    if two vectors have different sizes, then returned vector has the smaller size */

    int sum_size = vec1.size() <= vec2.size()? vec1.size():vec2.size();
    vector <float> ret_vec(sum_size);
    for (int i=0; i< sum_size; ++i) {
        ret_vec[i] = vec1[i]+vec2[i];
    }
    return ret_vec;
}

vector<float> sum_vecs(vector< vector<float> > vecs)
{
    /* return summation of a list of vectors.
    if vectors have different sizes, then returned vector has the smallest size */

    if (vecs.size() == 0) {
        vector <float> tmp_vec;
        return tmp_vec;
    }

    int sum_size = vecs[0].size();
    for (int i_vec=1;i_vec<vecs.size();++i_vec) {
        if (sum_size > vecs[i_vec].size()) {
            sum_size = vecs[i_vec].size();
        }
    }

    vector<float> ret_vec(sum_size);
    for (int i=0; i< sum_size; ++i) {
        for (int i_vec=0;i_vec<vecs.size();++i_vec) {
            ret_vec[i] += vecs[i_vec][i];
        }
    }
    return ret_vec;
}

void print_vec(const vector<float> vec, const string fmt)
{
    for (int i=0; i<vec.size(); ++i) {
        printf(fmt.c_str(), vec[i]);
        if (i < vec.size()-1) printf(" ");
    }
    printf("\n");
}

string vec2string(const vector<float> vec, const string fmt)
{
    char* c_str;
    string ret_str, tmp_str;
    int len = atoi(fmt.substr(1,1).c_str())+1;
    c_str = (char*)malloc(len * sizeof(char));

    for (int i=0; i<vec.size(); ++i) {
        sprintf(c_str, fmt.c_str(), vec[i]);
        tmp_str = c_str;
        ret_str += tmp_str;
        if (i < vec.size()-1) ret_str += " ";
    }

    free(c_str);
    return ret_str;
}

string vec2string(const vector<int> vec, const string fmt)
{
    char* c_str;
    string ret_str, tmp_str;
    int len = atoi(fmt.substr(1,1).c_str())+1;
    c_str = (char*)malloc(len * sizeof(char));

    for (int i=0; i<vec.size(); ++i) {
        sprintf(c_str, fmt.c_str(), vec[i]);
        tmp_str = c_str;
        ret_str += tmp_str;
        if (i < vec.size()-1) ret_str += " ";
    }

    free(c_str);
    return ret_str;
}

vector <float> get_stats(vector<float> list, vector<float> weight, float bin_size) {
    vector <float> stats;
    if (list.size() != weight.size()) return stats;

    int n_stat;
    float min_value = min(list);
    float max_value = max(list);

    n_stat = (int) (ceil(max_value/bin_size) - floor(min_value/bin_size));
    stats.resize(n_stat*2);

    for (int i=0; i<n_stat; ++i) {
        stats[i] = bin_size * (floor(min_value/bin_size) + (float)i + 0.5);
        stats[n_stat+i] = 0;
    }

    for (int k=0; k<list.size(); ++k) {
        stats[n_stat + (int) (floor(list[k]/bin_size) - floor(min_value/bin_size))] += weight[k];
    }

    return stats;
}

vector <float> get_stats(vector<float> list, float bin_size) {
    vector<float>stats;
    vector<float>weight;
    weight.resize(list.size());
    for (int i=0; i<weight.size();++i) {
        weight[i] = 1.;
    }
    stats = get_stats(list, weight, bin_size);
    return stats;
}

vector < vector <int> > get_bins(vector<int> list, vector<float> values, float bin_size) {
    return get_bins(list, values, bin_size, 0.);
}

vector < vector <int> > get_bins(vector<int> list, vector<float> values, float bin_size, float zero_point) {

    vector < vector <int> > bins;
    if (list.size() != values.size()) return bins;

    int n_stat;
    float min_value = min(values);
    float max_value = max(values);

    n_stat = (int) (ceil((max_value-zero_point)/bin_size) - floor((min_value-zero_point)/bin_size));
    bins.resize(n_stat);

    for (int k=0; k<values.size(); ++k) {
        int k_bin = (int) (floor((values[k]-zero_point)/bin_size) - floor((min_value-zero_point)/bin_size));
        bins[k_bin].push_back(list[k]);
    }
    return bins;
}

float get_dihedral(Vector4<float> v0, Vector4<float> v1, Vector4<float> v2, Vector4<float> v3) {
    /* a simpler version based on the formula from http://en.wikipedia.org/wiki/Dihedral_angle, tested to give the same result */
    float angle;
    Vector4<float> b1 = v1 - v0;
    Vector4<float> b2 = v2 - v1;
    Vector4<float> b3 = v3 - v2;

    angle = atan2(b2.Length3() * b1.DotProduct3(CrossProduct(b2,b3)), CrossProduct(b1,b2).DotProduct3(CrossProduct(b2,b3)));
    return angle;
}

