#include "mcce.hpp"
#include <string>
#include <iostream>

int main() {
    string uniqID("HAN+1A0001_001");
    Conf conf;
    conf.initialize_by_uniqID(uniqID);
    Res res(conf);
    cout << res;
    return 0;
}
