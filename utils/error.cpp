//
// Created by m4zz31 on 26/10/21.
//

#include <stdlib.h>
#include <iostream>
#include <string>

using namespace std;

__attribute__((unused)) [[ noreturn ]] void error_report(string s){
    cout << s << endl;
    exit(EXIT_FAILURE);
}