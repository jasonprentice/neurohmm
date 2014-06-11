//
//  paramsStruct.cpp
//  
//
//  Created by Jason Prentice on 6/11/14.
//
//


#include "paramsStruct.h"

using namespace std;

// paramsStruct
paramsStruct::paramsStruct() : nfields(0) {}

int paramsStruct::get_nfields() {return nfields;}

void paramsStruct::addField(string name, myMatrix<double>& value) {
    nfields++;
    fields[name] = &value;
    
    fieldNames.clear();
    for (map<string, myMatrix<double>* >::iterator it = fields.begin(); it!=fields.end(); ++it) {
        fieldNames.push_back((it->first).data());
    }
    
    return;
}

const char** paramsStruct::fieldNamesArray() {
    return fieldNames.data();
}

vector<double>* paramsStruct::getFieldData(int index) {
    return fields[fieldNames[index]]->data();
}

int paramsStruct::getFieldN(int index) {
    if (index >= nfields || index < 0) {
        cerr << "Index out of range." << endl;
        return 0;
    } else {
        return fields[fieldNames[index]]->get_N();
    }
}

int paramsStruct::getFieldM(int index) {
    if (index >= nfields || index < 0) {
        cerr << "Index out of range." << endl;
        return 0;
    } else {
        return fields[fieldNames[index]]->get_M();
    }
}

const char* paramsStruct::getFieldName(int index) {
    if (index >= nfields || index < 0) {
        cerr << "Index out of range." << endl;
        return nullptr;
    } else {
        return fieldNames[index];
    }
}

