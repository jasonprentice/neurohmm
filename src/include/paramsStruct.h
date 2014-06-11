//
//  paramsStruct.h
//  
//
//  Created by Jason Prentice on 6/11/14.
//
//

#ifndef ____paramsStruct__
#define ____paramsStruct__

#include <vector>
#include <map>

#include "myMatrix.h"

// *********************** paramsStruct ****************************
class paramsStruct {            // Variable-sized structure holding double matrices
public:
    paramsStruct();
    int                     get_nfields();
    void                    addField( std::string, myMatrix<double>& );
    const char**            fieldNamesArray();
    std::vector<double>*         getFieldData( int );
    int                     getFieldN( int );
    int                     getFieldM( int );
    const char*             getFieldName( int );
    
private:
    int                     nfields;
    std::map<std::string,
    myMatrix<double>* > fields;
    std::vector<const char*> fieldNames;
    
};
// *****************************************************************
#endif /* defined(____paramsStruct__) */
