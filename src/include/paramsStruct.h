//
//  paramsStruct.h
//  
//
//  Created by Jason Prentice on 6/11/14.
//
//

#ifndef ____paramsStruct__
#define ____paramsStruct__

#include <set>
#include <string>
#include <map>


class paramsStruct {
public:
    typedef std::set<std::string>::const_iterator key_iter_t;
    
    template <typename T>
    void             set( std::string, T );
    template <typename T>
    T&               get( std::string );
    template <typename T>
    const T&         get( std::string ) const;
    void             erase( std::string );
    int              get_nfields() const;
    
    key_iter_t       begin() const;
    key_iter_t       end() const;
    
private:
    typedef std::set<std::string> keylist_t;
    
    struct ElemBase {
        virtual ~ElemBase() {};
        virtual ElemBase* clone() const = 0;

    };
    template <typename T>
    struct Elem : public ElemBase {
        T item;
        
        Elem (T _item) : item (_item) {};
        inline ElemBase* clone() const {
            return new Elem<T> (item); };
    };
    class ElemPr {
    public:
        ElemPr() = default;
        ElemPr( const ElemPr& other ) : pr (other.pr->clone()) {};
        ElemPr( ElemPr&& other ) : ElemPr() { swap(*this,other); };
        ElemPr& operator=( ElemPr other ) {
            swap(*this,other);
            return *this;
        };
        ~ElemPr() { delete pr; };
        friend void swap(ElemPr& one, ElemPr& two ) {
            std::swap(one.pr, two.pr);
            return;
        };
        
        template <typename T>
        static inline ElemPr make( T _item ) {
            return ElemPr( new Elem<T> (_item) ); };
        
        template <typename T>
        inline T& get() { return (dynamic_cast<Elem<T>* >(pr))->item; };

        template <typename T>
        inline const T& get() const { return (dynamic_cast<Elem<T>*>(pr))->item; };
        
    private:
        ElemBase* pr = nullptr;
        
        ElemPr( ElemBase* _pr ) : pr(_pr) {};
    };
    
    int                               nfields = 0;
    std::map<std::string, ElemPr>     fields;
    keylist_t                         field_names;
    
};

template <typename T>
inline T& paramsStruct::get( std::string name ) {
    return fields.at(name).get<T>();
}

template <typename T>
inline const T& paramsStruct::get( std::string name ) const {
    return fields.at(name).get<T>();
}

inline int paramsStruct::get_nfields() const {
    return nfields;
}

inline paramsStruct::key_iter_t paramsStruct::begin() const {
    return field_names.begin();
}
inline paramsStruct::key_iter_t paramsStruct::end() const {
    return field_names.end();
}

template <typename T>
void paramsStruct::set( std::string name, T item) {
    ElemPr item_pr = ElemPr::make(item);
    auto ins = fields.insert(std::make_pair(name, item_pr));
    if (ins.second) {
        nfields++;
        field_names.insert(name);
    } else {
        fields[name] = item_pr;
    }
    return;
}

void paramsStruct::erase( std::string name ) {
    fields.erase(name);
    field_names.erase(name);
    nfields--;
    return;
}


#endif /* defined(____paramsStruct__) */
