#ifndef INCLUDE_registryClass
#define INCLUDE_registryClass

#include "../../main.H"


class registryClass
{

  const double bytes_to_megabytes=8;
//  enum data_types(int,double,std::string);
 int _ids;
protected:

// eine ID fuer jede Routine, gezaehlt nach Reihenfolge der Registrierung
// retrieve ID for name given
  std::map<int,std::string> _registered_routines;
  std::list<std::pair<int,std::string> > _list_of_routines;
//  struct struct_vars_allocated
//  {
      //pointer to allocated variable
      // want to retrieve pointer to variable U using
      // Reg.find("U")
//      std::map<int,std::pair<int,T> > variables_allocated;
//  };
////
public:

  registryClass(){this->_ids=0;};
  ~registryClass();
//  inline void registrate(const std::string name,std::unique_ptr<T> dptr) {this->_ids++;this->_registered_routines.insert(std::pair<int,std::string>(this->_ids,name));this->_fetch_information(dpointer)};
  inline std::map<int,std::string> write_routines() { return _registered_routines; };
  inline std::string actual_position() {return "not yet implemented!";};

  //  void push_list_of_routines(int routine_id,bool check_in)
//  {
//    //element of the routine which should be checked in
//    std::map<char,int>::iterator it=registered_routines.find(routine_id);
//
//    if(check_in==true)
//     {
//    //register a new routine
//     list_of_routines.push_back(it->first,it->second);
//     }
//    else
//     {
//    //delete routine from stack
//    list_of_routines.pop_back();
//    }
//
//  }
//  void actual_position_in_code()
//  {
//      std::cout << list_of_routines << std::endl;
//  }
//  void register_variable(data_types *pointer_t)
//  {
//     it = std::map<int,std::pair<int,data_types> >::iterator=variables_allocated.end();
//
//     it->first=size(*pointer_to_variable)*size(type(*pointer_to_variable));
//
//  }
//
//  inline int return_total_allocated_space()
//  {
//      int gesamtgroesse;
//      it = std::map<int,std::pair<int,data_types>>::iterator=variables_allocated.begin();
//
//      gesamtgroesse=0;
//      for(it=variables_allocated.begin();it=variables_allocated.end();it++)
//      {
//      gesamtgroesse+=it->first;
//      }
//      gesamtgroesse*=bytes_to_megabytes;
//
//      std::cout << "gesamt wurden " << gesamtgroesse " MB allokiert " << std::endl;
//
//  }
//
};


#endif
