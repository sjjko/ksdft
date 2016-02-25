#include "main.H"
#include "registry.H"

void registry::actual_position_in_code()
  {
      std::cout << "actual position in code:"
      std::cout << list_of_routines << std::endl;
      std::cout << "currently in line " << __LINE__ << " in file " << __FILE__ << std::endl;
  }
void registry::register_variable(data_types *pointer_t)
  {
     it = std::map<int,std::pair<int,data_types>>::iterator=variables_allocated.end();

     it->first=size(*pointer_to_variable)*size(type(*pointer_to_variable));

  }
int registry::return_total_allocated_space()
  {
      int gesamtgroesse;
      it = std::map<int,std::pair<int,data_types>>::iterator=variables_allocated.begin();

      gesamtgroesse=0;
      for(it=variables_allocated.begin();it=variables_allocated.end();it++)
      {
      gesamtgroesse+=it->first;
      }
      gesamtgroesse*=bytes_to_megabytes;

      std::cout << "gesamt wurden " << gesamtgroesse " MB allokiert " << std::endl;
      return gesamtgroesse;

  }
void registry::register(const int id, const std::string name)
{
  registered_routines.insert(std::pair<int,std::string>(id,name);
}

std::map<int,std::string> registry::write_routines()
{
    return registered_routines;
}

};
