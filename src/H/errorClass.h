#ifndef INCLUDE_ERRORCLASS
#define INCLUDE_ERRORCLASS

#include "registryClass.H"
//#include "../../main.H"
//class registryClass;


class errorClass
{

//private:

  int _error_identifier;
  registryClass *_pointer_Registry;

public:

  errorClass(registryClass *Ri)
  {
      this->_pointer_Registry=Ri;
  }
  ~errorClass(){};
  inline void assert(bool criterion, std::string text) {if(criterion==false)
  {
      std::cout << " Error in routine " << this->_pointer_Registry->actual_position() << " error message: " << text << std::endl;
      }
      }

};


#endif
