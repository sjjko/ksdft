  //(const operatorStruct Op,const paramStruct Pa,const arma::mat dr,const arma::cx_mat Sfi,const arma::mat Xi,arma::mat G2i)
  cout << "start assembling osci potential!" << endl;
  double omega=2.;
  arma::mat V=0.5*pow(omega,2)*(dr%dr);
  arma::mat Vdual=real(*Op.Jd*(*Op.O*(*Op.J*V)));