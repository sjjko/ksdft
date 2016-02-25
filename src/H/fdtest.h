verbosity("fdtest: an energy test routine",2,__FILE__,__LINE__);


verbosity("fdtest: define the oscillator potential",2,__FILE__,__LINE__);
double omega=1;
arma::mat Vosci(dr.n_rows,1);
Vosci=0.5*pow(omega,2)*(dr%dr);

gpL.plotMatrix3Slices("fdtest_Vosci",Vosci,S);

verbosity("fdtest: the oscillator potential has dimension "+std::to_string(Vosci.n_rows),2,__FILE__,__LINE__);
verbosity("fdtest: define the dual potential Vdual",2,__FILE__,__LINE__);
//arma::mat VDualOsci=getVdual(Op,Vosci);

verbosity("fdtest: performed for "+std::to_string(Pa.number_of_wavefunctions)+" wavefunctions!",2,__FILE__,__LINE__);

arma::cx_mat Wfd(prodS,Pa.number_of_wavefunctions);
Wfd=arma::randu<arma::cx_mat>(prodS,Pa.number_of_wavefunctions);
verbosity("fdtest: compute the energy first",2,__FILE__,__LINE__);
double E0=getE(Op,Pa,Wfd,Vosci);
verbosity("fdtest: E="+std::to_string(E0),2,__FILE__,__LINE__);
verbosity("fdtest: now the gradient",2,__FILE__,__LINE__);
arma::cx_mat g0=getgrad(Op,Pa,Wfd,Vosci);

gpL.plotMatrix3Slices("fdtest_grad_0",real(g0.col(0)),S);
//gpL.plotMatrix3Slices("fdtest_grad_1",real(g0.col(1)),S);

//cout << as_scalar(arma::max(g0))<<endl;
//cout << as_scalar(arma::min(g0))<<endl;
//cout << as_scalar(arma::mean(g0))<<endl;

verbosity("fdtest: now the wf differential",2,__FILE__,__LINE__);
arma::cx_mat dWfd(prodS,Pa.number_of_wavefunctions);
dWfd=arma::randu<arma::cx_mat>(prodS,Pa.number_of_wavefunctions);

double delta;
int i=1;
int j=0;
Col<double> energies(11),dEnergies(11);
for(delta=pow(10,i),i=1,j=0;delta>pow(10,-10);i--,delta=pow(10,i),j++)
{
    std::cout<<"fdtest: j:" << j << " delta: " << delta << " i " << i << std::endl;
    dEnergies(j)=2.*real(arma::trace(g0.t()*delta*dWfd));
    arma::cx_mat Wfdt=Wfd+delta*dWfd;
    energies(j)=getE(Op,Pa,Wfdt,Vosci);
}
for(int i=0; i<energies.n_elem;i++)
{
    std::cout<<"fdtest: (E(j)-E0)/DE(j) for i "<<i<<" is "<< (energies(i)-E0)/dEnergies(i) << std::endl;
}

verbosity("fdtest: finished testing",2,__FILE__,__LINE__);
