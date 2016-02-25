#ifndef POISSON_H_INCLUDED
#define POISSON_H_INCLUDED

latX.subSection(" solve the poisson equation ");

/*
Mat<double> mat_center_of_cell(ones(prodS,3));
verbosity("compute vec_center_of_cell2");
for(int j=0;j<prodS;j++)
    {
     mat_center_of_cell.row(j)=(arma::sum(R,1)/2).t();
    }*/

//vec_center_of_cell.col(1)=(arma::sum(R,0)/2)[1];
//cout << arma::sum((r-mat_center_of_cell)%(r-mat_center_of_cell),1) << endl;

//verbosity("compute dr");

latX.newLine(" $dr = \\sqrt{ \\Sigma_{i} (r-\\langle r \\rangle)^{2} $ }");


Mat<double> dr(arma::sqrt(arma::sum((r-mat_center_of_cell)%(r-mat_center_of_cell),1)));

cout << "sizedr" << sum(dr) << endl;

const double sigma1=0.75;
const double sigma2=0.5;

verbosity("compute g1,g2",2,__FILE__,__LINE__);

Col<double> g1=arma::exp(-(dr%dr)/(2.*pow(sigma1,2)))/sqrt(pow(2.*arma::datum::pi*sigma1*sigma1,3));

cout << sum(g1);

Col<double> g2=arma::exp(-(dr%dr)/(2.*pow(sigma2,2)))/sqrt(pow(2.*arma::datum::pi*sigma2*sigma2,3));

verbosity("compute n",2,__FILE__,__LINE__);

arma::cx_mat n(g1.n_elem,Pa.number_of_wavefunctions);
n.fill(0);
for(int i=0;i<n.n_cols;i++){n.col(i)=cx_vec(g2-g1,arma::zeros<arma::vec>(g2.n_elem));} //arma::fill::zeros);}
//for(int i=0;i<n.n_cols;i++){n.col(i).set_real(g1-g2);}

//%# Check norms and integral (should be near 1 and 0, respectively)
cout << "Normalization check on g1:" << sum(g1)*det(R)/prod(S) << endl;
cout << "Normalization check on g2:" << sum(g2)*det(R)/prod(S) << endl;
cout << "Total charge check:" << sum(n)*det(R)/prod(S) << endl;
verbosity("Normalization check on g1:",2,__FILE__,__LINE__);

arma::cx_mat phi=*Op.I*(*Op.Li*(-4.*arma::datum::pi*(*Op.O*(*Op.J*n))));

cout << arma::accu(phi)<< endl;
arma::cx_mat intermediateM(-4.*arma::datum::pi*(*Op.O*(*Op.J*n)));
cout << arma::accu(intermediateM.rows(0,12)) << endl;
//exit(-1);
/*arma::cx_mat toaccu((J*phi).t()%(Oo*(J*n)));
std::complex<double> cxvart=arma::accu(toaccu);
double num=real(cxvart);*/
double Unum=0.5*real(arma::accu((*Op.J*phi)%(*Op.O*(*Op.J*n))));
double Uanal=((1./sigma1+1./sigma2)/2.-std::sqrt(2.)/sqrt(pow(sigma1,2)+pow(sigma2,2)))/sqrt(arma::datum::pi);
//fprintf("Numeric, analytic Coulomb energy: %20.16f,%20.16f\n",Unum,Uanal);
std::cout << "numeric coulomb energy: " << Unum << " analytic coulomb energy " << Uanal << endl;

#endif // POISSON_H_INCLUDED
