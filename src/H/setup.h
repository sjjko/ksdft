#ifndef SETUP_H_INCLUDED
#define SETUP_H_INCLUDED

#include "main.h"

using namespace arma;

//Col<double> S(3);
//S[0]=20;S[1]=25;S[2]=16;
//S[0]=20;S[1]=25;S[2]=30;
//Pa.S=S;

lString="$\\mathbf{S}: 3 \\times 1$ number of sample points along each lattice vector ";
latX.newLine(lString);

lString="$ \\Bigg( \\begin{array}{ccc}"+std::to_string(Pa.S[0])+" \& 0 \& 0 \\\\ ";
lString+=" 0 \& "+ std::to_string(Pa.S[1]) +" \& 0  \\\\ ";
lString+=" 0 \& 0 \& " + std::to_string(Pa.S[2]) + "\\end{array} \\Bigg) $ ";
latX.newLine(lString);

//std::stringstream buffer;
//buffer << S << std::endl;
//latX.newLine(buffer.str());

//Col<int> Si(3);
//Si[0]=(int)S[0];
//Si[1]=(int)S[1];
//Si[2]=(int)S[2];
//double prodS=prod(S);
Col<double> ms(Pa.prodS);
for(unsigned int i=0;i<ms.n_elem;i++)
{ ms[i]=i; }

Col<double> m1(Pa.prodS);
Col<double> m2(Pa.prodS);
Col<double> m3(Pa.prodS);

for(unsigned int i=0;i<ms.n_elem;i++)
{
m1[i]=std::fmod(ms(i),Pa.S[0]);
m2[i]=std::floor(std::fmod(ms(i)/Pa.S[0],Pa.S[1]));
m3[i]=std::floor(std::fmod(ms(i)/(Pa.S[0]*Pa.S[1]),Pa.S[2]));
}


Mat<double> M(Pa.prodS,3);
M.fill(0);
M.col(0)=m1;
M.col(1)=m2;
M.col(2)=m3;

Col<double> n1(Pa.prodS);
Col<double> n2(Pa.prodS);
Col<double> n3(Pa.prodS);


for(unsigned int j=0;j<m1.n_rows;j++)
    {
        //value = (expression) ? (if true) : (if false);
        n1[j] = (m1[j]<=(Pa.S[0]/2)) ? (m1[j]) : (-Pa.S[0]+m1[j]);
        n2[j] = (m2[j]<=(Pa.S[1]/2)) ? (m2[j]) : (-Pa.S[1]+m2[j]);
        n3[j] = (m3[j]<=(Pa.S[2]/2)) ? (m3[j]) : (-Pa.S[2]+m3[j]);

        }

Mat<double> N(Pa.prodS,3);
N.fill(0);
N.col(0)=n1;
N.col(1)=n2;
N.col(2)=n3;


//Mat<double> R(3,3);

//R.eye();
//R*=6;

lString=" $\\mathbf{R}: 3 \\times 3$ matrix, each of whose columns is a lattice vector $ R_{k},\\mathbf{R}=\[R_{1},R_{2}, R_{3}\] $ ";
latX.newLine(lString);
lString="$ \\Bigg( \\begin{array}{ccc}"+std::to_string(Pa.R[0,0])+" \& 0 \& 0 \\\\ ";
lString+=" 0 \& "+ std::to_string(Pa.R[1,1]) +" \& 0  \\\\ ";
lString+=" 0 \& 0 \& " + std::to_string(Pa.R[2,2]) + "\\end{array} \\Bigg) $ ";
latX.newLine(lString);

//cout << R << endl;

Mat<double> diagS = diagmat(Pa.S);
Mat<double> idiagS = inv((mat)diagS);

latX.emptyLine();
lString=" $\\mathbf{r}: S_{k} \\times 3$ matrix, matrix, each of whose rows contains the real-space coordinates of a sample point"
" standard order $\\mathbf{r} ≡ \\mathbf{M} Diag(\\vec(S))^{-1}\\mathbf{R}^{T}$  ";
latX.newLine(lString);

mat r = M*idiagS*Pa.R.t();

Mat<double> mat_center_of_cell(ones(Pa.prodS,3));

verbosity(Pa,"compute center of cell",2,__FILE__,__LINE__);

for(int j=0;j<Pa.prodS;j++)
    {
     mat_center_of_cell.row(j)=(arma::sum(Pa.R,1)/2).t();
    }

mat dr(arma::sqrt(arma::sum((r-mat_center_of_cell)%(r-mat_center_of_cell),1)));

Mat<double> G;

latX.emptyLine();
latX.newLine(" The G operator: ");
latX.newLine(" components of reciprocal lattice vector: ");
latX.newLine(" $\\Pi_{k} S_{k} \\times 3 matrix$ ");
latX.newLine(" $ \\mathbf{G} = 2 \\pi \\mathbf{N} \\mathbf{R}^{-1}$ ");

G=2*arma::datum::pi*N*inv(Pa.R);

latX.emptyLine();
latX.newLine(" The G2 operator: ");
latX.newLine(" components of square magnitude of the vector in the corresponding row of G: ");
latX.newLine(" $\\Pi_{k} S_{k} \\times 1 matrix$ ");
latX.newLine(" $ G2 = G^{2}$ ");

arma::Mat<double> G2;
G2=arma::sum(G%G,1);



for(int i=0;i<diagS.n_rows;i++)
{
double remainder = diagS(i,i)-2.*floor(diagS(i,i)/2.);
if(remainder!=0.0)
{
cout << "in dimension " << (i+1) << " we have an uneven number of gridpoints - error!" << endl;
return -1;
}
}

verbosity(Pa,"compute compressed G2",2,__FILE__,__LINE__);

mat eS=Pa.S/2+0.5;

verbosity(Pa,"Mreduced first",2,__FILE__,__LINE__);

mat Mreduced = abs(M-arma::ones(M.n_rows,1)*eS.t());

verbosity(Pa,"find edges",2,__FILE__,__LINE__);

uvec edges=find(any(Mreduced<1,1));

verbosity(Pa,"get minimal element",2,__FILE__,__LINE__);

//find maximum element for reducing elements in G2
double G2mx=min(G2(edges));

verbosity(Pa,"get active indices",2,__FILE__,__LINE__);

uvec activeIndices=find(G2<G2mx/4.);
Pa.activeIndices=activeIndices;
Pa.activeIndicesPtr=&activeIndices;

Pa.numberOfActiveIndices=activeIndices.n_elem;
//store compressed G2 in G2 compressed!

verbosity(Pa,"retrieve G2 compressed",2,__FILE__,__LINE__);

mat G2comp=G2(activeIndices);



verbosity(Pa,"now compute the structure factor",2,__FILE__,__LINE__);

verbosity(Pa,"Compute structure factor using G: "+std::to_string(G.n_rows)+" x "+std::to_string(G.n_cols),2,__FILE__,__LINE__);
verbosity(Pa,"Compute structure factor using X: "+std::to_string(X.n_rows)+" x "+std::to_string(X.n_cols),2,__FILE__,__LINE__);

arma::mat tmpmat(-G*X);

verbosity(Pa,"Now compute exponent of temporary matrix: "+std::to_string(tmpmat.n_rows)+" x "+std::to_string(tmpmat.n_cols),2,__FILE__,__LINE__);

arma::cx_mat exponent(tmpmat.n_rows,tmpmat.n_cols,fill::zeros);
exponent.set_imag(tmpmat);

latX.emptyLine();
latX.newLine(" The Sf operator: ");
latX.newLine(" structure factor associated with the wave vector of the corresponding column of G: ");
latX.newLine(" $ \\Pi_{k} S_{k} \\times 1$ column vector (one column matrix) ");
latX.newLine(" $ Sf(\\vec(G)) = \\Sigma_{I} \\exp(-i \\vec(G) \\vec(X)_{I}) $ ");

verbosity(Pa,"Sf is a sum over all atoms I of exp(G*X_{I})!",2,__FILE__,__LINE__);
verbosity(Pa,"Sf is a column vector Sx1!",2,__FILE__,__LINE__);

arma::cx_mat Sf(sum(exp(exponent),1));

verbosity(Pa,"setup finished",2,__FILE__,__LINE__);



#endif // SETUP_H_INCLUDED
