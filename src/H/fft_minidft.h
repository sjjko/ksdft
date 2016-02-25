#ifndef FFT_MINIDFT_H_INCLUDED
#define FFT_MINIDFT_H_INCLUDED

#include <armadillo>


using namespace std;
using namespace arma;

arma::cx_vec func_fftn(arma::cx_vec vec_Inp,arma::Col<int> S,const string forward_or_backward);
arma::cx_vec func_dothefft(arma::cx_vec vec_Inp,arma::Col<double> S,double Faktor,fftw_direction dir);

inline arma::cx_vec func_fftn(arma::cx_vec vec_Inp,arma::Col<double> S,const string forward_or_backward)
{
fftw_direction dir;
double Faktor;
arma::cx_vec output(vec_Inp);

if(forward_or_backward=="forward")
{
    dir=FFTW_FORWARD;
    Faktor=1;
}
else
{
    dir=FFTW_BACKWARD;
    Faktor=prod(S);
}

output=func_dothefft(vec_Inp,S,Faktor,dir);

return output;
}

inline arma::cx_vec func_dothefft(arma::cx_vec vec_Inp,arma::Col<double> S,double Faktor,fftw_direction dir)
{
fftwnd_plan plan;

//cout << "fft of " << S[0] << " " << S[1] << " " << S[2] << endl;

plan=fftw3d_create_plan(S[2],S[1],S[0],dir,FFTW_IN_PLACE);
//column-major ordering in Armadillo, raw-major in fftw
// drehe die Dimensionen um

//reshape und recast - Transformation in fftw Format
//http://blog.joey-dumont.ca/2013/03/interfacing-armadillo-and-fftw.html

//reshape matrix to cube!
//Input als Vektor
//kopiere diesen Vektor in Cube - dann reshape

cx_cube cube_Inp(vec_Inp.n_elem,1,1);
cube_Inp.slice(0)=vec_Inp;
reshape(cube_Inp,S[0],S[1],S[2]);
// = reshape(vec_Inp,S[0])

//Dimensionen in umgekehrter Reihenfolge fuer fftw!
//reshape(cube_Inp,S[2],S[1],S[0]);

fftw_complex* fftw_cube_in = reinterpret_cast<fftw_complex*> (cube_Inp.memptr());
fftw_complex out[(int)S[0]][(int)S[1]][(int)S[2]];
//Adressierung ueber Pointer
//fftwnd_one(plan,fftw_cube_in,NULL);

//cout << "now fftw" << endl;
fftwnd_one(plan,fftw_cube_in,&out[0][0][0]);

fftwnd_destroy_plan(plan);

//double *d data=Inp.memptr();

//wieder in Vektor transformieren!
reshape(cube_Inp,prod(S),1,1);
//dann umwandeln - nur ein Schnitt des umgewandelten Arrays enthaelt Daten
cx_vec output_vec(cube_Inp.slice(0).col(0));

//return Faktor*output_vec;
return output_vec;

}


#endif // FFT_MINIDFT_H_INCLUDED
