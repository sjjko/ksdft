

latX.subSection("Solve the Schroedinger equation using an oscillator potential");
latX.subsubSection("define an oscillator potential");
latX.newLine(" $V(\\mathbf{r})={1 \\over 2} \\omega^{2} d\\mathbf{r}^{2}$ ");
verbosity("schroedinger: solve the schroedinger equation",2,__FILE__,__LINE__);
verbosity("schroedinger: define the oscillator potential",2,__FILE__,__LINE__);
#define SCHROEDINGER_PREPROCESSING
#ifdef SCHROEDINGER_PREPROCESSING
double omega=2.;
arma::mat Vosci(Pa.prodS,1);
Vosci=0.5*pow(omega,2)*(dr%dr);

string imgName=gpL.plotAMatrixSlice("schroedinger_Vosci",Vosci,Pa.S,0);
latX.newLine("A parabolic potential is assumed with the following shape:");
latX.insertImage(imgName,"parabolic ion potential for which Schroedinger equation is solved.",Pa.caseName);

verbosity("schroedinger: the oscillator potential has dimension "+std::to_string(Vosci.n_rows),2,__FILE__,__LINE__);
verbosity("schroedinger: initialize wavefunction W",2,__FILE__,__LINE__);
latX.subsubSection("Preparation of the wavefunction");
//arma::cx_mat Wm(prodS,Pa.number_of_wavefunctions);
latX.newLine(" Initialize wavefunctions randomĺy : ");
latX.newLine(" $W_{m}:$ complex matrix with $ \\Sigma S_{k} \\times N_{s} $ elements in r space ");

//Wm=arma::randu<arma::cx_mat>(prodS,Pa.number_of_wavefunctions);
//std::shared_ptr<arma::cx_mat> W(&Wm);
std::shared_ptr<arma::cx_mat> W(new arma::cx_mat(arma::randn<arma::cx_mat>(Pa.prodS,Pa.number_of_wavefunctions)));

#endif
#define SCHROEDINGER_DO_SD
#ifdef SCHROEDINGER_DO_SD
verbosity("schroedinger: create optimizer sd",2,__FILE__,__LINE__);
sdOptimizer sd(Op,Pa,Vosci,&latX);

sd.myLatexClass->commentMyFunction();
verbosity("schroedinger: prepare sd optimizer",2,__FILE__,__LINE__);
sd.setup(W);
verbosity("schroedinger: compute density and plot",2,__FILE__,__LINE__);

//check if W is normalized W=Y=W*U^⁻1/2
verbosity("schroedinger: check if W is normalized W=Y=W*U^⁻1/2",2,__FILE__,__LINE__);

cx_mat Ui=Uinvers(Op,Pa,*W.get());
mat density(computeDensityFromWavefuncs(Op,Pa,*W.get(),Ui));
gpL.plotMatrix3Slices("schroedinger_initialDensity",density,Pa.S);

latX.newLine("The following image shows the initial density computed using the initial state of random wavefunction.");

latX.insertImage("schroedinger_initialDensity_xy","density after initialization of the wavefunctions",Pa.caseName);

sd.myLatexClass->commentMyFunction();
verbosity("====================================",2,__FILE__,__LINE__);
verbosity("schroedinger: optimize using sd",2,__FILE__,__LINE__);
verbosity("====================================",2,__FILE__,__LINE__);
//int Nstep=250;
sd.optimize(W);

pccgOptimizer pcg(Op,Pa,Vosci,G2,&latX);
pcg.orthogonalizeWfunc(W);
pcg.optimize(W);

sd.myLatexClass->commentMyFunction();
verbosity("now get the wavefunction",2,__FILE__,__LINE__);
// Psi as a pointer to the wavefunction
std::shared_ptr<arma::cx_mat> Psi(new arma::cx_mat(*W));
// Epsilon points to the eigenvalues (real vector)
//std::shared_ptr<arma::cx_vec> Epsilon(new arma::cx_vec(W->n_rows));
std::shared_ptr<arma::mat> Epsilon(new arma::mat(Pa.number_of_wavefunctions,1));

#endif // SCHROEDINGER_DO_SD
#define SCHROEDINGER_CALC_PSI
#ifdef SCHROEDINGER_CALC_PSI
if(getPsi(Op,Pa,*W,Vosci,Psi,Epsilon))
{
    cout << "succesfully extracted eigenvalues and eigenvectors from the final solution W!" << std::endl;
     cout << "schroedinger: returnEpsilon " << Epsilon->col(0) << endl;

}
latX.newLine("The electron states have the following energies:");

#endif
#define SCHROEDINGER_PLOT_RESULT
#ifdef SCHROEDINGER_PLOT_RESULT
arma::mat dat(Psi->n_rows,1);
string name;
for(unsigned int st=0; st<Psi->n_cols; st++)
{
     cout << "schroedinger: returnEpsilon " << Epsilon->row(0) << endl;

    //arma::cx_mat AEpsilon = cx_mat(Epsilon.get()->zeros(),Epsilon.get()->zeros());
    std::cout << "=== State " << st << ", has energy = " << Epsilon->row(st) << " === " << std::endl;
    latX.newLine("=== state " + std::to_string(st) + " has energy " + std::to_string(as_scalar(Epsilon->row(st)))+ "hartree ===");

    verbosity("now we get the square of the wavefunction for output",2,__FILE__,__LINE__);
    dat=real(pow(*Op.I*Psi->col(st),2));
    verbosity("print three slices of result to as ppm / gnuplot",2,__FILE__,__LINE__);
    #ifdef PLOT_RESULT
    for(int k=0; k<3; k++)
    {
        string plane;
        if(k==0) {plane="yz";}
        else if(k==1) {plane="xz";}
        else if(k==2) {plane="xy";}
        verbosity("now take a slice of data at plane "+plane,2,__FILE__,__LINE__);
        arma::mat sl;
        verbosity("take a slice "+plane,2,__FILE__,__LINE__);
        sl=myFunctions::slice(dat,Pa.S,(int) Pa.S(k)/2.,k);
        #ifdef PLOT_PPM
        verbosity("get name for ppm file "+plane,2,__FILE__,__LINE__);
        name="psi"+std::to_string(st)+"d_m_"+plane; //std::to_string(k)
        verbosity("write ppt file "+name,2,__FILE__,__LINE__);
        ppm(name,sl*0.3,sl,sl,&latX);
        #endif
        #ifdef PLOT_GNUPLOT
        if(st==0) latX.newLine("The physical eingenstates of the hamiltonian Psi are shown in the next figure:");
        if(k==0)
        {
            name="schroedinger_psi_"+std::to_string(st);
            string imageName=gpL.plotAMatrixSlice(name,dat,Pa.S,k);
            string caption="result from solution of schroedinger equation: slice through wave function of electron " + std::to_string(st) + " plane "+plane;
            latX.insertImage(imageName,caption,Pa.caseName);
        }
        #endif
    }
    #endif
}
#endif

Ui=Uinvers(Op,Pa,*W.get());
density=computeDensityFromWavefuncs(Op,Pa,*W.get(),Ui);
imgName=gpL.plotAMatrixSlice("schroedinger_density_final",density,Pa.S,0);
latX.newLine("Finally a density distribution as shown in the next image is computed.");
latX.insertImage(imgName,"final result of Schroedinger equation: density",Pa.caseName);
