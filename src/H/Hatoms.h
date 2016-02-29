latX.subSection("compute electronic configuration of hydrogen atom(s)!");

verbosity("compute electronic energy of hydrogen atom(s)!",2,__FILE__,__LINE__);
verbosity("first compute the ionic potential from all the (static) ion cores [with core electrons]",2,__FILE__,__LINE__);

latX.subsubSection("first prepare the ion potential");

latX.newLine(" Fourier transform of a single ion potential: ");
latX.newLine(" $V_{ps}(\\vec{k})=-4 \\pi Z G2^-{1}$ ");

arma::mat Vps(-4.*PI*Pa.Z/G2);
Vps[0]=0;
verbosity("Hatom: use the structure factor to distribute the single potential to all the ionic sites!",2,__FILE__,__LINE__);
latX.newLine(" Use structure function to get distributed ion potential: ");
latX.newLine(" $V_{ion}(\\vec{k})=S_{f} V_{ps}$ ");
latX.newLine(" Transform this global potential back to r space : ");
latX.newLine(" $V_{dual}(\\vec{x})=FFT^{-1}(V_{ion}(\\vec{k}))=J(V_{ion}(\\vec{k})$ ");
latX.newLine(" $V_{dual}:$ vector with $S_{k}$ elements in r space ");

arma::mat Vdual(Vps);
Vdual=real(*Op.J*(Vps%Sf));
verbosity("Hatom: now plot the Vdual potential",2,__FILE__,__LINE__);

#ifdef PLOT_GNUPLOT
string name="hatoms_Vdual";
string imageName=gpL.plotAMatrixSlice(name,Vdual,Pa.S,0);
string caption="Hatoms: The initial ion potential as computed using the structure function - "+gpL.getPlaneName()+" plane.";
latX.insertImage(imageName,caption);
#endif

latX.subsubSection("Preparation of the wavefunction");
latX.newLine(" Initialize wavefunctions randomĺy : ");
latX.newLine(" $W_{m}:$ complex matrix with $ \\Sigma S_{k} \\times N_{s} $ elements in r space ");
verbosity("use random generator to setup random wavefunctions!",2,__FILE__,__LINE__);
//Wm=arma::randu<arma::cx_mat>(prodS,Pa.number_of_wavefunctions);
//std::shared_ptr<arma::cx_mat> W(&Wm);
//arma::cx_mat Wm(prodS,Pa.number_of_wavefunctions);
//std::shared_ptr<arma::cx_mat> W = new arma::randu<arma::cx_mat>(prodS,Pa.number_of_wavefunctions);
std::shared_ptr<arma::cx_mat> W(new arma::cx_mat
(arma::randu<arma::cx_mat>(Pa.prodS,Pa.number_of_wavefunctions)));

verbosity("create optimizer sd",2,__FILE__,__LINE__);

latX.subsubSection("Do a first energy minimisation using steepest descent method.");

sdOptimizer sd(Op,Pa,Vps,&latX);
verbosity("push the generated latex comment to output",2,__FILE__,__LINE__);
sd.myLatexClass->commentMyFunction();
verbosity("prepare sd optimizer",2,__FILE__,__LINE__);
sd.setup(W);
sd.myLatexClass->commentMyFunction();
verbosity("====================================",0,__FILE__,__LINE__);
verbosity("optimize using sd",0,__FILE__,__LINE__);
verbosity("====================================",0,__FILE__,__LINE__);
verbosity("20 optimization steps with steepest descent",2,__FILE__,__LINE__);
sd.optimize(W);
sd.myLatexClass->commentMyFunction();
verbosity("restart as orthonormal function",2,__FILE__,__LINE__);
//*W=W->t()*(*Op.O*(*W));
*W=*W*inv(sqrt((W->t())*(*Op.O*(*W)))); // Restart as orthonormal functions
#ifdef DO_PCCG
    verbosity("create optimizer pccg",2,__FILE__,__LINE__);
    pccgOptimizer pccg(Op,Pa,Vdual,G2,&latX);
    verbosity("prepare pccg optimizer",2,__FILE__,__LINE__);
    pccg.setup(W);
    pccg.orthogonalizeWfunc(W);
    verbosity("====================================",2,__FILE__,__LINE__);
    verbosity("optimize using pccg",2,__FILE__,__LINE__);
    verbosity("====================================",2,__FILE__,__LINE__);
    pccg.optimize(W);
#endif
verbosity("now get the wavefunction",2,__FILE__,__LINE__);
// Psi as a pointer to the wavefunction
std::shared_ptr<arma::cx_mat> Psi(new arma::cx_mat(*W));
// Epsilon points to the eigenvalues (real vector)
//std::shared_ptr<arma::cx_vec> Epsilon(new arma::cx_vec(W->n_rows));
std::shared_ptr<arma::mat> Epsilon(new arma::mat(W->n_rows,1));

#define HATOMS_CALC_PSI
#ifdef HATOMS_CALC_PSI
if(getPsi(Op,Pa,*W,Vps,Psi,Epsilon))
{
    cout << "succesfully extracted eigenvalues and eigenvectors from the final solution W!" << std::endl;
     cout << "schroedinger: returnEpsilon " << Epsilon->col(0) << endl;

}
latX.newLine("The electron states have the following energies:");

#endif
#define HATOMS_PLOT_RESULT
#ifdef HATOMS_PLOT_RESULT
arma::mat dat(Psi->n_rows,1);
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
            name="HATOMS_psi_"+std::to_string(st);
            string imageName=gpL.plotAMatrixSlice(name,dat,Pa.S,k);
            string caption="result from solution of schroedinger equation: slice through wave function of electron " + std::to_string(st) + " plane "+plane;
            latX.insertImage(imageName,caption);
        }
        #endif
    }
    #endif
}
#endif

cx_mat Ui=Uinvers(Op,Pa,*W.get());
mat density=computeDensityFromWavefuncs(Op,Pa,*W.get(),Ui);
string imgName=gpL.plotAMatrixSlice("HATOMS_density_final",density,Pa.S,0);







//if(getPsi(Op,Pa,*W,Vdual,Psi,Epsilon))
//{
//    cout << "succesfully extracted eigenvalues and eigenvectors from the final solution W!" << std::endl;
//}
//verbosity("now we get the square of the wavefunction for output",2,__FILE__,__LINE__);
//arma::mat dat(Psi->n_rows,Psi->n_cols);
//for(unsigned int st=0; st<dat.n_cols; st++)
//{
//    arma::cx_mat AEpsilon = cx_mat(Epsilon.get()->zeros(),Epsilon.get()->zeros());
//    std::cout << "=== State " << st << ", Energy = " << real(AEpsilon(st)) << " === " << std::endl;
//    dat.col(st)=real(pow(*Op.I*Psi->col(st),2));
//    verbosity("print result to file",2,__FILE__,__LINE__);
//    for(int k=0; k<3; k++)
//    {
//        string plane;
//        if(k==0) {plane="yz";}
//        else if(k==1) {plane="xz";}
//        else if(k==2) {plane="xy";}
//        verbosity("now take a slice of data at plane "+plane,2,__FILE__,__LINE__);
//        arma::mat sl;
//        verbosity("take a slice "+plane,2,__FILE__,__LINE__);
//        sl=myFunctions::slice(dat,Pa.S,(int) S(k)/2.,k);
//        verbosity("get name for ppm file "+plane,2,__FILE__,__LINE__);
//        string name="psi"+std::to_string(st)+"d_m_"+plane; //std::to_string(k)
//        verbosity("write ppt file "+name,2,__FILE__,__LINE__);
//        ppm(name,sl*0.3,sl,sl,&latX);
//    }
//}
