#include "nonadiabaticPropagator.hpp"
#include "utilities.hpp"
#include "constants.hpp"

#define Ith(v,i)    NV_Ith_S(v,i-1)       /* Ith numbers components 1..NEQ */

nonadiabaticPropagator::nonadiabaticPropagator(inputParameters &IP) :
  propagatorBase(IP),
  firstRun_(true),
  t0_(0.0),
  time_(t0_),
  tFinal_(IP.max_time_),
  noutputs_(IP.n_outputs_),
  pulses_(IP.pulses_),
  atol_(IP.atol_),
  rtol_(IP.rtol_)
{
  dt_ = (tFinal_ - t0_) / noutputs_;
  initializeCVODE();
  initializeOutputs(IP);
  if ( (molecule_->Us_ != nullptr) && (molecule_->invUs_ != nullptr) )
    transformObservables();
}

void nonadiabaticPropagator::initializeCVODE()
{
  for (int ii = 0; ii < basisSets_->size(); ii++)
  {
    int nEq = pow(basisSets_->at(ii)->size(),2)*2;
    atols_.push_back(N_VNew_Serial(nEq));
    if (check_flag((void *)atols_.back(), (char *)"N_VNew_Serial", 0))
      throw std::runtime_error("Error allocating space for atol vector");
    ys_.push_back(N_VNew_Serial(nEq));
    if (check_flag((void *)ys_.back(), (char *)"N_VNew_Serial", 0))
      throw std::runtime_error("Error allocating space for y vector");
    // cvode_managers_.push_back(CVodeCreate(CV_BDF, CV_NEWTON));
    cvode_managers_.push_back(CVodeCreate(CV_ADAMS, CV_FUNCTIONAL));
    if (check_flag((void *)cvode_managers_.back(), (char *)"CVodeCreate", 0))
      throw std::runtime_error("Error allocating cvode object");
  }
}

void nonadiabaticPropagator::initializeOutputs(inputParameters &IP)
{
  // Setup Output stream
  output_file_name_ = IP.molecule_name_ + "_nonad.txt";
  out_file_.open(output_file_name_);

  // Initialize observable objects
  if (IP.output_cos3D_) observables_.push_back(std::make_shared<obsCosTheta3D>(basisSets_,fieldFreeHamiltonians_));
  if (IP.output_energy_) observables_.push_back(std::make_shared<obsEnergy>(basisSets_,fieldFreeHamiltonians_));
  if (IP.output_cos3DAlt_) observables_.push_back(std::make_shared<obsCosThetaAlt>(basisSets_,fieldFreeHamiltonians_));
  if (IP.output_J_) observables_.push_back(std::make_shared<obsJ>(basisSets_,fieldFreeHamiltonians_));
  if (IP.output_K_) observables_.push_back(std::make_shared<obsK>(basisSets_,fieldFreeHamiltonians_));
  if (IP.output_M_) observables_.push_back(std::make_shared<obsM>(basisSets_,fieldFreeHamiltonians_));


  // Print data identifiers
  out_file_ << "#Time (ps)";
  for (auto obs : observables_)
    out_file_ << "\t" << obs->id_tag_;
  out_file_ << std::endl;
}

void nonadiabaticPropagator::printOutputs()
{
  out_file_ << time_/CONSTANTS::AUperFS/1000.0 << "\t";
  for (auto obs : observables_)
    out_file_ << "\t" << obs->evaluate_(densities_);
  out_file_ << std::endl;
}

void nonadiabaticPropagator::run()
{
  for (int ii = 0; ii < noutputs_; ii++)
    step();
}

int nonadiabaticPropagator::evalRHS(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  nonadiabaticPropagator* obj = (nonadiabaticPropagator*)(user_data);
  int ii = obj->index_flag_;
  int N  = obj->basisSets_->at(ii)->size();
  double efield = 0.0;
  for (auto &p : obj->pulses_)
    efield += p.evaluate(t);
  std::shared_ptr<matrixComp> totalH;
  if (efield > 1.0e-10)
    totalH = std::make_shared<matrixComp>(*(obj->fieldFreeHamiltonians_->at(ii)) + *(obj->intHamiltonians_->at(ii))*efield);
  else
    totalH = obj->fieldFreeHamiltonians_->at(ii);
  zgemm3m_("N", "N", N, N, N, cplx(0.0,-1.0), totalH->data(),                                        N, reinterpret_cast<const cplx *>(N_VGetArrayPointer(y)), N, cplx(0.0), reinterpret_cast<cplx *>(N_VGetArrayPointer(ydot)), N);
  zgemm3m_("N", "N", N, N, N, cplx(0.0,1.0),  reinterpret_cast<const cplx *>(N_VGetArrayPointer(y)), N, totalH->data(),                                        N, cplx(1.0), reinterpret_cast<cplx *>(N_VGetArrayPointer(ydot)), N);
  // if dissipation needed, add terms here
  return(0);
}

void nonadiabaticPropagator::step()
{
/// Step using the boost diff eq libraries
#if 0
  copy_n( reinterpret_cast<double*>(P_.data()),n_*n_*2, x_.data());

//odeint control parameters
  double tInit  = time_;
  double tFinal = time_+t;
  double step   = t/1000.0;

  boost::function<void (const vector<double> &, vector<double> &, const double)> RHS (boost::bind(&DMPropagator::rhs,this,_1,_2,_3));
  boost::numeric::odeint::integrate( RHS , x_, tInit , tFinal , step);
  copy_n( x_.data(), x_.size(), (double *)P_.data() );
  time_ += t;
#endif

/// Step using CVODE libraries
  int flag;
  double tInit  = time_;
  double tFinal = time_ + dt_;
  realtype t;
  printOutputs();
  for (int ii = 0; ii < basisSets_->size(); ii++)
  {
    // if the population is too small, skip this set
    if (densities_->at(ii)->trace().real() < 1.0e-4)
      continue;
    int nEq = pow(basisSets_->at(ii)->size(),2)*2;
    flag    = CVodeSetUserData(cvode_managers_.at(ii),(void*)this);
    std::fill_n(&Ith(atols_.at(ii),1), nEq, atol_);
    std::copy_n(reinterpret_cast<double*>(densities_->at(ii)->data()), nEq, N_VGetArrayPointer(ys_.at(ii)) );

    if (firstRun_)
    {
      flag = CVodeInit(cvode_managers_.at(ii), &nonadiabaticPropagator::evalRHS, tInit, ys_.at(ii));
      if (check_flag(&flag, (char *)"CVodeInit", 1))
        throw std::runtime_error("Error initializing CVODE");
    }
    else
    {
      flag = CVodeReInit(cvode_managers_.at(ii), tInit, ys_.at(ii));
      if (check_flag(&flag, (char *)"CVodeReInit", 1))
        throw std::runtime_error("Error reinitializing CVODE");
    }

    flag = CVodeSVtolerances(cvode_managers_.at(ii), rtol_, atols_.at(ii));
    if (check_flag(&flag, (char *)"CVodeSVtolerances", 1))
      throw std::runtime_error("Error during CVODE step with CVodeSVtolerances function.");

    // flag = CVLapackBand(cvode_mem_, nEq,0,10);
    flag = CVBand(cvode_managers_.at(ii),nEq,2,10);
    // flag = CVDense(cvode_mem_, nEq);
    if (check_flag(&flag,(char *)"CVSolver", 1))
      throw std::runtime_error("Error during CVODE step with solver function.");

    index_flag_ = ii;
    flag = CVode(cvode_managers_.at(ii), tFinal, ys_.at(ii), &t, CV_NORMAL);
    if (check_flag(&flag, (char *)"CVode", 1))
      throw std::runtime_error("Error during CVODE step with CVode function.");

    std::copy_n(N_VGetArrayPointer(ys_.at(ii)), nEq, reinterpret_cast<double *>(densities_->at(ii)->data()));
  }
  if (firstRun_) firstRun_ = false;
  time_ = tFinal;
}



/**************************
  CVODE Helper functions
**************************/

int check_flag(void *flagvalue, char *funcname, int opt)
{
  int *errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
      funcname);
    return(1); }

  /* Check if flag < 0 */
  else if (opt == 1) {
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
        funcname, *errflag);
      return(1); }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && flagvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
      funcname);
    return(1); }

  return(0);
}

void PrintFinalStats(void *cvode_mem)
{
  long int nst, nfe, nsetups, nje, nfeLS, nni, ncfn, netf, nge;
  int flag;

  flag = CVodeGetNumSteps(cvode_mem, &nst);
  check_flag(&flag, (char *)"CVodeGetNumSteps", 1);
  flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);
  check_flag(&flag, (char *)"CVodeGetNumRhsEvals", 1);
  flag = CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
  check_flag(&flag, (char *)"CVodeGetNumLinSolvSetups", 1);
  flag = CVodeGetNumErrTestFails(cvode_mem, &netf);
  check_flag(&flag, (char *)"CVodeGetNumErrTestFails", 1);
  flag = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
  check_flag(&flag, (char *)"CVodeGetNumNonlinSolvIters", 1);
  flag = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
  check_flag(&flag, (char *)"CVodeGetNumNonlinSolvConvFails", 1);

  flag = CVDlsGetNumJacEvals(cvode_mem, &nje);
  check_flag(&flag, (char *)"CVDlsGetNumJacEvals", 1);
  flag = CVDlsGetNumRhsEvals(cvode_mem, &nfeLS);
  check_flag(&flag, (char *)"CVDlsGetNumRhsEvals", 1);

  flag = CVodeGetNumGEvals(cvode_mem, &nge);
  check_flag(&flag, (char *)"CVodeGetNumGEvals", 1);

  printf("\nFinal Statistics:\n");
  printf("nst = %-6ld nfe  = %-6ld nsetups = %-6ld nfeLS = %-6ld nje = %ld\n",
   nst, nfe, nsetups, nfeLS, nje);
  printf("nni = %-6ld ncfn = %-6ld netf = %-6ld nge = %ld\n \n",
   nni, ncfn, netf, nge);
}
