#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

// SLC, Nov 9, 2011

#ifndef _CGOPTIMIZE_H_
#define _CGOPTIMIZE_H_

#define CGOPT_LINE_SEARCH_SECANT 0
#define CGOPT_LINE_SEARCH_BISECTION 1
#define CGOPT_LINE_SEARCH_REGULAR_FALSI 2

#include "parstream.H"
#include "REAL.H"
#include "Misc.H"
#include <cmath>
#include "NamespaceHeader.H"
#include "CH_assert.H"

/** @file */ 
inline Real sqr(const Real& a)
{
  return a*a;
}


/// Specify operations required to carry out gradient based optimization
template <class X>
class ObjectiveWithGradient {

public:

  ///virtual destructor
   virtual ~ObjectiveWithGradient(){};

  /** 
      Compute dot product a.b
      @return a.b
  */
  virtual Real dotProduct(X& a , const X& b) = 0; 


  /**
     Compute a main objective (fm), a penalty function (fp) and the 
     gradient (g) of f + p with respect to x

     @param fm (output) main objective
     @param fp (output) penalty function
     @param g  (output) gradient
     @param x  (input)  state vector
     @param iter (input) iteration number 
     @param inner  (input)  is this an 'inner' iteration?
  */
  virtual void computeObjectiveAndGradient(Real& fm, Real& fp ,X& g, const X& x , int iter, bool inner) = 0;

  /// Copy x to y
   virtual void assign(X& y, const X& x) = 0;

  /// Set y = y + s * x
   virtual void incr(X& y, const X& x, Real s) = 0;

  /// resize a / allocate storage to conform with b
   virtual void create (X &a, const X &b) = 0;

  /// free storage 
   virtual void free (X &a) = 0;

  /// set a = a * s
   virtual void scale(X &a, const Real s) = 0;

  /// apply precoditioner K, so a = Kb
   virtual void preCond(X &a, const X &b) = 0;

  /// carry out whatver operations are required prior to a restart
   virtual void restart() = 0;

  /// number of degrees-of-freedom (maximum number of iterations before restart)
   virtual int nDoF(const X& x) = 0;

  /// save a state that would allow CG to resume from iteration iter.
  virtual void saveCGState(int iter, const X& x, const X& r, const X& s, const X& d) = 0;
   /// read a state allowing CG to resume. Should return true if a state was read
  virtual bool readCGState(int& iter, bool& fresh_start, X& x,  X& r, X& s, X& d) = 0;

};

/// bisection method: find y = x + s*d such that |f'(y).f'(x)| << f'(x).f'(x)
/**
   x - initial point
   r - input - -f'(x) (steepest descent) 
   d - input - (preconditioned) search direction / steepest descent ( d = M^-1 r)
   y - output
   w - workspace
   a_fx - f(x) (where f' is the gradient of scalar f > 0)

 */
template < class X>
void bisectionLineSearch
(ObjectiveWithGradient<X>& a_F, X& x, X& r, X& d, X& y, X& w, Real a_fx, 
 int a_maxIter, Real a_initialStep, Real a_stepMaxGrow, Real a_tol)
{

  Real a, b, c, s; // brackets, mid-point, increment along d
  Real fa, fb, fc; // f'(x + s d).d @ s = a,b,c
  Real Fm, Fp; // misfit and penalty norms (F = Fm + Fp)
  Real f0, s0; // initial values

  int iter(0);
  a = 0.0;
  fa = a_F.dotProduct(r,d);
  f0 = fa;
  s0 = a_initialStep;
  s = s0;
  b = a;
  
  
  Real fb_prev = 0.0;
  Real b_prev;
  Real mult = 2.0;
  bool doubling = true;
  //walk the line by doubling until fb < 0 (since fa > 0) 
  do {
	b_prev = b;
	if (iter == 0) b += s0;
	if ((fb > fb_prev) && (iter > 1)) {
		b = b_prev;
		doubling = false;
		iter -= 2;
		pout() << " ... initial bisection: too far, switching to single steps " << std::endl;
	}
	else b *= mult;
	fb_prev = fb;
    //b += s0;
    a_F.assign(y, x);
    a_F.incr(y, d, b);
    a_F.computeObjectiveAndGradient(Fm,Fp,w,y,iter,true);
    fb = -a_F.dotProduct(w,d);

    pout() << " ... initial bisection iter " << iter
	   << " f'(x + sd).f'(x) = " << fb
	   << " f'(x).f'(x) = " << f0
	   << " s = " << b << " f(x + sd) = " << Fm + Fp << std::endl;

    iter++;
  } while ( (2*iter < a_maxIter) && (fb > 0.0) && (s > a_tol * s0 ) && (doubling));
  
  //walk the line in single steps until fb < 0 (since fa > 0) 
	Real step = (iter > 1) ? b/10 : s0;
  if (fb > 0.0) do {
	b_prev = b;
	if (iter == 0) step = s0/10.0;
	if ((iter > 1) && (fb > fb_prev)) { step *= -0.5; pout() << " ... initial reg falsi: changing search direction " << std::endl;}
    fb_prev = fb;
    b += step;
    a_F.assign(y, x);
    a_F.incr(y, d, b);
    a_F.computeObjectiveAndGradient(Fm,Fp,w,y,iter,true);
    fb = -a_F.dotProduct(w,d);

    pout() << " ... bisection: searching for a bracket: iter " << iter
	   << " f'(x + bd).d = " << fb
	   << " f'(x).d = " << f0
	   << " b = " << b << " f(x + bd) = " << Fm + Fp << std::endl;
    
    if (fb < fa)
      {
	if (fb > 0.0)
	  {
	    // might as well move the lower limit 
	    fa = fb; a = b;	
	  }
      }
    else
      {
	// if d really is a descent direction, a reduction should help
        // need to start again?
	b = 0.0;
	a = 0.0;;
	s *= 0.5;
      }  
    
    iter++;
  } while ( (2*iter < a_maxIter) && (fb > 0.0) && (s > a_tol * s0 )  && (!doubling));
	
  a = std::min(b,b_prev); 
  b = std::max(b,b_prev); 
  fa = (b<b_prev)?fb:fb_prev;
  if (iter < a_maxIter) iter = 0;
  
  do {
    c = 0.5*(a + b); //midpoint
    a_F.assign(y, x);
    a_F.incr(y, d, c);
    a_F.computeObjectiveAndGradient(Fm,Fp,w,y,iter,true);

    fc = -a_F.dotProduct(w,d);
    
    pout() << " ... bisection iter " << iter 
	   << "  f'(x + sd).d = " << fc
	   << " f'(x).d = " << f0
	   << " s = " << c << " f(x + sd) = " << Fm + Fp << std::endl;
    
      if (std::signbit(fc) == std::signbit(fa))
      {
	a = c;
	fa = fc;
      }
    else
      {
	b = c;
	fb = fc;
      }
    iter++;
    
  } while ( (iter < a_maxIter) && (Abs(fc) > sqr(a_tol) * Abs(f0)) && (Abs(b-a) > a_tol * s0));
  a_F.assign(x, y);

}



/// regular falsi method: find y = x + s*d such that |f'(y).f'(x)| << f'(x).f'(x)
/**
   x - initial point
   r - input - -f'(x) (steepest descent) 
   d - input - (preconditioned) search direction / steepest descent ( d = M^-1 r)
   y - output
   p - workspace
   a_fx - f(x) (where f' is the gradient of scalar f > 0)

 */
template < class X>
void regularFalsiLineSearch
(ObjectiveWithGradient<X>& a_F, X& x, X& r, X& d, X& y, X& p, Real a_fx, 
 int a_maxIter, Real a_initialStep, Real a_stepMaxGrow, Real a_tol)
{

  Real a, fa, b, fb, c, fc, fm, fp, f0, s0;
  int iter(0);
  a = 0.0;
  fa = a_F.dotProduct(r,d);
  f0 = fa;
  s0 = a_initialStep;
  b = a;
  Real fb_prev = 0.0;
  Real b_prev;
  Real mult = 2.0;
  bool doubling = true;
  //walk the line by doubling until fb < 0 (since fa > 0) 
  do {
	b_prev = b;
	if (iter == 0) b += s0;
	if ((fb > fb_prev) && (iter > 1)) {
		b = b_prev;
		doubling = false;
		iter -= 2;
		pout() << " ... initial reg falsi: too far, switching to single steps " << std::endl;
	}
	else b *= mult;
	fb_prev = fb;
    //b += s0;
    a_F.assign(y, x);
    a_F.incr(y, d, b);
    a_F.computeObjectiveAndGradient(fm,fp,p,y,iter,true);
    fb = -a_F.dotProduct(p,d);

    pout() << " ... initial reg falsi iter " << iter
	   << " f'(x + sd).f'(x) = " << fb
	   << " f'(x).f'(x) = " << f0
	   << " s = " << b << " f(x + sd) = " << fm + fp << std::endl;

    iter++;
  } while ( (iter < a_maxIter) && (fb > 0.0) && (doubling));
  
  //walk the line in single steps until fb < 0 (since fa > 0) 
	Real step = (iter > 1) ? b/10 : s0;
  if (fb > 0.0) do {
	b_prev = b;
	if (iter == 0) step = s0/10.0;
	if ((iter > 1) && (fb > fb_prev)) { step *= -0.5; pout() << " ... initial reg falsi: changing search direction " << std::endl;}
    fb_prev = fb;
    b += step;
	a_F.assign(y, x);
    a_F.incr(y, d, b);
    a_F.computeObjectiveAndGradient(fm,fp,p,y,iter,true);
    fb = -a_F.dotProduct(p,d);
	
	// Check the gradient sign. If negative, search backwards

    pout() << " ... initial reg falsi iter " << iter
	   << " f'(x + sd).f'(x) = " << fb
	   << " f'(x).f'(x) = " << f0
	   << " s = " << b << " f(x + sd) = " << fm + fp << std::endl;

    iter++;
  } while ( (iter < a_maxIter) && (fb > 0.0) && (!doubling) );
	
  a = std::min(b,b_prev); 
  b = std::max(b,b_prev); 
  fa = (b<b_prev)?fb:fb_prev;
  if (iter < a_maxIter) iter = 0;

  do {
    //c = 0.5*(a + b); //midpoint
	c = a - (fa*(b-a)) / (fb - fa); // linearly interpolated root
    a_F.assign(y, x);
    a_F.incr(y, d, c);
    a_F.computeObjectiveAndGradient(fm,fp,p,y,iter,true);

    fc = -a_F.dotProduct(p,d);
    
    pout() << " ... reg falsi iter " << iter 
	   << "  f'(x + sd).f'(x) = " << fc
	   << " f'(x).f'(x) = " << f0
	   << " s = " << c << " f(x + sd) = " << fm + fp << " a = " << a << " b = " << b << std::endl;
    
      if (std::signbit(fc) == std::signbit(fa))
      {
	a = c;
	fa = fc;
      }
    else
      {
	b = c;
	fb = fc;
      }
    iter++;
    
  } while ( (iter < a_maxIter) && (Abs(fc) > sqr(a_tol) * Abs(f0)) && (Abs(b-a) > a_tol * s0));
  a_F.assign(x, y);

}



/// Secant method: find y = x + s*d such that |f'(y).f'(x)| << f'(x).f'(x)
/**
   x - initial point
   r - input - ||f'(x)||^2
   d - input - search direction
   y - output
   p - workspace
   a_fx - f(x) (where f' is the gradient of scalar f > 0)

 */
template < class X>
void secantLineSearch
(ObjectiveWithGradient<X>& a_F, X& x, X& r, X& d, X& y, X& p, Real a_fx, 
 int a_maxIter, Real a_initialStep, Real a_stepMaxGrow, Real a_tol)
{

  //secant line search
  Real deltaD = a_F.dotProduct(r,d);
  Real eta = -deltaD;
  Real initialStep = a_initialStep;
  Real alpha = -initialStep;
  a_F.assign(y,x);
  a_F.incr(y,d,initialStep);
  Real fms,fps;
  a_F.computeObjectiveAndGradient(fms,fps,p,y,0,true);
  Real etaPrev = a_F.dotProduct(p,d);
  
  pout() << " ... initial secant ||f'(x+sd).d||^2 = " 
	 << -etaPrev << " s = " << initialStep << " f(x+sd) = " << fms << std::endl;

  int iter = 0;
  while ( (fms + fps <  a_fx) && (-etaPrev > -eta) && iter < a_maxIter)
    {
      //we have a problem, because although we have a descent direction in the  objective
      //function f , we have an ascent direction for f'(x+sd).d. For now, just walk along d till
      // we do have descent in both, then start the secant method offset by alpha + a_initialStep
      pout() << " d is a descent direction for f(x+sd) but not for f'(x+sd).d  " 
	     << " so walking the line (remind SLC to look up a better method) " << std::endl;
      alpha -= initialStep;
      a_F.incr(y,d,initialStep);
      a_F.computeObjectiveAndGradient(fms,fps,p,y,0,true);
      eta = etaPrev;
      etaPrev = a_F.dotProduct(p,d);
      pout() << " ... initial secant f'(x+ sd).d = " 
	     << -etaPrev << " s = " << alpha << " f(x+ sd) = " << a_fx << std::endl;
      iter++;
    }
  
  iter = 0;
  while ( (fms + fps >  a_fx) && iter < a_maxIter )
    {
      // here we have an ascent direction the objective function. Try reducing the step.
      // this can of course be a happy failure
      pout() << " d is a not descent direction for f(x+sd), reducing s ";
      a_F.incr(y,d,-0.5*initialStep);
      initialStep *= 0.5;
      alpha = - initialStep;
      a_F.computeObjectiveAndGradient(fms,fps,p,y,0,true);
      eta = etaPrev;
      etaPrev = a_F.dotProduct(p,d);
      pout() << " ... initial secant f'(x+ sd).d = " 
	     << -etaPrev << " s = " << initialStep << " f(x+ sd) = " << a_fx << std::endl;
      iter++;
    }
  
  int j = 0;
  Real q = -(alpha + a_initialStep); // normally alpha = - a_initialStep, but sometimes we needed to walk the line
  if (q > 0.0)
    a_F.incr(x,d,q);
  
  alpha = -a_initialStep;
  
  do {
    
    if (j > 0)
      {
	a_F.computeObjectiveAndGradient(fms,fps,r,x,j,true);
	eta = a_F.dotProduct(r,d);
      }
    
    pout() <<  " ... secant iteration j = " << j << " f'(x[0]+qd).d = " << -eta
	   << " q = " << q 
	   << std::endl;
	
    Real absAlphaPrev = std::abs(alpha);
    alpha = alpha * eta / (etaPrev-eta);
    
    //limit the rate at which alpha grows
    if ( (a_stepMaxGrow > 0.0) && (std::abs(alpha) > a_stepMaxGrow * absAlphaPrev) )
      {
	alpha = ((alpha > 0.0)?(1.0):(-1.0)) * a_stepMaxGrow * absAlphaPrev;
      }
    
    q += alpha;
    pout() << " ... secant step a = " << alpha << "  q = " << q  << std::endl;
    
    a_F.incr(x,d,alpha);
    etaPrev = eta;
    j++;
  } while (j < a_maxIter  && (sqr(alpha) * deltaD > sqr(a_tol)));
  
  
}



/**
   Nonlinear Conjugate-gradient method, used to mimimise F(x) w.r.t x 
*/
template <class X>
int CGOptimize(ObjectiveWithGradient<X>& a_F, X& a_x, 
	       int  a_maxIter,
	       Real a_tol,
	       Real a_hang,
	       int a_lineSearchMethod,
	       Real a_lineSearchInitialStep,
	       Real a_lineSearchStepMaxGrow,
	       int  a_lineSearchMaxIter,
	       Real a_lineSearchTol,
	       int a_iter = 0)
{

  //pout() << "CGOptimize -- a_iter = " << a_iter << endl;
  //Preconditioned Nonlinear Conjugate Gradients 
  //with Secant line OR bisection search and Polak-Ribiere updates
  //from J. R. Shewchuk, 1994 
  //"An Introduction to the Conjugate Gradient Method
  // Without the Agonizing Pain"

  X& x = a_x;
  X r,s,d,y,p;
  a_F.create(r,x);
  a_F.create(s,x);
  a_F.create(d,x);
  a_F.create(y,x);
  a_F.create(p,x);

  Real fm,fp, fmOld; // just for retreiving the the misfit norm and penalty norm
  int iter(a_iter);

  //continue CG, or start afresh
  bool fresh_start(false);
  bool readCG = a_F.readCGState(iter, fresh_start, x, r, s, d);
  if (readCG && !fresh_start)
    {
      // Need to recompute fm,fp. We don't need to recompute r (= - grad F),
      // but worth knowing whether we *would* re-compute it
      a_F.computeObjectiveAndGradient(fm,fp,p,x,iter,false);
      // check p = -r  
      a_F.incr(p,r,1.0);
      Real pdp = a_F.dotProduct(p,p);
      Real rdr = a_F.dotProduct(r,r);
      pout() << "CGOptimize read state iteration " << iter 
	     << " (r_read-r_compute).(r_read-r_compute) / r_read.r_read = " << pdp/(rdr + 1.0e-10)
	     << std::endl;
    }
  else
    {
      // fresh start, either because we read nothing, or we read x but specified fresh_start
      a_F.computeObjectiveAndGradient(fm,fp,r,x,iter,false);
      a_F.scale(r,-1.0);
      a_F.preCond(s,r);
      a_F.assign(d,s);
    }
  
  
  fmOld = fm;

  Real deltaNew = a_F.dotProduct(r,s);
  Real deltaZero = deltaNew;
  int iter0(iter);
  int k(0);
  bool hang(false);
  // if we have a_iter > a_maxIter, we have presumably restarted
  int niter = a_maxIter;
  if ( (iter > a_maxIter) && (a_maxIter > 0))  niter += (iter/a_maxIter) * a_maxIter;
  
  while ( (iter < niter)
	  && (deltaNew > sqr(a_tol) * deltaZero) 
	  && (fm > 1.0e-16) 
	  && (!hang))
    {

      if (iter > iter0)
	{
	  a_F.saveCGState(iter, x, r, s, d);
	}

      pout() << "CGOptimize iteration " << iter 
	     << " ||f'(x)||^2 = " << deltaNew 
	     << " ||fm(x)||^2 = " << fm 
	     << " ||fp(x)||^2 = " << fp 
	     << " ||fm(x)||^2 + ||fp(x)||^ = " << fm + fp;

      if (iter > iter0)
	{
	  pout() << " ||fm(x)||^2/||fm_old(x)||^2 = " << fm/fmOld;
	  hang = (fm > a_hang * fmOld);
	}
      pout() << std::endl;

      fmOld = fm;
  

      //find the (approximate) root of f'(y).f'(x)
      if (a_lineSearchMethod == CGOPT_LINE_SEARCH_BISECTION)
	{
	  bisectionLineSearch(a_F, x, r, d, y, p, fm+fp, a_lineSearchMaxIter, 
			      a_lineSearchInitialStep, a_lineSearchStepMaxGrow, a_lineSearchTol);
	}
      else if (a_lineSearchMethod == CGOPT_LINE_SEARCH_REGULAR_FALSI)
	{
	  regularFalsiLineSearch(a_F, x, r, d, y, p, fm+fp, a_lineSearchMaxIter, 
			      a_lineSearchInitialStep, a_lineSearchStepMaxGrow, a_lineSearchTol);
	}
      else
	{
	  // we had secant as the only option since 2011, so make it default
	  secantLineSearch(a_F, x, r, d, y, p, fm+fp,  a_lineSearchMaxIter, 
			   a_lineSearchInitialStep, a_lineSearchStepMaxGrow, a_lineSearchTol);
	}
	
      a_F.computeObjectiveAndGradient(fm,fp,r,x,iter+1,false);
      a_F.scale(r,-1.0);
      Real deltaOld = deltaNew;
      Real deltaMid = a_F.dotProduct(r,s);
      a_F.preCond(s,r);
      deltaNew = a_F.dotProduct(r,s);
      Real beta = (deltaNew - deltaMid)/deltaOld;
      k++;
      if (k == a_F.nDoF(x) || beta <= 0.0 )
	{
	  pout() << "CGOptimize restart k = " << k 
		 << " beta = " << beta << std::endl;
	  
	  a_F.assign(d,s);
	  a_F.restart();
	  k = 0;
	}
      else
	{
	  a_F.scale(d,beta);
	  a_F.incr(d,s,1.0);
	}
      iter++;
    }
  
  pout() << "CGOptimize iteration " << iter 
	 << " ||f'(x)||^2 = " << deltaNew 
	 << " ||fm(x)||^2 = " << fm 
	 << " ||fp(x)||^2 = " << fp   
	 << " ||fm(x)||^2 + ||fp(x)||^ = " << fm + fp
	 << " ||fm(x)||^2/||fm_old(x)||^2 = " << fm/fmOld 
	 << std::endl;
  

  a_F.free(r);
  a_F.free(s);
  a_F.free(d);
  a_F.free(y);
  a_F.free(p);

  return iter;

}


#include "NamespaceFooter.H"
#endif /*_CGOPTIM_H_*/
