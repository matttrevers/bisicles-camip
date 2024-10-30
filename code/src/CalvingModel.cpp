#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "CalvingModel.H"
#include "MaskedCalvingModel.H"
#include "CrevasseCalvingModel.H"
#include "IceConstants.H"
#include "AmrIce.H"
#include "ParmParse.H"
#include <fstream>
#include "computeSum.H"
#include "computeNorm.H"
#include "AMRPoissonOpF_F.H"
#include "NamespaceHeader.H"
#include <limits>
#include <time.h>

/// a default implementation
/**
   most models provide a criterion, rather than a rate.
 */
void
CalvingModel::getCalvingRate(LevelData<FArrayBox>& a_calvingRate, const AmrIce& a_amrIce,int a_level)
{
  DataIterator dit = a_calvingRate.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      a_calvingRate[dit].setVal(0.0);
    }
}

void
CalvingModel::getWaterDepth(LevelData<FArrayBox>& a_waterDepth, const AmrIce& a_amrIce,int a_level)
{
  DataIterator dit = a_waterDepth.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      a_waterDepth[dit].setVal(0.0);
    }
}

void
VariableRateCalvingModel::getCalvingRate(LevelData<FArrayBox>& a_calvingRate, const AmrIce& a_amrIce,int a_level)
{
  m_calvingRate->evaluate(a_calvingRate, a_amrIce, a_level, 0.0);
}


void 
DeglaciationCalvingModelA::applyCriterion
(LevelData<FArrayBox>& a_thickness, 
 LevelData<FArrayBox>& a_calvedIce,
 LevelData<FArrayBox>& a_addedIce,
 LevelData<FArrayBox>& a_removedIce, 
 LevelData<FArrayBox>& a_iceFrac, 
 const AmrIce& a_amrIce,
 int a_level,
 Stage a_stage)
{
  const LevelSigmaCS& levelCoords = *a_amrIce.geometry(a_level);
  for (DataIterator dit(levelCoords.grids()); dit.ok(); ++dit)
    {
      const BaseFab<int>& mask = levelCoords.getFloatingMask()[dit];
      FArrayBox& thck = a_thickness[dit];
      FArrayBox& calved = a_calvedIce[dit];
      FArrayBox& added = a_addedIce[dit];
      FArrayBox& removed = a_removedIce[dit];
      Box b = thck.box();

      for (BoxIterator bit(b); bit.ok(); ++bit)
	{
	  const IntVect& iv = bit();
	  Real prevThck = thck(iv);
	  if (mask(iv) == OPENSEAMASKVAL)
	    {
	      thck(iv) = 0.0;
	    }
	  else if (mask(iv) == OPENLANDMASKVAL)
	    {
	      thck(iv) = 0.0;
	    }
	  else
	    {
	      thck(iv) = std::max(thck(iv),m_minThickness);
	    }
	      
	  // Record gain/loss of ice
	  if (calved.box().contains(iv))
	    {
	      updateCalvedIce(thck(iv),prevThck,mask(iv),added(iv),calved(iv),removed(iv));
	    }
	}
    }
}

void DomainEdgeCalvingModel::applyCriterion
(LevelData<FArrayBox>& a_thickness,
 LevelData<FArrayBox>& a_calvedIce,
 LevelData<FArrayBox>& a_addedIce,
 LevelData<FArrayBox>& a_removedIce,  
 LevelData<FArrayBox>& a_iceFrac, 
 const AmrIce& a_amrIce,
 int a_level,
 Stage a_stage)
{

  const LevelSigmaCS& levelCoords = *a_amrIce.geometry(a_level);
  const DisjointBoxLayout& grids = levelCoords.grids();
  const ProblemDomain domain = grids.physDomain();
  const LevelData<BaseFab<int> >& levelMask = levelCoords.getFloatingMask();
  const IntVect ghost = a_thickness.ghostVect();
  //const LevelData<FArrayBox>& vt  = *a_amrIce.viscousTensor(a_level);
  DataIterator dit = grids.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      //const Box& gridBox = grids[dit];
      for (int dir=0; dir<SpaceDim; dir++)
	{
	  if (!domain.isPeriodic(dir))
	    {

	      if (m_frontLo[dir] > 0)
		{
		  Box loBox = adjCellLo(domain,dir,ghost[dir]);
                  // (DFM 5-25-15) grow in transverse direction
                  // to ensure that we don't wind up with corner
                  // cells with ice in them
                  IntVect transverseVect = ghost;
                  transverseVect[dir] = 0;
                  loBox.grow(transverseVect);
		  loBox &= a_thickness[dit].box();
		  for (BoxIterator bit(loBox); bit.ok(); ++bit)
		    {
		      const IntVect& iv = bit();
		      const IntVect ip = iv + BASISV(dir);
		      //if (levelMask[dit](ip) != GROUNDEDMASKVAL)
		      Real prevThck = a_thickness[dit](iv);
		      a_thickness[dit](iv) = 0.0;
		      if (a_iceFrac[dit].box().contains(iv))
			a_iceFrac[dit](iv) = 0.0;
		      // Record gain/loss of ice
		      if (a_calvedIce[dit].box().contains(iv))
			{
			  updateCalvedIce(a_thickness[dit](iv),prevThck,levelMask[dit](iv),
					  a_addedIce[dit](iv),a_calvedIce[dit](iv),a_removedIce[dit](iv));
			}

		    }
		}
	      
	      if (m_frontHi[dir] > 0)
		{
		  Box hiBox = adjCellHi(domain,dir,ghost[dir]);
                  // (DFM 5-25-15) grow in transverse direction
                  // to ensure that we don't wind up with corner
                  // cells with ice in them
                  IntVect transverseVect = ghost;
                  transverseVect[dir] = 0;
                  hiBox.grow(transverseVect);
		  hiBox &= a_thickness[dit].box();
		  for (BoxIterator bit(hiBox); bit.ok(); ++bit)
		    {
		      const IntVect& iv = bit();
		      const IntVect ip = iv - BASISV(dir);
		      //if (levelMask[dit](ip) != GROUNDEDMASKVAL)
		      Real prevThck = a_thickness[dit](iv);
		      a_thickness[dit](iv) = 0.0;
		      if (a_iceFrac[dit].box().contains(iv))
			a_iceFrac[dit](iv) = 0.0;
		      // Record gain/loss of ice
		      if (a_calvedIce[dit].box().contains(iv))
			{
			  updateCalvedIce(a_thickness[dit](iv),prevThck,levelMask[dit](iv),
				      a_addedIce[dit](iv),a_calvedIce[dit](iv),a_removedIce[dit](iv));
			}

		    }
		} 
	    } // end if (!domain.isPeriodic(dir))
	} // end loop over dirs
      
      const BaseFab<int>& mask = levelMask[dit];
      FArrayBox& thck = a_thickness[dit];
      FArrayBox& calved = a_calvedIce[dit];
      FArrayBox& added = a_addedIce[dit];
      FArrayBox& removed = a_removedIce[dit];
      const Box& b = grids[dit];
      for (BoxIterator bit(b); bit.ok(); ++bit)
	{
	  const IntVect& iv = bit();
	  Real prevThck = thck(iv);
	  if (m_preserveSea && mask(iv) == OPENSEAMASKVAL)
	    {
	      thck(iv) = 0.0;
	    }
	  else if (m_preserveLand && mask(iv) == OPENLANDMASKVAL)
	    {
	      thck(iv) = 0.0;
	    }
	  thck(iv) = std::max(thck(iv),0.0);

	  // Record gain/loss of ice
	  if (calved.box().contains(iv))
	    {
	      updateCalvedIce(thck(iv),prevThck,mask(iv),added(iv),calved(iv),removed(iv));
	    }

	}

    } // end loop over boxes

}

void ProximityCalvingModel::applyCriterion
(LevelData<FArrayBox>& a_thickness,
 LevelData<FArrayBox>& a_calvedIce,
 LevelData<FArrayBox>& a_addedIce,
 LevelData<FArrayBox>& a_removedIce,  
 LevelData<FArrayBox>& a_iceFrac, 
 const AmrIce& a_amrIce,
 int a_level,
 Stage a_stage)
{

  Real time = a_amrIce.time();
  bool calvingActive = (time >= m_startTime && time < m_endTime);
  calvingActive = false;
  pout() << " time = " << time 
	 << " m_startTime = " <<  m_startTime
	 << " m_endTime = " <<  m_endTime
	 << "calvingActive = " << calvingActive
	 << std::endl;
  if (true || calvingActive)
    {
      const LevelSigmaCS& levelCoords = *a_amrIce.geometry(a_level);
      const LevelData<FArrayBox>& proximity = *a_amrIce.groundingLineProximity(a_level);
      const LevelData<FArrayBox>& velocity = *a_amrIce.velocity(a_level);
      for (DataIterator dit(levelCoords.grids()); dit.ok(); ++dit)
	{
	  const BaseFab<int>& mask = levelCoords.getFloatingMask()[dit];
	  FArrayBox& thck = a_thickness[dit];
	  FArrayBox& calved = a_calvedIce[dit];
	  FArrayBox& added = a_addedIce[dit];
	  FArrayBox& removed = a_removedIce[dit];
	  const FArrayBox& prox = proximity[dit];
	  const FArrayBox& vel = velocity[dit];
	  Box b = thck.box();b &= prox.box();
	  for (BoxIterator bit(b); bit.ok(); ++bit)
	    {
	      const IntVect& iv = bit();
	      Real prevThck = thck(iv);
	      Real vmod = std::sqrt(vel(iv,0)*vel(iv,0) + vel(iv,1)*vel(iv,1));
	      if (prox(iv) < m_proximity && calvingActive && vmod > m_velocity)
		{
		  //thck(iv) *= 0.5; thck(iv) = max(thck(iv),10.0);
		  thck(iv) = 0.0;
		}
	      if (mask(iv) == OPENSEAMASKVAL)
		{
		   thck(iv) = 0.0;
		}
	      if (mask(iv) == FLOATINGMASKVAL)
		{
		  thck(iv) = max(thck(iv),1.0);
		}

	      // Record gain/loss of ice
	      if (calved.box().contains(iv))
		{
		  updateCalvedIce(thck(iv),prevThck,mask(iv),added(iv),calved(iv),removed(iv));
		}

	    }
	}
    }
}



CalvingModel* CalvingModel::parseCalvingModel(const char* a_prefix)
{

  CalvingModel* ptr = NULL;
  std::string type = "";
  ParmParse pp(a_prefix);
  pp.query("type",type);
  
  if (type == "NoCalvingModel")
    {
      ptr = new NoCalvingModel;
    }
  else if (type == "DomainEdgeCalvingModel")
    {
      Vector<int> frontLo(2,false); 
      pp.getarr("front_lo",frontLo,0,frontLo.size());
      Vector<int> frontHi(2,false);
      pp.getarr("front_hi",frontHi,0,frontHi.size());
      bool preserveSea = false;
      pp.query("preserveSea",preserveSea);
      bool preserveLand = false;
      pp.query("preserveLand",preserveLand);
      ptr = new DomainEdgeCalvingModel
	(frontLo, frontHi,preserveSea,preserveLand);
    }
  else if (type == "FixedFrontCalvingModel")
    {
      Real minThickness = 0.0;
      pp.get("min_thickness", minThickness );
      ptr = new DeglaciationCalvingModelA
	(0.0,  1.0e+10, minThickness, -1.2345678e+300, 1.2345678e+300);
    }
  else if (type == "DeglaciationCalvingModelA")
    {  
      Real minThickness = 0.0;
      pp.get("min_thickness", minThickness );
      Real calvingThickness = 0.0;
      pp.get("calving_thickness", calvingThickness );
      Real calvingDepth = 0.0;
      pp.get("calving_depth", calvingDepth );
      Real startTime = -1.2345678e+300;
      pp.query("start_time",  startTime);
      Real endTime = 1.2345678e+300;
      pp.query("end_time",  endTime);
      ptr = new DeglaciationCalvingModelA
	(calvingThickness,  calvingDepth, minThickness, startTime, endTime); 
    }
  else if (type == "DeglaciationCalvingModelB")
    {  
      Real minThickness = 0.0;
      pp.get("min_thickness", minThickness );
      Real calvingThickness = 0.0;
      pp.get("calving_thickness", calvingThickness );
      Real calvingDepth = 0.0;
      pp.query("calving_depth", calvingDepth );
      Real startTime = -1.2345678e+300;
      pp.query("start_time",  startTime);
      Real endTime = 1.2345678e+300;
      pp.query("end_time",  endTime);
      ptr = new DeglaciationCalvingModelB
	(calvingThickness,  calvingDepth, minThickness, startTime, endTime); 
    }
  else if (type == "ProximityCalvingModel")
    {
      Real proximity = 0.0;
      pp.get("proximity", proximity );
      Real velocity = 0.0;
      pp.query("velocity", velocity );
      Real startTime = -1.2345678e+300;
      pp.get("startTime",  startTime);
      Real endTime = 1.2345678e+300;
      pp.get("endTime",  endTime);
      ptr = new ProximityCalvingModel(proximity,velocity, startTime, endTime);
    }
  else if (type == "FlotationCalvingModel")
    {
      Vector<int> frontLo(2,false); 
      pp.getarr("front_lo",frontLo,0,frontLo.size());
      Vector<int> frontHi(2,false);
      pp.getarr("front_hi",frontHi,0,frontHi.size());
      bool preserveSea = false;
      pp.query("preserveSea",preserveSea);
      bool preserveLand = false;
      pp.query("preserveLand",preserveLand);
      ptr = new FlotationCalvingModel
	(frontLo, frontHi,preserveSea,preserveLand);
    }
  else if (type == "BennCalvingModel")
    {
      ptr = new BennCalvingModel(pp);
    }
  else if (type == "VanDerVeenCalvingModel")
    {
      ptr = new VdVCalvingModel(pp);
    }
  else if (type == "ThicknessCalvingModel")
    {  
      Real minThickness = 0.0;
      pp.get("min_thickness", minThickness );
      Real calvingThickness = 0.0;
      pp.get("calving_thickness", calvingThickness );
      Real calvingDepth = 0.0;
      pp.query("calving_depth", calvingDepth );
      Real startTime = -1.2345678e+300;
      pp.query("start_time",  startTime);
      Real endTime = 1.2345678e+300;
      pp.query("end_time",  endTime);
      ptr = new ThicknessCalvingModel
	(calvingThickness,  calvingDepth, minThickness, startTime, endTime); 
    }
  else if (type == "CliffCollapseCalvingModel")
    {  
      Real maxCliffHeight = 100.0;
      pp.get("max_cliff_height", maxCliffHeight);
      Real recessionRate = 0.0;
      pp.get("recession_rate", recessionRate );
      Real startTime = -1.2345678e+300;
      pp.query("start_time",  startTime);
      Real endTime = 1.2345678e+300;
      pp.query("end_time",  endTime);
      ptr = new CliffCollapseCalvingModel
	(maxCliffHeight, recessionRate, startTime, endTime); 
    }
  else if (type == "MaxiumumExtentCalvingModel")  
    {
      Real startTime = -1.2345678e+300;
      pp.query("start_time",  startTime);
      Real endTime = 1.2345678e+300;
      pp.query("end_time",  endTime);

      Vector<Real> vect(SpaceDim,0.0);

      pp.getarr("lowLoc",vect,0,SpaceDim);
      RealVect lowLoc(D_DECL(vect[0], vect[1],vect[2]));      

      pp.getarr("highLoc",vect,0,SpaceDim);
      RealVect highLoc(D_DECL(vect[0], vect[1],vect[2]));      

      MaximumExtentCalvingModel* Ptr = new MaximumExtentCalvingModel(highLoc,
                                                                     lowLoc,
                                                                     startTime,
                                                                     endTime);
      ptr = static_cast<CalvingModel*>(Ptr);

    }
  else if (type == "MaskedCalvingModel")
    {
      Real minThickness = 0.0;
      pp.get("min_thickness", minThickness );

      // masked calving model uses a surfaceFlux as a mask
      std::string mask_prefix(a_prefix);
      mask_prefix += ".mask";
      SurfaceFlux* mask_ptr = SurfaceFlux::parse(mask_prefix.c_str());

      MaskedCalvingModel* Ptr = new MaskedCalvingModel(mask_ptr, minThickness);

      ptr = static_cast<CalvingModel*>(Ptr);

      // MaskedCalvingModel makes a copy of the mask, so clean up here
      if (mask_ptr != NULL)
        {
          delete mask_ptr;
        }
    }
  else if (type == "CompositeCalvingModel")
    {
      int nElements;
      pp.get("nElements",nElements);
     
      std::string elementPrefix(a_prefix);
      elementPrefix += ".element";

      Vector<CalvingModel*> elements(nElements);
      for (int i = 0; i < nElements; i++)
        {
          std::string prefix(elementPrefix);
          char s[32];
          sprintf(s,"%i",i);
          prefix += s;
          ParmParse pe(prefix.c_str());
          elements[i] = parseCalvingModel(prefix.c_str());
          CH_assert(elements[i] != NULL);
        }
      CompositeCalvingModel* compositePtr = new CompositeCalvingModel(elements);
      ptr = static_cast<CalvingModel*>(compositePtr);
    }
  
  else if (type == "VariableRateCalvingModel")
    {
      ptr = new VariableRateCalvingModel(pp);
    }

   else if (type == "RateProportionalToSpeedCalvingModel")
    {
      ptr = new RateProportionalToSpeedCalvingModel(pp);
    }
  

  return ptr;
}


void 
DeglaciationCalvingModelB::applyCriterion
(LevelData<FArrayBox>& a_thickness,
 LevelData<FArrayBox>& a_calvedIce,
 LevelData<FArrayBox>& a_addedIce,
 LevelData<FArrayBox>& a_removedIce,  
 LevelData<FArrayBox>& a_iceFrac, 
 const AmrIce& a_amrIce,
 int a_level,
 Stage a_stage)
{
  
  const LevelSigmaCS& levelCoords = *a_amrIce.geometry(a_level);
  
  for (DataIterator dit(levelCoords.grids()); dit.ok(); ++dit)
    {
      const BaseFab<int>& mask = levelCoords.getFloatingMask()[dit];
      FArrayBox& thck = a_thickness[dit];
      FArrayBox& calved = a_calvedIce[dit];
      FArrayBox& added = a_addedIce[dit];
      FArrayBox& removed = a_removedIce[dit];
      Box b = thck.box();
      
      for (BoxIterator bit(b); bit.ok(); ++bit)
	{
	  const IntVect& iv = bit();
	  Real prevThck = thck(iv);
	  if (mask(iv) == OPENSEAMASKVAL)
	    {
	      thck(iv) = 0.0;
	    }
	  else if (mask(iv) == OPENLANDMASKVAL)
	    {
	      thck(iv) = 0.0;
	    }
          else if ((mask(iv) == FLOATINGMASKVAL) && (thck(iv) < m_calvingThickness))
            {
	      thck(iv) = m_minThickness;              
            }
	  else
	    {
	      thck(iv) = std::max(thck(iv),m_minThickness);
	    }

	  // Record gain/loss of ice
	  if (calved.box().contains(iv))
	    {
	      updateCalvedIce(thck(iv),prevThck,mask(iv),added(iv),calved(iv),removed(iv));
	    }
	}
    }
}


void 
ThicknessCalvingModel::applyCriterion
(LevelData<FArrayBox>& a_thickness,
 LevelData<FArrayBox>& a_calvedIce,
 LevelData<FArrayBox>& a_addedIce,
 LevelData<FArrayBox>& a_removedIce,  
 LevelData<FArrayBox>& a_iceFrac, 
 const AmrIce& a_amrIce,
 int a_level,
 Stage a_stage)
{
  
  const LevelSigmaCS& levelCoords = *a_amrIce.geometry(a_level);
  for (DataIterator dit(levelCoords.grids()); dit.ok(); ++dit)
    {
      const BaseFab<int>& mask = levelCoords.getFloatingMask()[dit];
      FArrayBox& iceFrac = a_iceFrac[dit];
      FArrayBox& thck = a_thickness[dit];
      FArrayBox& calved = a_calvedIce[dit];
      FArrayBox& added = a_addedIce[dit];
      FArrayBox& removed = a_removedIce[dit];
      FArrayBox effectiveThickness(thck.box(), 1);
      effectiveThickness.copy(thck);
      Box b = thck.box();
      b &= iceFrac.box();
      
      for (BoxIterator bit(b); bit.ok(); ++bit)
	{
	  const IntVect& iv = bit();          
          // if iceFrac > 0, then rescale effectiveThickness
          // by dividing by iceFrac value, which gives "actual" thickness
          // in the partial cell. Probably eventually want to move this to 
          // fortran
	  Real prevThck = thck(iv);
          if (iceFrac(iv,0) > 0.0)
            {
              effectiveThickness(iv,0) /= iceFrac(iv,0);
            }
            
          if (mask(iv) == OPENLANDMASKVAL)
	    {
	      thck(iv) = 0.0;
	    }
          // allow ice to spread into open sea regions too, if appropriate
          else if (((mask(iv) == FLOATINGMASKVAL) || (mask(iv) == OPENSEAMASKVAL))
                   && (effectiveThickness(iv) < m_calvingThickness))
            {
              // note that we're setting thck here, not effectiveThickness, 
              // which is a temporary
              // also set the iceFrac to zero in these cells
	      thck(iv) = m_minThickness; 
              iceFrac(iv,0) = 0.0;
            }
	  else
	    {
	      thck(iv) = std::max(thck(iv),m_minThickness);
	    }

	  // Record gain/loss of ice
	  if (calved.box().contains(iv))
	    {
	      updateCalvedIce(thck(iv),prevThck,mask(iv),added(iv),calved(iv),removed(iv));
	    }

	}
    }
}



  
//alter the thickness field at the end of a time step
void
MaximumExtentCalvingModel::applyCriterion(LevelData<FArrayBox>& a_thickness,
					  LevelData<FArrayBox>& a_calvedIce,
					  LevelData<FArrayBox>& a_addedIce,
					  LevelData<FArrayBox>& a_removedIce,  
					  LevelData<FArrayBox>& a_iceFrac, 
					  const AmrIce& a_amrIce,
					  int a_level,
					  Stage a_stage)
{

  const LevelSigmaCS& levelCoords = *a_amrIce.geometry(a_level);
  const Real dx = a_amrIce.amrDx()[a_level];
  
  for (DataIterator dit(levelCoords.grids()); dit.ok(); ++dit)
    {
      const BaseFab<int>& mask = levelCoords.getFloatingMask()[dit];
      FArrayBox& thck = a_thickness[dit];
      FArrayBox& calved = a_calvedIce[dit];
      FArrayBox& added = a_addedIce[dit];
      FArrayBox& removed = a_removedIce[dit];
      Box b = thck.box();
      
      for (BoxIterator bit(b); bit.ok(); ++bit)
	{
	  const IntVect& iv = bit();
	  Real prevThck = thck(iv);
          // compute location of cell center
          RealVect loc(iv);          
          loc += 0.5*RealVect::Unit;
          loc *= dx;
          
          // check high and low extents
          if ((mask(iv) == FLOATINGMASKVAL) || (mask(iv) == OPENSEAMASKVAL))
            {
              if (loc[0] <= m_lowLoc[0])
                {
                  thck(iv) = 0.0;
                }
              else if (loc[1] <= m_lowLoc[1])
                {
                  thck(iv) = 0.0;
                }
              else if (loc[0] > m_highLoc[0]) 
                {
                  thck(iv) = 0.0;
                }
              else if (loc[1] > m_highLoc[1])
                {
                  thck(iv) = 0.0;
                }
            } // end if floating or opensea

	  // Record gain/loss of ice
	  if (calved.box().contains(iv))
	    {
	      updateCalvedIce(thck(iv),prevThck,mask(iv),added(iv),calved(iv),removed(iv));
	    }

        } // end loop over cells
  
    }

}




//alter the thickness field at the end of a time step
void 
CompositeCalvingModel::applyCriterion(LevelData<FArrayBox>& a_thickness, 
				      LevelData<FArrayBox>& a_calvedIce,
				      LevelData<FArrayBox>& a_addedIce,
				      LevelData<FArrayBox>& a_removedIce, 
				      LevelData<FArrayBox>& a_iceFrac, 
				      const AmrIce& a_amrIce,
				      int a_level,
				      Stage a_stage)
{
  for (int n=0; n<m_vectModels.size(); n++)
    {
      m_vectModels[n]->applyCriterion( a_thickness, a_calvedIce, a_addedIce, a_removedIce, a_iceFrac,a_amrIce, a_level, a_stage);
    }
}

  
CompositeCalvingModel::~CompositeCalvingModel()
{
  for (int n=0; n<m_vectModels.size(); n++)
    {
      delete m_vectModels[n];
      m_vectModels[n] = NULL;
    }
}

void FlotationCalvingModel::applyCriterion
(LevelData<FArrayBox>& a_thickness,
 LevelData<FArrayBox>& a_calvedIce,
 LevelData<FArrayBox>& a_addedIce,
 LevelData<FArrayBox>& a_removedIce,  
 LevelData<FArrayBox>& a_iceFrac, 
 const AmrIce& a_amrIce,
 int a_level,
 Stage a_stage)
{

  m_domainEdgeCalvingModel.applyCriterion( a_thickness, a_calvedIce, a_addedIce, a_removedIce, a_iceFrac,a_amrIce, a_level, a_stage);
  const LevelSigmaCS& levelCoords = *a_amrIce.geometry(a_level);
  for (DataIterator dit(levelCoords.grids()); dit.ok(); ++dit)
    {
      FArrayBox& thck = a_thickness[dit];
      FArrayBox& calved = a_calvedIce[dit];
      FArrayBox& added = a_addedIce[dit];
      FArrayBox& removed = a_removedIce[dit];
      const BaseFab<int>& mask = levelCoords.getFloatingMask()[dit];
      const Box& b = levelCoords.grids()[dit];
      for (BoxIterator bit(b); bit.ok(); ++bit)
	{
	  const IntVect& iv = bit();
	  Real prevThck = thck(iv);
	  if (mask(iv) == FLOATINGMASKVAL)
	    {
	      thck(iv) = 0.0; 
	    }

	  // Record gain/loss of ice
	  if (calved.box().contains(iv))
	    {
	      updateCalvedIce(thck(iv),prevThck,mask(iv),added(iv),calved(iv),removed(iv));
	    }

	}
    }
}

void
CalvingModel::updateCalvedIce(const Real& a_thck, const Real a_prevThck, const int a_mask, Real& a_added, Real& a_calved, Real& a_removed)
{

  if (a_thck > a_prevThck)
    {
      a_added += (a_prevThck-a_thck);
    }
  else 
    {
      if ((a_mask == OPENSEAMASKVAL) || (a_mask == FLOATINGMASKVAL))
	{
	  a_calved += (a_prevThck-a_thck);
	}
      else
	{
	  a_removed += (a_prevThck-a_thck);
	}
    } 

}


VariableRateCalvingModel::VariableRateCalvingModel(ParmParse& a_pp)
{
      Real startTime = -1.2345678e+300;
      a_pp.query("start_time",  startTime);
      Real endTime = 1.2345678e+300;
      a_pp.query("end_time",  endTime);
 
      Vector<int> frontLo(2,false); 
      a_pp.getarr("front_lo",frontLo,0,frontLo.size());
      Vector<int> frontHi(2,false);
      a_pp.getarr("front_hi",frontHi,0,frontHi.size());
      bool preserveSea = false;
      a_pp.query("preserveSea",preserveSea);
      bool preserveLand = false;
      a_pp.query("preserveLand",preserveLand);

      m_domainEdgeCalvingModel = new DomainEdgeCalvingModel(frontLo,frontHi,preserveSea,preserveLand);

      std::string prefix (a_pp.prefix());
      m_calvingRate = SurfaceFlux::parse( (prefix + ".CalvingRate").c_str());

}

void VariableRateCalvingModel::applyCriterion
(LevelData<FArrayBox>& a_thickness,
 LevelData<FArrayBox>& a_calvedIce,
 LevelData<FArrayBox>& a_addedIce,
 LevelData<FArrayBox>& a_removedIce,  
 LevelData<FArrayBox>& a_iceFrac, 
 const AmrIce& a_amrIce,
 int a_level,
 Stage a_stage)
{

  (*m_domainEdgeCalvingModel).applyCriterion( a_thickness, a_calvedIce, a_addedIce, a_removedIce, a_iceFrac,a_amrIce, a_level, a_stage);

  const LevelSigmaCS& levelCoords = *a_amrIce.geometry(a_level);
  for (DataIterator dit(levelCoords.grids()); dit.ok(); ++dit)
    {
      FArrayBox& thck = a_thickness[dit];
      FArrayBox& calved = a_calvedIce[dit];
      FArrayBox& added = a_addedIce[dit];
      FArrayBox& removed = a_removedIce[dit];
      FArrayBox& frac = a_iceFrac[dit];
      const BaseFab<int>& mask = levelCoords.getFloatingMask()[dit];
      Real frac_eps = 1.0e-6;
      const Box& b = levelCoords.grids()[dit];
      for (BoxIterator bit(b); bit.ok(); ++bit)
	{
	  const IntVect& iv = bit();
	  Real prevThck = thck(iv);
	  
	  if (frac(iv) < frac_eps)
	    {
	      thck(iv)=0.0;
	    }

	  // Record gain/loss of ice
	  if (calved.box().contains(iv))
	    {
	      updateCalvedIce(thck(iv),prevThck,mask(iv),added(iv),calved(iv),removed(iv));
	    }

	}
    }
}



CalvingModel* VariableRateCalvingModel::new_CalvingModel()
  {
    VariableRateCalvingModel* ptr = new VariableRateCalvingModel(*this);
    ptr->m_startTime = m_startTime;
    ptr->m_endTime = m_endTime;
    ptr->m_calvingRate = m_calvingRate->new_surfaceFlux();
    ptr->m_domainEdgeCalvingModel = new DomainEdgeCalvingModel(*m_domainEdgeCalvingModel);
    return ptr; 
  }

VariableRateCalvingModel::~VariableRateCalvingModel()
{

  if (m_domainEdgeCalvingModel != NULL)
    {
      delete m_domainEdgeCalvingModel;
      m_domainEdgeCalvingModel = NULL;
    }

  if (m_calvingRate != NULL)
    {
      delete m_calvingRate;
      m_calvingRate = NULL;
    }

}


void 
CliffCollapseCalvingModel::applyCriterion(LevelData<FArrayBox>& a_thickness,
					  LevelData<FArrayBox>& a_calvedIce,
					  LevelData<FArrayBox>& a_addedIce,
					  LevelData<FArrayBox>& a_removedIce,  
					  LevelData<FArrayBox>& a_iceFrac, 
					  const AmrIce& a_amrIce,
					  int a_level,
					  Stage a_stage)
{

  // only do this at the end of a timestep
  // (since that's the only time a time-integrated recession rate makes any sense)
  if (a_stage == PostThicknessAdvection)
    {
      Real dt = a_amrIce.dt();
      Real dx = a_amrIce.amrDx()[a_level];
      
      const LevelSigmaCS& levelCoords = *a_amrIce.geometry(a_level);
      const LevelData<FArrayBox>& surfaceHeight = levelCoords.getSurfaceHeight();
      
      for (DataIterator dit(levelCoords.grids()); dit.ok(); ++dit)
	{
	  const BaseFab<int>& mask = levelCoords.getFloatingMask()[dit];
	  FArrayBox& iceFrac = a_iceFrac[dit];
	  FArrayBox& thck = a_thickness[dit];
	  FArrayBox& calved = a_calvedIce[dit];
	  FArrayBox& added = a_addedIce[dit];
	  FArrayBox& removed = a_removedIce[dit];
	  const FArrayBox& surface = surfaceHeight[dit];
	  FArrayBox effectiveSurface(surface.box(),1);
	  effectiveSurface.copy(surface);
	  FArrayBox effectiveThickness(thck.box(), 1);
	  effectiveThickness.copy(thck);
	  Box b = thck.box();
	  b &= iceFrac.box();
	  
	  // keep track of which cells we've already done in order to avoid double-counting
	  BaseFab<int> alreadyDone(b,1);
	  alreadyDone.setVal(0);
	  
	  Real phiNew;
	  
	  for (BoxIterator bit(b); bit.ok(); ++bit)
	    {
	      const IntVect& iv = bit();          
	      // if iceFrac > 0, then rescale effectiveThickness
	      // by dividing by iceFrac value, which gives "actual" thickness
	      // in the partial cell. Probably eventually want to move this to 
	      // fortran
	      // also compute "effective surface", which is the upper surface height based
	      // on the effective thickness rather than the cell-averaged thickness
	      // Probably eventually want to move this to fortran
	      Real prevThck = thck(iv);
	      if (iceFrac(iv,0) > 0.0)
		{
		  effectiveThickness(iv,0) /= iceFrac(iv,0);
		  effectiveSurface(iv,0) += effectiveThickness(iv,0) - thck(iv,0);	      
		}
	      
	      // if ice is grounded, look at neighbors to see if there are any empty neighbors, then look at
	      // surface differences.
	      if (mask(iv) == GROUNDEDMASKVAL) 
		{
		  // loop over directions
		  for (int dir=0; dir<SpaceDim; dir++)
		    {
		      IntVect shiftVect = BASISV(dir);
		      IntVect ivp = iv + shiftVect;
		      IntVect ivm = iv - shiftVect;
		      
		      // look in both high and low directions at once
		      if (((mask(ivp,0) != GROUNDEDMASKVAL) && (mask(ivp,0) != OPENLANDMASKVAL) &&  ((effectiveSurface(iv,0) - effectiveSurface(ivp,0)) > m_maxCliffThickness)) ||
			  ((mask(ivm,0) != GROUNDEDMASKVAL) && (mask(ivm,0) != OPENLANDMASKVAL) && ((effectiveSurface(iv,0) - effectiveSurface(ivm,0)) > m_maxCliffThickness)))
			{
			  // we have a cliff!  only adjust this cell if we haven't already
			  if (alreadyDone(iv,0) == 0)
			    {
			      alreadyDone(iv,0) = 1;
			      phiNew = iceFrac(iv,0) - m_recessionRate*dt/dx;
			      // don't go below zero
			      phiNew = Max(phiNew, 0.0);
			      
			      // note that we're setting thck here, not effectiveThickness, 
			      // which is a temporary
			      // also modify the iceMask to zero in these cells
			      
			      thck(iv,0) = thck(iv,0)*phiNew/iceFrac(iv,0);
			      iceFrac(iv,0) = phiNew;
			    } // end we haven't already done this one
			} // end if we have a cliff
		    } // end loop over directions
		} // end if this cell is grounded (no floating cliffs)
	      
	      
	      // Record gain/loss of ice
	      if (calved.box().contains(iv))
		{
		  updateCalvedIce(thck(iv),prevThck,mask(iv),added(iv),calved(iv),removed(iv));
		}
	      
	    } // end loop over cells in this box
	} // end loop over grids on this level	  
    } // end if we're at the post-advection stage
}


RateProportionalToSpeedCalvingModel::RateProportionalToSpeedCalvingModel(ParmParse& a_pp)
{
      Real startTime = -1.2345678e+300;
      a_pp.query("start_time",  startTime);
      Real endTime = 1.2345678e+300;
      a_pp.query("end_time",  endTime);
 
      Vector<int> frontLo(2,false); 
      a_pp.getarr("front_lo",frontLo,0,frontLo.size());
      Vector<int> frontHi(2,false);
      a_pp.getarr("front_hi",frontHi,0,frontHi.size());
      bool preserveSea = false;
      a_pp.query("preserveSea",preserveSea);
      bool preserveLand = false;
      a_pp.query("preserveLand",preserveLand);

      m_domainEdgeCalvingModel = new DomainEdgeCalvingModel(frontLo,frontHi,preserveSea,preserveLand);

      std::string prefix (a_pp.prefix());
      m_proportion = SurfaceFlux::parse( (prefix + ".proportion").c_str());

      m_normalize_vel_to_flowline_terminus_vel = false;
      a_pp.query("normalize_vel_to_flowline_terminus_vel", m_normalize_vel_to_flowline_terminus_vel);
      if (m_normalize_vel_to_flowline_terminus_vel)
	{
	  m_flowline_terminus_speed_update_interval = 0.25;
	  m_time_last = -1.0;
	  m_level_last = -2;
	  a_pp.query("flowline_terminus_speed_update_interval", m_flowline_terminus_speed_update_interval);
	  m_flowline_terminus_speed_last_update = -1.2345678e+300;
	  m_max_flowline_terminus_speed = 2.0e+4;
	  m_min_flowline_terminus_speed = 5.0e+3;
	  m_uT = 0.0;
	  a_pp.query("max_flowline_terminus_speed",m_max_flowline_terminus_speed);
	  a_pp.query("min_flowline_terminus_speed",m_min_flowline_terminus_speed);
      m_measure_along_flowline = false;
      a_pp.query("measure_along_flowline", m_measure_along_flowline);
	  m_update_every_timestep = false;	  
      a_pp.query("update_every_timestep", m_update_every_timestep);
	  
	  //need to read a space separated file of flowline co-ordinates
	  // format needs to be x\sy\n
	  std::string file_name;
	  a_pp.get("flowline_file",file_name);
	  std::ifstream infile(file_name);
	  Real x,y;
	  while (infile >> x >> y)
	    {
	      flowline_xy.push_back(std::make_pair(x,y));
	    }
	  
      // MJT - new module to calculate required calving rate to hit a prescribed cfl
      m_prescribe_cfl = true;	  
      a_pp.query("prescribe_cfl", m_prescribe_cfl);
	  if (m_prescribe_cfl) {
	    // Get the file containing cfl times and positions
	    std::string cfl_file;
	    a_pp.get("cfl_file",cfl_file);
	    std::ifstream cflfile(cfl_file);
	    Real t,cfl;
    	while (cflfile >> t >> cfl)
	      {
	        calvingFront_time_locn.push_back(std::make_pair(t,cfl));
	      }	    
	  }
	  
      // MJT - use a prescribed calving rate
      m_prescribe_calving_rate = false;	  
      a_pp.query("prescribe_calving_rate", m_prescribe_calving_rate);
	  if (m_prescribe_calving_rate) {
	    // Get the file containing cfl times and positions
	    std::string calving_rate_file;
	    a_pp.get("calving_rate_file",calving_rate_file);
	    std::ifstream ratefile(calving_rate_file);
	    Real t,rate;
    	while (ratefile >> t >> rate)
	      {
	        calvingRate_time_rate.push_back(std::make_pair(t,rate));
	      }	    
	  }
	}
}

void RateProportionalToSpeedCalvingModel::applyCriterion
(LevelData<FArrayBox>& a_thickness,
 LevelData<FArrayBox>& a_calvedIce,
 LevelData<FArrayBox>& a_addedIce,
 LevelData<FArrayBox>& a_removedIce,  
 LevelData<FArrayBox>& a_iceFrac, 
 const AmrIce& a_amrIce,
 int a_level,
 Stage a_stage)
{

  (*m_domainEdgeCalvingModel).applyCriterion( a_thickness, a_calvedIce, a_addedIce, a_removedIce, a_iceFrac,a_amrIce, a_level, a_stage);

  const LevelSigmaCS& levelCoords = *a_amrIce.geometry(a_level);
  for (DataIterator dit(levelCoords.grids()); dit.ok(); ++dit)
    {
      FArrayBox& thck = a_thickness[dit];
      FArrayBox& calved = a_calvedIce[dit];
      FArrayBox& added = a_addedIce[dit];
      FArrayBox& removed = a_removedIce[dit];
      FArrayBox& frac = a_iceFrac[dit];
      const BaseFab<int>& mask = levelCoords.getFloatingMask()[dit];
      Real frac_eps = 1.0e-6;
      const Box& b = levelCoords.grids()[dit];
      for (BoxIterator bit(b); bit.ok(); ++bit)
	{
	  const IntVect& iv = bit();
	  Real prevThck = thck(iv);
	  
	  if (frac(iv) < frac_eps)
	    {
	      thck(iv)=0.0;
	    }

	  // Record gain/loss of ice
	  if (calved.box().contains(iv))
	    {
	      updateCalvedIce(thck(iv),prevThck,mask(iv),added(iv),calved(iv),removed(iv));
	    }

	}
    }
	
  a_thickness.exchange();
  a_iceFrac.exchange();

  /*if (a_stage == CalvingModel::PostVelocitySolve)
    {
      // now seems like a good time to update the terminus velocity
      if (m_normalize_vel_to_flowline_terminus_vel)
	{
	  const Real& t = a_amrIce.time();
	  if (m_update_every_timestep)
	    {
		  m_flowline_terminus_speed_update_interval = a_amrIce.dt();
		}  
	  if ( ( t - m_flowline_terminus_speed_last_update) >= m_flowline_terminus_speed_update_interval )
	    {
	      m_flowline_terminus_speed_last_update = t;
		  if (m_prescribe_cfl)
		    {
			  //pair<Real,Real> alongFlowTerminusSpeed_output =  alongFlowTerminusSpeed(a_amrIce);
			  //m_flowline_terminus_speed = alongFlowTerminusSpeed_output.first;
			  //m_flowline_terminus_distance = alongFlowTerminusSpeed_output.second;		
			}
		else
	  	  {
		    if (m_measure_along_flowline)
		      {
	            //pair<Real,Real> alongFlowTerminusSpeed_output =  alongFlowTerminusSpeed(a_amrIce);
			    //m_flowline_terminus_speed = alongFlowTerminusSpeed_output.first;
			    //m_flowline_terminus_distance = alongFlowTerminusSpeed_output.second;
	  		  }
		    else
		      {
	            //m_flowline_terminus_speed = flowlineTerminusSpeed(a_amrIce);
			  }
		  }
	    }
	}
    }*/
}



CalvingModel* RateProportionalToSpeedCalvingModel::new_CalvingModel()
  {
    RateProportionalToSpeedCalvingModel* ptr = new RateProportionalToSpeedCalvingModel(*this);
    ptr->m_startTime = m_startTime;
    ptr->m_endTime = m_endTime;
    ptr->m_proportion = m_proportion->new_surfaceFlux();
    ptr->m_domainEdgeCalvingModel = new DomainEdgeCalvingModel(*m_domainEdgeCalvingModel);
    return ptr; 
  }

RateProportionalToSpeedCalvingModel::~RateProportionalToSpeedCalvingModel()
{

  if (m_domainEdgeCalvingModel != NULL)
    {
      delete m_domainEdgeCalvingModel;
      m_domainEdgeCalvingModel = NULL;
    }

  if (m_proportion != NULL)
    {
      delete m_proportion;
      m_proportion = NULL;
    }

}



Real
RateProportionalToSpeedCalvingModel::flowlineTerminusSpeed(const AmrIce& a_amrIce)
{

  /// we want to find the speed at the point where the flowline intersects the calving front.
  /// MJT did this in a python function
  /// by stepping through the flowline points and choosing the velocity
  /// from the first point where h > 1.0e-2

  /// However, this is quite fiddly in AMR / parallel so 
  ///  SLC has decided to compute 
  /// v_front = integral ( v * mesh_dirac_delta (x-xf,y-yf) , dx*dy)
  /// mesh_dirac_delta(x,y) is the (normalised) W(x,y)*|laplacian(ice_frac)*dx * h|^2

  /// The obvious prescription for W(x,y) is to set
  /// W(x,y) = 1 where the cell centred on x,y contains a flowline point.
  /// this method might break down if the flowline points are sparser than the mesh
  /// More robust (but slower) is to set W(x,y) = exp ( - |(x-x_f, y-y_f)|^2 / scale^2 )
  /// and normalize
  
  /// MJT - I've updated this so that the weight function W(x,y) is based on distance
  /// to the nearest line segment, not distance to the nearest point.

  CH_TIME("RateProportionalToSpeedCalvingModel::flowlineTerminusSpeed");
  
  Real scale = 1.0e+2;
  //Real scale = a_amrIce.dx(0)[0];
  Real max_speed = 0.0;
  Vector<LevelData<FArrayBox>* > mesh_dirac_delta;
  Vector<LevelData<FArrayBox>* > speed_weighted;
  double inf = std::numeric_limits<double>::infinity();
  double max = std::numeric_limits<double>::max();

  //compute mesh_dirac_delta. Will need to normalize later
  //compute speed at the same time
  for (int lev = 0; lev <= a_amrIce.finestLevel(); lev++)
    {
      mesh_dirac_delta.push_back(new LevelData<FArrayBox>( a_amrIce.grids(lev), 1, IntVect::Zero));
      speed_weighted.push_back(new LevelData<FArrayBox>( a_amrIce.grids(lev), 1, IntVect::Zero));
      
      const RealVect& dx = a_amrIce.dx(lev);
      for (DataIterator dit( a_amrIce.grids(lev) ); dit.ok(); ++dit)
	{
	  const FArrayBox& f = (*a_amrIce.iceFrac(lev))[dit];
	  const FArrayBox& h = (*a_amrIce.geometry(lev)).getH()[dit];
	  const Box& box = a_amrIce.grids(lev)[dit];

	  // laplacian(ice_frac) (undivided, absolute)
	  FArrayBox lap_f(box,1);
	  lap_f.setVal(0.0);
	  Real alpha = 0;
      Real beta = 1.0;
	  Real bogusDx = 1.0;
	  FORT_OPERATORLAP(CHF_FRA(lap_f),
                           CHF_FRA(f),
                           CHF_BOX(box),
                           CHF_CONST_REAL(bogusDx),
                           CHF_CONST_REAL(alpha),
                           CHF_CONST_REAL(beta));
	  
	  lap_f *= lap_f; // keep it positive
	  lap_f *= h; // don't want any values of h from zero thickness cells.
	    
	  // W(x,y) 
	  FArrayBox flowlineweight(box,1);
	  flowlineweight.setVal(0.0);
	  Real sumflowlineweight = 0.0;
	  Real sumlapf = 0.0;
	  for ( BoxIterator bit(box); bit.ok(); ++bit)
	    {
	      const IntVect& iv = bit();
	    /*  for (int k = 0; k < flowline_xy.size(); k++)
		{
		  Real xiv = ( Real(iv[0]) + 0.5 )*dx[0];
		  Real yiv = ( Real(iv[1]) + 0.5 )*dx[1];
		  Real rx = (xiv - flowline_xy[k].first)/scale;
		  Real ry = (yiv - flowline_xy[k].second)/scale;
		  Real rsq = rx*rx + ry*ry;
		  flowlineweight(iv) += std::exp(-rsq);
		  sumflowlineweight += flowlineweight(iv);
		  if (lap_f(iv) > max)
		    { lap_f(iv) = 1.0; }
		  sumlapf += lap_f(iv);
		}*/
		
	      for (int k = 0; k < flowline_xy.size() - 1; k++)
		{
		  Real xiv = ( Real(iv[0]) + 0.5 )*dx[0];
		  Real yiv = ( Real(iv[1]) + 0.5 )*dx[1];
		  Real x0 = flowline_xy[k].first;
		  Real x1 = flowline_xy[k].first;
		  Real y0 = flowline_xy[k+1].second;
		  Real y1 = flowline_xy[k+1].second;
		  Real px = x1 - x0;
		  Real py = y1 - y0;
		  Real norm = px*px + py*py;
		  Real ua = ((xiv - x0) * px + (yiv - y0) * py) / norm;
		  ua = std::max(0.0,std::min(1.0,ua));
		  Real xl = x0 + ua*px;
		  Real yl = y0 + ua*py;
		  Real dist = std::sqrt((xiv-xl)*(xiv-xl) + (yiv-yl)*(yiv-yl)); 
		  flowlineweight(iv) = std::max(flowlineweight(iv),std::exp(-(dist*dist)/(scale*scale)));
		  if (lap_f(iv) > max) { 
		  lap_f(iv) = 1.0; 
		  pout() << "RateProportionalToSpeedCalvingModel::Infinity found" << std::endl;
		  }
		  sumlapf += lap_f(iv);
		}
	    }

	  // compute mesh_dirac_delta. Will need to normalize later.
	  (*mesh_dirac_delta[lev])[dit].copy(flowlineweight);
	  (*mesh_dirac_delta[lev])[dit] *= lap_f;
	  
      //pout() << "RateProportionalToSpeedCalvingModel::sumflowlineweight = " << sumflowlineweight << ", sumlapf = " << sumlapf << std::endl;

	  /// compute speed everywhere.
	  FArrayBox& uu = (*speed_weighted[lev])[dit];
	  uu.setVal(0.0);
	  const FArrayBox& u = (*a_amrIce.velocity(lev))[dit];
	  for ( BoxIterator bit(box); bit.ok(); ++bit)
	    {
	      const IntVect& iv = bit();
	      Real t = 0.0;
	      for (int dir = 0; dir < SpaceDim; dir++)
		{
		  t += u(iv,dir)*u(iv,dir);
		}
	      
	      uu(iv) = std::sqrt(t);
	    }
	  uu *=  (*mesh_dirac_delta[lev])[dit];
	}
    }
  
  Real denom = computeSum(mesh_dirac_delta,  a_amrIce.refRatios(), a_amrIce.dx(0)[0],  Interval(0,0), 0);
  Real num = computeSum(speed_weighted,  a_amrIce.refRatios(), a_amrIce.dx(0)[0],  Interval(0,0), 0);
  Real result = std::min(num/denom, m_max_flowline_terminus_speed);
  //pout() << "RateProportionalToSpeedCalvingModel::flowlineTerminusSpeed = " << result << std::endl;
  //clean up
  for (int lev = 0; lev <= a_amrIce.finestLevel(); lev++)
    {
      delete ( mesh_dirac_delta[lev] );  mesh_dirac_delta[lev] = NULL;
      delete ( speed_weighted[lev] );  speed_weighted[lev] = NULL;
    }

  return result;
}

std::pair<Real, Real>
RateProportionalToSpeedCalvingModel::alongFlowTerminusSpeed(const AmrIce& a_amrIce)
{

  /// CalvingModel.measure_along_flowline = true
  /// This essentially does the same as flowlineTerminusSpeed, but measures along the flowline
  /// like MJT's original Python function.
  
  /// For each flowline point, create a function that's 0 everywhere and 1 for the cell 
  /// containing the point. Then multiply by the ice fraction until we find the first cell
  /// containing ice.
  
  /// Update (27/11/2020) - Also return the distance of the cfl along the flowline

  CH_TIME("RateProportionalToSpeedCalvingModel::alongFlowTerminusSpeed");
  
  Real result = 1.0;  
  
  Real distanceAlongFlow = 0.0;
  
  for (int k = 0; k < flowline_xy.size(); k++)
    {  
	
	  // Calculate the along flowline distance
	  if (k > 0) {
	    Real delx = flowline_xy[k].first - flowline_xy[k-1].first;
	    Real dely = flowline_xy[k].second - flowline_xy[k-1].second;
		distanceAlongFlow += std::sqrt(delx*delx + dely*dely);
	  }
  
      //Real scale = 1.0e+2;
      Vector<LevelData<FArrayBox>* > mesh_dirac_delta;
      Vector<LevelData<FArrayBox>* > speed_weighted;
      Vector<LevelData<FArrayBox>* > mesh_frac;

      //compute mesh_dirac_delta. Will need to normalize later
      //compute speed at the same time
      for (int lev = 0; lev <= a_amrIce.finestLevel(); lev++)
        {
          mesh_dirac_delta.push_back(new LevelData<FArrayBox>( a_amrIce.grids(lev), 1, IntVect::Zero));
          speed_weighted.push_back(new LevelData<FArrayBox>( a_amrIce.grids(lev), 1, IntVect::Zero));
          mesh_frac.push_back(new LevelData<FArrayBox>( a_amrIce.grids(lev), 1, IntVect::Zero));
      
          const RealVect& dx = a_amrIce.dx(lev);
		  
	      (*mesh_dirac_delta[lev]).exchange();
	      (*mesh_frac[lev]).exchange();
	      (*speed_weighted[lev]).exchange();		  
		  
          for (DataIterator dit( a_amrIce.grids(lev) ); dit.ok(); ++dit)
        	{
        	  const FArrayBox& f = (*a_amrIce.iceFrac(lev))[dit];
        	  //const FArrayBox& h = (*a_amrIce.geometry(lev)).getH()[dit];
        	  const Box& box = a_amrIce.grids(lev)[dit];
	          FArrayBox ice_frac(box,1);
	          ice_frac.setVal(0.0);
	          ice_frac += f;

        	  // W(x,y) 
        	  FArrayBox flowlineweight(box,1);
        	  flowlineweight.setVal(0.0);
              Real scale = dx[0] / 0.01;
        	  for ( BoxIterator bit(box); bit.ok(); ++bit)
	            {
        	      const IntVect& iv = bit();
        		  Real xiv = ( Real(iv[0]) + 0.5 )*dx[0];
		          Real yiv = ( Real(iv[1]) + 0.5 )*dx[1];
		          Real rx = (xiv - flowline_xy[k].first)/scale;
		          Real ry = (yiv - flowline_xy[k].second)/scale;
		          Real rsq = rx*rx + ry*ry;
		          flowlineweight(iv) += std::exp(-rsq);
	            }

	          // compute mesh_dirac_delta. Will need to normalize later.
	          (*mesh_dirac_delta[lev])[dit].copy(flowlineweight);
	          (*mesh_frac[lev])[dit].copy(ice_frac);
	          (*mesh_frac[lev])[dit].copy(ice_frac);
	        }
	      (*mesh_dirac_delta[lev]).exchange();
	      (*mesh_frac[lev]).exchange();
	      (*speed_weighted[lev]).exchange();
        } // end level loop
	
      // Set the weight function to be 1 for the max, 0 everywhere else
      Real maxWeight = computeMax(mesh_dirac_delta, a_amrIce.refRatios(), Interval(0,0), 0);
      for (int lev = 0; lev <= a_amrIce.finestLevel(); lev++)
        {      
          for (DataIterator dit( a_amrIce.grids(lev) ); dit.ok(); ++dit)
	        {
	          const Box& box = a_amrIce.grids(lev)[dit];
	          for ( BoxIterator bit(box); bit.ok(); ++bit)
	            {
                  const IntVect& iv = bit();
		          if ( (*mesh_dirac_delta[lev])[dit](iv) < maxWeight)
		            {
                      (*mesh_dirac_delta[lev])[dit](iv) = 0.0;
		            }
		          else
		            {
		              (*mesh_dirac_delta[lev])[dit](iv) = 1.0;
		            }
		        }
	        }
	      (*mesh_dirac_delta[lev]).exchange();
	      (*mesh_frac[lev]).exchange();
	      (*speed_weighted[lev]).exchange();
	    }
	
      maxWeight = computeMax(mesh_dirac_delta, a_amrIce.refRatios(), Interval(0,0), 0);
      Real sumWeight = computeSum(mesh_dirac_delta, a_amrIce.refRatios(), a_amrIce.dx(0)[0], Interval(0,0), 0);
	  
	  // Multiply the ice fraction by the weight function
      for (int lev = 0; lev <= a_amrIce.finestLevel(); lev++)
        {      
	      for (DataIterator dit( a_amrIce.grids(lev) ); dit.ok(); ++dit)
	        {
	          FArrayBox& mf = (*mesh_frac[lev])[dit];
	          const FArrayBox& f = (*a_amrIce.iceFrac(lev))[dit];
	          mf.setVal(0.0);
	          mf += (*mesh_dirac_delta[lev])[dit];
	          mf *= f;
	        }
	      (*mesh_dirac_delta[lev]).exchange();
	      (*mesh_frac[lev]).exchange();
	      (*speed_weighted[lev]).exchange();
	    }
	  // Calculate the sum of ice_frac(x,y)*w(x,y)
      Real sumMeshFrac = computeSum(mesh_frac,  a_amrIce.refRatios(), a_amrIce.dx(0)[0],  Interval(0,0), 0);
	  
	  // Break out of the flowline loop the first time the sum > 0, after calculating the velocity
      //pout() << "RateProportionalToSpeedCalvingModel::flowlineIteration = " << k << ", distanceAlongFlow = " << distanceAlongFlow << ", sumWeight = " << sumWeight << ", sumMeshFrac = " << sumMeshFrac << ", alongFlowTerminusSpeed = " << result << std::endl;	
	  if ((sumMeshFrac/sumWeight) > 1.0e-3)
	    {
          for (int lev = 0; lev <= a_amrIce.finestLevel(); lev++)
            {      
              for (DataIterator dit( a_amrIce.grids(lev) ); dit.ok(); ++dit)
	            {
	              // compute speed everywhere.
            	  FArrayBox& uu = (*speed_weighted[lev])[dit];
            	  uu.setVal(0.0);
	              const FArrayBox& u = (*a_amrIce.velocity(lev))[dit];
	              const Box& box = a_amrIce.grids(lev)[dit];
	              for ( BoxIterator bit(box); bit.ok(); ++bit)
	                {
	                  const IntVect& iv = bit();
	                  Real t = 0.0;
	                  for (int dir = 0; dir < SpaceDim; dir++)
		                {
		                  t += u(iv,dir)*u(iv,dir);
		                }
	      
            	      uu(iv) = std::sqrt(t);
	                }
	              uu *=  (*mesh_dirac_delta[lev])[dit];
	            }
	      (*mesh_dirac_delta[lev]).exchange();
	      (*mesh_frac[lev]).exchange();
	      (*speed_weighted[lev]).exchange();
	        }
          Real denom = computeSum(mesh_dirac_delta,  a_amrIce.refRatios(), a_amrIce.dx(0)[0],  Interval(0,0), 0);
          Real num = computeSum(speed_weighted,  a_amrIce.refRatios(), a_amrIce.dx(0)[0],  Interval(0,0), 0);
          result = std::min(num/denom, m_max_flowline_terminus_speed);
	      distanceAlongFlow += (1.0-(sumMeshFrac/sumWeight)) * a_amrIce.dx(a_amrIce.finestLevel())[0];
          pout() << "RateProportionalToSpeedCalvingModel::flowlineIteration = " << k << ", sumWeight = " << sumWeight << ", sumMeshFrac = " << sumMeshFrac << ", alongFlowTerminusSpeed = " << result << std::endl;	
	      pout() << "RateProportionalToSpeedCalvingModel::alongFlowTerminusSpeed: flowlineIteration = " << k << ", distanceAlongFlow = " << distanceAlongFlow << std::endl;	
          //clean up
          for (int lev = 0; lev <= a_amrIce.finestLevel(); lev++)
            {
	      (*mesh_dirac_delta[lev]).exchange();
	      (*mesh_frac[lev]).exchange();
	      (*speed_weighted[lev]).exchange();
              delete ( mesh_dirac_delta[lev] );  mesh_dirac_delta[lev] = NULL;
              delete ( speed_weighted[lev] );  speed_weighted[lev] = NULL;
              delete ( mesh_frac[lev] );  mesh_frac[lev] = NULL;
            }		  
		  break;
		}

      //clean up
      for (int lev = 0; lev <= a_amrIce.finestLevel(); lev++)
        {
          delete ( mesh_dirac_delta[lev] );  mesh_dirac_delta[lev] = NULL;
          delete ( speed_weighted[lev] );  speed_weighted[lev] = NULL;
          delete ( mesh_frac[lev] );  mesh_frac[lev] = NULL;
        }
    } // End of the flowline loop

  /*struct return_values {
    Real terminus_vel, terminus_distance;
  };
  return return_values {result,distanceAlongFlow};*/
  std::pair <Real,Real> return_values (result,distanceAlongFlow);
  //return result;
  return return_values;
}

void
RateProportionalToSpeedCalvingModel::getCalvingRate
(LevelData<FArrayBox>& a_calvingRate, const AmrIce& a_amrIce,int a_level)
{
  //srand((unsigned) time(NULL));
  m_proportion->evaluate(a_calvingRate, a_amrIce, a_level, 0.0);

  // Calculate the factor
  Real factor = 1.0;
  const Real& t = a_amrIce.time();
  if (m_update_every_timestep)
    {
	  m_flowline_terminus_speed_update_interval = a_amrIce.dt();
	}  
  if ( ( t - m_flowline_terminus_speed_last_update) >= m_flowline_terminus_speed_update_interval )
    {
      m_flowline_terminus_speed_last_update = t;
  if (m_time_last != a_amrIce.time()) {
    m_time_last = a_amrIce.time();
    //pout() << "RateProportionalToSpeedCalvingModel::flowlineTerminusSpeed = " << m_flowline_terminus_speed << ", flowlineTerminusDistance = " << m_flowline_terminus_distance << std::endl;
    if (m_prescribe_cfl) {
      pair<Real,Real> alongFlowTerminusSpeed_output =  alongFlowTerminusSpeed(a_amrIce);
  	  m_flowline_terminus_speed = alongFlowTerminusSpeed_output.first;
	  m_flowline_terminus_distance = alongFlowTerminusSpeed_output.second;		
    }
    else {
	  if (m_measure_along_flowline) {
	    pair<Real,Real> alongFlowTerminusSpeed_output =  alongFlowTerminusSpeed(a_amrIce);
	    m_flowline_terminus_speed = alongFlowTerminusSpeed_output.first;
	    m_flowline_terminus_distance = alongFlowTerminusSpeed_output.second;
	  }
	  else {
	    m_flowline_terminus_speed = flowlineTerminusSpeed(a_amrIce);
	  }
    }
    pout() << "RateProportionalToSpeedCalvingModel::time = " << a_amrIce.time() << ", flowlineTerminusSpeed = " << m_flowline_terminus_speed << ", flowlineTerminusDistance = " << m_flowline_terminus_distance << std::endl;
  }
  //m_velocity[a_level]->exchange(); // MJT
  }
  
  if (m_prescribe_cfl) {
  // uc(x,y) = u(x,y)/uT * [uT + (cfl' - cfl(t))/(t'-t)]
  // cfl' and t' are the prescribed cfl and time
  // uT is the terminus flowline speed
  
  // Get t, cfl(t), cfl' and t' and uT
  Real t_now = a_amrIce.time();
  Real cfl_now = m_flowline_terminus_distance;
  Real t_next = 0.0;
  Real cfl_next = 0.0;
  Real uT = std::max(m_flowline_terminus_speed,1.0);
  uT = m_flowline_terminus_speed;
  
  for (int k = 0; k < calvingFront_time_locn.size(); k++)
    { 
	  if (t_now < calvingFront_time_locn[k].first) {
	    t_next = calvingFront_time_locn[k].first;
	    cfl_next = calvingFront_time_locn[k].second;
		break;
	  }
    }
  Real uC = (cfl_next - cfl_now)/(t_next - t_now) + uT;
  factor = uC / (uT + 1.0);
  //if (uT <= 1.0) {factor = 1.0;}			// For the first timestep, otherwise we get massive initial rate
  if (uT <= 1.0) 
    {
	  if (t_now == m_startTime) {factor = 1.0;}			// For the first timestep, otherwise we get massive initial rate
	  else {											// Otherwise, the flowline terminus speed function has messed up. Set it to whatever was measured at the last timestep
	    uT = m_uT;
	    factor = uC / (uT + 1.0);
	  }
	}
  if (factor <= 0.0) {factor = 0.0001;}		
  pout() << "RateProportionalToSpeedCalvingModel::t=" << t_now << ", t'=" << t_next << ", cfl(t)=" << cfl_now << ", cfl'=" << cfl_next << ", u_T=" << uT << ", u_C=" << uC << ", factor=" << factor << std::endl;	
  }
  else {
  if (m_normalize_vel_to_flowline_terminus_vel)
    {
      Real uT = m_flowline_terminus_speed;
      //pout() << "RateProportionalToSpeedCalvingModel::getCalvingRate level = " << a_level << " v = " << v << std::endl;
	  Real uC = 1.0;
      Real t_now = a_amrIce.time();
	  if (m_prescribe_calving_rate) {
        Real t_last = 0.0;
        Real t_next = 0.0;
        Real rate_last = 0.0;
        Real rate_next = 0.0;  
        for (int k = 0; k < calvingRate_time_rate.size(); k++) {
      	  if (t_now < calvingRate_time_rate[k].first) {
            t_last = calvingRate_time_rate[k-1].first;
            t_next = calvingRate_time_rate[k].first;
      	    rate_last = calvingRate_time_rate[k-1].second;
      	    rate_next = calvingRate_time_rate[k].second;
        	break;
          }
        }
		uC = rate_last + (rate_next-rate_last)*(t_now-t_last)/(t_next-t_last);
		pout() << "RateProportionalToSpeedCalvingModel::t=" << t_now << ", t'=" << t_next << ", uC=" << uC << ", rate_next'=" << rate_next << ", u_T=" << uT << ", u_C=" << uC << ", factor=" << uC / (uT + 1.0) << std::endl;	
 	  
	  }
      factor = uC / (uT + 1.0);
	  
	  
	  if (uT <= 1.0) 
		{
		  if (t_now == m_startTime) {factor = 1.0;}			// For the first timestep, otherwise we get massive initial rate
		  else {											// Otherwise, the flowline terminus speed function has messed up. Set it to whatever was measured at the last timestep
			uT = m_uT;
			factor = uC / (uT + 1.0);
		  }
		}
	  if (factor <= 0.0) {factor = 0.0001;}	
	  
    }
  }
  
  m_uT = m_flowline_terminus_speed;
  
   // set uc(x,y) = u(x,y) and multiply by the factor
  const LevelData<FArrayBox>& vel = *a_amrIce.velocity(a_level);
  bool outputOnce = false;
  //m_prescribe_cfl = true;	// Seems to be a memory leak here somewhere!
  for (DataIterator dit(vel.dataIterator()); dit.ok(); ++dit)
    {
      Box b = a_calvingRate[dit].box();
      for (BoxIterator bit(b); bit.ok(); ++bit) {
        const IntVect& iv = bit();
        Real usq = 0.0;
        for (int dir = 0; dir < SpaceDim; dir++) {
          usq += vel[dit](iv,dir)*vel[dit](iv,dir);
        }
        if (m_prescribe_cfl) {    
          a_calvingRate[dit](iv) = factor * (1.0 + std::sqrt(usq));
        }
        else {
          if (m_normalize_vel_to_flowline_terminus_vel) {
            if (m_prescribe_calving_rate) {
              a_calvingRate[dit](iv) = factor * (1.0 + std::sqrt(usq));
            }
            else {
              a_calvingRate[dit](iv) *= factor * (1.0 + std::sqrt(usq));
            }
          }
          else {
            a_calvingRate[dit](iv) *= (0.001 + std::sqrt(usq));
            //pout() << "RateProportionalToSpeedCalvingModel::time = " << a_amrIce.time() << ", uC = " << a_calvingRate[dit](iv) << ", u = " << std::sqrt(usq) << ", iv = " << iv << std::endl;
          }
        }
      }
    }
	//a_calvingRate.exchange();
	
  int dbg = 0;dbg++; 
}

#include "NamespaceFooter.H"
