/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 ******************************************************************************/

#pragma once

#include <string>
#include <vector>
#include <map>

namespace Manta
{

struct PbClassData;
class FluidSolver;



//! Base class for all classes exposed to Python
class PbClass
{
public:
	PbClass(FluidSolver* parent, const std::string& name="")
	    : mParent(parent)
        , mName(name)
    {
    }

	PbClass(const PbClass& a)
        : mParent(a.mParent)
        , mName("_unnamed")
    {
    }

	virtual ~PbClass()
    {
    }
	
	// basic property setter/getters
	void setName(const std::string& name)   { mName = name; }
	std::string getName() const             { return mName; }

	FluidSolver* getParent() const          { return mParent; }
	void setParent(FluidSolver* v)          { mParent = v; }
	void checkParent()
    {
	    if (getParent() == nullptr) errMsg("New class " + mName + ": no parent given -- specify using parent=xxx !");
    }
	
protected:
	FluidSolver* mParent;
	std::string  mName;
};

} // namespace        
