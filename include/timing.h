/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Plugin timing
 *
 ******************************************************************************/

#pragma once

#include "manta.h"
#include <map>

namespace Manta
{ 

class TimingData
{
private:
	TimingData();

public:
	static TimingData& instance() { static TimingData a; return a; }

	void print();
	void saveMean(const std::string& filename);
	void start(FluidSolver* parent, const std::string& name);
	void stop(FluidSolver* parent, const std::string& name);

protected:
	void step();

	struct TimingSet
    {
		TimingSet()
            : num(0)
            , updated(false) 
        {
            cur.clear();
            total.clear();
        }

		MuTime cur, total;
		int num;
		bool updated;
		std::string solver;
	};

	bool updated;

	int num;
	MuTime mPluginTimer;
	std::string mLastPlugin;
	std::map<std::string, std::vector<TimingSet> > mData;
};

// Python interface
class Timings : public PbClass
{
public:
	Timings() 
    : PbClass(0)
    {}
	
	void display() { TimingData::instance().print(); }
};

}
