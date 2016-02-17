/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2013 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * GNU General Public License (GPL) 
 * http://www.gnu.org/licenses
 *
 * Particle data functionality
 *
 ******************************************************************************/

#include <fstream>
#include  <cstring>
#include  <limits>
#if NO_ZLIB!=1
#include <zlib.h>
#endif
#include "particle.h"
#include "levelset.h"
#include "fileio.h"

using namespace std;

namespace Manta
{

//******************************************************************************
// ParticleBase members
//******************************************************************************
ParticleBase::ParticleBase(FluidSolver* parent)
    : PbClass(parent)
    , mAllowCompress(true)
    , mFreePdata(false)
{
}

ParticleBase::~ParticleBase()
{
	// make sure data fields now parent system is deleted
	for (int i = 0; i < (int)mPartData.size(); ++i)
		mPartData[i]->setParticleSys(NULL);
	
	if(mFreePdata)
    {
		for(int i=0; i<(int)mPartData.size(); ++i)
			delete mPartData[i];
	}
	
}

std::string ParticleBase::infoString() const
{ 
	return "ParticleSystem " + mName + " <no info>"; 
}

void ParticleBase::cloneParticleData(ParticleBase* nm)
{
	// clone additional data , and make sure the copied particle system deletes it
	nm->mFreePdata = true;

	for(int i=0; i<(int)mPartData.size(); ++i)
    {
		ParticleDataBase* pdata = mPartData[i]->clone();
		nm->registerPdata(pdata);
	} 
}

void ParticleBase::deregister(ParticleDataBase* pdata)
{
	bool done = false;
	// remove pointer from particle data list
	for(int i=0; i<(int)mPartData.size(); ++i)
    {
		if(mPartData[i] == pdata)
        {
			if(i<(int)mPartData.size()-1)
				mPartData[i] = mPartData[mPartData.size()-1];
			mPartData.pop_back();
			done = true;
		}
	}

	if(!done) errMsg("Invalid pointer given, not registered!");
}

void ParticleBase::registerPdata(ParticleDataBase* pdata)
{
	pdata->setParticleSys(this);
	mPartData.push_back(pdata);

	if (pdata->getType() == ParticleDataBase::TypeReal)
    {
		ParticleDataImpl<Real>* pd = dynamic_cast< ParticleDataImpl<Real>* >(pdata);
		if(!pd) errMsg("Invalid pdata object posing as real!");
		this->registerPdataReal(pd);
	} else if (pdata->getType() == ParticleDataBase::TypeInt)
    {
		ParticleDataImpl<int>* pd = dynamic_cast< ParticleDataImpl<int>* >(pdata);
		if(!pd) errMsg("Invalid pdata object posing as int!");
		this->registerPdataInt(pd);
	} else if (pdata->getType() == ParticleDataBase::TypeVec3)
    {
		ParticleDataImpl<Vec3>* pd = dynamic_cast< ParticleDataImpl<Vec3>* >(pdata);
		if (!pd) errMsg("Invalid pdata object posing as vec3!");
		this->registerPdataVec3(pd);
	}
}

void ParticleBase::registerPdataReal(ParticleDataImpl<Real>* pd)
{
    mPdataReal.push_back(pd);
}

void ParticleBase::registerPdataVec3(ParticleDataImpl<Vec3>* pd)
{
    mPdataVec3.push_back(pd);
}

void ParticleBase::registerPdataInt(ParticleDataImpl<int >* pd)
{
    mPdataInt.push_back(pd);
}

void ParticleBase::addAllPdata()
{
	for (int i=0; i<(int)mPartData.size(); ++i)
    {
		mPartData[i]->addEntry();
	} 
}
 


//******************************************************************************
// BasicParticleSystem members
//******************************************************************************
BasicParticleSystem::BasicParticleSystem(FluidSolver* parent)
    : ParticleSystem<BasicParticleData>(parent)
{
	this->mAllowCompress = false;
}

// file io
void BasicParticleSystem::writeParticlesText(string name)
{
	ofstream ofs(name.c_str());
	if (!ofs.good())
		errMsg("can't open file!");
	ofs << this->size()<<", pdata: "<< mPartData.size()<<" ("<<mPdataInt.size()<<","<<mPdataReal.size()<<","<<mPdataVec3.size()<<") \n";
	
    for(int i=0; i<this->size(); ++i)
    {
		ofs << i<<": "<< this->getPos(i) <<" , "<< this->getStatus(i) <<". "; 
		for(int pd=0; pd<(int)mPdataInt.size() ; ++pd) ofs << mPdataInt [pd]->get(i)<<" ";
		for(int pd=0; pd<(int)mPdataReal.size(); ++pd) ofs << mPdataReal[pd]->get(i)<<" ";
		for(int pd=0; pd<(int)mPdataVec3.size(); ++pd) ofs << mPdataVec3[pd]->get(i)<<" ";
		ofs << "\n"; 
	}

	ofs.close();
}

void BasicParticleSystem::writeParticlesRawPositionsGz(string name)
{
#	if NO_ZLIB!=1
	gzFile gzf = gzopen(name.c_str(), "wb1");
	if (!gzf) errMsg("can't open file "<<name);

	for (int i=0; i<this->size(); ++i)
    {
		Vector3D<float> p = toVec3f( this->getPos(i) );
		gzwrite(gzf, &p, sizeof(float)*3);
	}

	gzclose(gzf);
#	else
	cout << "file format not supported without zlib" << endl;
#	endif
}

void BasicParticleSystem::writeParticlesRawVelocityGz(string name)
{
#	if NO_ZLIB!=1
	gzFile gzf = gzopen(name.c_str(), "wb1");
	if (!gzf) errMsg("can't open file "<<name);

	if(mPdataVec3.size() < 1) errMsg("no vec3 particle data channel found!");

	// note , assuming particle data vec3 0 is velocity! make optional...
	for (int i=0; i<this->size(); ++i)
    {		
		Vector3D<float> p = toVec3f( mPdataVec3[0]->get(i) );
		gzwrite(gzf, &p, sizeof(float)*3);
	}

	gzclose(gzf);
#	else
	cout << "file format not supported without zlib" << endl;
#	endif
}

void BasicParticleSystem::load(string name)
{
	if (name.find_last_of('.') == string::npos)
		errMsg("file '" + name + "' does not have an extension");
	string ext = name.substr(name.find_last_of('.'));
	if ( ext == ".uni") 
		readParticlesUni(name, this );
	else if ( ext == ".raw") // raw = uni for now
		readParticlesUni(name, this );
	else 
		errMsg("particle '" + name +"' filetype not supported for loading");
}

void BasicParticleSystem::save(string name)
{
	if (name.find_last_of('.') == string::npos)
		errMsg("file '" + name + "' does not have an extension");
	string ext = name.substr(name.find_last_of('.'));
	if (ext == ".txt") 
		this->writeParticlesText(name);
	else if (ext == ".uni") 
		writeParticlesUni(name, this);
	else if (ext == ".raw") // raw = uni for now
		writeParticlesUni(name, this);
	// raw data formats, very basic for simple data transfer to other programs
	else if (ext == ".posgz") 
		this->writeParticlesRawPositionsGz(name);
	else if (ext == ".velgz") 
		this->writeParticlesRawVelocityGz(name);
	else
		errMsg("particle '" + name +"' filetype not supported for saving");
}

void BasicParticleSystem::printParts(int start, int stop, bool printIndex)
{
	std::ostringstream sstr;
	int s = (start>0 ? start : 0                 );
	int e = (stop>0  ? stop  : (int)mData.size() );
	s = Manta::clamp(s, 0, (int)mData.size());
	e = Manta::clamp(e, 0, (int)mData.size());

	for(int i=s; i<e; ++i)
    {
		if(printIndex) sstr << i<<": ";
		sstr<<mData[i].pos<<" "<<mData[i].flag<<"\n";
	} 

	debMsg( sstr.str() , 1 );
}



//******************************************************************************
// ParticleData members
//******************************************************************************
ParticleDataBase::ParticleDataBase(FluidSolver* parent)
    : PbClass(parent)
    , mpParticleSys(NULL)
{
}

ParticleDataBase::~ParticleDataBase()
{
	// notify parent of deletion 
	if(mpParticleSys)
		mpParticleSys->deregister(this);
}



//******************************************************************************
// ParticleData implementation
//******************************************************************************
template<class T>
ParticleDataImpl<T>::ParticleDataImpl(ParticleBase* particle_system)
    : ParticleDataBase(particle_system->getParent())
    , mpGridSource(nullptr)
    , mGridSourceMAC(false)
{
    this->setParticleSys(particle_system);
    mpParticleSys->registerPdata(this);
}

template<class T>
ParticleDataImpl<T>::ParticleDataImpl(FluidSolver* parent)
    : ParticleDataBase(parent)
    , mpGridSource(NULL)
    , mGridSourceMAC(false)
{
}

template<class T>
ParticleDataImpl<T>::ParticleDataImpl(FluidSolver* parent, ParticleDataImpl<T>* other)
    : ParticleDataBase(parent)
    , mpGridSource(NULL)
    , mGridSourceMAC(false)
{
	this->mData = other->mData;
}

template<class T>
ParticleDataImpl<T>::~ParticleDataImpl()
{
}

// set contents to zero, as for a grid
template<class T>
void ParticleDataImpl<T>::clear()
{
	for(int i=0; i<(int)mData.size(); ++i) mData[i] = 0.;
}

template<class T>
int ParticleDataImpl<T>::getSizeSlow() const
{
	return mData.size();
}

template<class T>
void ParticleDataImpl<T>::addEntry()
{
	// add zero'ed entry
	T tmp = T(0.);
	// for debugging, force init:
	//tmp = T(0.02 * mData.size()); // increasing
	//tmp = T(1.); // constant 1
	return mData.push_back(tmp);
}

template<class T>
void ParticleDataImpl<T>::resize(int s)
{
	mData.resize(s);
}

template<class T>
void ParticleDataImpl<T>::copyValueSlow(int from, int to)
{
	this->copyValue(from,to);
}

template<class T>
ParticleDataBase* ParticleDataImpl<T>::clone()
{
	ParticleDataImpl<T>* npd = new ParticleDataImpl<T>(getParent(), this);
	return npd;
}

template<class T>
void ParticleDataImpl<T>::setSource(Grid<T>* grid, bool isMAC)
{
	mpGridSource = grid;
	mGridSourceMAC = isMAC;
	if(isMAC) assertMsg( dynamic_cast<MACGrid*>(grid) != NULL , "Given grid is not a valid MAC grid");
}

template<class T>
void ParticleDataImpl<T>::initNewValue(int idx, Vec3 pos)
{
	if(!mpGridSource)
		mData[idx] = 0; 
	else
    {
		mData[idx] = mpGridSource->getInterpolated(pos);
	}
}

// special handling needed for velocities
template<>
void ParticleDataImpl<Vec3>::initNewValue(int idx, Vec3 pos)
{
	if(!mpGridSource)
		mData[idx] = 0;
	else {
		if(!mGridSourceMAC)
			mData[idx] = mpGridSource->getInterpolated(pos);
		else
			mData[idx] = ((MACGrid*)mpGridSource)->getInterpolated(pos);
	}
}

template<typename T>
void ParticleDataImpl<T>::load(string name)
{
	if (name.find_last_of('.') == string::npos)
		errMsg("file '" + name + "' does not have an extension");
	string ext = name.substr(name.find_last_of('.'));
	if ( ext == ".uni") 
		readPdataUni<T>(name, this);
	else if ( ext == ".raw") // raw = uni for now 
		readPdataUni<T>(name, this);
	else 
		errMsg("particle data '" + name +"' filetype not supported for loading");
}

template<typename T>
void ParticleDataImpl<T>::save(string name)
{
	if (name.find_last_of('.') == string::npos)
		errMsg("file '" + name + "' does not have an extension");
	string ext = name.substr(name.find_last_of('.'));
	if (ext == ".uni") 
		writePdataUni<T>(name, this);
	else if (ext == ".raw") // raw = uni for now
		writePdataUni<T>(name, this);
	else
		errMsg("particle data '" + name +"' filetype not supported for saving");
}

// specializations
template<>
ParticleDataBase::PdataType ParticleDataImpl<Real>::getType() const
{
	return ParticleDataBase::TypeReal;
}

template<>
ParticleDataBase::PdataType ParticleDataImpl<int>::getType() const
{
	return ParticleDataBase::TypeInt;
}

template<>
ParticleDataBase::PdataType ParticleDataImpl<Vec3>::getType() const
{
	return ParticleDataBase::TypeVec3;
}

// note, we need a flag value for functions such as advection
// ideally, this value should never be modified
int ParticleIndexData::flag = 0; 
Vec3 ParticleIndexData::pos = Vec3(0.,0.,0.); 

template <class T>
struct knSetPdataConst : public KernelBase
{
    knSetPdataConst(ParticleDataImpl<T>& pdata, T value)
        : KernelBase(pdata.size())
        , pdata(pdata)
        , value(value)
    {
        run();
    }
    
    inline void op(int idx, ParticleDataImpl<T>& pdata, T value )
    {
        pdata[idx] = value;
    }
    
    inline ParticleDataImpl<T>& getArg0() { return pdata; }
    typedef ParticleDataImpl<T> type0;
    inline T& getArg1() { return value; }
    typedef T type1;
    
    void run()
    {
        const int _sz = size;
#pragma omp parallel
        {
            this->threadId = omp_get_thread_num();
            this->threadNum = omp_get_num_threads();
#pragma omp for
            for (int i=0; i < _sz; i++)
                op(i,pdata,value);
        }
    }
    
    ParticleDataImpl<T>& pdata;
    T value;
};

template <class T, class S>
struct knPdataSet : public KernelBase
{
    knPdataSet(ParticleDataImpl<T>& me, const ParticleDataImpl<S>& other)
        : KernelBase(me.size())
        , me(me)
        , other(other)
    {
        run();
    }
    
    inline void op(int idx, ParticleDataImpl<T>& me, const ParticleDataImpl<S>& other)
    {
        me[idx] += other[idx];
    }
    
    inline ParticleDataImpl<T>& getArg0() { return me; }
    typedef ParticleDataImpl<T> type0;
    inline const ParticleDataImpl<S>& getArg1() { return other; }
    typedef ParticleDataImpl<S> type1;
    
    void run()
    {
        const int _sz = size;
#pragma omp parallel
        {
            this->threadId = omp_get_thread_num();
            this->threadNum = omp_get_num_threads();
#pragma omp for
            for (int i=0; i < _sz; i++)
                op(i,me,other);
        }
    }
    
    ParticleDataImpl<T>& me;
    const ParticleDataImpl<S>& other;
};

template <class T, class S>
struct knPdataAdd : public KernelBase
{
    knPdataAdd(ParticleDataImpl<T>& me, const ParticleDataImpl<S>& other)
        : KernelBase(me.size())
        , me(me)
        , other(other)
    {
        run();
    }
    
    inline void op(int idx, ParticleDataImpl<T>& me, const ParticleDataImpl<S>& other)
    {
        me[idx] += other[idx];
    }
    
    inline ParticleDataImpl<T>& getArg0() { return me; }
    typedef ParticleDataImpl<T> type0;
    inline const ParticleDataImpl<S>& getArg1() { return other; }
    typedef ParticleDataImpl<S> type1;
    
    void run()
    {
        const int _sz = size;
#pragma omp parallel
        {
            this->threadId = omp_get_thread_num();
            this->threadNum = omp_get_num_threads();
#pragma omp for
            for (int i=0; i < _sz; i++)
                op(i,me,other);
        }
    }
    
    ParticleDataImpl<T>& me;
    const ParticleDataImpl<S>& other;
};

template <class T, class S>
struct knPdataSub : public KernelBase
{
    knPdataSub(ParticleDataImpl<T>& me, const ParticleDataImpl<S>& other)
        : KernelBase(me.size())
        , me(me)
        , other(other)
    {
        run();
    }
    
    inline void op(int idx, ParticleDataImpl<T>& me, const ParticleDataImpl<S>& other)
    {
        me[idx] -= other[idx];
    }
    
    inline ParticleDataImpl<T>& getArg0() { return me; }
    typedef ParticleDataImpl<T> type0;
    inline const ParticleDataImpl<S>& getArg1() { return other; }
    typedef ParticleDataImpl<S> type1;
    
    void run()
    {
        const int _sz = size;
#pragma omp parallel
        {
            this->threadId = omp_get_thread_num();
            this->threadNum = omp_get_num_threads();
#pragma omp for
            for (int i=0; i < _sz; i++)
                op(i,me,other);
        }
    }
    
    ParticleDataImpl<T>& me;
    const ParticleDataImpl<S>& other;
};

template <class T, class S>
struct knPdataMult : public KernelBase
{
    knPdataMult(ParticleDataImpl<T>& me, const ParticleDataImpl<S>& other)
        : KernelBase(me.size())
        , me(me)
        , other(other)
    {
        run();
    }
    
    inline void op(int idx, ParticleDataImpl<T>& me, const ParticleDataImpl<S>& other)
    {
        me[idx] *= other[idx];
    }
    
    inline ParticleDataImpl<T>& getArg0() { return me; }
    typedef ParticleDataImpl<T> type0;
    inline const ParticleDataImpl<S>& getArg1() { return other; }
    typedef ParticleDataImpl<S> type1;
    
    void run()
    {
        const int _sz = size;
#pragma omp parallel
        {
            this->threadId = omp_get_thread_num();
            this->threadNum = omp_get_num_threads();
#pragma omp for
            for (int i=0; i < _sz; i++)
                op(i,me,other);
        }
    }
    
    ParticleDataImpl<T>& me;
    const ParticleDataImpl<S>& other;
};

template <class T, class S>
struct knPdataDiv : public KernelBase
{
    knPdataDiv(ParticleDataImpl<T>& me, const ParticleDataImpl<S>& other)
        : KernelBase(me.size())
        , me(me)
        , other(other)
    {
        run();
    }
    
    inline void op(int idx, ParticleDataImpl<T>& me, const ParticleDataImpl<S>& other)
    {
        me[idx] /= other[idx];
    }
    
    inline ParticleDataImpl<T>& getArg0() { return me; }
    typedef ParticleDataImpl<T> type0;
    inline const ParticleDataImpl<S>& getArg1() { return other; }
    typedef ParticleDataImpl<S> type1;
    
    void run()
    {
        const int _sz = size;
#pragma omp parallel
        {
            this->threadId = omp_get_thread_num();
            this->threadNum = omp_get_num_threads();
#pragma omp for
            for (int i=0; i < _sz; i++)
                op(i,me,other);
        }
    }
    
    ParticleDataImpl<T>& me;
    const ParticleDataImpl<S>& other;
};

template <class T, class S>
struct knPdataSetScalar : public KernelBase
{
    knPdataSetScalar(ParticleDataImpl<T>& me, const S& other)
        : KernelBase(me.size())
        , me(me)
        , other(other)
    {
        run();
    }
    
    inline void op(int idx, ParticleDataImpl<T>& me, const S& other)
    {
        me[idx]  = other;
    }
    
    inline ParticleDataImpl<T>& getArg0() { return me; }
    typedef ParticleDataImpl<T> type0;
    inline const S& getArg1() { return other; }
    typedef S type1;
    
    void run()
    {
        const int _sz = size;
#pragma omp parallel
        {
            this->threadId = omp_get_thread_num();
            this->threadNum = omp_get_num_threads();
#pragma omp for
            for (int i=0; i < _sz; i++)
                op(i,me,other);
        }
    }
    
    ParticleDataImpl<T>& me;
    const S& other;
};

template <class T, class S>
struct knPdataAddScalar : public KernelBase
{
    knPdataAddScalar(ParticleDataImpl<T>& me, const S& other)
        : KernelBase(me.size())
        , me(me)
        , other(other)
    {
        run();
    }
    
    inline void op(int idx, ParticleDataImpl<T>& me, const S& other)
    {
        me[idx] += other;
    }
    
    inline ParticleDataImpl<T>& getArg0() { return me; }
    typedef ParticleDataImpl<T> type0;
    inline const S& getArg1() { return other; }
    typedef S type1;
    
    void run()
    {
        const int _sz = size;
#pragma omp parallel
        {
            this->threadId = omp_get_thread_num();
            this->threadNum = omp_get_num_threads();
#pragma omp for
            for (int i=0; i < _sz; i++)
                op(i,me,other);
        }
    }
    
    ParticleDataImpl<T>& me;
    const S& other;
};

template <class T, class S>
struct knPdataMultScalar : public KernelBase
{
    knPdataMultScalar(ParticleDataImpl<T>& me, const S& other)
        : KernelBase(me.size())
        , me(me)
        , other(other)
    {
        run();
    }
    
    inline void op(int idx, ParticleDataImpl<T>& me, const S& other)
    {
        me[idx] *= other;
    }
    
    inline ParticleDataImpl<T>& getArg0() { return me; }
    typedef ParticleDataImpl<T> type0;
    inline const S& getArg1() { return other; }
    typedef S type1;
    
    void run()
    {
        const int _sz = size;
#pragma omp parallel
        {
            this->threadId = omp_get_thread_num();
            this->threadNum = omp_get_num_threads();
#pragma omp for
            for (int i=0; i < _sz; i++)
                op(i,me,other);
        }
    }
    
    ParticleDataImpl<T>& me;
    const S& other;
};

template <class T, class S>
struct knPdataScaledAdd : public KernelBase
{
    knPdataScaledAdd(ParticleDataImpl<T>& me, const ParticleDataImpl<T>& other, const S& factor)
        : KernelBase(me.size())
        , me(me)
        , other(other)
        , factor(factor)
    {
        run();
    }
    
    inline void op(int idx, ParticleDataImpl<T>& me, const ParticleDataImpl<T>& other, const S& factor)
    {
        me[idx] += factor * other[idx];
    }
    
    inline ParticleDataImpl<T>& getArg0() { return me; }
    typedef ParticleDataImpl<T> type0;
    inline const ParticleDataImpl<T>& getArg1() { return other; }
    typedef ParticleDataImpl<T> type1;
    inline const S& getArg2() { return factor; }
    typedef S type2;
    
    void run()
    {
        const int _sz = size;
#pragma omp parallel
        {
            this->threadId = omp_get_thread_num();
            this->threadNum = omp_get_num_threads();
#pragma omp for
            for (int i=0; i < _sz; i++)
                op(i,me,other,factor);
        }
    }
    
    ParticleDataImpl<T>& me;
    const ParticleDataImpl<T>& other;
    const S& factor;
};

template <class T>
struct knPdataSafeDiv : public KernelBase
{
    knPdataSafeDiv(ParticleDataImpl<T>& me, const ParticleDataImpl<T>& other)
        : KernelBase(me.size())
        , me(me)
        , other(other)
    {
        run();
    }
    
    inline void op(int idx, ParticleDataImpl<T>& me, const ParticleDataImpl<T>& other)
    {
        me[idx] = safeDivide(me[idx], other[idx]);
    }
    
    inline ParticleDataImpl<T>& getArg0() { return me; }
    typedef ParticleDataImpl<T> type0;
    inline const ParticleDataImpl<T>& getArg1() { return other; }
    typedef ParticleDataImpl<T> type1;
    
    void run()
    {
        const int _sz = size;
#pragma omp parallel
        {
            this->threadId = omp_get_thread_num();
            this->threadNum = omp_get_num_threads();
#pragma omp for
            for (int i=0; i < _sz; i++)
                op(i,me,other);
        }
    }
    
    ParticleDataImpl<T>& me;
    const ParticleDataImpl<T>& other;
};

template <class T>
struct knPdataSetConst : public KernelBase
{
    knPdataSetConst(ParticleDataImpl<T>& pdata, T value)
        : KernelBase(pdata.size())
        , pdata(pdata)
        , value(value)
    {
        run();
    }
    
    inline void op(int idx, ParticleDataImpl<T>& pdata, T value)
    {
        pdata[idx] = value;
    }
    
    inline ParticleDataImpl<T>& getArg0() { return pdata; }
    typedef ParticleDataImpl<T> type0;
    inline T& getArg1() { return value; }
    typedef T type1;
    
    void run()
    {
        const int _sz = size;
#pragma omp parallel
        {
            this->threadId = omp_get_thread_num();
            this->threadNum = omp_get_num_threads();
#pragma omp for
            for (int i=0; i < _sz; i++)
                op(i,pdata,value);
        }
    }
    
    ParticleDataImpl<T>& pdata;
    T value;
};

template <class T>
struct knPdataClamp : public KernelBase
{
    knPdataClamp(ParticleDataImpl<T>& me, T min, T max)
        : KernelBase(me.size())
        , me(me)
        , min(min)
        , max(max)
    {
        run();
    }
    
    inline void op(int idx, ParticleDataImpl<T>& me, T min, T max)
    {
        me[idx] = clamp( me[idx], min, max);
    }
    
    inline ParticleDataImpl<T>& getArg0() { return me; }
    typedef ParticleDataImpl<T> type0;
    inline T& getArg1() { return min; }
    typedef T type1;
    inline T& getArg2() { return max; }
    typedef T type2;
    
    void run()
    {
        const int _sz = size;
#pragma omp parallel
        {
            this->threadId = omp_get_thread_num();
            this->threadNum = omp_get_num_threads();
#pragma omp for
            for (int i=0; i < _sz; i++)
                op(i,me,min,max);
        }
    }
    
    ParticleDataImpl<T>& me;
    T min;
    T max;
};

// python operators
template<typename T>
ParticleDataImpl<T>& ParticleDataImpl<T>::copyFrom(const ParticleDataImpl<T>& a)
{
	assertMsg (a.mData.size() == mData.size() , "different pdata size "<<a.mData.size()<<" vs "<<this->mData.size() );
	memcpy(&mData[0], &a.mData[0], sizeof(T) * mData.size());
	return *this; 
}

template<typename T>
void ParticleDataImpl<T>::setConst(T s)
{
	knPdataSetScalar<T,T>op(*this, s);
}

template<typename T>
void ParticleDataImpl<T>::add(const ParticleDataImpl<T>& a)
{
	knPdataAdd<T,T>op(*this, a);
}

template<typename T>
void ParticleDataImpl<T>::sub(const ParticleDataImpl<T>& a)
{
	knPdataSub<T,T>op(*this, a);
}

template<typename T>
void ParticleDataImpl<T>::addConst(T s)
{
	knPdataAddScalar<T,T>op(*this, s);
}

template<typename T>
void ParticleDataImpl<T>::addScaled(const ParticleDataImpl<T>& a, const T& factor)
{
	knPdataScaledAdd<T,T>op(*this, a, factor);
}

template<typename T>
void ParticleDataImpl<T>::mult(const ParticleDataImpl<T>& a)
{
	knPdataMult<T,T>op(*this, a);
}

template<typename T>
void ParticleDataImpl<T>::multConst(T s)
{
	knPdataMultScalar<T,T>op(*this, s);
}

template<typename T>
void ParticleDataImpl<T>::clamp(Real min, Real max)
{
	knPdataClamp<T>op(*this, min,max);
}

template<typename T>
struct CompPdata_Min : public KernelBase
{
    CompPdata_Min(const ParticleDataImpl<T>& val)
        : KernelBase(val.size())
        , val(val)
        , minVal(std::numeric_limits<Real>::max())
    {
        run();
    }
    
    inline void op(int idx, const ParticleDataImpl<T>& val, Real& minVal)
    {
	    if (val[idx] < minVal)
		    minVal = val[idx];
    }
    
    inline operator Real () { return minVal; }
    inline Real  & getRet() { return minVal; }
    inline const ParticleDataImpl<T>& getArg0() { return val; }
    typedef ParticleDataImpl<T> type0;
    
    void run()
    {
        const int _sz = size;
#pragma omp parallel
        {
            this->threadId = omp_get_thread_num();
            this->threadNum = omp_get_num_threads();
            Real minVal = std::numeric_limits<Real>::max();
#pragma omp for nowait
            for (int i=0; i < _sz; i++)
                op(i,val,minVal); 
#pragma omp critical
            {
                this->minVal = min(minVal, this->minVal);
            }
        }
    }
    
    const ParticleDataImpl<T>& val;
    Real minVal;
};

template<typename T>
struct CompPdata_Max : public KernelBase
{
    CompPdata_Max(const ParticleDataImpl<T>& val)
        : KernelBase(val.size())
        , val(val)
        , maxVal(-std::numeric_limits<Real>::max())
    {
        run();
    }
    
    inline void op(int idx, const ParticleDataImpl<T>& val, Real& maxVal)
    {
	    if (val[idx] > maxVal)
		    maxVal = val[idx];
    }
    
    inline operator Real () { return maxVal; }
    inline Real  & getRet() { return maxVal; }
    inline const ParticleDataImpl<T>& getArg0() { return val; }
    typedef ParticleDataImpl<T> type0;
    
    void run()
    {
        const int _sz = size;
#pragma omp parallel
        {
            this->threadId = omp_get_thread_num();
            this->threadNum = omp_get_num_threads();
            Real maxVal = -std::numeric_limits<Real>::max();
#pragma omp for nowait
            for (int i=0; i < _sz; i++)
                op(i,val,maxVal); 
#pragma omp critical
            {
                this->maxVal = max(maxVal, this->maxVal);
            }
        }
    }
    
    const ParticleDataImpl<T>& val;
    Real maxVal;
};

template<typename T>
Real ParticleDataImpl<T>::getMinValue()
{
	return sqrt(CompPdata_Min<T>(*this));
}

template<typename T>
Real ParticleDataImpl<T>::getMaxAbsValue()
{
	Real amin = CompPdata_Min<T> (*this);
	Real amax = CompPdata_Max<T> (*this);
	return max(fabs(amin), fabs(amax));
}

template<typename T>
Real ParticleDataImpl<T>::getMaxValue()
{
	return sqrt(CompPdata_Max<T>(*this));
} 

template<typename T>
void ParticleDataImpl<T>::printPdata(int start, int stop, bool printIndex)
{
	std::ostringstream sstr;
	int s = (start>0 ? start : 0                 );
	int e = (stop>0  ? stop  : (int)mData.size() );
	s = Manta::clamp(s, 0, (int)mData.size());
	e = Manta::clamp(e, 0, (int)mData.size());

	for(int i=s; i<e; ++i)
    {
		if(printIndex) sstr << i<<": ";
		sstr<<mData[i]<<" "<<"\n";
	}

	debMsg( sstr.str() , 1 );
}

// specials for vec3
struct CompPdata_MinVec3 : public KernelBase
{
    CompPdata_MinVec3(const ParticleDataImpl<Vec3>& val)
        : KernelBase(val.size())
        , val(val)
        , minVal(-std::numeric_limits<Real>::max())
    {
        run();
    }
    
    inline void op(int idx, const ParticleDataImpl<Vec3>& val, Real& minVal)
    {
	    const Real s = normSquare(val[idx]);
	    if (s < minVal)
		    minVal = s;
    }
    
    inline operator Real () { return minVal; }
    inline Real  & getRet() { return minVal; }
    inline const ParticleDataImpl<Vec3>& getArg0() { return val; }
    typedef ParticleDataImpl<Vec3> type0;
    
    void run()
    {
        const int _sz = size;
#pragma omp parallel
        {
            this->threadId = omp_get_thread_num();
            this->threadNum = omp_get_num_threads();
            Real minVal = -std::numeric_limits<Real>::max();
#pragma omp for nowait
            for (int i=0; i < _sz; i++)
                op(i,val,minVal); 
#pragma omp critical
            {
                this->minVal = min(minVal, this->minVal);
            }
        }
    }
    
    const ParticleDataImpl<Vec3>& val;
    Real minVal;
};

struct CompPdata_MaxVec3 : public KernelBase
{
    CompPdata_MaxVec3(const ParticleDataImpl<Vec3>& val)
        : KernelBase(val.size())
        , val(val)
        , maxVal(-std::numeric_limits<Real>::min())
    {
        run();
    }
    
    inline void op(int idx, const ParticleDataImpl<Vec3>& val, Real& maxVal)
    {
	    const Real s = normSquare(val[idx]);
	    if (s > maxVal)
		    maxVal = s;
    }
    
    inline operator Real () { return maxVal; }
    inline Real  & getRet() { return maxVal; }
    inline const ParticleDataImpl<Vec3>& getArg0() { return val; }
    typedef ParticleDataImpl<Vec3> type0;
    
    void run()
    {
        const int _sz = size;
#pragma omp parallel
        {
            this->threadId = omp_get_thread_num();
            this->threadNum = omp_get_num_threads();
            Real maxVal = -std::numeric_limits<Real>::min();
#pragma omp for nowait
            for (int i=0; i < _sz; i++) op(i,val,maxVal);
#pragma omp critical
            {
                this->maxVal = max(maxVal, this->maxVal);
            }
        }
    }
    
    const ParticleDataImpl<Vec3>& val;
    Real maxVal;
};

template<>
Real ParticleDataImpl<Vec3>::getMinValue()
{
	return sqrt(CompPdata_MinVec3 (*this));
}

template<>
Real ParticleDataImpl<Vec3>::getMaxAbsValue()
{
	Real amin = CompPdata_MinVec3 (*this);
	Real amax = CompPdata_MaxVec3 (*this);
	return max( fabs(amin), fabs(amax));
}

template<>
Real ParticleDataImpl<Vec3>::getMaxValue()
{
	return sqrt(CompPdata_MaxVec3 (*this));
}



// explicit instantiation
template class ParticleDataImpl<int>;
template class ParticleDataImpl<Real>;
template class ParticleDataImpl<Vec3>;

} // namespace
