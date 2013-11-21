#ifndef DATABASE_H
#define DATABASE_H

#ifndef DONT_USE_NETCDF

#include "../Resources/std_include.h"
#include "../Resources/Exception.h"
#include "../Resources/Resources.h"
//#include "NcHelper.h"
#include <netcdfcpp.h>

namespace LiuJamming
{

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
//////////////////////////////   CLASS DEFINITION  //////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

class CDatabase
{
public:
	string filename;
	NcFile File;
	const int Mode;

	//To help keep a common api, subclasses should implement the following methods.
	//
	//	int GetNumRecs() const;
	//
	//	void Write(OBJ const &o, int rec=-1);
	//	
	//	void ReadFirst(OBJ &o);
	//	void ReadLast(OBJ &o);
	//	void Read(OBJ &o, int rec);
	//
	//	Here, OBJ is whatever needs to be written/read. It can be a class (as in 
	//	the examples) or multiple parameters (i.e. void Write(int const &i1, int const &i2, int rec=-1);).
	//	In the write method, the rec parameter defaults to -1, in which case the
	//	data is written to a new record at the end of the file. 

	CDatabase(string fn="temp.nc", NcFile::FileMode mode=NcFile::ReadOnly);
};

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
//////////////////////////////   IMPLEMENTATION   ///////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

CDatabase::CDatabase(string fn, NcFile::FileMode mode)
	: filename(fn),
	  Mode(mode),
	  File(fn.c_str(), mode)
{};


}

#include "StaticDB.h"
#include "StdDataDB.h"

#endif //DONT_USE_NETCDF

#endif //DATABASE_H
