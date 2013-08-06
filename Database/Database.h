#ifndef DATABASE_H
#define DATABASE_H

#include "../Resources/std_include.h"
#include "../Resources/Exception.h"
#include "../Resources/Resources.h"
#include "NcHelper.h"

namespace LiuJamming
{
using namespace netCDF;
using namespace netCDF::exceptions;

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
	enum FileMode{_read=NcFile::read, _write=NcFile::write, _replace=NcFile::replace, _newFile=NcFile::newFile};

	string filename;
	NcFile File;
	const int Mode;

	CDatabase(string fn="temp.nc", NcFile::FileMode mode=NcFile::read, NcFile::FileFormat format=NcFile::nc4);


};

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
//////////////////////////////   IMPLEMENTATION   ///////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

CDatabase::CDatabase(string fn, NcFile::FileMode mode, NcFile::FileFormat format)
	: filename(fn),
	  Mode(mode),
	  //File(fn, mode)
	  File(fn, mode, format)
{};

/*
CDatabase::CDatabase(string fn, NcFile::FileMode mode)
	: filename(fn)
{
	Open(mode);
};
*/



}

#endif //DATABASE_H
