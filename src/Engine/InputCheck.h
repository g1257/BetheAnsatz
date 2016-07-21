/*
Copyright (c) 2015, UT-Battelle, LLC

BetheAnsatz, Version 0.1

This file is part of BetheAnsatz.
BetheAnsatz is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
BetheAnsatz is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with BetheAnsatz. If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef BETHE_INPUT_CHECK_H
#define BETHE_INPUT_CHECK_H
#include <vector>
#include <stdexcept>
#include "Vector.h"

namespace BetheAnsatz {

class InputCheck {

public:

	bool check(const PsimagLite::String&,
	           const PsimagLite::Vector<PsimagLite::String>::Type&,
	           SizeType) const
	{
		return false;
	}

	void check(const PsimagLite::String&, const PsimagLite::String&, SizeType)
	{}

	bool checkSimpleLabel(const PsimagLite::String&,
	                      SizeType) const
	{
		// FIXME: needs implementation
		return true;
	}

	void usageMain(const PsimagLite::String& name) const
	{
		std::cerr<<"USAGE: "<<name<<"\n";
	}

	void checkForThreads(SizeType nthreads) const
	{
#ifndef USE_PTHREADS
		if (nthreads==1) return;

		PsimagLite::String message1(__FILE__);
		message1 += " FATAL: You are requesting nthreads>0 but you ";
		message1 += "did not compile with USE_PTHREADS enabled\n";
		message1 += " Either set Threads=1 in the input file (you won't ";
		message1 += "have threads though) or\n";
		message1 += " add -DUSE_PTHREADS to the CPP_FLAGS in your Makefile ";
		message1 += "and recompile\n";
		throw PsimagLite::RuntimeError(message1.c_str());
#endif
	}

private:

	bool checkForVector(const PsimagLite::Vector<PsimagLite::String>::Type& vec) const
	{
		if (vec.size() == 0) return false;
		SizeType n = atoi(vec[0].c_str());
		return (vec.size() == n+1);
	}

	bool checkForMatrix(const PsimagLite::Vector<PsimagLite::String>::Type& vec) const
	{
		if (vec.size() < 2) return false;
		SizeType row = atoi(vec[0].c_str());
		SizeType col = atoi(vec[1].c_str());
		SizeType n = row*col;
		return (vec.size() == n+2);
	}

	bool error1(const PsimagLite::String& message,SizeType line) const
	{
		PsimagLite::String s(__FILE__);
		s += " : Input error for label " + message + " near line " + ttos(line) + "\n";
		throw PsimagLite::RuntimeError(s.c_str());

	}

}; // class InputCheck
} // namespace BetheAnsatz

/*@}*/
#endif // BETHE_INPUT_CHECK_H

