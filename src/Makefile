# Modify Config.make instead of this file
# BetheAnsatz by G.A.

include Config.make
CPPFLAGS += -I../../PsimagLite -I../../PsimagLite/src -IEngine
all: hubbard 

hubbard.o: hubbard.cpp  Makefile Config.make
	$(CXX) $(CPPFLAGS) -c hubbard.cpp

hubbard: hubbard.o
	$(CXX) -o  hubbard hubbard.o $(LDFLAGS)
	$(STRIP_COMMAND) hubbard

../../PsimagLite/lib/libpsimaglite.a:
	$(MAKE) -f Makefile -C ../../PsimagLite/lib/

Makefile.dep: hubbard.cpp 
	$(CXX) $(CPPFLAGS) -MM hubbard.cpp  > Makefile.dep

clean: Makefile.dep
	rm -f core* hubbard *.o *.dep  gitrev

include Makefile.dep

