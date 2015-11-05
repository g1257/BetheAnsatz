# BetheAnsatz by G.A.

include Config.make
CPPFLAGS += -I../../PsimagLite -I../../PsimagLite/src 
all: hubbard heisenberg tJ

hubbard.o: hubbard.cpp  Makefile Config.make
	$(CXX) $(CPPFLAGS) -c hubbard.cpp

heisenberg.o: heisenberg.cpp  Makefile Config.make
	$(CXX) $(CPPFLAGS) -c heisenberg.cpp

tJ.o: tJ.cpp Makefile Config.make
	$(CXX) $(CPPFLAGS) -c tJ.cpp

hubbard: hubbard.o
	$(CXX) -o  hubbard hubbard.o $(LDFLAGS)
	$(STRIP_COMMAND) hubbard

tJ: tJ.o
	$(CXX) -o  tJ tJ.o $(LDFLAGS)
	$(STRIP_COMMAND) tJ

heisenberg: heisenberg.o
	$(CXX) -o  heisenberg heisenberg.o $(LDFLAGS)
	$(STRIP_COMMAND) heisenberg

../../PsimagLite/lib/libpsimaglite.a:
	$(MAKE) -f Makefile -C ../../PsimagLite/lib/

Makefile.dep: hubbard.cpp  heisenberg.cpp
	$(CXX) $(CPPFLAGS) -MM hubbard.cpp heisenberg.cpp > Makefile.dep

clean: Makefile.dep
	rm -f core* hubbard *.o *.dep  gitrev

include Makefile.dep

