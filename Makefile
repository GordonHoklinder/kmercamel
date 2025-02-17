.PHONY: all clean test cpptest converttest verify quick-verify

CXX=         g++
CXXFLAGS=    -g -Wall -Wno-unused-function -std=c++17 -O2 -I/opt/homebrew/include
LDFLAGS=     -lz -lglpk -L/opt/homebrew/lib
LARGEFLAGS=  -DLARGE_KMERS
ELARGEFLAGS=  -DEXTRA_LARGE_KMERS
SRC=         src
TESTS=       tests
GTEST=       $(TESTS)/googletest/googletest
DATA=        data


all: kmercamel

test: cpptest verify

mask-verify:
	./verify.py --k 13 --superstring_path $(DATA)/global-k13c.fa $(DATA)/spneumoniae.fa
	./verify.py --k 31 --superstring_path $(DATA)/global-k31c.fa $(DATA)/spneumoniae.fa
	./verify.py --k 63 --superstring_path $(DATA)/global-k63c.fa $(DATA)/spneumoniae.fa
	./verify.py --k 127 --superstring_path $(DATA)/global-k127c.fa $(DATA)/spneumoniae.fa


verify: verify.py kmercamel
	./verify.py $(DATA)/spneumoniae.fa

quick-verify: verify.py kmercamel
	./verify.py --quick $(DATA)/spneumoniae.fa

cpptest: kmercameltest kmercameltest-large kmercameltest-extra-large
	./kmercameltest
	./kmercameltest-large
	./kmercameltest-extra-large

kmercamel: $(SRC)/main.cpp $(SRC)/$(wildcard *.cpp *.h *.hpp) src/version.h
	./create-version.sh
	$(CXX) $(CXXFLAGS) $(SRC)/main.cpp -o $@ $(LDFLAGS)
	cp kmercamel  🐫 || true

kmercameltest: $(TESTS)/unittest.cpp gtest-all.o $(SRC)/$(wildcard *.cpp *.h *.hpp) $(TESTS)/$(wildcard *.cpp *.h *.hpp)
	$(CXX) $(CXXFLAGS) -isystem $(GTEST)/include -I $(GTEST)/include $(TESTS)/unittest.cpp gtest-all.o -pthread -o $@ $(LDFLAGS)

kmercameltest-large: $(TESTS)/unittest.cpp gtest-all.o $(SRC)/$(wildcard *.cpp *.h *.hpp) $(TESTS)/$(wildcard *.cpp *.h *.hpp)
	$(CXX) $(CXXFLAGS) -isystem $(GTEST)/include -I $(GTEST)/include $(TESTS)/unittest.cpp gtest-all.o -pthread -o $@ $(LDFLAGS) $(LARGEFLAGS)

kmercameltest-extra-large: $(TESTS)/unittest.cpp gtest-all.o $(SRC)/$(wildcard *.cpp *.h *.hpp) $(TESTS)/$(wildcard *.cpp *.h *.hpp)
	$(CXX) $(CXXFLAGS) -isystem $(GTEST)/include -I $(GTEST)/include $(TESTS)/unittest.cpp gtest-all.o -pthread -o $@ $(LDFLAGS) $(ELARGEFLAGS)

gtest-all.o: $(GTEST)/src/gtest-all.cc $(wildcard *.cpp *.h *.hpp)
	$(CXX) $(CXXFLAGS) -isystem $(GTEST)/include -I $(GTEST)/include -I $(GTEST) -DGTEST_CREATE_SHARED_LIBRARY=1 -c -pthread $(GTEST)/src/gtest-all.cc -o $@

src/version.h: src/version
	./create-version.sh


clean:
	rm -f kmercamel
	rm -f 🐫 || true
	rm -f kmercameltest
	rm -f kmercameltest-large
	rm -r -f ./bin
	rm -f gtest-all.o
	rm -f src/version.h
