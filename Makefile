# derived from hal2vg/Makefile
rootPath = ./

all : panpatch

libbdsgPath = deps/libbdsg-easy
libbdsgLibs = ${libbdsgPath}/lib/libbdsg.a ${libbdsgPath}/lib/libhandlegraph.a ${libbdsgPath}/lib/libsdsl.a ${libbdsgPath}/lib/libdivsufsort.a ${libbdsgPath}/lib/libdivsufsort64.a

# Multithreading with OpenMP.
PARALLEL_FLAGS=-fopenmp -pthread
LIBS=-L$(CWD)

ifeq ($(shell uname -s),Darwin)
    # Our compiler might be clang that lacks -fopenmp support.
    # Sniff that
    ifeq ($(strip $(shell $(MY_CXX) -fopenmp /dev/null -o/dev/null 2>&1 | grep fopenmp | wc -l)), 1)
        # The compiler complained about fopenmp instead of its nonsense input file.
        # We need to use the hard way of getting OpenMP not bundled with the compiler.

        # The compiler only needs to do the preprocessing
        PARALLEL_FLAGS = -Xpreprocessor -fopenmp -pthread

        ifeq ($(shell if [ -d /opt/local/lib/libomp ];then echo 1;else echo 0;fi), 1)
            # Use /opt/local/lib/libomp if present, because Macports installs libomp there.
            # Brew is supposed to put it somewhere the compiler can find it by default.
            LIBS += -L/opt/local/lib/libomp
            # And we need to find the includes. Homebrew puts them in the normal place
            # but Macports hides them in "libomp"
            PARALLEL_FLAGS += -I/opt/local/include/libomp
        endif

        # We also need to link it
        LIBS += -lomp

    endif
endif

CXX ?= g++
ifeq (${CXXFLAGS},)
CXXFLAGS = 
endif
#CXXFLAGS += -O0 -fno-inline -fno-omit-frame-pointer -fsanitize=address
CXXFLAGS += -O3
CXXFLAGS += -Werror=return-type -std=c++14 -ggdb -g -MMD -MP $(PARALLEL_FLAGS)

CXXFLAGS += -I deps/libbdsg-easy/include 

static:
	CFLAGS="$${CFLAGS} -static" \
	CXXFLAGS="$${CXXFLAGS} -static" \
	${MAKE} all

check-static: static
	if [ $(shell ls panpatch | xargs ldd 2>& 1 | grep "not a dynamic" | wc -l) = $(shell ls panpatch | wc -l) ] ; then\
		echo "ldd verified that all files in bin/ are static";\
	else\
		echo "ldd found dynamic linked binary in bin/";\
		exit 1;\
	fi

cleanFast : 
	rm -f panpatch *.o

clean :
	rm -f panpatch ${libbdsgPath}/lib/libbdsg.a 
	cd deps/libbdsg-easy && ${MAKE} clean

panpatch_main.o : panpatch_main.cpp ${libbdsgPath}/lib/libbdsg.a 
	${CXX} ${CXXFLAGS} -I . panpatch_main.cpp -c


${libbdsgPath}/lib/libbdsg.a:
	cd deps/libbdsg-easy && ${MAKE}


panpatch : panpatch_main.o  ${libbdsgPath}/lib/libbdsg.a 
	${CXX} ${CXXFLAGS} panpatch_main.o  ${libbdsgLibs}  -lcurl -lm -lz -llzma -lbz2 -ldeflate -fopenmp -pthread -o panpatch

all : panpatch

test : panpatch
	cd test && python3 panpatchTest.py
