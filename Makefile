CFLAGS_WARN = -Wall -Wextra
CFLAGS = $(CFLAGS_WARN) -O3 -g -pthread -march=native -mtune=native -fomit-frame-pointer

LIBS = -lm $(HEXL_DIR)/build/hexl/lib/libhexl.a -lstdc++

libmpfr = -lmpfr
LIBS += $(libmpfr)

libgmp = -lgmp
LIBS += $(libgmp)

.PHONY: default all
default: lib
all: lib-all

LAZER_DIR = lazer

THIRD_PARTY_DIR = $(LAZER_DIR)/third_party
HEXL_SUBDIR = hexl-development
HEXL_DIR = $(THIRD_PARTY_DIR)/$(HEXL_SUBDIR)
HEXL_ZIP = $(HEXL_DIR).zip

$(HEXL_DIR): $(HEXL_ZIP)
	cd $(THIRD_PARTY_DIR) && unzip $(HEXL_SUBDIR).zip
	cd $(HEXL_DIR) && cmake -S . -B build -DHEXL_BENCHMARK=OFF -DHEXL_TESTING=OFF
	cd $(HEXL_DIR) && cmake --build build


#### lib lazer
LAZER_INC = \
 $(LAZER_DIR)/src/aes256ctr.h \
 $(LAZER_DIR)/src/brandom.h \
 $(LAZER_DIR)/src/dom.h \
 $(LAZER_DIR)/src/grandom.h \
 $(LAZER_DIR)/src/intvec.h \
 $(LAZER_DIR)/src/lazer-in1.h \
 $(LAZER_DIR)/src/lazer-in2.h \
 $(LAZER_DIR)/src/lnp-tbox.h \
 $(LAZER_DIR)/src/memory.h \
 $(LAZER_DIR)/src/mont.h \
 $(LAZER_DIR)/src/poly.h \
 $(LAZER_DIR)/src/rng.h \
 $(LAZER_DIR)/src/shake128.h \
 $(LAZER_DIR)/src/stopwatch.h \
 $(LAZER_DIR)/src/urandom.h \

LAZER_SRC = \
 $(LAZER_DIR)/src/lazer.c \
 $(LAZER_DIR)/src/abdlop.c \
 $(LAZER_DIR)/src/aes256ctr.c \
 $(LAZER_DIR)/src/aes256ctr-amd64.c \
 $(LAZER_DIR)/src/brandom.c \
 $(LAZER_DIR)/src/bytes.c \
 $(LAZER_DIR)/src/coder.c \
 $(LAZER_DIR)/src/dcompress.c \
 $(LAZER_DIR)/src/dump.c \
 $(LAZER_DIR)/src/grandom.c \
 $(LAZER_DIR)/src/int.c \
 $(LAZER_DIR)/src/intmat.c \
 $(LAZER_DIR)/src/intvec.c \
 $(LAZER_DIR)/src/lin-proofs.c \
 $(LAZER_DIR)/src/lnp.c \
 $(LAZER_DIR)/src/lnp-quad.c \
 $(LAZER_DIR)/src/lnp-quad-eval.c \
 $(LAZER_DIR)/src/lnp-quad-many.c \
 $(LAZER_DIR)/src/lnp-tbox.c \
 $(LAZER_DIR)/src/memory.c \
 $(LAZER_DIR)/src/poly.c \
 $(LAZER_DIR)/src/polymat.c \
 $(LAZER_DIR)/src/polyring.c \
 $(LAZER_DIR)/src/polyvec.c \
 $(LAZER_DIR)/src/quad.c \
 $(LAZER_DIR)/src/rejection.c \
 $(LAZER_DIR)/src/rng.c \
 $(LAZER_DIR)/src/shake128.c \
 $(LAZER_DIR)/src/spolymat.c \
 $(LAZER_DIR)/src/spolyvec.c \
 $(LAZER_DIR)/src/stopwatch.c \
 $(LAZER_DIR)/src/urandom.c \
 $(LAZER_DIR)/src/version.c \

#### lib modified-lazer
LIB_DIR = src

LIBSOURCES = \
 $(LIB_DIR)/lazer.c \
 $(LIB_DIR)/tools.c \
 $(LIB_DIR)/encoding.c \
 $(LIB_DIR)/modified-abdlop.c \
 $(LIB_DIR)/modified-quad.c \
 $(LIB_DIR)/modified-quad-eval.c \
 $(LIB_DIR)/modified-quad-many.c \

TESTS = \
 tests/modified-abdlop-test \
 tests/modified-quad-test \
 tests/modified-quad-eval-test \
 tests/modified-quad-many-test


.PHONY: lib lib-all lib-static lib-shared lib-static-all lib-shared-all
lib-all: lazer.h lib-static-all lib-shared-all
lib: lazer.h lib-static lib-shared

lib-shared-all: lazer.h liblazer.so
lib-shared: lazer.h liblazer.so
lib-static-all: lazer.h liblazer.a
lib-static: lazer.h liblazer.a

liblazer.a: src/lazer_static.o src/hexl_static.o
	ar rcs liblazer.a src/lazer_static.o src/hexl_static.o

liblazer.so: src/lazer_shared.o src/hexl_shared.o
	$(CC) $(CPPFLAGS) $(CFLAGS) -I. -shared -o liblazer.so src/lazer_shared.o src/hexl_shared.o

src/lazer_static.o: $(LIBSOURCES) $(LAZER_INC) $(LAZER_SRC) lazer.h src/lazer.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -I. -I$(LAZER_DIR) -I$(LAZER_DIR)/src -c -o src/lazer_static.o src/lazer.c

src/lazer_shared.o: $(LIBSOURCES) $(LAZER_INC) $(LAZER_SRC) lazer.h src/lazer.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -I. -I$(LAZER_DIR) -I$(LAZER_DIR)/src -c -fPIC -o src/lazer_shared.o src/lazer.c

src/hexl_static.o: $(LAZER_DIR)/src/hexl.h $(HEXL_DIR)
	$(CXX) $(CPPFLAGS) $(CFLAGS) -I$(LAZER_DIR) -I$(LAZER_DIR)/src -I$(HEXL_DIR)/hexl/include -c -o src/hexl_static.o $(LAZER_DIR)/src/hexl.cpp

src/hexl_shared.o: $(LAZER_DIR)/src/hexl.h $(HEXL_DIR)
	$(CXX) $(CPPFLAGS) $(CFLAGS) -I$(LAZER_DIR) -I$(LAZER_DIR)/src -I$(HEXL_DIR)/hexl/include -c -fPIC -o src/hexl_shared.o $(LAZER_DIR)/src/hexl.cpp

src/lazer.c:
	{ head -n33 $(LAZER_DIR)/src/lazer.c; cat src/lazer-modified.c; tail -n+34 $(LAZER_DIR)/src/lazer.c; } | sed 's/blindsig.c/stopwatch.h/' > src/lazer.c

lazer.h: $(LAZER_DIR)/src/lazer-in1.h $(LAZER_DIR)/src/lazer-in2.h src/lazer-modified.h $(LAZER_DIR)/src/moduli.h $(LAZER_DIR)/config.h
	cat $(LAZER_DIR)/src/lazer-in1.h > tmp.h
	echo "" >> tmp.h

	echo "#ifndef LAZER_CONFIG_H" >> tmp.h
	echo "#define LAZER_CONFIG_H" >> tmp.h
	echo "" >> tmp.h
	cat $(LAZER_DIR)/config.h >> tmp.h
	echo "" >> tmp.h
	echo "#endif" >> tmp.h
	echo "" >> tmp.h

	cat $(LAZER_DIR)/src/lazer-in2.h >> tmp.h
	head --lines=-2 tmp.h > lazer.h
	rm tmp.h

	echo "" >> lazer.h

	cat src/lazer-modified.h >> lazer.h
	echo "" >> lazer.h
	
	cat lazer/src/moduli.h >> lazer.h	
	
TESTDEPS = $(LAZER_DIR)/tests/test.h tests/test.o lazer.h liblazer.a
TESTLIBS = tests/test.o liblazer.a $(LIBS)

.PHONY: params

params: params-setup params-abdlop params-quad params-many

params-setup: 
	cd scripts && $(MAKE) $@ $<
	

params-abdlop: params-setup
	cd scripts && $(MAKE) $@ $<

params-quad: params-setup
	cd scripts && $(MAKE) $@ $<

params-many: params-setup
	cd scripts && $(MAKE) $@ $<

.PHONY: check
check: $(TESTS)
	cd tests && ./run-modified-tests

tests/test.o: $(LAZER_DIR)/tests/test.c $(LAZER_DIR)/tests/test.h lazer.h liblazer.a
	$(CC) $(CPPFLAGS) $(CFLAGS) -I$(LAZER_DIR) -I. -c -o $@ $<

tests/lazer-test: $(LAZER_DIR)/tests/lazer-test.c $(TESTDEPS)
	$(CC) $(CPPFLAGS) $(CFLAGS) -I$(LAZER_DIR) -I. -o $@ $< $(TESTLIBS)

tests/modified-abdlop-test: tests/modified-abdlop-test.c $(TESTDEPS) tests/modified-abdlop-params1.h tests/modified-abdlop-params2.h tests/modified-abdlop-params3.h tests/modified-abdlop-params4.h tests/modified-abdlop-params5.h
	$(CC) $(CPPFLAGS) $(CFLAGS) -I. -I$(LAZER_DIR) -I$(LAZER_DIR)/tests -o $@ $< $(TESTLIBS)

tests/modified-quad-test: tests/modified-quad-test.c $(TESTDEPS) tests/modified-quad-params1.h tests/modified-quad-params2.h tests/modified-quad-params3.h tests/modified-quad-params4.h tests/modified-quad-params5.h
	$(CC) $(CPPFLAGS) $(CFLAGS) -I. -I$(LAZER_DIR) -I$(LAZER_DIR)/tests -o $@ $< $(TESTLIBS)

tests/modified-quad-eval-test: tests/modified-quad-eval-test.c $(TESTDEPS) tests/modified-quad-eval-params1.h tests/modified-quad-eval-params2.h
	$(CC) $(CPPFLAGS) $(CFLAGS) -I. -I$(LAZER_DIR) -I$(LAZER_DIR)/tests -o $@ $< $(TESTLIBS)

tests/modified-quad-many-test: tests/modified-quad-many-test.c $(TESTDEPS) tests/modified-quad-params1.h tests/modified-quad-params2.h tests/modified-quad-params3.h tests/modified-quad-params4.h tests/modified-quad-params5.h
	$(CC) $(CPPFLAGS) $(CFLAGS) -I. -I$(LAZER_DIR) -I$(LAZER_DIR)/tests -o $@ $< $(TESTLIBS)

.PHONY: clean
clean:
	rm -f lazer.h src/lazer.c liblazer.a liblazer.so 
	cd scripts && rm -f moduli.sage.py modified-abdlop-codegen.sage.py modified-quad-codegen.sage.py modified-quad-eval-codegen.sage.py
	cd src && rm -f *.o
	cd lazer && $(MAKE) clean
	cd $(THIRD_PARTY_DIR) && rm -rf $(HEXL_SUBDIR)
	cd tests && rm -f *.o *.dSYM && cd .. && rm -f $(TESTS) && rm -f sage-test.sage.py  

