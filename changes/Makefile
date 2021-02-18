CC = gcc -fopenmp
FFTW3_INC = -I/gpfs/home/ajamieson/local/include
FFTW3_LIB = -L/gpfs/home/ajamieson/local/lib
GL_INCL = -I/gpfs/software/gsl-2.40-gnu/include/gsl
GSL_LIBS = -L/gpfs/software/gsl-2.40-gnu/lib
LIB = $(GSL_LIBS) $(FFTW3_LIB) -lfftw3 -lfftw3_threads -lm -lgsl -lgslcblas
INCLUDE = $(FFTW3_INC) $(GSL_INCL) 
CFLAGS = -O2 $(INCLUDE)

TARGET = ic2lpt
BIN = /gpfs/home/ajamieson/local/bin
ic2lpt = write_power.o get_rhohat.o populate_rhohat.o power_interp.o bsrch.o ranc_ctr.o enforce_hermitian.o greens.o get_l2.o displace_2lpt.o writegadget.o ic2lpt.o
all: $(TARGET)

$(TARGET): $(ic2lpt)
	$(CC) -o $(TARGET) $(ic2lpt) $(LIB)

.PHONY: install
install:
	cp $(TARGET) $(BIN)

.PHONY: uninstall
uninstall:
	rm -f $(BIN)/$(TARGET)

.PHONY: tidy
tidy:
	rm -f *.o

.PHONY: clean
clean:
	rm -f $(TARGET) *.o
