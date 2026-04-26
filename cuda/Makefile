## UniversalFFT — CUDA / HIP Makefile
##
## Usage:
##   make           # build CUDA binary (requires nvcc + cuFFT)
##   make HIP=1     # build HIP binary  (requires hipcc + hipFFT)
##   make test      # build + run demo

ifdef HIP
  CC       = hipcc
  LDFLAGS  = -lhipfft
  TARGET   = ufft_demo_hip
  DEFS     = -DUFFT_USE_HIP
else
  CC       = nvcc
  LDFLAGS  = -lcufft
  TARGET   = ufft_demo_cu
  DEFS     =
endif

CFLAGS   = -O2

.PHONY: all test clean

all: $(TARGET)

$(TARGET): demo.cu universalfft.cu universalfft.cuh
	$(CC) $(CFLAGS) $(DEFS) -o $@ demo.cu universalfft.cu $(LDFLAGS)

test: $(TARGET)
	mkdir -p ../tests/data
	./$(TARGET) ../tests/data

clean:
	rm -f ufft_demo_cu ufft_demo_hip *.csv
