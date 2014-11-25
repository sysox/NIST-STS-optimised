NIST-STS-optimised
==================

Currently 2 main versions:
  Version 1 - 30x faster than original used, Look Up tables, fast histogram computation, Berlekamp Massey alg. optimised
              using word - word (32bit or 64bit ) operations 
  Version 2 - improved 2x Berlekamp Massey - computes just linear complexity of the sequence, not costruction of  LSFR.
            - incorporated FFTW for Spectral test - better for larger factor of sequence size 
              for sequence size of the form  with small factors n = 2^k or n = 10^k ... is slightly better original FFT from NIST 2.2.1 
              for larger factors of n (103, 1007, ...)  FFTW is significantly better 
