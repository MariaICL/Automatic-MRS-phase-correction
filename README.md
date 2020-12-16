# Automatic-MRS-phase-correction

Automatic zero and first order phase order correction for MRS, implemented in Matlab

USE:

[out,outw,degree, N_flat]=Final_auto_phase_correction(in,inw);


DESCRIPTION:
Perform automatic zero and first order phase correction

INPUTS:
in     = water suppressed input data in matlab structure format.
inw    = water unsuppressed input data in matlab structure format.

OUTPUTS:
out    = Water suppressed output following phase correction  
outw   = Water unsuppressed output following phase correction
degree   = Phase correction in degrees following zero order phase correction
