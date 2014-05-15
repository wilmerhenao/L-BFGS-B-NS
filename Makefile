FC = gfortran

#FFLAGS = -O -Wall -fbounds-check -g -Wno-uninitialized 
FFLAGS = -g -pg

DRIVER1_90 = DriverRosenbrockNonSmooth.f90
DRIVER2_90 = Driveryurirosen_ns1.f90
DRIVER3_90 = Driveryurirosen_ns2.f90
DRIVER4_90 = driver3.f90
DRIVER5_90 = chained_CB3v1.f90
DRIVER6_90 = chained_CB3v2.f90
DRIVER7_90 = chained_crescent1.f90
DRIVER8_90 = chained_crescent2.f90
DRIVER9_90 = chained_LQ.f90
DRIVER10_90 = chained_mifflin2.f90
DRIVER11_90 = gen_brownfunc2.f90
DRIVER12_90 = gen_maxhilb.f90
DRIVER13_90 = gen_maxq.f90
DRIVER14_90 = nactfaces.f90
DRIVER15_90 = Driveryurirosen_ns1_p.f90
DRIVER16_90 = DriverRosenbrockp.f90
DRIVER18_90 = DriverRosenbrockprintingp.f90
DRIVER19_90 = gen_maxqPrinting.f90
DRIVER20_90 = mifflin2Printing.f90

LBFGS = lbfgsb.f
LBFGSB  = lbfgsbNS.f90
LBFGSBSTRONG = lbfgsbnomessagesStrong.f90
LINPACK = linpack.f
BLAS    = blas.f
TIMER   = timer.f
LAPACKFLAG = -llapack
BLASFLAG = -lblas

#all :  lbfgsb_90_1 lbfgsb_90_2 lbfgsb_90_3 lbfgsb_90_4 lbfgsb_90_5 lbfgsb_90_6 lbfgsb_90_7 lbfgsb_90_8 lbfgsb_90_9 lbfgsb_90_10 lbfgsb_90_11 lbfgsb_90_12 lbfgsb_90_13 lbfgsb_90_14 lbfgsb_90_15 lbfgsb_90_16
all : lbfgsb_90_16

lbfgsb_90_1 : $(DRIVER1_90) $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER)
	$(FC) $(FFLAGS) $(DRIVER1_90) $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER) $(LAPACKFLAG) $(BLASFLAG) -o rosenbrocknonsmooth

lbfgsb_90_2 : $(DRIVER2_90) $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER)
	$(FC) $(FFLAGS) $(DRIVER2_90) $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER) $(LAPACKFLAG) $(BLASFLAG) -o yurirosen_ns1

lbfgsb_90_3 : $(DRIVER3_90) $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER)
	$(FC) $(FFLAGS) $(DRIVER3_90) $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER) $(LAPACKFLAG) $(BLASFLAG) -o yurirosen_ns2

lbfgsb_90_4 : $(DRIVER4_90) $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER)
	$(FC) $(FFLAGS) $(DRIVER4_90) $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER) $(LAPACKFLAG) $(BLASFLAG) -o driver3

lbfgsb_90_5 : $(DRIVER5_90) $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER)
	$(FC) $(FFLAGS) $(DRIVER5_90) $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER) $(LAPACKFLAG) $(BLASFLAG) -o chained_CB3v1

lbfgsb_90_6 : $(DRIVER6_90) $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER)
	$(FC) $(FFLAGS) $(DRIVER6_90) $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER) $(LAPACKFLAG) $(BLASFLAG) -o chained_CB3v2

lbfgsb_90_7 : $(DRIVER7_90) $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER)
	$(FC) $(FFLAGS) $(DRIVER7_90) $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER) $(LAPACKFLAG) $(BLASFLAG) -o chained_crescent1

lbfgsb_90_8 : $(DRIVER8_90) $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER)
	$(FC) $(FFLAGS) $(DRIVER8_90) $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER) $(LAPACKFLAG) $(BLASFLAG) -o chained_crescent2

lbfgsb_90_9 : $(DRIVER9_90) $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER)
	$(FC) $(FFLAGS) $(DRIVER9_90) $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER) $(LAPACKFLAG) $(BLASFLAG) -o chained_LQ

lbfgsb_90_10 : $(DRIVER10_90) $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER)
	$(FC) $(FFLAGS) $(DRIVER10_90) $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER) $(LAPACKFLAG) $(BLASFLAG) -o chained_mifflin2

lbfgsb_90_11 : $(DRIVER11_90) $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER)
	$(FC) $(FFLAGS) $(DRIVER11_90) $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER) $(LAPACKFLAG) $(BLASFLAG) -o gen_brownfunc2

lbfgsb_90_12 : $(DRIVER12_90) $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER)
	$(FC) $(FFLAGS) $(DRIVER12_90) $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER) $(LAPACKFLAG) $(BLASFLAG) -o gen_maxhilb

lbfgsb_90_13 : $(DRIVER13_90) $(LBFGS) $(LINPACK) $(BLAS) $(TIMER)
	$(FC) $(FFLAGS) $(DRIVER13_90) $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER) $(LAPACKFLAG) $(BLASFLAG) -o gen_maxq

lbfgsb_90_14 : $(DRIVER14_90) $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER)
	$(FC) $(FFLAGS) $(DRIVER14_90) $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER) $(LAPACKFLAG) $(BLASFLAG) -o nactfaces

lbfgsb_90_15 : $(DRIVER15_90) $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER)
	$(FC) $(FFLAGS) $(DRIVER15_90) $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER) $(LAPACKFLAG) $(BLASFLAG) -o yurirosenp

lbfgsb_90_16 : $(DRIVER16_90) $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER)
	$(FC) $(FFLAGS) $(DRIVER16_90) $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER) $(LAPACKFLAG) $(BLASFLAG) -o rosenbrockp

lbfgsb_90_17 : $(DRIVER16_90) $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER)
	$(FC) $(FFLAGS) $(DRIVER16_90) $(LBFGSBSTRONG) $(LINPACK) $(BLAS) $(TIMER) $(LAPACKFLAG) $(BLASFLAG) -o rosenbrockStrongp

lbfgsb_90_18 : $(DRIVER18_90) $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER)
	$(FC) $(FFLAGS) $(DRIVER18_90) $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER) $(LAPACKFLAG) $(BLASFLAG) -o rosenbrockprintingp

lbfgsb_90_19 : $(DRIVER19_90) $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER)
	$(FC) $(FFLAGS) $(DRIVER19_90) $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER) $(LAPACKFLAG) $(BLASFLAG) -o maxq

lbfgsb_90_20 : $(DRIVER20_90) $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER)
	$(FC) $(FFLAGS) $(DRIVER20_90) $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER) $(LAPACKFLAG) $(BLASFLAG) -o mifflin2
