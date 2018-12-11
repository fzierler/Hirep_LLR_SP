TOPDIR = ..
MKDIR = $(TOPDIR)/Make

EXES = hmc hmc_forces hmc_wf

hmc_OBJS = hmc_utils.o ../Spectrum/meson_measurements.o ../GaugeFix/gaugefix.o

hmc_wf_OBJS = hmc_utils.o ../Spectrum/meson_measurements.o ../GaugeFix/gaugefix.o  ../WilsonFlow/wilsonflow.o

hmc_forces_OBJS =  hmc_utils.o


LIBS += $(TOPDIR)/LibHR/libhr.a


include $(MKDIR)/MkRules

