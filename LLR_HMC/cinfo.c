static char MACROS[     145+1] = "-DBC_T_PERIODIC -DBC_X_PERIODIC -DBC_Y_PERIODIC -DBC_Z_PERIODIC -DNDEBUG -DIO_FLUSH -DREPR_FUNDAMENTAL -DREPR_NAME=@REPR_FUNDAMENTAL@ -DGAUGE_SUN";

static char CI_mkflags[    1046+1] = "NG = 2\nREPR = REPR_FUNDAMENTAL\n#REPR = REPR_SYMMETRIC\n#REPR = REPR_ANTISYMMETRIC\n#REPR = REPR_ADJOINT\n\n#CHOICES ARE GAUGE_SUN AND GAUGE_SON \nGAUGE_GROUP = GAUGE_SUN\n\n#T => PERIODIC, ANTIPERIODIC, OPEN, THETA\n#X => PERIODIC, ANTIPERIODIC, THETA\n#Y => PERIODIC, ANTIPERIODIC, THETA\n#Z => PERIODIC, ANTIPERIODIC, THETA\n\n#MACRO += -DBC_T_THETA\nMACRO += -DBC_T_PERIODIC\nMACRO += -DBC_X_PERIODIC\nMACRO += -DBC_Y_PERIODIC\nMACRO += -DBC_Z_PERIODIC\n\n#MACRO += -DBC_XYZ_TWISTED\n\n#MACRO += -DHALFBG_SF\n#MACRO += -DBASIC_SF\n#MACRO += -DROTATED_SF\n\n\n#MACRO += -DUPDATE_EO\n#MACRO += -DWITH_MPI\n#MACRO += -DWITH_UMBRELLA\n#MACRO += -DWITH_QUATERNIONS\nMACRO += -DNDEBUG\n#MACRO += -DCHECK_SPINOR_MATCHING\n#MACRO += -DMPI_TIMING\nMACRO += -DIO_FLUSH\n#MACRO += -DUNROLL_GROUP_REPRESENT\n#MACRO += -DTIMING\n#MACRO += -DTIMING_WITH_BARRIERS\n#MACRO += -DAMALLOC_MEASURE\n\n#Compiler\nCC = /opt/MPI/bin/mpicc\n#CC = clang\nCFLAGS = -Wall -std=c99 -O3 -fomit-frame-pointer -msse -msse2 -msse3\n#CFLAGS = -Wall -std=c99 -g\nLDFLAGS = \n";

static char CI_cpuinfo[      13+1] = "No CPU info\n";

static char CI_linux[      17+1] = "No VERSION info\n";

static char CI_gcc[     425+1] = "Configured with: --prefix=/Applications/Xcode.app/Contents/Developer/usr --with-gxx-include-dir=/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.14.sdk/usr/include/c++/4.2.1\nApple LLVM version 10.0.0 (clang-1000.11.45.2)\nTarget: x86_64-apple-darwin17.7.0\nThread model: posix\nInstalledDir: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin\n";

static char CI_svninfo[     514+1] = "Path: /Users/djl/Physics/Plymouth/SU3/HiRep-LLR\nWorking Copy Root Path: /Users/djl/Physics/Plymouth/SU3/HiRep-LLR\nURL: svn+ssh://llr_hmc@py14240.swan.ac.uk/home/llr_hmc/SVN/HiRep-LLR\nRelative URL: ^/\nRepository Root: svn+ssh://llr_hmc@py14240.swan.ac.uk/home/llr_hmc/SVN/HiRep-LLR\nRepository UUID: 29680204-7f03-47a5-b968-235bbe8675bf\nRevision: 16\nNode Kind: directory\nSchedule: normal\nLast Changed Author: llr_hmc\nLast Changed Rev: 16\nLast Changed Date: 2017-04-26 16:00:11 +0100 (Wed, 26 Apr 2017)\n\n";

static char CI_svnstatus[     884+1] = "M       ../HMC/hmc.c\n!       ../HiRep.xcodeproj\n!       ../HiRep.xcodeproj/project.pbxproj\nM       ../Include/communications.h\nM       ../Include/global.h\nM       ../Include/update.h\n!       ../LLR_HMC/input_file\n!       ../LLR_HMC/input_file_llr_hmc_Nc_2\nM       ../LLR_HMC/llr_hmc.c\nM       ../LibHR/Geometry/communications.c\nM       ../LibHR/IO/read_action.c\nM       ../LibHR/Update/force0.c\nM       ../LibHR/Update/force_hmc.c\nM       ../LibHR/Update/integrators.c\nM       ../LibHR/Update/llr.c\nM       ../LibHR/Update/local_action.c\nM       ../LibHR/Update/mon_llr_gauge.c\nM       ../LibHR/Update/mon_llr_hmc.c\nM       ../LibHR/Update/mon_pg.c\nM       ../LibHR/Update/update_mt.c\nM       ../LibHR/Utils/background_field.c\nM       ../LibHR/Utils/boundary_conditions.c\nM       ../LibHR/Utils/eva_deflation.c\nM       ../Make/MkFlags\nM       ../Make/MkRules\n";

static int CI_svnrevision = 16;

#include <stdio.h>

static void print_compiling_info() {
  printf("MkFlags\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
  printf("%s",CI_mkflags);
  printf("\n\n");
  printf("svn info\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
  printf("%s",CI_svninfo);
  printf("\n\n");
  printf("svn st -q\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
  printf("%s",CI_svnstatus);
  printf("\n\n");
  printf("/proc/cpuinfo\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
  printf("%s",CI_cpuinfo);
  printf("\n\n");
  printf("/proc/version\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
  printf("%s",CI_linux);
  printf("\n\n");
  printf("gcc -v\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
  printf("%s",CI_gcc);
}
