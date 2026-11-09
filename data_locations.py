"""Where the raw data are on your computer. Edit this file before the first run.

Only the two preparation stages read raw files; the Results read the prepared
tables they make. Each line can also be set by the environment variable named
on it, which takes precedence over the value written here.

Where to get each one (all free, public-use files from the US National Center
for Health Statistics):
  NHANES 1999-2018       https://wwwn.cdc.gov/nchs/nhanes/  (the .XPT files of each
                         cycle, one folder per cycle: 1999-2000/, 2001-2002/, ...;
                         1_prepare_nhanes/config/*_sources.csv lists every file)
  Linked mortality       https://www.cdc.gov/nchs/data-linkage/mortality-public.htm
                         (the 2019 public-use files, NHANES_<cycle>_MORT_2019_PUBLIC.dat)
  NHANES III (1988-94)   https://wwwn.cdc.gov/nchs/nhanes/nhanes3/ (adult, exam and lab
                         .dat files with their .sas layouts, and NHANES_III_MORT_2019_PUBLIC.dat)
  HEI-2020 scoring       the CRAN archives heiscore_0.1.4.tar.gz and heiscore.data_0.0.1.tar.gz
                         (https://cran.r-project.org/src/contrib/Archive/)
"""
import os

NHANES_RAW = os.environ.get("NHANES_RAW_ROOT") or \
    r"EDIT_ME/folder_with_NHANES_1999_2018_files"
NHANES_MORTALITY = os.environ.get("NHANES_MORTALITY_ROOT") or \
    r"EDIT_ME/folder_with_NHANES_linked_mortality_files"
HEI_PACKAGES = os.environ.get("NHANES_HEI_PACKAGES") or \
    r"EDIT_ME/folder_with_heiscore_CRAN_archives"
NHANES3_RAW = os.environ.get("NHANES3_RAW_ROOT") or \
    r"EDIT_ME/folder_with_NHANES_III_files"
