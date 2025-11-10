################################################################################
# Script for loading the species vital rates
#
# Written by Young Jun Lee
# Oct 2024

sp1_data <- read_csv(here("StageDependentComp","result_data","wp2","sp1_vr.csv"),
                     col_types = cols(Plot = col_character()))

sp2_data <- read_csv(here("StageDependentComp","result_data","wp2","sp2_vr.csv"),
                     col_types = cols(Plot = col_character()))

spp_stage_struct <- read_csv(here("StageDependentComp","result_data","wp2","sp_stage_struct.csv"))

sp1_thres_len <- spp_stage_struct$Stage.boundary..cm[1]
sp2_thres_len <- spp_stage_struct$Stage.boundary..cm[2]
