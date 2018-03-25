*****************************************
*************** README ******************
*****************************************

This folder contains a way for the parralel computation of a Validation Data.

Sequence of operations:

1. split an existing ExperimentData.mat using: breakExData.m
2. run: `create_validation_data.sh $result_folder_from_above`
3. assemble the different validation data splited in (2), it was splited to 1 target in each experimentData.mat, using: assemble_validationData.m .
