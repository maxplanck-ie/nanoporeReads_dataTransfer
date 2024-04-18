#!/bin/bash

INSTALLER='mamba'

for YAML in src/npr/rules_dorado/envs/*.yaml
do
    echo "Creating conda environment for $YAML"
    $INSTALLER env create -f $YAML
done

echo "Done!"