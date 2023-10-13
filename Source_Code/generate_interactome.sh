#!/bin/bash
./updat_PDB.sh
python3 -m analysis
python3 -m PPI_extraction
python3 -m main