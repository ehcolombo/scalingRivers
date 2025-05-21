# scalingRivers
Generation of OCNs matching a real-world counterpart.
getBasin_metrics is a simple example where the OCN has area and number of nodes matching those from a basin.
Basins are identified by their row in basins_info. Different trials can be generated. The script needs OCNet and other dependencies.

Usage:

Rscript basin_id trial

basin_id goes from 1 to 1139; trial is integer connected to the seed.
When finished the script print the metrics and save a RDS file.
