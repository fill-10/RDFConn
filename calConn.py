from class_connectivity import connectivity

conn = connectivity('7', 1.25)
# name of atom, cutoff value

conn.loadfile('step72M.lammpstrj')
# input file

conn.findcluster()

conn.caldistro('clustersize.out')
# output file

# output of atoms to file is not completed
# need to decide what kind of format
# may use lammpstrj as the format, or the pdb
# to check the clusters, run this script in the interactive mode.
# data is in conn.clusters
