The batchfile in this directory, batch.txt, will attempt to reconcile
the following five gene trees with the species tree "speciestree":

 - b.b
 - doesnotexist
 - b.nb
 - nb.b
 - nb.nb

The gene tree "doesnotexist" does not exist, so when Notung tries to
read this file, it will report an error and continue.

Example command:
java -jar ../../Notung-2.9.jar -b batch.txt --reconcile --prune --outputdir testing --speciestag prefix 

The non-binary species tree contains two sets of species: {a, b, c, d}
and {w, x, y, z).  It has been constructed so that when it is pruned
to contain only species from one set, it is binary.  If it is pruned
to contain a mixture of species from both sets, it may be non-binary.

The names of the gene trees indicate the type of reconciliation that
Notung will attempt:

  - b.b : The binary gene tree b.b will prune the species tree to a binary
    	  species tree, and will reconcile using the standard method.

  - b.nb : The binary gene tree b.nb will prune the species tree to a non-binary
    	   species tree, and reconcile using our new non-binary species trees
    	   reconciliation method.

  - nb.b : The non-binary gene tree nb.b will prune the species tree to a binary
    	   species tree, and will reconcile using our resolve-and-collapse
    	   method.

  - nb.nb : The non-binary gene tree nb.nb will prune the species tree to a
    	    non-binary species tree, and will report an error when it is
    	    reconciled.

