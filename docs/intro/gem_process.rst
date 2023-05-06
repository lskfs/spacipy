.. _`gem_process`:

========================================
Handle gem files for downstream analysis
========================================

.. note:: 
    GEM is a tabular structure format specified for saving stereo-seq data, which generally contains columns of x, y, geneID, MIDCount and decodes the feature information captured by each DNB on the stereo-seq chip. Additional columns counld be included for DNB/feature annotation. A demo description can be found at `stereopy site <https://stereopy.readthedocs.io/en/latest/Tutorials/IO.html#GEM>`_.

.. code-block:: python3

    import spacipy.Gem as Gem

Read in example gem data, which holds mouse brain stereo-seq data from `here <https://www.stomics.tech/demoData>`_.

.. code-block:: python3

    gem = Gem.readin('example.gem.gz')

preliminary understanding of your data
======================================

GEM saves spatial transcriptomics data, which can be easily understanding by seeing what it looks like.

.. code-block:: python3

    gem.plot(color='nCount', bin_size=20)

here you may focus on two points:
    1. whether the spatial distributuion of the data is well decoded
    2. how many features were captured in each square bin or DNB

extract data where tissue covered
=================================

Not all DNBs on the stereo-seq chip were covered by the tissue. An effective way for deleting non-related data is to use a mask, which is often a two dimensional numpy ndarray with tissue regions labeled as non-zero integer (see `here <https://numpy.org/doc/stable/reference/maskedarray.generic.html#what-is-a-masked-array>`_). Once prepared a ready-to-use mask matrix, you can process as following:

.. code-block:: python3

    import numpy as np
    mask = np.load('./data/mask_matrix.npy')
    gem = gem.mask(matrix=mask)

now the gem data will be

.. code-block:: python3

    gem.plot(color='nCount', bin_size=20)

aggregate the data into cell/bin unit
=====================================

After getting the effective data set, the next step is to determine how to group the nanometer-resolution DNB array into a analysis-ready unit, here always in two way:
    1. group into square bin with equal width and height
    2. group into putative cells

1. group into square bin
************************

Square binning strategy do NOT consider any other informations like nucleus location and molecular homogeny. The only parameter in this step is *bin_size*, which specify how many DNBs in width and height should be grouped as a single bin unit. People can change this value based on the physical area of each bin and the captured features in each bin unit. 

.. code-block:: python3

    binning_gem = gem.binning(bin_size=20)

2. group into putative cells
****************************

Currently stereo-seq comprehend ssDNA stainning protocal, which stains the nucleus of cells on the same section. With the location of cell nucleus decoded, labeled 



