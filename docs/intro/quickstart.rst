.. _intro-tutorial:

===========
Quick Start
===========

Installion
==========

See :doc:`install`.

Read in the gem file as Gem object
==================================

Import this module as:

>>> import spacipy.Gem as Gem

Examples:

>>> gem = Gem.readin('example.gem.gz')
>>> gem.head()

Convert gem object to anndata
=============================

>>> adata = gem.to_anndata(bin_size=50)
>>> adata

Plotting the gem object
=======================

>>> gem.plot(bin_size=50)

View casting mask array as Stoarr object
========================================

Import module and get example data

>>> import spacipy.Stoarr as Stoarr
>>> import numpy as np
>>> mask = np.load('cell_mask.txt', dtype='int')

Examples:

>>> mask = mask.view(Stoarr)
>>> mask

Slice the mask with nearest contour box
=======================================

>>> sliced_mask = mask.contour_slice()
>>> sliced_mask.shape

Imaging the stoarr object
=========================

>>> annotation = pd.DataFrame()
>>> annotation.head()

>>> sliced_mask.plot(annotation=annotation)





