.. _topics-index:

.. include:: _links.rst

.. role:: small

.. role:: smaller

=====================
Spacipy documentation
=====================

Spacipy is collection of in-house python modules, targeted at the stereo-seq, that aims to make basic data analysis, modeling and visualization easier. It builds on the well-known NumPy, Pandas and MatPlotLib packages.

    * Quickly obtain data
    * Easily convert format
    * Create well formatted plots
    * Change coordinates effortlessly


.. note:: 
    You may need to be familiar with some image processing like crop, segment to use some modules.

Quick start
===========

.. toctree:: 
	:caption: Quick start
	:hidden: 
	:maxdepth: 1

	intro/install
	intro/basicusage

* :doc:`/intro/install`
* :doc:`/intro/basicusage`

Tutorial
========

.. toctree::
	:caption: Tutorial
	:hidden: 
	:maxdepth: 1
	
	intro/gem_process
	intro/stoarr_process
	intro/mip_plot

* :doc:`/intro/gem_process`
* :doc:`/intro/stoarr_process`
* :doc:`/intro/mip_plot`

API
===

.. toctree::
	:caption: API
	:hidden: 
	:maxdepth: 1

	gem
	stoarr
	coords3d
	mipplot

* :doc:`/gem`: *object constructed based on pandas.DataFrame with adapted method for stereo-seq GEM file.*
* :doc:`/stoarr`: *numpy.ndarray view casting with adapted method for spatial matrix*
* :doc:`/mipplot`: *MIP plot*

