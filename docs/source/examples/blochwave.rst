==========
Blochwave
==========

The documentation for the formulation can be found elsewhere. These are simple examples
on how to run the library.

.. note::

   If you want to try them out quickly without installing the library, these examples
   can be tested directly on Jupyter notebooks here.

Basic usage
*******************
.. code-block:: python

   from blochwave import bloch
   b0 = bloch.Bloch('diamond',path='',u=[0,0,1],Nmax=3,Smax=0.1,opts='svt')

   b0.show_beam_vs_thickness()
   b0.convert2tiff()
