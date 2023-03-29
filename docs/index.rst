.. libquickbam documentation master file, created by
   sphinx-quickstart on Thu Aug 26 17:24:10 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


Introduction to libquickbam
===========================

libquickbam offers stunningly fast reading of indexed BAM files, taking advantage
of parallel data processing. Together with modern, high throughput hardware
such as NVME SSDs or Lustre, libquickbam enables gigabytes per second read
speed, and allows for existing bioinformatics ananlysis algorithms to be
accelerated significantly, or a new generation of analysis software to be
created, for rapid analysis turnaround.

In addition to sequence read iteration, libquickbam also offers a multi-file
pileup engine to enable parallel piling up. An example of utilizing this
function is the ``snp-pileup-quickbam`` example program.

.. toctree::
   :maxdepth: 2
   :caption: Getting started

   getting_started

.. toctree::
   :maxdepth: 2
   :caption: Parallel processing

   parallel

.. toctree::
   :maxdepth: 2
   :caption: API Reference

   api



..
  Indices and tables

..
  ==================

..
  * :ref:`genindex`

..
  * :ref:`modindex`

..
  * :ref:`search`
