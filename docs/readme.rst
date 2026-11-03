Gaussian process Optimized Photometric Regression of Extragalactic Archival Ultraviolet-infrared eXplosions

A package for multi-dimensional Gaussian Process Regression of extragalactic astronomical transient light curves, enabling a full characterization of their spectral energy distribution evolution with time.

Background
----------

Time-domain astronomy is entering a golden age of discovery, powered by wide-field surveys such as Rubin Observatory's Legacy Survey of Space and Time. The wealth of data that will be produced over the coming years offers time-domain astronomers the opportunity to conduct population-level studies of different types of explosive extragalactic transients, many for the first time. However, our knowledge of the physics underpinning these explosions is lagging, and parsing the deluge of data in real time to identify interesting transients requires detailed knowledge of their time evolution.

``GOPREAUX`` addresses these problems by producing multidimensional template light curve and spectral surfaces of different classes of transients for the first time. ``GOPREAUX`` models aggregated archival data, spanning the ultraviolet to the infrared, using Gaussian Process Regression across both phase and wavelength simultaneously. The models produce "template" time-evolving spectral energy distribution surfaces that can be used for a variety of use cases—such as performing physical parameter inference, generating machine learning training sets of different transients at arbitrary phase and redshift, or identifying rare and unusual transients in real time.

Installation and Setup
----------------------

Dependencies are managed using ``poetry``. The recommended installation is to create a new Python environment for this repository and install ``poetry`` within that environment. Poetry can then be used to install the dependencies::

    conda create -n gopreaux python=3.10
    conda activate gopreaux
    pip install poetry
    poetry install

Finally, to make use of Milky Way extinction correction, the dust map files of Schlafly and Finkbeiner (2011) must be fetched and saved locally. To do so, open a Python shell and run:

.. code-block:: python

    import dustmaps.sfd
    dustmaps.sfd.fetch()

That's it! ``gopreaux`` should now be ready to use.

Code Example
----------------------

.. code-block:: python

    from caat import SN, SNModel

    # Load the GP model
    model = SNModel(
        surface="SESNe_SNIIb_GP_model.fits"
    )

    # Predict a light curve from -20 to 45 days at 5000 Angstroms. This returns the
    # phases (days), the magnitudes relative to the light curve peak, and their uncertainties
    phases, mags, mag_errs = model.predict_lightcurve(-20.0, 45.0, 5000, show=False)

    # Compare the prediction to real photometry from a supernova
    gkg = SN(name="SN2016gkg")
    model.compare_lightcurve_with_photometry(gkg, filt="V")

The ``predict_lightcurve``, ``predict_sed`` and ``predict_photometry_points`` methods all return
magnitudes relative to the light curve peak (so a value of -2 is two magnitudes fainter than peak),
with uncertainties in magnitudes.

More examples can be found as Jupyter notebooks in :doc:`Tutorials </tutorials>`.

Citations
----------------------

If your work makes use of the `GOPREAUX` software or data reductions, please cite 
our paper which is available on the arXiv:


.. code-block::

    @ARTICLE{2026arXiv260403372P,
        author = {{Pellegrino}, C. and {Pritchard}, T.~A. and {Modjaz}, M. and {Crawford}, A. and {Khakpash}, S. and {Bianco}, F.},
            title = "{GOPREAUX I: Open-source Code and Data to Model Multi-wavelength Emission of Extragalactic Transients using Gaussian Processes}",
        journal = {arXiv e-prints},
        keywords = {Instrumentation and Methods for Astrophysics, High Energy Astrophysical Phenomena},
            year = 2026,
            month = apr,
            eid = {arXiv:2604.03372},
            pages = {arXiv:2604.03372},
            doi = {10.48550/arXiv.2604.03372},
    archivePrefix = {arXiv},
        eprint = {2604.03372},
    primaryClass = {astro-ph.IM},
        adsurl = {https://ui.adsabs.harvard.edu/abs/2026arXiv260403372P},
        adsnote = {Provided by the SAO/NASA Astrophysics Data System}
    }

as well as the sources of any data used in your analysis.