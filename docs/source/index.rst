Multiphase Test-Bench - A slim CFD framework for multiphase flows
=================================================================

..
    Moin Moin!

    **Lumache** (/lu'make/) is a Python library for cooks and food lovers
    that creates recipes mixing random ingredients.
    It pulls data from the `Open Food Facts database <https://world.openfoodfacts.org/>`_
    and offers a *simple* and *intuitive* API.

    Check out the :doc:`usage` section for further information, including
    how to :ref:`installation` the project.

.. note::

   This project is under active development.

Contents
--------

.. toctree::
   
   usage
   api
   examples

```{image} _static/logo.png
---
width: 675px
align: center
---
```
**poliastro** is an open source ([MIT](https://opensource.org/licenses/MIT)) pure Python library
for interactive Astrodynamics and Orbital Mechanics,
with a focus on ease of use, speed, and quick visualization.
It provides a simple and intuitive {ref}`API <api-reference>`,
and handles physical quantities with units.

View the [source code](https://github.com/poliastro/poliastro) of poliastro!

Some of its awesome features are:

- Analytical and numerical orbit propagation
- Conversion between position and velocity vectors and classical
  orbital elements
- Coordinate frame transformations
- Hohmann and bielliptic maneuvers computation
- Trajectory plotting
- Initial orbit determination (Lambert problem)
- Planetary ephemerides (using [SPICE kernels](https://naif.jpl.nasa.gov/naif/data.html) via [Astropy](https://www.astropy.org/))
- Computation of Near-Earth Objects (NEOs)

And more to come!

poliastro is developed by an open, international community. Release
announcements and general discussion take place on our [mailing
list](https://groups.io/g/poliastro-dev) and
[chat](http://chat.poliastro.space/).

..
    ```{eval-rst}
    .. raw:: html

        <div class="classictemplate template" style="display: block;">
        <style type="text/css">
          #groupsio_embed_signup input { border:1px solid #999; -webkit-appearance:none; }
          #groupsio_embed_signup .email { display:block; padding:8px 0; margin:0 4% 10px 0; text-indent:5px; width:58%; min-width:130px; }
          #groupsio_embed_signup { background:#fff; clear:left; font:14px Helvetica,Arial,sans-serif; }
          #groupsio_embed_signup .button {
              width:25%; margin:0 0 10px 0; min-width:90px;
              background-image: linear-gradient(to bottom,#337ab7 0,#265a88 100%);
              background-repeat: repeat-x;
              border-color: #245580;
              text-shadow: 0 -1px 0 rgba(0,0,0,.2);
              box-shadow: inset 0 1px 0 rgba(255,255,255,.15),0 1px 1px rgba(0,0,0,.075);
              padding: 5px 10px;
              font-size: 12px;
              line-height: 1.5;
              border-radius: 3px;
              color: #fff;
              background-color: #337ab7;
              display: inline-block;
              margin-bottom: 0;
              font-weight: 400;
              text-align: center;
              white-space: nowrap;
              vertical-align: middle;
              position: relative;
              top: -15px;
          }
        </style>
        <div id="groupsio_embed_signup">
        <form action="https://groups.io/g/poliastro-dev/signup?u=3341734695881558152" method="post" id="groupsio-embedded-subscribe-form" name="groupsio-embedded-subscribe-form" target="_blank">
            <div id="groupsio_embed_signup_scroll">
              <input value="" name="email" class="email" id="email" placeholder="email address" required="" type="email">

            <div style="position: absolute; left: -5000px;" aria-hidden="true"><input name="b_3341734695881558152" tabindex="-1" value="" type="text"></div>
            <div id="templatearchives"><p></p></div>
            <input value="Subscribe to mailing list! 🚀" name="subscribe" id="groupsio-embedded-subscribe" class="button" type="submit">
          </div>
        </form>
        </div>
        </div>
    ```
..
    ```{figure} _static/molniya.png
    ---
    align: right
    width: 300px
    alt: Molniya Orbit
    ---
    Plot of a [Molniya orbit](https://en.wikipedia.org/wiki/Molniya_orbit) around the Earth <br />
    ({math}`a = 26600~\mathrm{km}, e = 0.75, i = 63.4\mathrm{^\circ}`).
    ```

The [source code](https://github.com/poliastro/poliastro), [issue
tracker](https://github.com/poliastro/poliastro/issues) and
[wiki](https://github.com/poliastro/poliastro/wiki/) are hosted on
GitHub, and all contributions and feedback are more than welcome. You
can test poliastro in your browser using [Binder](https://mybinder.org/), a cloud Jupyter
notebook server:

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/poliastro/poliastro/main?labpath=index.ipynb)

See [benchmarks](https://benchmarks.poliastro.space/) for the
performance analysis of poliastro.

poliastro works on the recent Python versions and is released under the
MIT license, allowing commercial use of the library.

```python
from poliastro.examples import molniya
molniya.plot()
```
------------------------------------------------------------------------

..
    ## Contents

    ```{toctree}
    ---
    maxdepth: 2
    caption: Tutorials
    ---
    installation
    quickstart
    ```

    ```{toctree}
    ---
    maxdepth: 2
    caption: How-to guides & Examples
    ---
    gallery
    contributing
    ```

    ```{toctree}
    ---
    maxdepth: 2
    caption: Reference
    ---
    api
    bibliography
    changelog
    ```

    ```{toctree}
    ---
    maxdepth: 2
    caption: Background
    ---
    history
    related
    background
    ```


