# Primerize (NA_Thermo)

<img src="https://primerize.stanford.edu/site_media/images/logo_primerize.png" alt="Primerize Logo" width="200" align="right">

**Primerize** (previously named **NA_thermo**), is an archive of *Python* and *MATLAB* scripts for primer design and nucleic acid thermodynamic scripts developed by the [Das Lab](https://daslab.stanford.edu/) at Stanford University for high-throughput RNA synthesis and design.

The algorithm designs *forward* (sense strand) and *reverse* (anti-sense strand) primers that minimize the total length, and therefore the total synthesis cost, of the oligonucleotides. Although developed independently, **Primerize** is a special case of the general ‘*Gapped Oligo Design*’ algorithm, optimizing the mispriming score and sequence span instead of *T<sub>m<sub>*.

An online user-friendly GUI is available as the [**Primerize Server**](https://primerize.stanford.edu/).

## Installation

To install **Primerize**, simply:
```bash
cd path/to/primerize/repo
python setup.py install
```

For system-wide installation, you must have permissions and use with `sudo`.

**Primerize** requires the following *Python* packages as dependencies, all of which can be installed through [`pip`](https://pip.pypa.io/).
```json
matplotlib >= 1.5.0
numba == 0.22.1
numpy >= 1.10.1
xlwt >= 1.0.0
```

Note that the [`numba`](http://numba.pydata.org/) is used for its [`@jit`](http://numba.pydata.org/numba-doc/0.22.1/user/jit.html) decorator on loop optimization. `numba` requires [`llvm`](http://llvm.org/), which can be installed through [`apt-get`](https://help.ubuntu.com/lts/serverguide/apt-get.html) on *Linux* or [`brew`](http://brew.sh/) on Mac *OSX*.


## Usage

For simple Primer Design tasks, follow this example:

```python
import primerize

prm_1d = primerize.Primerize_1D()
job_1d = prm_1d.design('TTCTAATACGACTCACTATA...AAAAAGAAACAACAACAACAAC', MIN_TM=60.0, NUM_PRIMERS=None, MIN_LENGTH=15, MAX_LENGTH=60, prefix='P4P6_2HP')
if job_1d.is_success:
	print job_1d

prm_2d = primerize.Primerize_2D()
job_2d = prm_2d.design('TTCTAATACGACTCACTATA...AAAAAGAAACAACAACAACAAC', offset=-51, which_muts=range(102, 261 + 1), which_libs=[1], prefix='P4P6_2HP')
if job_2d.is_success:
	print job_2d
	job_2d.save('table')
	job_2d.save('image')
	job_2d.save('construct')
```

For advanced users, the returned `Design_1D` and `Design_2D` result classes offer methods for `get()`, `save()` and `echo()`:

```python
MIN_TM = job_1d.get('MIN_TM')
print job_1d.get('MISPRIME')
print job_1d.echo('WARNING')
if job_1d.is_success:
	job_1d.save(path='result/', name='Primer')

LIB = job_2d.get('which_libs')
N_PRIMER = job_2d.get('N_PRIMER')
print job_2d.get('CONSTRUCT')
if job_2d.is_success:
	print job_2d.echo()
	job_2d.save('assembly', path='result/', name='Lib')
```

Besides `design()`, the `Primerize_1D` and `Primerize_2D` worker classes offer methods for `get()`, `set()`, and `reset()`:

```python
COL_SIZE = prm_1d.get('COL_SIZE')
prm_1d.set('MIN_LENGTH', 30)
prm_1d.reset()
```

There are also `Assembly` and `Plate_96Well` helper classes. For more details, please refer to the **Documentation**.


#### MATLAB Code (Deprecated)

Instructions on MATLAB usage is available at old [README.md](https://github.com/DasLab/Primerize/blob/master/deprecated_MATLAB/README.md).

## Documentation

Documentation is available at https://primerize.stanford.edu/docs/.

## License

Copyright &copy; of **Primerize** source code is described in [LICENSE.md](https://github.com/DasLab/Primerize/blob/master/LICENSE.md).

## Reference

>Tian, S., *et al.* (**2015**)<br/>
>[Primerize: Automated Primer Assembly for Transcribing Interesting RNAs.](http://nar.oxfordjournals.org/content/43/W1/W522.full)<br/>
>*Nucleic Acid Research* **43 (W1)**: W522-W526.


>Thachuk, C., and Condon, A. (**2007**)<br/>
>[On the Design of Oligos for Gene Synthesis.](http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=4375554)<br/>
>*Proceedings of the 7th IEEE International Conference on Bioinformatics and Bioengineering* **2007**: 123-130.

<br/>
by [**t47**](http://t47.io/), *Jan 2016*.

