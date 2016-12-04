# Primerize (NA_Thermo)

<img src="https://primerize.stanford.edu/site_media/images/logo_primerize.png" alt="Primerize Logo" width="200" align="right">

**Primerize** (previously named **NA_thermo**), is an archive of *Python* and *MATLAB* scripts for primer design and nucleic acid thermodynamic scripts developed by the [Das Lab](https://daslab.stanford.edu/) at Stanford University for high-throughput RNA synthesis and design.

The algorithm designs *forward* (sense strand) and *reverse* (anti-sense strand) primers that minimize the total length, and therefore the total synthesis cost, of the oligonucleotides. Although developed independently, **Primerize** is a special case of the general ‘*Gapped Oligo Design*’ algorithm, optimizing the mispriming score and sequence span instead of *T<sub>m<sub>*.

An online user-friendly GUI is available as the [**Primerize Server**](https://primerize.stanford.edu/).

## Installation

To install **Primerize**, simply:
```bash
cd path/to/Primerize/
python setup.py install
```

For system-wide installation, you must have permissions and use with `sudo`.

**Primerize** requires the following *Python* packages as dependencies, all of which can be installed through [`pip`](https://pip.pypa.io/).
```json
matplotlib >= 1.5.0
numpy >= 1.10.1
xlwt >= 1.0.0
```

#### Loop Optimization with `numba` _(Optional)_

To speed up **Primerize** code, we take advantage of [`@jit`](http://numba.pydata.org/numba-doc/0.23.1/user/jit.html) decorator of [`numba`](http://numba.pydata.org/) on loop optimization. **This is totally optional.** Enabling such feature may speed up the run for up to _10x_.

`numba` requires [`llvm`](http://llvm.org/), which can be installed through [`apt-get`](https://help.ubuntu.com/lts/serverguide/apt-get.html) on *Linux* or [`brew`](http://brew.sh/) on Mac *OSX*. It also requires `llvmlite`, which can be installed through `pip`. The compatibility between `numba`, `llvmlite`, and `llvm` needs to pay special attention to. The below specified `numba` and `llvmlite` versions have been tested to work with `llvm 3.6.2` on *Linux* machines.

```json
llvmlite == 0.8.0
numba == 0.23.1
```

Or the following works with `llvm 3.7.1`.

```json
llvmlite == 0.12.1
numba == 0.27.0
```

Newer version combinations may work, but we haven't test since.

#### Test

To test if **Primerize** is functioning properly in local installation, run the *unit test* scripts:

```bash
cd path/to/Primerize/tests/
python -m unittest discover
```

All test cases should pass.


## Usage

For simple Primer Design tasks, follow this example:

```python
import primerize

prm_1d = primerize.Primerize_1D
job_1d = prm_1d.design('TTCTAATACGACTCACTATAGGCCAAAGGCGUCGAGUAGACGCCAACAACGGAAUUGCGGGAAAGGGGUCAACAGCCGUUCAGUACCAAGUCUCAGGGGAAACUUUGAGAUGGCCUUGCAAAGGGUAUGGUAAUAAGCUGACGGACAUGGUCCUAACCACGCAGCCAAGUCCUAAGUCAACAGAUCUUCUGUUGAUAUGGAUGCAGUUCAAAACCAAACCGUCAGCGAGUAGCUGACAAAAAGAAACAACAACAACAAC', MIN_TM=60.0, NUM_PRIMERS=None, MIN_LENGTH=15, MAX_LENGTH=60, prefix='P4P6_2HP')
if job_1d.is_success:
	print job_1d

job_2d = primerize.Primerize_2D.design(job_1d, offset=-51, which_muts=range(102, 261 + 1), which_lib=1)
if job_2d.is_success:
	print job_2d
	job_2d.save()

job_3d = primerize.Primerize_3D.design(job_1d, offset=-51, structures=['...........................((((((.....))))))...........((((((...((((((.....(((.((((.(((..(((((((((....)))))))))..((.......))....)))......)))))))....))))))..)).))))((...((((...(((((((((...)))))))))..))))...)).............((((((.....))))))......................'], N_mutations=1, which_lib=1, is_single=True, is_fillWT=True)
if job_3d.is_success:
    print job_3d
    job_3d.save()
```

For advanced users, the returned `Design_Single` and `Design_Plate` result classes (refer to the [**Documentation**](https://daslab.github.io/Primerize/primerize.wrapper)) offer methods for `get()`, `save()` and `echo()`:

```python
MIN_TM = job_1d.get('MIN_TM')
print job_1d.get('MISPRIME')
print job_1d.echo('WARNING')
if job_1d.is_success:
	job_1d.save(path='result/', name='Primer')

LIB = job_2d.get('which_lib')
N_PRIMER = job_2d.get('N_PRIMER')
print job_2d.get('CONSTRUCT')
if job_2d.is_success:
	print job_2d.echo('region')
	job_2d.save('assembly', path='result/', name='Lib')

N_PLATE = job_3d.get('N_PLATE')
if job_3d.is_success:
    print job_3d.echo('plate')
    print repr(job_3d)
```

Besides `design()`, the `Primerize_1D`, `Primerize_2D`, and `Primerize_3D` factory instances offer methods for `get()`, `set()`, and `reset()`:

```python
COL_SIZE = prm_1d.get('COL_SIZE')
prm_1d.set('MIN_LENGTH', 30)
prm_1d.reset()
```

There are also `Assembly`, `Mutation`, `Construct_List`, and `Plate_96Well` helper classes. For more details, please refer to the [**Documentation**](https://daslab.github.io/Primerize/primerize.util).

Additionally, you can specify your customized list of mutations through the `Primerize_Custom` factory instance:

```python
mut_list = primerize.util.Construct_List()
mut_list.push(['T120C'])
job_cm = primerize.Primerize_Custom.design(job_1d, offset=-51, mut_list=mut_list)
```

More importantly, you can use the `merge()` method of `Construct_List` to combine multiple results into one:

```python
mut_list = job_2d.get('CONSTRUCT')
mut_list.merge(job_3d.get('CONSTRUCT'))
job_cm = primerize.Primerize_Custom.design(job_1d, offset=-51, mut_list=mut_list)
if job_cm.is_success:
    print job_2d
    job_2d.save()
```


#### MATLAB Code _(Deprecated)_

Instructions on *MATLAB* usage is available at old [README.md](https://github.com/DasLab/Primerize/blob/master/MATLAB/README.md). Please note that *MATLAB* code is no longer actively under development or fully maintained.

## Documentation

Code Documentation is available at https://daslab.github.io/Primerize/ or https://primerize.stanford.edu/docs/.

Experimental Protocol is available at https://primerize.stanford.edu/protocol/.

## License

Copyright &copy; of **Primerize** _Source Code_ is described in [LICENSE.md](https://github.com/DasLab/Primerize/blob/master/LICENSE.md).

## Reference

>Tian, S., *et al.* (**2015**)<br/>
>[Primerize: Automated Primer Assembly for Transcribing Interesting RNAs.](http://nar.oxfordjournals.org/content/43/W1/W522.full)<br/>
>*Nucleic Acid Research* **43 (W1)**: W522-W526.


>Thachuk, C., and Condon, A. (**2007**)<br/>
>[On the Design of Oligos for Gene Synthesis.](http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=4375554)<br/>
>*Proceedings of the 7th IEEE International Conference on Bioinformatics and Bioengineering* **2007**: 123-130.

<br/>
Developed by **Das lab**, _Leland Stanford Junior University_.
<br/>
README by [**t47**](http://t47.io/), *January 2016*.
