from setuptools import setup, find_packages

from primerize.__init__ import __version__

setup(
    name='primerize',
    description='PCR Assembly Primer Design',
    keywords='primerize PCR assembly misprime',
    version=__version__,

    author='Siqi Tian, Rhiju Das',
    author_email='rhiju@stanford.edu',

    url='https://github.com/DasLab/Primerize/',
    license='https://primerize.stanford.edu/license',

    packages=find_packages(),
    install_requires=[
        'numpy >= 1.10.1',
        'numba >= 0.22.1',
        'matplotlib >= 1.5.0',
        'xlwt >= 1.0.0'
    ],
    classifiers=(
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Developers',
        'Natural Language :: English',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7'
    )
)
