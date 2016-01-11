from setuptools import setup


setup(
    name='primerize',
    description='PCR Assembly Primer Design',
    keywords='primerize PCR assembly misprime',
    version='1.0.0',

    author='Siqi Tian, Rhiju Das',
    author_email='rhiju@stanford.edu',

    url='https://github.com/DasLab/Primerize/',
    license='https://primerize.stanford.edu/license',

    packages=['primerize'],
    install_requires=[
        'numpy >= 1.10.1',
        'numba >= 0.22.1',
        'matplotlib >= 1.5.0',
        'xlwt >= 1.0.0'
    ],
    classifiers=(
        'Intended Audience :: Developers',
        'Natural Language :: English',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7'
    )
)