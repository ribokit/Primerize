Example Scripts
------------------

For simple Primer Design tasks, follow this example:

.. code-block:: python

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

For advanced users, the returned ``Design_Single`` and ``Design_Plate`` result classes (please refer to the `Documentation <../primerize.wrapper>`_) offer methods for ``get()``, ``save()`` and ``echo()``:

.. code-block:: python

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

Besides ``design()``, the ``Primerize_1D``, ``Primerize_2D``, and ``Primerize_3D`` factory instances offer methods for ``get()``, ``set()``, and ``reset()``:

.. code-block:: python

    COL_SIZE = prm_1d.get('COL_SIZE')
    prm_1d.set('MIN_LENGTH', 30)
    prm_1d.reset()

There are also ``Assembly``, ``Mutation``, ``Construct_List``, and ``Plate_96Well`` helper classes. For more details, please refer to the `Documentation <../primerize.util>`_.


Additionally, you can specify your customized list of mutations through the ``Primerize_Custom`` factory instance:

.. code-block:: python

    mut_list = primerize.util.Construct_List()
    mut_list.push(['T120C'])
    job_cm = primerize.Primerize_Custom.design(job_1d, offset=-51, mut_list=mut_list)

More importantly, you can use the ``merge()`` method of ``Construct_List`` to combine multiple results into one:

.. code-block:: python

    mut_list = job_2d.get('CONSTRUCT')
    mut_list.merge(job_3d.get('CONSTRUCT'))
    job_cm = primerize.Primerize_Custom.design(job_1d, offset=-51, mut_list=mut_list)
    if job_cm.is_success:
        print job_2d
        job_2d.save()
