--------------------
Import and Export
--------------------

Details of modes solved can be exported as CSV file. Assuming we have a list
of modes stored in ``modes``, they can be saved to a file ``filename = example.csv``
using::

    lpmodes.modes_to_csv(modes, filename)
    
The modes can also be reloaded by::

    modes = lpmodes.csv_to_modes(filename)
    
By default, a header line is included in the CSV, pass ``header = False`` as
a keyword argument to not include a header (and also when loading in from a file
with no header).

       

