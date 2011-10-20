#!/usr/bin/env python

# from Jeet Sukumaran's code
# https://gist.github.com/1015874

import os
import sys
import dendropy

src_filepath = sys.argv[1]
full_data = dendropy.DnaCharacterMatrix.get_from_path(src_filepath, 'nexus')
items = full_data.character_subsets.items()
basename = os.path.splitext(os.path.basename(src_filepath))[0]
for idx, (char_subset_name, char_subset) in enumerate(items):
    sys.stderr.write("%d/%d: %s\n" % (idx+1, len(items), char_subset_name))
    subset_data = full_data.export_character_subset(char_subset)
    subset_data.write_to_path(basename + '_' + char_subset_name + ".nex", 'nexus')
