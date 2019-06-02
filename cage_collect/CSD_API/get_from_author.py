#!/usr/bin/env python2
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to search for and collect CIFs using a list of authors.

Author: Andrew Tarzia

Date Created: 1 Mar 2019

"""
import ccdc.search
import sys
sys.path.insert(0, '/home/atarzia/thesource/')
import CSD_f


def write_entry(file, author, number, DOI, CSD, solvent, disorder):
    '''Write entry to CIF DB file that contains all names and references for a
    structure.

    '''
    with open(file, 'a') as f:
        f.write(author+','+number+','+DOI+','+CSD+','+solvent+','+disorder+'\n')


def write_REFCODES(file, CSD):
    '''Write REFCODE to file.

    '''
    with open(file, 'a') as f:
        f.write(CSD+'\n')


def main():
    if (not len(sys.argv) == 4):
        print """
    Usage: get_from_author.py author_file cage_type output_prefix
        author_file (str) - file with list of authors
        cage_type (str) - organic if organic cages, metal if is_organometallic
            organic: sets is_organometallic is False
            metal: sets is_organometallic is True
            anything else: passes this test
        output_prefix (str) - prefix of .txt and .gcd file to output
        """
        sys.exit()
    else:
        author_file = sys.argv[1]
        cage_type = sys.argv[2]
        output_prefix = sys.argv[3]

    out_txt = output_prefix+'.txt'
    out_gcd = output_prefix+'.gcd'

    # files = []
    authors = []
    DOIs = []
    CSD = []
    for line in open(author_file, 'r'):
        authors.append(line.rstrip())

    with open(out_txt, 'w') as f:
        f.write('author,number,DOI,CSD,solvent,disorder\n')

    with open(out_gcd, 'w') as f:
        f.write('')

    count = 0
    count_no = 0
    idents = []
    for i, author in enumerate(authors):
        # break at '-----'
        if '-----' in author:
            break
        count_no += 1
        query = ccdc.search.TextNumericSearch()
        query.add_author(author)
        hits = query.search(database='CSD')
        print author+': '+str(len(hits))
        if len(hits) == 0:
            print(author)
        for hit in hits:
            author_list = [i.strip() for i in hit.entry.publication.authors.split(',')]
            # skip polymeric structures
            if hit.entry.chemical_name is not None:
                if 'catena' in hit.entry.chemical_name:
                    continue
            if hit.entry.is_polymeric is True:
                continue
            # skip if structure is powder study
            if hit.entry.is_powder_study is True:
                continue
            if cage_type == 'organic':
                # skip structures that are NOT purely organic
                if hit.entry.is_organometallic is True:
                    continue
            elif cage_type == 'metal':
                # skip structures that are purely organic
                if hit.entry.is_organometallic is False:
                    continue
            else:
                # do not skip any
                pass
            # note structures with solvent
            solvent = 'n'
            if hit.entry.chemical_name is not None:
                if len(hit.entry.chemical_name.split(' ')) > 1:
                    solvent = 'y'
            # note structures with disorder
            disorder = 'n'
            if hit.entry.has_disorder is True:
                disorder = 'y'
            crystal = hit.crystal
            # write REFCODE to file
            if hit.identifier not in idents:
                idents.append(hit.identifier)
                write_entry(out_txt, author, str(hit.entry.ccdc_number),
                            hit.entry.doi, hit.identifier, solvent,
                            disorder)
                write_REFCODES(out_gcd, hit.identifier)
                count += 1

    print str(count)+' cifs found from '+str(count_no)+' authors'


if __name__ == "__main__":
    main()
