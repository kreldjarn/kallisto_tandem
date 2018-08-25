#!/usr/bin/env python3
import sys
import csv
import argparse
import pathlib

# Utility functions
# =================
def _populate_tpm(rows):
    denominators = {}
    # 1st pass to calculate the TPM denominators
    for row in rows:
        lengths = {k[11:]: v for k,v in row.items() if k.startswith('eff_length')}
        counts = {k[11:]: v for k,v in row.items() if k.startswith('est_counts')}
        for k in [*lengths]:
            try:
                denominators[k] += float(counts[k]) / float(lengths[k])
            except KeyError as e:
                denominators[k] = float(counts[k]) / float(lengths[k])

        # TODO:
        # Ask Lior whether this is the correct way to calculate the aggregate TPM
        row['aggregate'] = sum(map(float, counts.values())) / sum(map(float, lengths.values()))
        try:
            denominators['aggregate'] += row['aggregate']
        except KeyError as e:
            denominators['aggregate'] = row['aggregate']

    # 2nd pass to calculate the actual TPM
    for row in rows:
        lengths = {k[11:]: v for k,v in row.items() if k.startswith('eff_length')}
        counts = {k[11:]: v for k,v in row.items() if k.startswith('est_counts')}
        for k in [*lengths]:
            if len(k) > 0:
                key = f'tpm_{k}'
            else:
                key = 'tpm'
            row[key] = float(counts[k]) / float(lengths[k]) * 10e6 / denominators[k]
        row['tpm_aggregate'] = row['aggregate'] * 10e6 / denominators['aggregate']
        del row['aggregate']

        if len(lengths) == 1:
            del row['tpm_aggregate']

def _collapse(classes, cls, variant=None):
    # Collapses the input matrix into separate matrices for each
    # class and for each class-variant combination
    if not variant:
        matrix = {}
        for var, rows in classes[cls].items():
            if len(classes[cls].items()) > 1:
                for row in rows:
                    if row['target_id'] not in matrix:
                        matrix[row['target_id']] = {}
                    partial_row = {f'{k}_{var}': v for k, v in row.items()}
                    del partial_row[f'target_id_{var}']
                    matrix[row['target_id']] = {**matrix[row['target_id']], **partial_row}
            else:
                matrix[row['target_id']] = {**row}
                del matrix[row['target_id']]['target_id']
        collapsed = [{'target_id':key, **row} for key, row in matrix.items()]
    else:
        collapsed = [{**row} for row in classes[cls][variant]]
    return collapsed

# Public functions
# ================
def postprocess_abundance(filename, output):
    # Creates separate abundance tsv files for each class and each
    # class-variant combination. Renormalizes the TPM for for each
    # subset of the sample, respectively
    print('Parsing kallisto output into classes and variants')
    classes = {}
    with open(filename, 'r') as fh:
        tsvreader = csv.DictReader(fh, delimiter='\t')
        for row in tsvreader:
            # Parsing kallisto output into classes and variants
            target_id = row['target_id']
            v_idx = target_id.rfind('|')
            c_idx = target_id.rfind('|', 0, v_idx)
            cls = target_id[c_idx + 1: v_idx]
            var = target_id[v_idx + 1:]
            gene = target_id[:c_idx]
            obj = {
                'target_id': gene,
                'length': row['length'],
                'eff_length': row['eff_length'],
                'est_counts': row['est_counts']
            }
            if cls in classes:
                if var in classes[cls]:
                    classes[cls][var].append(obj)
                else:
                    classes[cls][var] = [obj]
            else:
                classes[cls] = {var: [obj]}

    if filename.endswith('.tsv'):
        fn = filename[:-4]
    else:
        fn = filename

    matrices = {}
    for cls, variants in classes.items():
        matrices[f'{fn}_{cls}'] = _collapse(classes, cls)
        for var in variants.keys():
            matrices[f'{fn}_{cls}_{var}'] = _collapse(classes, cls, var)

    for key, rows in matrices.items():
        print(f'Renormalizing TPM for {key}...')
        _populate_tpm(rows)

    print('Writing results to disk...')
    pathlib.Path(output).mkdir(parents=True, exist_ok=True)
    for key, matrix in matrices.items():
        with open(f'{output}/{key}', 'w') as fh:
            tsvwriter = csv.DictWriter(fh, fieldnames=matrix[0].keys(), delimiter='\t')
            tsvwriter.writeheader()
            tsvwriter.writerows(matrix)

def concatenate_fasta(filename, output):
    # Concatenates the fasta files denoted in filename. Identifies
    # the entries from each fasta file with the corresponding class
    # and variant, also denoted in filename. Writes the resulting
    # fasta file to output
    with open(filename, 'r') as fh:
        args = [{'path': l[0],
                 'cls': l[1],
                 'var': l[2]} for l in map(lambda l: l.split(), fh.read().rstrip('\n').split('\n'))]
    of = open(output, 'w')
    for f in args:
        with open(f['path'], 'r') as fh:
            for line in fh:
                if line.startswith('>'):
                    line = line.rstrip('\n')
                    of.write(f'{line}|{f["cls"]}|{f["var"]}\n')
                else:
                    of.write(line)
    of.close()

def driver():
    parser = argparse.ArgumentParser()
    parser.add_argument('action', choices=['concatenate', 'renormalize'], help='Concatenates one or more fasta files to be used to create a joint kallisto index. Renormalize splits up the data from a tandem kallisto run into its constituent classes and variants, renormalizes TPM w.r.t. each submatrix.')
    parser.add_argument('-i', help='Path to the input file to be processed.', metavar='input_file', required=True)
    parser.add_argument('-o', help='Path to the output file to which results are to be written.', metavar='output_directory', required=True)
    parser.add_argument('--version', action='version', version='%(prog)s v0.1')
    args = parser.parse_args()

    {
        'renormalize': postprocess_abundance,
        'concatenate': concatenate_fasta
    }.get(args.action)(filename=args.i, output=args.o)

if __name__ == '__main__':
    driver()
