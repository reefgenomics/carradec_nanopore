import os
import subprocess
from collections import defaultdict
from hume_core_functions import *
from pathlib import Path

"""The ITS2 amplifications have been done with the LaJeunesse primers. The fastq files are located here:
/home/humebc/projects/carradec_nanopore/all_reads
We should create a sequence collection from each of the fastq.gz files
convert them to fasta files
then run a mothur analysis on them that will do a PCR using the laj primers
we should then run a blast analysis on them to get only symbiodiniaceae sequences
we should then seperate them by clade and then run MED.
If MED isn't working then I think it should be quite simple to come up with an algrythm that approximates the concepts
behind MED. I.e. find the consensus sequence, then find the next most numerate nucleotide etc etc. and compare these
to the preMED SymPortal results.
"""

def do_analysis(already_moved_files=False):
    list_of_fastq_file_paths = get_list_of_fastq_files(already_moved_files)

    process_on_sample_by_sample_basis(list_of_fastq_file_paths)


def get_list_of_fastq_files(already_moved_files):
    if not already_moved_files:
        list_of_fastq_file_paths = move_fastq_gz_files_to_new_subdirectories_and_decompress()
    else:
        list_of_fastq_file_paths = get_list_of_fastq_files_already_moved()
    return list_of_fastq_file_paths


def get_list_of_fastq_files_already_moved():
    list_of_fastq_file_paths = []
    new_directory_head = os.path.abspath((os.path.join(Path(__file__).parents[1], 'working_directory')))
    for directory_path_outer in return_list_of_directory_paths_in_directory(new_directory_head):
        for directory_path_inner in return_list_of_directory_paths_in_directory(directory_path_outer):
            temp_list_of_fastq_files = [file_path for file_path in
                                        return_list_of_file_paths_in_directory(directory_path_inner) if
                                        file_path.endswith('fastq')]
            if len(temp_list_of_fastq_files) != 1:
                raise RuntimeError('There appear to be more than 1 fastq files in the directory')
            else:
                list_of_fastq_file_paths.append(temp_list_of_fastq_files[0])
    return list_of_fastq_file_paths


def process_on_sample_by_sample_basis(list_of_fastq_file_paths):
    """This is the top function for the actual processing of the samples
    - Create sequence collections from each of the fastq
    """
    for fastq_to_process_path in list_of_fastq_file_paths:

        new_sequence_collection_from_fastq_file = SequenceCollection(
            name=fastq_to_process_path.split('/')[-1].replace('.fastq', ''), path_to_file=fastq_to_process_path
        )
        mothur_analysis = MothurAnalysis.from_sequence_collection(
            sequence_collection=new_sequence_collection_from_fastq_file,
            pcr_analysis_name='lajeunesse', pcr_fwd_primer_mismatch=6, pcr_rev_primer_mismatch=6, num_processors=20
        )

        # conduct a primer pcr
        mothur_analysis.execute_pcr(do_reverse_pcr_as_well=True)
        apples = 'asdf'
        #todo do a size screening very roughly, maybe 100bp (write a new method for this)

        #todo do a symbiodinium blast analysis (write a class for this)


def move_fastq_gz_files_to_new_subdirectories_and_decompress():
    sequence_file_directory = '/home/humebc/projects/carradec_nanopore/all_reads'
    list_of_fastq_gz_file_paths = get_list_of_fasta_gz_files(sequence_file_directory)
    wkd_head_path = setup_directory_system_to_work_within()
    list_of_fastq_file_paths = move_and_extract_each_fastq_gz_to_corresponding_subdirectory(
        wkd_head_path, list_of_fastq_gz_file_paths)
    return list_of_fastq_file_paths


def move_and_extract_each_fastq_gz_to_corresponding_subdirectory(wkd_head_path, list_of_fastq_gz_file_paths):
    list_of_fastq_file_paths = []
    for fastq_gz_path in list_of_fastq_gz_file_paths:
        target_path = move_fastq_gz_file_to_wkd(fastq_gz_path, wkd_head_path)
        extract_fastq_gz_in_new_directory(target_path)
        list_of_fastq_file_paths.append(target_path)
    return list_of_fastq_file_paths



def extract_fastq_gz_in_new_directory(target_path):
    subprocess.run(['gunzip', f'{target_path}'])


def move_fastq_gz_file_to_wkd(fastq_gz_path, wkd_head_path):
    file_name = fastq_gz_path.split('/')[-1]
    site = file_name.replace('fastq.gz', '').split('_')[0]
    sample = file_name.replace('.fastq.gz', '').split('_')[1]
    target_path = os.path.join(wkd_head_path, site, sample, file_name)
    subprocess.run(['cp', f'{fastq_gz_path}', f'{target_path}'],
                   stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return target_path

def setup_directory_system_to_work_within():
    """ There are three sites and 12 samples within each site. As such we will create a directory system that reflects
    this. Its location will be relative to the sequence_file_directory"""

    new_directory_head = create_directory_head()

    create_sub_directories_within_directory_head(new_directory_head)
    return new_directory_head


def create_directory_head():
    new_directory_head = os.path.abspath((os.path.join(Path(__file__).parents[1], 'working_directory')))
    os.makedirs(new_directory_head, exist_ok=True)
    return new_directory_head


def create_sub_directories_within_directory_head(new_directory_head):
    site_names = ('KBS1', 'KBS2', 'KBS3', 'KBS4', 'KBS4bis')
    for site_name in site_names:
        for j in range(1, 13):
            if j < 10:
                zero_str = '0'
            else:
                zero_str = ''
            os.makedirs(os.path.join(new_directory_head, f'{site_name}', f'BC{zero_str}{j}'), exist_ok=True)


def get_list_of_fasta_gz_files(sequence_file_directory):

    return [
        file_path for file_path in return_list_of_file_paths_in_directory(sequence_file_directory)
        if file_path.endswith('fastq.gz')
    ]

do_analysis(already_moved_files=True)





# def create_fasta_from_fastq_with_certain_primer_seq():
#     fastq_file = read_defined_file_to_list('/home/humebc/projects/carradec_nanopore/KBS1/BC01.fastq')
#     counter = 0
#     fasta_out = []
#
#     for i  in range(len(fastq_file)):
#         if fastq_file[i].startswith('@'):
#             get_name_from_fastq_def_line
#             if i < len(fastq_file) - 1:
#                 if 'AATGGCCTCCTGAACGTG' in fastq_file[i+1]:
#                     fasta_out.extend(['>seq_{}'.format(counter), fastq_file[i+1]])
#                     counter += 1
#
#     # write out to have a look at and find the
#     write_list_to_destination('/home/humebc/projects/carradec_nanopore/python_code/sample.fasta', fasta_out)


# create_fasta_from_fastq()

# Pull out the Symbiodinium sequences and have a look at how many sequences are returned for the number
# of mismatches allowed to the LaJeunesse sequences
# The fastqs contain coral 18S, symbiodinium ITS2 and bacterial 16S so we can start by simply
# running all of the sequences against a blast search and then only keep the symbiodinium sequences
# we still need to get rid of the primer sequences before analysis obviously but we don't really
# need to use this as a type of QC so we can have quite a high mismatch allowance I think

# do quality screening using seqtk fqchk set at value of 10
# then use seqtk to convert the resultant fastq into a fasta
# then read this in below.

# def pull_out_symbiodinium_seqs_from_fastq(required_symbiodinium_matches):
#     '''generate a fasta that only contains symbiodinium sequences. Return this as a list
#     do this by running a blast search against the nt database and checking for Symbiodinium hits'''
#
#
#     # sample_base_dir = '/home/humebc/projects/carradec_nanopore/KBS1'
#     # # Write out the hidden file that points to the ncbi database directory.
#     # ncbircFile = []
#     # # db_path = os.path.abspath(os.path.join(os.path.dirname(__file__), 'symbiodiniumDB'))
#     # db_path = '/home/humebc/phylogeneticSoftware/ncbi-blast-2.6.0+/ntdbdownload'
#     # ncbircFile.extend(["[BLAST]", "BLASTDB={}".format(db_path)])
#     # writeListToDestination("{}/.ncbirc".format(sample_base_dir), ncbircFile)
#     #
#     #
#     # # Create a fasta file from the fastq
#     # fasta_file = create_fasta_from_fastq('{}/BC01.fastq'.format(sample_base_dir))
#     # # write the fast_file
#     # path_to_input_fasta = '{}/BC01_all_seqs.fasta'.format(sample_base_dir)
#     # writeListToDestination(path_to_input_fasta, fasta_file)
#     # fasta_file_dict = createDictFromFasta(fasta_file)
#     #
#     # # Set up environment for running local blast
#     # blastOutputPath = '{}/BC01.blast.out'.format(sample_base_dir)
#     # outputFmt = "6 qseqid sseqid staxids evalue pident qcovs staxid stitle ssciname"
#     # # inputPath = r'{}/below_e_cutoff_seqs.fasta'.format(data_sub_data_dir)
#     # os.chdir('{}'.format(sample_base_dir))
#     #
#     # # Run local blast
#     # # completedProcess = subprocess.run([blastnPath, '-out', blastOutputPath, '-outfmt', outputFmt, '-query', inputPath, '-db', 'symbiodinium.fa', '-max_target_seqs', '1', '-num_threads', '1'])
#     # # completedProcess = subprocess.run(
#     # #     ['blastn', '-out', blastOutputPath, '-outfmt', outputFmt, '-query', path_to_input_fasta, '-db', 'nt',
#     # #      '-max_target_seqs', '10', '-num_threads', '20'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#     #
#     # # Read in blast output
#     # blast_output_file = readDefinedFileToList('{}/BC01.blast.out'.format(sample_base_dir))
#     #
#     # # create a dict that is the query name key and a list of subject return value
#     # blast_output_dict = defaultdict(list)
#     # for line in blast_output_file:
#     #     blast_output_dict[line.split('\t')[0]].append('\t'.join(line.split('\t')[1:]))
#     #
#     # verified_sequence_list = []
#     # for k, v in blast_output_dict.items():
#     #     sym_count = 0
#     #     for result_str in v:
#     #         if 'Symbiodinium' in result_str:
#     #             percentage_coverage = float(result_str.split('\t')[4])
#     #             if percentage_coverage > 90:
#     #                 sym_count += 1
#     #                 if sym_count == required_symbiodinium_matches:
#     #                     verified_sequence_list.append(k)
#     #                     break
#     #
#     # # We only need to proceed from here to make a new database if we have sequences that have been verified as
#     # # Symbiodinium
#     #
#     # # here we have a list of the Symbiodinium sequences that we can add to the reference db fasta
#     # new_fasta = []
#     # for seq_to_add in verified_sequence_list:
#     #     new_fasta.extend(['>{}'.format(seq_to_add), '{}'.format(fasta_file_dict[seq_to_add])])
#     #
#     #
#     # # Here we have the symbiodinium sequence fasta that we can work with
#     #
#     #
#     # path_to_input_fasta = '{}/BC01_symbiodinium.fasta'.format(sample_base_dir)
#     # writeListToDestination(path_to_input_fasta, new_fasta)
#     #
#     # # here we should use mothur to remove the primers
#     # # we will have to do this twice to make sure that we get the rev complement seqs as well
#     # # do it first, then get the remainder file and reverse complement all of the sequences in there
#     # # then run the same pcr.seqs on this rev comped file
#     # primerFwdSeq = 'GAATTGCAGAACTCCGTGAACC'  # Written 5'-->3'
#     # primerRevSeq = 'CGGGTTCWCTTGTYTGACTTCATGC'  # Written 5'-->3'
#     #
#     # oligoFile = [
#     #     r'#SYM_VAR_5.8S2',
#     #     'forward\t{0}'.format(primerFwdSeq),
#     #     r'#SYM_VAR_REV',
#     #     'reverse\t{0}'.format(primerRevSeq)
#     # ]
#     #
#     # path_to_oligo_file = '{}/oligo_file.oligo'.format(sample_base_dir)
#     # writeListToDestination(path_to_oligo_file, oligoFile)
#     #
#     # mBatchFile = [
#     #     'pcr.seqs(fasta={}, oligos={}, pdiffs=6, rdiffs=6)'.format(
#     #         path_to_input_fasta, path_to_oligo_file)
#     # ]
#     #
#     # path_to_m_batch_file = '{}/m_batch_file'.format(sample_base_dir)
#     # writeListToDestination(path_to_m_batch_file, mBatchFile)
#     # subprocess.run(['mothur', path_to_m_batch_file])
#     # for_comp_fasta_out_path = path_to_input_fasta.replace('.fasta', '.pcr.fasta')
#     # # find out what the output file will be called so that we can reverse complement and then pcr again.
#     # scrapped_seq_fasta_path = path_to_input_fasta.replace('.fasta', '.scrap.pcr.fasta')
#     # # need to remove the comments about matches etc.
#     # scrapped_fasta_clean = remove_annotations_from_fasta(scrapped_seq_fasta_path)
#     # # write out the scrapped fasta
#     # writeListToDestination(scrapped_seq_fasta_path, scrapped_fasta_clean)
#     # rev_complemented_fasta_path = scrapped_seq_fasta_path.replace('.fasta', '.rc.fasta')
#     # #reverse complement the fasta and then perorm the seqs.pcr again
#     # mBatchFile_two = [
#     #     'reverse.seqs(fasta={})'.format(scrapped_seq_fasta_path),
#     #     'pcr.seqs(fasta={}, oligos={}, pdiffs=6, rdiffs=6)'.format(rev_complemented_fasta_path, path_to_oligo_file)
#     # ]
#     #
#     # path_to_m_batch_file_two = '{}/m_batch_file_two'.format(sample_base_dir)
#     # writeListToDestination(path_to_m_batch_file_two, mBatchFile_two)
#     # subprocess.run(['mothur', path_to_m_batch_file_two])
#     # rev_comp_fasta_out_path = path_to_input_fasta.replace('.fasta', '.scrap.pcr.rc.pcr.fasta')
#     # # now add both of these fastas into one file
#     # for_comp_fasta = readDefinedFileToList(for_comp_fasta_out_path)
#     # rev_comp_fasta = readDefinedFileToList(rev_comp_fasta_out_path)
#     # primer_fasta = []
#     # primer_fasta += for_comp_fasta
#     # primer_fasta += rev_comp_fasta
#     # clean_primer_fasta = createNoSpaceFastaFile(primer_fasta)
#     # path_to_input_fasta_clade_blast = '{}/BC01_symbiodinium.pcr.rev.for.fasta'.format(sample_base_dir)
#     # blast_output_two = '{}/BC01.blast.two.out'.format(sample_base_dir)
#     # # we should clade separate these sequences
#     #
#     # writeListToDestination(path_to_input_fasta_clade_blast, clean_primer_fasta)
#     # fasta_dict = createDictFromFasta(clean_primer_fasta)
#     #
#     # ncbircFile = []
#     # db_path = '/home/humebc/phylogeneticSoftware/SymPortal_interim_250318/symbiodiniumDB'
#     # ncbircFile.extend(["[BLAST]", "BLASTDB={}".format(db_path)])
#     # outputFmt = "6 qseqid sseqid staxids evalue"
#     # writeListToDestination("{}/.ncbirc".format(sample_base_dir), ncbircFile)
#     #
#     # completedProcess = subprocess.run(
#     #     ['blastn', '-out', blast_output_two, '-outfmt', outputFmt, '-query', path_to_input_fasta_clade_blast, '-db', 'symClade_2_2.fa',
#     #      '-max_target_seqs', '1', '-num_threads', '1'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#     #
#     # clade_list = list('ABCDEFGHI')
#     # fasta_holder = [[] for clade in clade_list]
#     #
#     # # Read in blast output
#     # blast_output_file = readDefinedFileToList(blast_output_two)
#     #
#     # # Read in blast output
#     # blast_dict = {a.split('\t')[0]: a.split('\t')[1][-1] for a in blast_output_file}
#     #
#     # for k, v in blast_dict.items():
#     #     fasta_holder[clade_list.index(v)].extend(['>{}'.format(k), fasta_dict[k]])
#     #
#     # # now we just need to write out the cladal fasta
#     # for clade_fasta in fasta_holder:
#     #     if clade_fasta:
#     #         path_to_clade_fasta = '{0}/{1}/BC01_clade{1}.fasta'.format(sample_base_dir, clade_list[fasta_holder.index(clade_fasta)])
#     #         writeListToDestination(path_to_clade_fasta, clade_fasta)
#
#     # for now lets just have a look at the the single C alignment and see what MED can do with it.
#     # add dashes to the short seqs.
#     gaps_input = '/home/humebc/projects/carradec_nanopore/KBS1/C/cladeC.fasta'
#     MED_out_dir = '/home/humebc/projects/carradec_nanopore/KBS1/C'
#
#     completedProcess = subprocess.run([r'o-pad-with-gaps', r'{}'.format(gaps_input)], stdout=subprocess.PIPE,
#                                       stderr=subprocess.PIPE)
#
#     listOfFiles = []
#     for (dirpath, dirnames, filenames) in os.walk(MED_out_dir):
#         listOfFiles.extend(filenames)
#         break
#     for file in listOfFiles:
#         if 'PADDED' in file:
#             path_to_file = '{0}/{1}'.format(MED_out_dir, file)
#             break
#
#     # completedProcess = subprocess.run(
#     #     [r'decompose', '--skip-gexf-files', '--skip-gen-figures', '--skip-gen-html', '--skip-check-input', '-o',
#     #      MED_out_dir, path_to_file], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#
#     # lets unique
#
#
#     completedProcess = subprocess.run(
#         [r'decompose', '--skip-gexf-files', '--skip-gen-figures', '--skip-gen-html', '--skip-check-input', '-M', '1', '-o',
#          MED_out_dir, path_to_file])
#
#
#     # once we have the sequences clade separated we should write them out in a clade directory
#     # we should then have a look at them to see how they look.
#
#     return new_fasta




