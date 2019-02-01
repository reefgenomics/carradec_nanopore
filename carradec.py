import os
import subprocess
from collections import defaultdict
from hume_core_functions import *
from pathlib import Path
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt

""" So.. The idea of this analysis was to test out my own version of alignment decomposition.
    Sadly, I'm going to have to abandon it at this point. I've spent a lot of time on it but I think
    that at the end of the day, the sequencing is just too innacurate to get anything close to what the illumina
    data gives us. Instead, what I'll do is simply annotate the sequences that he currently has according to the
    same naming system that the SymPortal remote repository uses and then ship this back to him.
    
The ITS2 amplifications have been done with the LaJeunesse primers. The fastq files are located here:
/home/humebc/projects/carradec_nanopore/all_reads

The dcomposition was going to work as such: You do an alignment of the sequences you have,
then you work your way column by column through that alignment calculating each states score. The state is either
gap, A, C, G or T. The score is +1 for every occurence of that score in the column you are in. You work all the way
through the alignment doing this. Then, you set a threshold for the score that is required for a given state to be
recognised as non-error and to be kept. 
Then you work through the original sequences, one by one. For each sequences, go through each of the nucleotides.
At each position in the sequences, you look at what nucleotide is there (gap, A, C, T, or G) and you look up the 
corresponding score of the nucleotide from when you did the state scoring. This score will be specific for that 
nucleotide in that column. If you are above the threshold that has been set, you keep the nucleotide in that 
position. if not, you revert it to the highest scoring state in that column. Simples :). I was considering doing
this through several iterations, ie. set the threshold quite low, do the decomp, then realign, then do the decomp
again. Problem was that there was just soooo much error in the nanopore data and that each ITS2 seq variant is only
a single nucleotide apart that it really wasn't possible to get this to show good agreement with the Illumina
sequencing results. Sad times.
"""

# ### Code for doing the annotations according to the SymPortal db.
def annotate_quentin_samples():
    """ Quentin has a tab delim file that has each of the sequences that were called in each of the samples
    We will read this in and parse through the symportal master fasta to find matches. We will look for exact matches
    If we don't have an exact match then we can work with the seq name he already has.
    """
    path_to_tab_delim = '/home/humebc/projects/carradec_nanopore/ITS2_consensus.txt'
    path_to_symportal_master_fasta = '/home/humebc/projects/carradec_nanopore/master_fasta.fa'

    sp_dict = create_dict_from_fasta(read_defined_file_to_list(path_to_symportal_master_fasta))

    # make a pandas dataframe from the tab delim file
    cols = ['site', 'sample', 'q_type', 'count', 'samp_count', 'sequence']
    quentin_df = pd.read_csv(sep='\t', filepath_or_buffer=path_to_tab_delim, names=cols)
    quentin_df.index = range(len(quentin_df.index.values.tolist()))

    # get a unique set of sequences. We will use these to make a dictionary of what they should be called with the
    # full sequence as the key and the name of the sequence as the value.
    # this way we are not redoing searches for the same sequences just becuase they are in different samples
    unique_sequence_list = []
    for index_value in quentin_df.index.values.tolist():
        seq_in_question = quentin_df.iloc[index_value]['sequence']
        if seq_in_question not in unique_sequence_list:
            unique_sequence_list.append(seq_in_question)

    # here we have a list of the sequences we will need to run through.

    seq_to_check_count = 0
    match_dict = {}
    for seq_to_check in unique_sequence_list:
        print(f'Checking seq {seq_to_check_count} out of {len(unique_sequence_list)}')
        found = False
        for k, v in sp_dict.items():
            if seq_to_check in v or v in seq_to_check:
                # then we have a match and we should use k as the name for this sequences
                match_dict[seq_to_check] = k
                found = True
                break
        if not found:
            match_dict[seq_to_check] = 'no_match'
        seq_to_check_count += 1


    # here we have checked all of the sequences and we now have a dict
    # we can now go back through the original quentin dict and add a column that can be the sp_name
    quentin_df['sp_name'] = 'no_value'
    # reorder so that the sp_name is next to the q_type
    quentin_df = quentin_df[['site', 'sample', 'q_type', 'sp_name', 'count', 'samp_count', 'sequence']]
    for index_value in quentin_df.index.values.tolist():
        quentin_df.at[index_value, 'sp_name'] = match_dict[quentin_df.iloc[index_value]['sequence']]

    # here we have the new row of the df popualted and it is now time to
    # write it out as tab delim again and ship it back
    quentin_df.to_csv(
        path_or_buf='/home/humebc/projects/carradec_nanopore/ITS2_consensus_sp.txt', sep='\t', index=False
    )

    # baddabing.

annotate_quentin_samples()
# #### Code for trying to do the decomposition of the samples
def do_analysis_decomposition(already_moved_files=False):
    """The idea of this analysis was to test out my own version of alignment decomposition.
    Sadly, I'm going to have to abandon it at this point. I've spent a lot of time on it but I think
    that at the end of the day, the sequencing is just too innacurate to get anything close to what the illumina
    data gives us. Instead, what I'll do is simply annotate the sequences that he currently has according to the
    same naming system that the SymPortal remote repository uses."""
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

        # conduct a primer pcr with the lajeunesse primer set using mismatch at 6 fwd and rev
        mothur_analysis.execute_pcr(do_reverse_pcr_as_well=True)

        # run a rough size screening to get rid of really small reads
        mothur_analysis.execute_screen_seqs({'maxlength' : '1000', 'minlength' : '100', 'processors' : '20'})

        # todo do a symbiodinium blast analysis (write a class for this)
        blastn_analysis = BlastnAnalysis(
            input_file_path=mothur_analysis.fasta_path,
            output_file_path=os.path.join(os.path.dirname(mothur_analysis.fasta_path), 'blast.out'),
            db_path='/home/humebc/phylogeneticSoftware/SymPortal_framework/symbiodiniumDB/symClade.fa',
            db_name='symClade.fa', max_target_seqs=1, num_threads=20,
            output_format_string="6 qseqid sseqid staxids evalue pident qcovs staxid stitle ssciname",
            blastn_exec_path='blastn')

        blastn_analysis.execute_blastn()
        blast_output_file = read_defined_file_to_list(blastn_analysis.output_file_path)

        fasta_dictionary = create_dict_from_fasta(read_defined_file_to_list(mothur_analysis.fasta_path))

        clade_list = list('ABCDEFGHI')
        lists_of_clade_separated_seqs = [[] for _ in clade_list]

        name_of_sequences_already_added = []
        for blast_output_line in blast_output_file:
            # get the clade of the sequence in question
            clade_of_seq = blast_output_line.split('\t')[1][-1]
            # put the sequence in the respective list
            clade_index = clade_list.index(clade_of_seq)
            sequence_in_question_name = blast_output_line.split('\t')[0]
            # first check that we haven't already added the sequence
            if sequence_in_question_name not in name_of_sequences_already_added:
                name_of_sequences_already_added.append(sequence_in_question_name)
                lists_of_clade_separated_seqs[clade_index].extend([f'>{sequence_in_question_name}', f'{fasta_dictionary[sequence_in_question_name]}'])
            else:
                this = 'asdf'
        # here we have lists of fastas that are separated by clade

        # then work clade by clade and align the sequences
        for index_of_clade_separated_fasta, clade_separated_fasta_list in enumerate(lists_of_clade_separated_seqs):
            clade_in_question = clade_list[index_of_clade_separated_fasta]
            # only do it if we have more than 10 sequencs otherwise don't bother
            if len(clade_separated_fasta_list) > 20:
                # write out the fasta
                post_blast_fasta_path = os.path.join(os.path.dirname(mothur_analysis.fasta_path), f'symbiodinium_fasta_to_align_clade_{clade_in_question}.fa')
                write_list_to_destination(post_blast_fasta_path, clade_separated_fasta_list)


                # DEBUG check this fasta for name duplicates to see where it is happening
                list_of_names = []
                dup_list = []
                dup_count = 0
                for i in range(0, len(clade_separated_fasta_list), 2):
                    if clade_separated_fasta_list[i][1:] not in list_of_names:
                        list_of_names.append(clade_separated_fasta_list[i][1:])
                    else:
                        this = 'is'
                        dup_list.append(clade_separated_fasta_list[i][1:])
                        dup_count += 1

                output_path_for_aligned_fasta = os.path.join(os.path.dirname(mothur_analysis.fasta_path), f'symbiodinium_fasta_aligned_clade_{clade_in_question}.fa')

                # align the fasta
                # this takes some time to do well so let's check to see if the file already exists
                # if it does, skip doing the alignment again and simply move on to doing the decomposition
                if not os.path.exists(output_path_for_aligned_fasta):
                    mafft_align_fasta(
                        input_path=post_blast_fasta_path, output_path=output_path_for_aligned_fasta,
                        num_proc=20, method='linsi', iterations=1000)


                aligned_fasta_for_decomposition_interleaved = read_defined_file_to_list(output_path_for_aligned_fasta)
                aligned_fasta_for_decomposition = convert_interleaved_to_sequencial_fasta(aligned_fasta_for_decomposition_interleaved)
                fasta_to_decompose_as_df = fasta_to_pandas_df(aligned_fasta_for_decomposition)
                # this df is what we should pass into the method that will do the decomposition
                # TODO make the below case insensitive just incase we want to decompose non lowercase alignments
                # create a dataframe that will hold the results to this
                state_score_df = pd.DataFrame(index=list('-acgt'), columns=list(fasta_to_decompose_as_df))
                for i in list(fasta_to_decompose_as_df):
                    column_as_series = fasta_to_decompose_as_df[i]
                    count_of_unique_values = column_as_series.value_counts()
                    for index in state_score_df.index.values.tolist():
                        if index in count_of_unique_values.index:
                            state_score_df.loc[index, i] = count_of_unique_values[index]
                        else:
                            state_score_df.loc[index, i] = 0

                ## this is simply
                # # here we have all of the states processed for the alignment
                # apples = 'asdf'
                # # now we can go about trying to plot them. I guess just as a scatter on a single point or maybe as the function of the column position
                # ax1 = plt.subplot()
                # x_data = []
                # y_data = []
                # for i in list(state_score_df):
                #     x_data.extend([i for _ in range(len(state_score_df.index))])
                #     y_data.extend(state_score_df[i].values.tolist())
                #
                # ax1.hist(y_data, bins=400, color='blue', edgecolor='black')
                # ax1.set_ylim(0, 200)
                # plt.show()
                # ax1.scatter(x=x_data, y=y_data)

                apples = 'asdf'
                # the alignment path will be located at the given output path
                # we should have a look at this.
                # Ok so we can see that there is indeed a break in the data and we can work with a value of about 300
                # so now we can work through each of the sequeneces (we will skip doing the uniqueing
                # because we are here in the first place due to the high error rate there are likely to be any unique
                # sequenecs and so this will save us very little time.
                # we will now go sequence by sequence through the alignment df and for each of the nucleotides
                # we can simply do a look up in the state df and unless that is higher than the 300 value then
                # we will convert this back to the consensus value of this value (which will be the highest value in
                # the state matrix for that (column)

                # DEBUG there seem to be two of the same entry in the df.
                # lets check if that is here
                name_list = []
                dup_count_list = []
                dup_count = 0
                for seq_name in fasta_to_decompose_as_df.index.values.tolist():
                    if seq_name not in name_list:
                        name_list.append(seq_name)
                    else:
                        dup_count +=1
                        dup_count_list.append(seq_name)
                apples = 'asdf'



                # for each sequence in the df
                seq_number = len(fasta_to_decompose_as_df.index.values.tolist())
                seq_count = 0
                # TODO there should be a much faster way to do this.
                for sequence_name_index in fasta_to_decompose_as_df.index.values.tolist():
                    seq_count += 1
                    print(f'sequence {seq_count} our of {seq_number}')
                    seq_to_check_as_series = fasta_to_decompose_as_df.loc[sequence_name_index]
                    # for each nucleotide in the sequence
                    for i in list(seq_to_check_as_series.index.values):
                        # here we can look at a given nuclotide and look it up in the state table and see if it is above the threshold
                        current_nucleotide_at_pos_in_seq = seq_to_check_as_series[i]
                        current_state_score = state_score_df.loc[current_nucleotide_at_pos_in_seq, i]
                        if current_state_score < (seq_number*.025): #TODO make this value to be 2.5% of the total number of seqs
                            # then this will need changing to the highest scoring state of the column
                            #idxmax
                            high_scoring_nucleotide_for_column_in_question = state_score_df[i].idxmax()
                            fasta_to_decompose_as_df.loc[sequence_name_index, i] = high_scoring_nucleotide_for_column_in_question




                # for the fasta that we are working on it should now have been fully decomposed with the 300 cutoff.
                # now we should write out the df as a fasta file and take a look at it.
                decomposed_fasta_as_list = pandas_df_to_fasta(fasta_to_decompose_as_df)

                output_path_for_decomposed_fasta = os.path.join(os.path.dirname(output_path_for_aligned_fasta), f'decomposed_fasta_clade_{clade_in_question}_with_gaps.fa')
                print(output_path_for_decomposed_fasta)
                write_list_to_destination(output_path_for_decomposed_fasta, decomposed_fasta_as_list)

                # this will now need to be uniqued and then plotted up to have a look at it.
        # here we have a set of sequences that have found matches to some degree with the symClade.fa dictionary
        # we should now attempt to align them, within clade
        # crop them
        # and then run the algorythm thingy on them.
        apples = 'asdf'

def pause_point_decomposed_fasta():
    decomposed_fasta_path = '/home/humebc/projects/carradec_nanopore/working_directory/KBS4/BC12/decomposed_fasta_clade_C_with_gaps.fa'
    decomposed_fasta_without_gaps_as_list = remove_gaps_from_fasta(read_defined_file_to_list(decomposed_fasta_path))
    decomposed_fasta_path_without_gaps = '/home/humebc/projects/carradec_nanopore/working_directory/KBS4/BC12/decomposed_fasta_clade_C_without_gaps.fa'
    write_list_to_destination(decomposed_fasta_path_without_gaps, decomposed_fasta_without_gaps_as_list)
    apples = 'asdf'
    # here we should unique and then plot up
    # update the mothur analysis object and write the code for uniqueing
    # for time being just do it on the command line
    path_of_uniqued_decomposed_fasta = '/home/humebc/projects/carradec_nanopore/working_directory/KBS4/BC12/decomposed_fasta_clade_C_without_gaps.unique.fa'
    path_of_uniqued_decomposed_name = '/home/humebc/projects/carradec_nanopore/working_directory/KBS4/BC12/decomposed_fasta_clade_C_without_gaps.names'


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

# do_analysis(already_moved_files=True)
# pause_point_decomposed_fasta()
