import os
import subprocess
import sys

class MothurAnalysis:

    def __init__(
            self, sequence_collection=None,  input_dir=None, output_dir=None, name=None,
            fastq_gz_fwd_path=None, fastq_gz_rev_path=None,  fasta_path=None,
             name_file_path=None, mothur_execution_path='mothur', auto_convert_fastq_to_fasta=True,
            pcr_fwd_primer=None, pcr_rev_primer=None, pcr_oligo_file_path=None,
            pcr_fwd_primer_mismatch=2, pcr_rev_primer_mismatch=2, pcr_analysis_name=None, num_processors=10
            ):

        self.setup_core_attributes(auto_convert_fastq_to_fasta, fasta_path, fastq_gz_fwd_path, fastq_gz_rev_path,
                                   input_dir, mothur_execution_path, name, name_file_path, output_dir,
                                   sequence_collection, num_processors)


        self.setup_pcr_analysis_attributes(pcr_analysis_name, pcr_fwd_primer, pcr_fwd_primer_mismatch,
                                           pcr_oligo_file_path, pcr_rev_primer, pcr_rev_primer_mismatch)

    def setup_core_attributes(self, auto_convert_fastq_to_fasta, fasta_path, fastq_gz_fwd_path, fastq_gz_rev_path,
                              input_dir, mothur_execution_path, name, name_file_path, output_dir,
                              sequence_collection, num_processors):

        self.verify_that_is_either_sequence_collection_or_fastq_pair(fastq_gz_fwd_path, fastq_gz_rev_path,
                                                                 sequence_collection)

        if sequence_collection is not None:
            self.setup_sequence_collection_attribute(auto_convert_fastq_to_fasta, name, sequence_collection)
        elif sequence_collection is None:
            self.setup_fastq_attributes(fastq_gz_fwd_path, fastq_gz_rev_path)

        self.setup_remainder_of_core_attributes(fasta_path, input_dir, mothur_execution_path, name_file_path,
                                                output_dir, sequence_collection, num_processors)

    def setup_remainder_of_core_attributes(self, fasta_path, input_dir, mothur_execution_path, name_file_path,
                                           output_dir, sequence_collection, num_processors):
        self.exec_path = mothur_execution_path
        if input_dir is None:
            self.input_dir = os.path.dirname(sequence_collection.file_path)
        else:
            self.input_dir = input_dir
        if output_dir is None:
            self.output_dir = os.path.dirname(sequence_collection.file_path)
        else:
            self.output_dir = input_dir

        self.name_file_path = name_file_path
        self.mothur_batch_file_path = None
        self.processors = num_processors

    def setup_fastq_attributes(self, fastq_gz_fwd_path, fastq_gz_rev_path):
        self.fastq_gz_fwd_path = fastq_gz_fwd_path
        self.fastq_gz_rev_path = fastq_gz_rev_path
        self.sequence_collection = None
        self.fasta_path = None

    def setup_sequence_collection_attribute(self, auto_convert_fastq_to_fasta, name, sequence_collection):
        self.fastq_gz_fwd_path = None
        self.fastq_gz_rev_path = None
        if sequence_collection.file_type == 'fastq':
            self.convert_to_fasta_or_raise_value_error(auto_convert_fastq_to_fasta, sequence_collection)
        if name is None:
            self.name = sequence_collection.name
        else:
            self.name = name
        self.sequence_collection = sequence_collection
        self.fasta_path = self.sequence_collection.file_path

    def convert_to_fasta_or_raise_value_error(self, auto_convert_fastq_to_fasta, sequence_collection):
        if auto_convert_fastq_to_fasta:
            print('SequenceCollection must be of type fasta\n. Running SeqeunceCollection.convert_to_fasta.\n')
            sequence_collection.convert_to_fasta()
        else:
            ValueError('SequenceCollection must be of type fasta. You can use the SequenceCollection')

    def verify_that_is_either_sequence_collection_or_fastq_pair(self, fastq_gz_fwd_path, fastq_gz_rev_path,
                                                            sequence_collection):
        if sequence_collection and (fastq_gz_fwd_path or fastq_gz_rev_path):
            raise ValueError(
                'Please create a MothurAnalysis from either a sequence_collection OR a pair of fastq_gz files.\n'
                'MothurAnalysis.from_pair_of_fastq_gz_files or MothurAnalysis.from_sequence_collection')

    def setup_pcr_analysis_attributes(self, pcr_analysis_name, pcr_fwd_primer, pcr_fwd_primer_mismatch,
                                      pcr_oligo_file_path, pcr_rev_primer, pcr_rev_primer_mismatch):
        if pcr_analysis_name:
            if pcr_analysis_name.lower() in ['symvar', 'sym_var']:
                self.pcr_fwd_primer = 'GAATTGCAGAACTCCGTGAACC'
                self.rev_primer = 'GAATTGCAGAACTCCGTGAACC',

            elif pcr_analysis_name.lower() in ['laj', 'lajeunesse']:
                self.pcr_fwd_primer = 'GAATTGCAGAACTCCGTG'
                self.pcr_rev_primer = 'CGGGTTCWCTTGTYTGACTTCATGC'
            else:
                raise ValueError(
                    'pcr_analysis_name \'{}\' is not recognised.\nOptions are \'symvar\' or \'lajeunesse\'.'
                )
        else:
            self.pcr_fwd_primer = pcr_fwd_primer
            self.pcr_rev_primer = pcr_rev_primer
        self.pcr_fwd_primer_mismatch = pcr_fwd_primer_mismatch
        self.pcr_rev_primer_mismatch = pcr_rev_primer_mismatch
        self.pcr_oligo_file_path = pcr_oligo_file_path

    @classmethod
    def from_pair_of_fastq_gz_files(cls, name, fastq_gz_fwd_path, fastq_gz_rev_path,
                                    output_dir=None, mothur_execution_path='mothur', num_processors=10):
        return cls(name=name, sequence_collection=None, mothur_execution_path=mothur_execution_path,
                   input_dir=os.path.dirname(os.path.abspath(fastq_gz_fwd_path)), output_dir=output_dir,
                   fastq_gz_fwd_path=fastq_gz_fwd_path, fastq_gz_rev_path=fastq_gz_rev_path,
                   fasta_path=None, name_file_path=None, num_processors=num_processors)

    @classmethod
    def from_sequence_collection(cls, sequence_collection, name=None, input_dir=None,
                                 output_dir=None, mothur_execution_path='mothur',
                                 pcr_fwd_primer=None, pcr_rev_primer=None, pcr_oligo_file_path=None,
                                 pcr_fwd_primer_mismatch=2, pcr_rev_primer_mismatch=2, pcr_analysis_name=None, num_processors=10):
        return cls(
            name=name, sequence_collection=sequence_collection, input_dir=input_dir,
            output_dir=output_dir, mothur_execution_path=mothur_execution_path, pcr_fwd_primer=pcr_fwd_primer,
            pcr_rev_primer=pcr_rev_primer, pcr_oligo_file_path=pcr_oligo_file_path, pcr_fwd_primer_mismatch=pcr_fwd_primer_mismatch,
            pcr_rev_primer_mismatch=pcr_rev_primer_mismatch, pcr_analysis_name=pcr_analysis_name, num_processors=num_processors
        )


    def execute_pcr(self, do_reverse_pcr_as_well=False):
        """This will perform a mothur pcr.seqs analysis.
        if do_reverse_pcr__as_well is true then we will also reverse complement the fasta a perform the
        """

        self.pcr_validate_attributes_are_set()

        self.pcr_make_and_write_oligo_file_if_doesnt_exist()

        self.pcr_make_and_write_mothur_batch_file()

        completed_process = self.run_mothur_batch_file()

        fwd_output_scrapped_fasta_path, fwd_output_good_fasta_path = self.pcr_extract_good_and_scrap_output_paths(
            completed_process
        )

        remove_primer_mismatch_annotations_from_fasta(fwd_output_scrapped_fasta_path)
        remove_primer_mismatch_annotations_from_fasta(fwd_output_good_fasta_path)


        # then we should clean up the output_bad_fasta
        # then reverse complement it
        # then do a pcr on it again using the same oligo set as the first run
        # we should then get the output from that pcr and add it to the previous run
        if do_reverse_pcr_as_well:
            self.fasta_path = fwd_output_scrapped_fasta_path
            self.rev_comp_make_and_write_mothur_batch_file()
            completed_process = self.run_mothur_batch_file()
            self.fasta_path = self.rev_comp_extract_new_fasta_path(completed_process)
            self.pcr_make_and_write_mothur_batch_file()
            completed_process = self.run_mothur_batch_file()
            rev_output_good_fasta_path = self.pcr_extract_good_and_scrap_output_paths(completed_process)[1]
            remove_primer_mismatch_annotations_from_fasta(rev_output_good_fasta_path)
            self.make_new_fasta_path_for_fwd_rev_combined(rev_output_good_fasta_path)
            # now create a fasta that is the good fasta from both of the pcrs. this will become the new mothuranalysis fasta.

            combine_two_fasta_files(
                path_one=fwd_output_good_fasta_path,
                path_two=rev_output_good_fasta_path,
                path_for_combined=self.fasta_path
            )
        else:
            self.fasta_path = fwd_output_good_fasta_path

        self.update_sequence_collection_from_fasta_file()

    def make_new_fasta_path_for_fwd_rev_combined(self, rev_output_good_fasta_path):
        self.fasta_path = rev_output_good_fasta_path.replace('.scrap.pcr.rc.pcr', '.pcr.combined')

    def update_sequence_collection_from_fasta_file(self):
        self.sequence_collection.generate_sequence_collection(self.fasta_path)

    def rev_comp_make_and_write_mothur_batch_file(self):
        mothur_batch_file = self.make_rev_complement_mothur_batch_file()
        self.mothur_batch_file_path = os.path.join(self.input_dir, 'mothur_batch_file')
        write_list_to_destination(self.mothur_batch_file_path, mothur_batch_file)

    def make_rev_complement_mothur_batch_file(self):
        mothur_batch_file = [
            f'set.dir(input={self.input_dir})',
            f'set.dir(output={self.output_dir})',
            f'reverse.seqs(fasta={self.fasta_path})'
        ]
        return mothur_batch_file

    def rev_comp_extract_new_fasta_path(self, completed_process):
        stdout_string_as_list = completed_process.stdout.decode('utf-8').split('\n')
        for i in range(len(stdout_string_as_list)):
            print(stdout_string_as_list[i])
            if 'Output File Names' in stdout_string_as_list[i]:
                return stdout_string_as_list[i+1]

    def pcr_extract_good_and_scrap_output_paths(self, completed_process):
        stdout_string_as_list = completed_process.stdout.decode('utf-8').split('\n')
        for i in range(len(stdout_string_as_list)):
            print(stdout_string_as_list[i])
            if 'Output File Names' in stdout_string_as_list[i]:
                output_good_fasta_path = stdout_string_as_list[i + 1]
                output_scrapped_fasta_path = stdout_string_as_list[i + 3]
                return output_scrapped_fasta_path, output_good_fasta_path


    def run_mothur_batch_file(self):
        completed_process = subprocess.run(
            [self.exec_path, self.mothur_batch_file_path],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        return completed_process

    def pcr_make_and_write_mothur_batch_file(self):
        mothur_batch_file = self.pcr_make_mothur_batch_file()
        self.mothur_batch_file_path = os.path.join(self.input_dir, 'mothur_batch_file')
        write_list_to_destination(self.mothur_batch_file_path, mothur_batch_file)

    def pcr_make_mothur_batch_file(self):
        if self.name_file_path:
            mothur_batch_file = [
                f'set.dir(input={self.input_dir})',
                f'set.dir(output={self.output_dir})',
                f'pcr.seqs(fasta={self.fasta_path}, name={self.name_file_path}, oligos={self.pcr_oligo_file_path}, '
                f'pdiffs={self.pcr_fwd_primer_mismatch}, rdiffs={self.pcr_rev_primer_mismatch}, processors={self.processors})'
            ]

        else:
            mothur_batch_file = [
                f'set.dir(input={self.input_dir})',
                f'set.dir(output={self.output_dir})',
                f'pcr.seqs(fasta={self.fasta_path}, oligos={self.pcr_oligo_file_path}, '
                f'pdiffs={self.pcr_fwd_primer_mismatch}, rdiffs={self.pcr_rev_primer_mismatch}, processors={self.processors})'
            ]
        return mothur_batch_file

    def pcr_make_and_write_oligo_file_if_doesnt_exist(self):
        if self.pcr_oligo_file_path is None:
            oligo_file = [
                f'forward\t{self.pcr_fwd_primer}',
                f'reverse\t{self.pcr_rev_primer}'
            ]
            self.pcr_oligo_file_path = os.path.join(self.input_dir, 'oligo_file.oligo')
            write_list_to_destination(self.pcr_oligo_file_path, oligo_file)

    def pcr_validate_attributes_are_set(self):
        sys.stdout.write(f'\nValidating PCR attributes are set\n')
        if self.fasta_path is None:
            raise ValueError('fasta_path is None. A valid fasta_path is required to perform the pcr method.')
        if self.pcr_fwd_primer is None or self.pcr_rev_primer is None:
            if self.pcr_fwd_primer is None and self.pcr_rev_primer is None:
                raise ValueError('Please set fwd_primer and rev_primer: ')
            elif self.pcr_fwd_primer is None:
                raise ValueError('Please set fwd_primer.')
            elif self.pcr_rev_primer is None:
                raise ValueError('Please set fwd_primer.')
        sys.stdout.write(f'\nPCR attributes: OK\n')


class SequenceCollection:
    """ A sequence collection is a set of sequences either generated from a fastq file or from a fasta file.
    It cannot be created directly from binary files or from paired files. As such, to generate a SequenceCollection
    for example from a pair of fastaq.gz files, you would first have to run a mothur contig analysis and create
    the SeqeunceCollection from the resultant fasta file that is generated."""
    def __init__(self, name, path_to_file=None, auto_convert_to_fasta=True):
        self.name = name
        self.file_path = path_to_file
        # self.file_as_list = read_defined_file_to_list(self.file_path)
        self.file_type = self.infer_file_type()
        self.sequence_collection = self.generate_sequence_collection()
        if auto_convert_to_fasta:
            self.convert_to_fasta()


    def convert_to_fasta(self):
        self.file_path = self.write_out_as_fasta()
        self.file_type = 'fasta'

    def __len__(self):
        return(len(self.sequence_collection))

    def write_out_as_fasta(self, path_for_fasta_file = None):
        if self.file_type == 'fasta':
            print(f'SequenceCollection is already of type fasta and a fasta file already exists: {self.file_path}')
            return
        if self.file_type == 'fastq':
            if path_for_fasta_file is None:
                fasta_path = self.infer_fasta_path_from_current_fastq_path()
                write_list_to_destination(destination=fasta_path, list_to_write=self.as_fasta())
            else:
                fasta_path = path_for_fasta_file
                write_list_to_destination(destination=fasta_path, list_to_write=self.as_fasta())
            return fasta_path


    def infer_fasta_path_from_current_fastq_path(self):
        return self.file_path.replace('fastq', 'fasta')

    def generate_sequence_collection(self, alt_fasta_path=None):
        if self.file_type == 'fasta':
            return self.parse_fasta_file_and_extract_nucleotide_sequence_objects(alternative_fasta_file_path=alt_fasta_path)
        elif self.file_type == 'fastq':
            return self.parse_fastq_file_and_extract_nucleotide_sequence_objects()

    def parse_fasta_file_and_extract_nucleotide_sequence_objects(self, alternative_fasta_file_path=None):
        list_of_nuleotide_sequence_objects = []
        if alternative_fasta_file_path:
            self.file_path = alternative_fasta_file_path
        fasta_file = read_defined_file_to_list(self.file_path)
        for i in range(0, len(fasta_file), 2):
            list_of_nuleotide_sequence_objects.append(
                NucleotideSequence(sequence=fasta_file[i+1], name=fasta_file[i][1:])
            )
        return list_of_nuleotide_sequence_objects

    def parse_fastq_file_and_extract_nucleotide_sequence_objects(self):
        list_of_nuleotide_sequence_objects = []
        fastq_file_as_list = read_defined_file_to_list(self.file_path)
        for i in range(len(fastq_file_as_list)):
            if i < len(fastq_file_as_list) - 2:
                if self.is_fastq_defline(fastq_file_as_list, i):
                    self.create_new_nuc_seq_object_and_add_to_list(fastq_file_as_list, i, list_of_nuleotide_sequence_objects)
        return list_of_nuleotide_sequence_objects

    def is_fastq_defline(self, fastsq_file, index_value):
        if fastsq_file[index_value].startswith('@') and fastsq_file[index_value + 2][0] == '+':
            return True

    def create_new_nuc_seq_object_and_add_to_list(self, fastq_file_as_list, index_val, list_of_nuleotide_sequence_objects):
        name, sequence = self.get_single_fastq_info_from_fastq_file_by_index(fastq_file_as_list, index_val)
        list_of_nuleotide_sequence_objects.append(NucleotideSequence(sequence=sequence, name=name))

    def get_single_fastq_info_from_fastq_file_by_index(self, fastq_file_as_list, index_val):
        name = fastq_file_as_list[index_val][1:].split(' ')[0]
        sequence = fastq_file_as_list[index_val + 1]
        return name, sequence

    def infer_file_type(self):
        if 'fasta' in self.file_path:
            return 'fasta'
        elif 'fastq' in self.file_path:
            return 'fastq'
        else:
            raise ValueError('Input file used to create the SequenceCollection must be either fasta or fastq')

    def as_fasta(self):
        fasta_file = []
        for seq_obj in self.sequence_collection:
            fasta_file.extend([f'>{seq_obj.name}', f'{seq_obj.sequence}'])
        return fasta_file

class NucleotideSequence:
    def __init__(self, sequence, name=None):
        self.sequence = sequence
        self.length = len(sequence)
        self.name = name

def return_list_of_file_names_in_directory(directory_to_list):
    """
    return a list that contains the filenames found in the specified directory
    :param directory_to_list: the directory that the file names should be returned from
    :return: list of strings that are the file names found in the directory_to_list
    """
    list_of_file_names_in_directory = []
    for (dirpath, dirnames, filenames) in os.walk(directory_to_list):
        list_of_file_names_in_directory.extend(filenames)
        return list_of_file_names_in_directory

def return_list_of_file_paths_in_directory(directory_to_list):
    """
    return a list that contains the full paths of each of the files found in the specified directory
    :param directory_to_list: the directory that the file paths should be returned from
    :return: list of strings that are the file paths found in the directory_to_list
    """
    list_of_file_paths_in_directory = []
    for (dirpath, dirnames, filenames) in os.walk(directory_to_list):
        list_of_file_paths_in_directory.extend([os.path.join(directory_to_list, file_name) for file_name in filenames])
        return list_of_file_paths_in_directory

def return_list_of_directory_names_in_directory(directory_to_list):
    """
        return a list that contains the directory names found in the specified directory
        :param directory_to_list: the directory that the directory names should be returned from
        :return: list of strings that are the directory names found in the directory_to_list
        """
    list_of_directory_names_in_directory = []
    for (dirpath, dirnames, filenames) in os.walk(directory_to_list):
        list_of_directory_names_in_directory.extend(dirnames)
        return list_of_directory_names_in_directory


def return_list_of_directory_paths_in_directory(directory_to_list):
    """
        return a list that contains the full paths of each of the directories found in the specified directory
        :param directory_to_list: the directory that the directory paths should be returned from
        :return: list of strings that are the directory paths found in the directory_to_list
        """
    list_of_directory_paths_in_directory = []
    for (dirpath, dirnames, filenames) in os.walk(directory_to_list):
        list_of_directory_paths_in_directory.extend([os.path.join(directory_to_list, dir_name) for dir_name in dirnames])
        return list_of_directory_paths_in_directory

def read_defined_file_to_list(filename):
    with open(filename, mode='r') as reader:
        return [line.rstrip() for line in reader]


def read_defined_file_to_generator(filename):
    with open(filename, mode='r') as reader:
        return (line.rstrip() for line in reader)


def write_list_to_destination(destination, list_to_write):
    #print('Writing list to ' + destination)
    try:
        os.makedirs(os.path.dirname(destination))
    except FileExistsError:
        pass

    with open(destination, mode='w') as writer:
        i = 0
        while i < len(list_to_write):
            if i != len(list_to_write)-1:
                writer.write(list_to_write[i] + '\n')
            elif i == len(list_to_write)-1:
                writer.write(list_to_write[i])
            i += 1


def create_dict_from_fasta(fasta_list):
    temp_dict = {}
    i = 0
    while i < len(fasta_list):
        sequence = fasta_list[i][1:]
        temp_dict[sequence] = fasta_list[i + 1]
        i += 2
    return temp_dict


def remove_primer_mismatch_annotations_from_fasta(fasta_path):
    temp_fasta = []
    fasta_to_clean = read_defined_file_to_list(fasta_path)
    for i in range(len(fasta_to_clean) - 1):
        if fasta_to_clean[i]:
            if fasta_to_clean[i][0] == '>' and fasta_to_clean[i + 1]:
                if '|' in fasta_to_clean[i]:
                    temp_fasta.extend([fasta_to_clean[i].split('|')[0], fasta_to_clean[i + 1]])
                else:
                    temp_fasta.extend([fasta_to_clean[i].split('\t')[0], fasta_to_clean[i+1]])
    write_list_to_destination(fasta_path, temp_fasta)



def create_no_space_fasta_file(fasta_list):
    temp_list = []
    i = 0
    while i < len(fasta_list):
        temp_list.extend([fasta_list[i].split('\t')[0], fasta_list[i + 1]])
        i += 2
    return temp_list

def combine_two_fasta_files(path_one, path_two, path_for_combined):
    one_file_one = read_defined_file_to_list(path_one)
    one_file_two = read_defined_file_to_list(path_two)
    one_file_one.extend(one_file_two)
    write_list_to_destination(path_for_combined, one_file_one)