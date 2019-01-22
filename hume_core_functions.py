import os

class MothurAnalysis:

    def __init__(
            self, name, sequence_collection,  input_dir, output_dir,
            fastq_gz_fwd_path=None, fastq_gz_rev_path=None,  fasta_path=None,
             name_file_path=None, mothur_execution_string='mothur'):

        self.name = name
        self.sequence_collection = sequence_collection
        self.exec_str = mothur_execution_string
        self.input_dir = input_dir
        self.output_dir = output_dir
        self.fastq_gz_fwd_path = fastq_gz_fwd_path
        self.fastq_gz_rev_path = fastq_gz_rev_path
        self.fasta_path = fasta_path
        self.name_file_path = name_file_path

    @classmethod
    def from_pair_of_fastq_gz_files(cls, name, fastq_gz_fwd_path, fastq_gz_rev_path,
                                    output_dir=None, mothur_execution_string='mothur'):
        return cls(name=name, sequence_collection=None, mothur_execution_string=mothur_execution_string,
                   input_dir=os.path.dirname(os.path.abspath(fastq_gz_fwd_path)), output_dir=output_dir,
                   fastq_gz_fwd_path=fastq_gz_fwd_path, fastq_gz_rev_path=fastq_gz_rev_path,
                   fasta_path=None, name_file_path=None)

    @classmethod
    def from_sequence_collection(cls, name, sequence_collection, input_dir, output_dir,
                                 mothur_execution_string='mothur'):
        if sequence_collection.file_type == 'fasta':
            fasta_path = sequence_collection.file_path
        else:
            raise ValueError('SequenceCollection must be of type fasta')
        return cls(name=name, sequence_collection=sequence_collection, input_dir=input_dir, output_dir=output_dir)



class SequenceCollection:
    """ A sequence collection is a set of sequences either generated from a fastq file or from a fasta file.
    It cannot be created directly from binary files or from paired files. As such, to generate a SequenceCollection
    for example from a pair of fastaq.gz files, you would first have to run a mothur contig analysis and create
    the SeqeunceCollection from the resultant fasta file that is generated."""
    def __init__(self, name, path_to_file=None):
        self.name = name
        self.file_path = path_to_file
        # self.file_as_list = read_defined_file_to_list(self.file_path)
        self.file_type = self.infer_file_type()
        self.sequence_collection = self.generate_sequence_collection()

    def __len__(self):
        return(len(self.sequence_collection))

    def generate_sequence_collection(self):
        if self.file_type == 'fasta':
            return self.parse_fasta_file_and_extract_nucleotide_sequence_objects()
        elif self.file_type == 'fastq':
            return self.parse_fastq_file_and_extract_nucleotide_sequence_objects()

    def parse_fasta_file_and_extract_nucleotide_sequence_objects(self):
        list_of_nuleotide_sequence_objects = []
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
        if fastsq_file[index_value].startswith('@') and fastsq_file[index_value + 1][0] == '+':
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


def remove_annotations_from_fasta(scrapped_seq_fasta_path):
    temp_fasta = []
    fasta_to_clean = read_defined_file_to_list(scrapped_seq_fasta_path)
    for i in range(len(fasta_to_clean) - 1):
        if fasta_to_clean[i]:
            if fasta_to_clean[i][0] == '>' and fasta_to_clean[i + 1]:
                temp_fasta.extend([fasta_to_clean[i].split('|')[0], fasta_to_clean[i+1]])
    return temp_fasta


def create_no_space_fasta_file(fasta_list):
    temp_list = []
    i = 0
    while i < len(fasta_list):
        temp_list.extend([fasta_list[i].split('\t')[0], fasta_list[i + 1]])
        i += 2
    return temp_list
