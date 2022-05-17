import unittest
import os
import subprocess
import porechop.misc


class TestCustomAdapters(unittest.TestCase):
    """
    This test tests custom adapter function as well as the additional Parameters for cropping and filtering. The basis of the test is the TestTwoAdapter sets
    This test set contains a combination of two adapter sets: SQK-MAP006 and SQK-NSK007
    """
    def run_command(self, command):
        runner_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'porechop-runner.py')
        input_path = os.path.join(os.path.dirname(__file__), 'test_custom_adapters.fastq')
        input_path_fasta = os.path.join(os.path.dirname(__file__), 'test_custom_adapters.fasta')
        command = command.replace('porechop', runner_path)
        command = command.replace('INPUT_FASTQ', input_path)
        command = command.replace('INPUT_FASTA', input_path_fasta)
        output_name = 'TEMP_' + str(os.getpid())
        command = command.replace('OUTPUT', output_name)
        try:
            self.output_file = [x for x in command.split() if output_name in x][0]
        except IndexError:
            self.output_file = ''
        p = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        out, err = p.communicate()
        return out.decode(), err.decode()

    def tearDown(self):
        if os.path.isfile(self.output_file):
            os.remove(self.output_file)

    def load_trimmed_reads(self):
        trimmed_reads, read_type = porechop.misc.load_fasta_or_fastq(self.output_file)
        if read_type == 'FASTA':
            trimmed_reads = [(x[2], x[1], '') for x in trimmed_reads]
        else:  # FASTQ
            trimmed_reads = [(x[4], x[1], x[3]) for x in trimmed_reads]
        return trimmed_reads, read_type

    def test_one_custom_adapter(self):
        """
        Checks with only one of the adapters (SQK-MAP006 as "Custom Adapter 1") is provided
        """
        out, _ = self.run_command('porechop -i INPUT_FASTQ -o OUTPUT.fastq --custom_adapter GGTTGTTTCTGTTGGTGCTGATATTGCT GCAATATCAGCACCAACAGAAA "Custom Adapter 1"')
        self.assertTrue('Trimming adapters' in out)
        self.assertTrue('Custom Adapter 1_(start)' in out)
        self.assertTrue('Custom Adapter 1_(end)' in out)
        self.assertTrue('MAP006' not in out)
        self.assertTrue('NSK007' not in out)
        self.assertTrue('1 / 4 reads had adapters trimmed from their start (28 bp removed)' in out)
        self.assertTrue('0 / 4 reads had adapters trimmed from their end (0 bp removed)' in out)
        self.assertTrue('0 bp are removed due to cropping.' in out)
        self.assertTrue('0 sequences are reversed.' in out)
        self.assertTrue('0 sequences are discarded due not meeting filter requirements' in out)
        self.assertTrue('Splitting reads containing middle adapters' in out)
        self.assertTrue('1 / 4 reads were split based on middle adapters' in out)

        trimmed_reads, read_type = self.load_trimmed_reads()
        
        self.assertEqual(len(trimmed_reads), 5)

        read_1 = trimmed_reads[0]
        self.assertEqual(read_1[0], '1')
        self.assertEqual(len(read_1[1]), 10000)
        self.assertTrue(read_1[1].startswith('GCGGAAGTCA'))
        self.assertTrue(read_1[1].endswith('TCCAGCCCCA'))
        if read_type == 'FASTQ':
            self.assertEqual(len(read_1[2]), 10000)
            self.assertTrue(read_1[2].startswith("4,)4+311.'"))
            self.assertTrue(read_1[2].endswith("3-16('...,"))

        read_2 = trimmed_reads[1]
        self.assertEqual(read_2[0], '2')
        self.assertEqual(len(read_2[1]), 10972)
        self.assertTrue(read_2[1].startswith('GTTGCAGCAC'))
        self.assertTrue(read_2[1].endswith('AAAGAACTGG'))
        if read_type == 'FASTQ':
            self.assertEqual(len(read_2[2]), 10972)
            self.assertTrue(read_2[2].startswith("#$&)*-++$+"))
            self.assertTrue(read_2[2].endswith("(*'%)%.2+)"))

        read_3 = trimmed_reads[2]
        self.assertEqual(read_3[0], '3')
        self.assertEqual(len(read_3[1]), 8000)
        self.assertTrue(read_3[1].startswith('AATGTACTTC'))
        self.assertTrue(read_3[1].endswith('TGAACGAAGT'))
        if read_type == 'FASTQ':
            self.assertEqual(len(read_3[2]), 8000)
            self.assertTrue(read_3[2].startswith("&-),1)5500"))
            self.assertTrue(read_3[2].endswith("/-((+1-/%."))

        read_4_1 = trimmed_reads[3]
        read_4_2 = trimmed_reads[4]
        self.assertEqual(read_4_1[0], '4_1')
        self.assertEqual(read_4_2[0], '4_2')
        self.assertEqual(len(read_4_1[1]), 1990)
        self.assertEqual(len(read_4_2[1]), 4878)
        self.assertTrue(read_4_1[1].startswith('TCATGACAGT'))
        self.assertTrue(read_4_1[1].endswith('CGAAGTAGGG'))
        self.assertTrue(read_4_2[1].startswith('ACGTGCTCGA'))
        self.assertTrue(read_4_2[1].endswith('AGTTCTCGTG'))
        if read_type == 'FASTQ':
            self.assertEqual(len(read_4_1[2]), 1990)
            self.assertEqual(len(read_4_2[2]), 4878)
            self.assertTrue(read_4_1[2].startswith("45.2*25-',"))
            self.assertTrue(read_4_1[2].endswith("$**+,,&$%3"))
            self.assertTrue(read_4_2[2].startswith("]213438::5"))
            self.assertTrue(read_4_2[2].endswith("+2+,853322"))

    def test_two_custom_adapter(self):
        """
        Checks with two custom  adapters (SQK-MAP006 as "Custom Adapter 1" and "SQK-NSK007" as "Custom Adapter 2") is provided
        """
        out, _ = self.run_command('porechop -i INPUT_FASTQ -o OUTPUT.fastq --custom_adapter GGTTGTTTCTGTTGGTGCTGATATTGCT GCAATATCAGCACCAACAGAAA "Custom Adapter 1" --custom_adapter AATGTACTTCGTTCAGTTACGTATTGCT GCAATACGTAACTGAACGAAGT "Custom Adapter 2"')

        self.assertTrue('Trimming adapters' in out)
        self.assertTrue('Custom Adapter 1_(start)' in out)
        self.assertTrue('Custom Adapter 1_(end)' in out)
        self.assertTrue('Custom Adapter 2_(start)' in out)
        self.assertTrue('Custom Adapter 2_(end)' in out)
        self.assertTrue('MAP006' not in out)
        self.assertTrue('NSK007' not in out)
        self.assertTrue('2 / 4 reads had adapters trimmed from their start (56 bp removed)' in out)
        self.assertTrue('1 / 4 reads had adapters trimmed from their end (22 bp removed)' in out)
        self.assertTrue('0 bp are removed due to cropping.' in out)
        self.assertTrue('0 sequences are reversed.' in out)
        self.assertTrue('0 sequences are discarded due not meeting filter requirements' in out)
        self.assertTrue('Splitting reads containing middle adapters' in out)
        self.assertTrue('1 / 4 reads were split based on middle adapters' in out)
        
        trimmed_reads, read_type = self.load_trimmed_reads()

        self.assertEqual(len(trimmed_reads), 6)
        
        read_1 = trimmed_reads[0]
        self.assertEqual(read_1[0], '1')
        self.assertEqual(len(read_1[1]), 10000)
        self.assertTrue(read_1[1].startswith('GCGGAAGTCA'))
        self.assertTrue(read_1[1].endswith('TCCAGCCCCA'))
        if read_type == 'FASTQ':
            self.assertEqual(len(read_1[2]), 10000)
            self.assertTrue(read_1[2].startswith("4,)4+311.'"))
            self.assertTrue(read_1[2].endswith("3-16('...,"))

        read_2 = trimmed_reads[1]
        self.assertEqual(read_2[0], '2')
        self.assertEqual(len(read_2[1]), 10972)
        self.assertTrue(read_2[1].startswith('GTTGCAGCAC'))
        self.assertTrue(read_2[1].endswith('AAAGAACTGG'))
        if read_type == 'FASTQ':
            self.assertEqual(len(read_2[2]), 10972)
            self.assertTrue(read_2[2].startswith("#$&)*-++$+"))
            self.assertTrue(read_2[2].endswith("(*'%)%.2+)"))

        read_3 = trimmed_reads[2]
        self.assertEqual(read_3[0], '3')
        self.assertEqual(len(read_3[1]), 7950)
        self.assertTrue(read_3[1].startswith('CCTATTTCGC'))
        self.assertTrue(read_3[1].endswith('GAGAGTGTGA'))
        if read_type == 'FASTQ':
            self.assertEqual(len(read_3[2]), 7950)
            self.assertTrue(read_3[2].startswith(",.2:2.4%%%"))
            self.assertTrue(read_3[2].endswith(",)*,)+00++"))

        read_4_1 = trimmed_reads[3]
        read_4_2 = trimmed_reads[4]
        read_4_3 = trimmed_reads[5]
        self.assertEqual(read_4_1[0], '4_1')
        self.assertEqual(read_4_2[0], '4_2')
        self.assertEqual(read_4_3[0], '4_3')
        self.assertEqual(len(read_4_1[1]), 1990)
        self.assertEqual(len(read_4_2[1]), 2868)
        self.assertEqual(len(read_4_3[1]), 1878)
        self.assertTrue(read_4_1[1].startswith('TCATGACAGT'))
        self.assertTrue(read_4_1[1].endswith('CGAAGTAGGG'))
        self.assertTrue(read_4_2[1].startswith('ACGTGCTCGA'))
        self.assertTrue(read_4_2[1].endswith('CCTTCCTAGA'))
        self.assertTrue(read_4_3[1].startswith('GGCTCCCCCC'))
        self.assertTrue(read_4_3[1].endswith('AGTTCTCGTG'))
        if read_type == 'FASTQ':
            self.assertEqual(len(read_4_1[2]), 1990)
            self.assertEqual(len(read_4_2[2]), 2868)
            self.assertEqual(len(read_4_3[2]), 1878)
            self.assertTrue(read_4_1[2].startswith("45.2*25-',"))
            self.assertTrue(read_4_1[2].endswith("$**+,,&$%3"))
            self.assertTrue(read_4_2[2].startswith("]213438::5"))
            self.assertTrue(read_4_2[2].endswith(".3/4:23]/]"))
            self.assertTrue(read_4_3[2].startswith("0430/6202,"))
            self.assertTrue(read_4_3[2].endswith("+2+,853322"))


    def test_cropping_tail(self):
        out, _ = self.run_command('porechop -i INPUT_FASTQ -o OUTPUT.fastq --custom_adapter GGTTGTTTCTGTTGGTGCTGATATTGCT GCAATATCAGCACCAACAGAAA "Custom Adapter 1" --custom_adapter AATGTACTTCGTTCAGTTACGTATTGCT GCAATACGTAACTGAACGAAGT "Custom Adapter 2" --tail_crop 100')

        self.assertTrue('Trimming adapters' in out)
        self.assertTrue('2 / 4 reads had adapters trimmed from their start (56 bp removed)' in out)
        self.assertTrue('1 / 4 reads had adapters trimmed from their end (22 bp removed)' in out)
        self.assertTrue('400 bp are removed due to cropping.' in out)
        self.assertTrue('0 sequences are reversed.' in out)
        self.assertTrue('0 sequences are discarded due not meeting filter requirements' in out)
        self.assertTrue('1 / 4 reads were split based on middle adapters' in out)
        
        trimmed_reads, read_type = self.load_trimmed_reads()

        self.assertEqual(len(trimmed_reads), 6)
        
        read_1 = trimmed_reads[0]
        self.assertEqual(read_1[0], '1')
        self.assertEqual(len(read_1[1]), 9900)
        self.assertTrue(read_1[1].startswith('GCGGAAGTCA'))
        self.assertTrue(read_1[1].endswith('GAACCAGTGG'))
        if read_type == 'FASTQ':
            self.assertEqual(len(read_1[2]), 9900)
            self.assertTrue(read_1[2].startswith("4,)4+311.'"))
            self.assertTrue(read_1[2].endswith("]]]754/,('"))

        read_2 = trimmed_reads[1]
        self.assertEqual(read_2[0], '2')
        self.assertEqual(len(read_2[1]), 10872)
        self.assertTrue(read_2[1].startswith('GTTGCAGCAC'))
        self.assertTrue(read_2[1].endswith('TTGGTCCGAT'))
        if read_type == 'FASTQ':
            self.assertEqual(len(read_2[2]), 10872)
            self.assertTrue(read_2[2].startswith("#$&)*-++$+"))
            self.assertTrue(read_2[2].endswith("%#.(&(%*'("))

        read_3 = trimmed_reads[2]
        self.assertEqual(read_3[0], '3')
        self.assertEqual(len(read_3[1]), 7850)
        self.assertTrue(read_3[1].startswith('CCTATTTCGC'))
        self.assertTrue(read_3[1].endswith('GGGTTATTAC'))
        if read_type == 'FASTQ':
            self.assertEqual(len(read_3[2]), 7850)
            self.assertTrue(read_3[2].startswith(",.2:2.4%%%"))
            self.assertTrue(read_3[2].endswith("0'($&#%'()"))

        read_4_1 = trimmed_reads[3]
        read_4_2 = trimmed_reads[4]
        read_4_3 = trimmed_reads[5]
        self.assertEqual(read_4_1[0], '4_1')
        self.assertEqual(read_4_2[0], '4_2')
        self.assertEqual(read_4_3[0], '4_3')
        self.assertEqual(len(read_4_1[1]), 1890)
        self.assertEqual(len(read_4_2[1]), 2768)
        self.assertEqual(len(read_4_3[1]), 1778)
        self.assertTrue(read_4_1[1].startswith('TCATGACAGT'))
        self.assertTrue(read_4_1[1].endswith('CAGACGGAAG'))
        self.assertTrue(read_4_2[1].startswith('ACGTGCTCGA'))
        self.assertTrue(read_4_2[1].endswith('CACGATTACC'))
        self.assertTrue(read_4_3[1].startswith('GGCTCCCCCC'))
        self.assertTrue(read_4_3[1].endswith('TCACACACCA'))
        if read_type == 'FASTQ':
            self.assertEqual(len(read_4_1[2]), 1890)
            self.assertEqual(len(read_4_2[2]), 2768)
            self.assertEqual(len(read_4_3[2]), 1778)
            self.assertTrue(read_4_1[2].startswith("45.2*25-',"))
            self.assertTrue(read_4_1[2].endswith("00254/3/63"))
            self.assertTrue(read_4_2[2].startswith("]213438::5"))
            self.assertTrue(read_4_2[2].endswith(".*)3+*+05+"))
            self.assertTrue(read_4_3[2].startswith("0430/6202,"))
            self.assertTrue(read_4_3[2].endswith("102--4533/"))

    def test_cropping_head(self):
        out, _ = self.run_command('porechop -i INPUT_FASTQ -o OUTPUT.fastq --custom_adapter GGTTGTTTCTGTTGGTGCTGATATTGCT GCAATATCAGCACCAACAGAAA "Custom Adapter 1" --custom_adapter AATGTACTTCGTTCAGTTACGTATTGCT GCAATACGTAACTGAACGAAGT "Custom Adapter 2" --head_crop 50')

        self.assertTrue('Trimming adapters' in out)
        self.assertTrue('2 / 4 reads had adapters trimmed from their start (56 bp removed)' in out)
        self.assertTrue('1 / 4 reads had adapters trimmed from their end (22 bp removed)' in out)
        self.assertTrue('200 bp are removed due to cropping.' in out)
        self.assertTrue('0 sequences are reversed.' in out)
        self.assertTrue('0 sequences are discarded due not meeting filter requirements' in out)
        self.assertTrue('1 / 4 reads were split based on middle adapters' in out)
        
        trimmed_reads, read_type = self.load_trimmed_reads()

        self.assertEqual(len(trimmed_reads), 6)
        
        read_1 = trimmed_reads[0]
        self.assertEqual(read_1[0], '1')
        self.assertEqual(len(read_1[1]), 9950)
        self.assertTrue(read_1[1].startswith('CGGTAGGGCA'))
        self.assertTrue(read_1[1].endswith('TCCAGCCCCA'))
        if read_type == 'FASTQ':
            self.assertEqual(len(read_1[2]), 9950)
            self.assertTrue(read_1[2].startswith("//3/013.+)"))
            self.assertTrue(read_1[2].endswith("3-16('...,"))

        read_2 = trimmed_reads[1]
        self.assertEqual(read_2[0], '2')
        self.assertEqual(len(read_2[1]), 10922)
        self.assertTrue(read_2[1].startswith('AAAGGGAGAA'))
        self.assertTrue(read_2[1].endswith('AAAGAACTGG'))
        if read_type == 'FASTQ':
            self.assertEqual(len(read_2[2]), 10922)
            self.assertTrue(read_2[2].startswith("+20*&.-.)1"))
            self.assertTrue(read_2[2].endswith("(*'%)%.2+)"))

        read_3 = trimmed_reads[2]
        self.assertEqual(read_3[0], '3')
        self.assertEqual(len(read_3[1]), 7900)
        self.assertTrue(read_3[1].startswith('ACTCAATTGA'))
        self.assertTrue(read_3[1].endswith('GAGAGTGTGA'))
        if read_type == 'FASTQ':
            self.assertEqual(len(read_3[2]), 7900)
            self.assertTrue(read_3[2].startswith("/$2]1%).]6"))
            self.assertTrue(read_3[2].endswith(",)*,)+00++"))

        read_4_1 = trimmed_reads[3]
        read_4_2 = trimmed_reads[4]
        read_4_3 = trimmed_reads[5]
        self.assertEqual(read_4_1[0], '4_1')
        self.assertEqual(read_4_2[0], '4_2')
        self.assertEqual(read_4_3[0], '4_3')
        self.assertEqual(len(read_4_1[1]), 1940)
        self.assertEqual(len(read_4_2[1]), 2818)
        self.assertEqual(len(read_4_3[1]), 1828)
        self.assertTrue(read_4_1[1].startswith('CGCAAGTAGA'))
        self.assertTrue(read_4_1[1].endswith('CGAAGTAGGG'))
        self.assertTrue(read_4_2[1].startswith('AGCACCCATT'))
        self.assertTrue(read_4_2[1].endswith('CCTTCCTAGA'))
        self.assertTrue(read_4_3[1].startswith('TATACTCTTT'))
        self.assertTrue(read_4_3[1].endswith('AGTTCTCGTG'))
        if read_type == 'FASTQ':
            self.assertEqual(len(read_4_1[2]), 1940)
            self.assertEqual(len(read_4_2[2]), 2818)
            self.assertEqual(len(read_4_3[2]), 1828)
            self.assertTrue(read_4_1[2].startswith("%&//10//45"))
            self.assertTrue(read_4_1[2].endswith("$**+,,&$%3"))
            self.assertTrue(read_4_2[2].startswith("+,--'&$.]]"))
            self.assertTrue(read_4_2[2].endswith(".3/4:23]/]"))
            self.assertTrue(read_4_3[2].startswith(")()1-&-'(0"))
            self.assertTrue(read_4_3[2].endswith("+2+,853322"))


    def test_min_length(self):
        out, _ = self.run_command('porechop -i INPUT_FASTQ -o OUTPUT.fastq --custom_adapter GGTTGTTTCTGTTGGTGCTGATATTGCT GCAATATCAGCACCAACAGAAA "Custom Adapter 1" --custom_adapter AATGTACTTCGTTCAGTTACGTATTGCT GCAATACGTAACTGAACGAAGT "Custom Adapter 2" --min_length 2000')

        self.assertTrue('Trimming adapters' in out)
        self.assertTrue('Custom Adapter 1_(start)' in out)
        self.assertTrue('Custom Adapter 1_(end)' in out)
        self.assertTrue('Custom Adapter 2_(start)' in out)
        self.assertTrue('Custom Adapter 2_(end)' in out)
        self.assertTrue('MAP006' not in out)
        self.assertTrue('NSK007' not in out)
        self.assertTrue('2 / 4 reads had adapters trimmed from their start (56 bp removed)' in out)
        self.assertTrue('1 / 4 reads had adapters trimmed from their end (22 bp removed)' in out)
        self.assertTrue('0 bp are removed due to cropping.' in out)
        self.assertTrue('0 sequences are reversed.' in out)
        self.assertTrue('0 sequences are discarded due not meeting filter requirements' in out)
        self.assertTrue('Splitting reads containing middle adapters' in out)
        self.assertTrue('1 / 4 reads were split based on middle adapters' in out)
        
        trimmed_reads, read_type = self.load_trimmed_reads()

        self.assertEqual(len(trimmed_reads), 4)
        
        read_1 = trimmed_reads[0]
        self.assertEqual(read_1[0], '1')
        self.assertEqual(len(read_1[1]), 10000)
        self.assertTrue(read_1[1].startswith('GCGGAAGTCA'))
        self.assertTrue(read_1[1].endswith('TCCAGCCCCA'))
        if read_type == 'FASTQ':
            self.assertEqual(len(read_1[2]), 10000)
            self.assertTrue(read_1[2].startswith("4,)4+311.'"))
            self.assertTrue(read_1[2].endswith("3-16('...,"))

        read_2 = trimmed_reads[1]
        self.assertEqual(read_2[0], '2')
        self.assertEqual(len(read_2[1]), 10972)
        self.assertTrue(read_2[1].startswith('GTTGCAGCAC'))
        self.assertTrue(read_2[1].endswith('AAAGAACTGG'))
        if read_type == 'FASTQ':
            self.assertEqual(len(read_2[2]), 10972)
            self.assertTrue(read_2[2].startswith("#$&)*-++$+"))
            self.assertTrue(read_2[2].endswith("(*'%)%.2+)"))

        read_3 = trimmed_reads[2]
        self.assertEqual(read_3[0], '3')
        self.assertEqual(len(read_3[1]), 7950)
        self.assertTrue(read_3[1].startswith('CCTATTTCGC'))
        self.assertTrue(read_3[1].endswith('GAGAGTGTGA'))
        if read_type == 'FASTQ':
            self.assertEqual(len(read_3[2]), 7950)
            self.assertTrue(read_3[2].startswith(",.2:2.4%%%"))
            self.assertTrue(read_3[2].endswith(",)*,)+00++"))

        read_4_2 = trimmed_reads[3]
        self.assertEqual(read_4_2[0], '4_2')
        self.assertEqual(len(read_4_2[1]), 2868)
        self.assertTrue(read_4_2[1].startswith('ACGTGCTCGA'))
        self.assertTrue(read_4_2[1].endswith('CCTTCCTAGA'))
        if read_type == 'FASTQ':
            self.assertEqual(len(read_4_2[2]), 2868)
            self.assertTrue(read_4_2[2].startswith("]213438::5"))
            self.assertTrue(read_4_2[2].endswith(".3/4:23]/]"))

    def test_max_length(self):
        out, _ = self.run_command('porechop -i INPUT_FASTQ -o OUTPUT.fastq --custom_adapter GGTTGTTTCTGTTGGTGCTGATATTGCT GCAATATCAGCACCAACAGAAA "Custom Adapter 1" --custom_adapter AATGTACTTCGTTCAGTTACGTATTGCT GCAATACGTAACTGAACGAAGT "Custom Adapter 2" --max_length 10000')

        self.assertTrue('Trimming adapters' in out)
        self.assertTrue('Custom Adapter 1_(start)' in out)
        self.assertTrue('Custom Adapter 1_(end)' in out)
        self.assertTrue('Custom Adapter 2_(start)' in out)
        self.assertTrue('Custom Adapter 2_(end)' in out)
        self.assertTrue('MAP006' not in out)
        self.assertTrue('NSK007' not in out)
        self.assertTrue('2 / 4 reads had adapters trimmed from their start (56 bp removed)' in out)
        self.assertTrue('1 / 4 reads had adapters trimmed from their end (22 bp removed)' in out)
        self.assertTrue('0 bp are removed due to cropping.' in out)
        self.assertTrue('0 sequences are reversed.' in out)
        self.assertTrue('1 sequences are discarded due not meeting filter requirements' in out)
        self.assertTrue('Splitting reads containing middle adapters' in out)
        self.assertTrue('1 / 4 reads were split based on middle adapters' in out)
        
        trimmed_reads, read_type = self.load_trimmed_reads()

        self.assertEqual(len(trimmed_reads), 5)
        
        read_1 = trimmed_reads[0]
        self.assertEqual(read_1[0], '1')
        self.assertEqual(len(read_1[1]), 10000)
        self.assertTrue(read_1[1].startswith('GCGGAAGTCA'))
        self.assertTrue(read_1[1].endswith('TCCAGCCCCA'))
        if read_type == 'FASTQ':
            self.assertEqual(len(read_1[2]), 10000)
            self.assertTrue(read_1[2].startswith("4,)4+311.'"))
            self.assertTrue(read_1[2].endswith("3-16('...,"))

        read_3 = trimmed_reads[1]
        self.assertEqual(read_3[0], '3')
        self.assertEqual(len(read_3[1]), 7950)
        self.assertTrue(read_3[1].startswith('CCTATTTCGC'))
        self.assertTrue(read_3[1].endswith('GAGAGTGTGA'))
        if read_type == 'FASTQ':
            self.assertEqual(len(read_3[2]), 7950)
            self.assertTrue(read_3[2].startswith(",.2:2.4%%%"))
            self.assertTrue(read_3[2].endswith(",)*,)+00++"))

        read_4_1 = trimmed_reads[2]
        read_4_2 = trimmed_reads[3]
        read_4_3 = trimmed_reads[4]
        self.assertEqual(read_4_1[0], '4_1')
        self.assertEqual(read_4_2[0], '4_2')
        self.assertEqual(read_4_3[0], '4_3')
        self.assertEqual(len(read_4_1[1]), 1990)
        self.assertEqual(len(read_4_2[1]), 2868)
        self.assertEqual(len(read_4_3[1]), 1878)
        self.assertTrue(read_4_1[1].startswith('TCATGACAGT'))
        self.assertTrue(read_4_1[1].endswith('CGAAGTAGGG'))
        self.assertTrue(read_4_2[1].startswith('ACGTGCTCGA'))
        self.assertTrue(read_4_2[1].endswith('CCTTCCTAGA'))
        self.assertTrue(read_4_3[1].startswith('GGCTCCCCCC'))
        self.assertTrue(read_4_3[1].endswith('AGTTCTCGTG'))
        if read_type == 'FASTQ':
            self.assertEqual(len(read_4_1[2]), 1990)
            self.assertEqual(len(read_4_2[2]), 2868)
            self.assertEqual(len(read_4_3[2]), 1878)
            self.assertTrue(read_4_1[2].startswith("45.2*25-',"))
            self.assertTrue(read_4_1[2].endswith("$**+,,&$%3"))
            self.assertTrue(read_4_2[2].startswith("]213438::5"))
            self.assertTrue(read_4_2[2].endswith(".3/4:23]/]"))
            self.assertTrue(read_4_3[2].startswith("0430/6202,"))
            self.assertTrue(read_4_3[2].endswith("+2+,853322"))

    def test_trimmed_only(self):
        out, _ = self.run_command('porechop -i INPUT_FASTQ -o OUTPUT.fastq --custom_adapter GGTTGTTTCTGTTGGTGCTGATATTGCT GCAATATCAGCACCAACAGAAA "Custom Adapter 1" --custom_adapter AATGTACTTCGTTCAGTTACGTATTGCT GCAATACGTAACTGAACGAAGT "Custom Adapter 2" --trimmed_only')

        self.assertTrue('Trimming adapters' in out)
        self.assertTrue('Custom Adapter 1_(start)' in out)
        self.assertTrue('Custom Adapter 1_(end)' in out)
        self.assertTrue('Custom Adapter 2_(start)' in out)
        self.assertTrue('Custom Adapter 2_(end)' in out)
        self.assertTrue('MAP006' not in out)
        self.assertTrue('NSK007' not in out)
        self.assertTrue('2 / 4 reads had adapters trimmed from their start (56 bp removed)' in out)
        self.assertTrue('1 / 4 reads had adapters trimmed from their end (22 bp removed)' in out)
        self.assertTrue('0 bp are removed due to cropping.' in out)
        self.assertTrue('0 sequences are reversed.' in out)
        self.assertTrue('2 sequences are discarded due not meeting filter requirements' in out)
        self.assertTrue('Splitting reads containing middle adapters' in out)
        self.assertTrue('1 / 4 reads were split based on middle adapters' in out)
        
        trimmed_reads, read_type = self.load_trimmed_reads()

        self.assertEqual(len(trimmed_reads), 5)

        read_2 = trimmed_reads[0]
        self.assertEqual(read_2[0], '2')
        self.assertEqual(len(read_2[1]), 10972)
        self.assertTrue(read_2[1].startswith('GTTGCAGCAC'))
        self.assertTrue(read_2[1].endswith('AAAGAACTGG'))
        if read_type == 'FASTQ':
            self.assertEqual(len(read_2[2]), 10972)
            self.assertTrue(read_2[2].startswith("#$&)*-++$+"))
            self.assertTrue(read_2[2].endswith("(*'%)%.2+)"))

        read_3 = trimmed_reads[1]
        self.assertEqual(read_3[0], '3')
        self.assertEqual(len(read_3[1]), 7950)
        self.assertTrue(read_3[1].startswith('CCTATTTCGC'))
        self.assertTrue(read_3[1].endswith('GAGAGTGTGA'))
        if read_type == 'FASTQ':
            self.assertEqual(len(read_3[2]), 7950)
            self.assertTrue(read_3[2].startswith(",.2:2.4%%%"))
            self.assertTrue(read_3[2].endswith(",)*,)+00++"))

        read_4_1 = trimmed_reads[2]
        read_4_2 = trimmed_reads[3]
        read_4_3 = trimmed_reads[4]
        self.assertEqual(read_4_1[0], '4_1')
        self.assertEqual(read_4_2[0], '4_2')
        self.assertEqual(read_4_3[0], '4_3')
        self.assertEqual(len(read_4_1[1]), 1990)
        self.assertEqual(len(read_4_2[1]), 2868)
        self.assertEqual(len(read_4_3[1]), 1878)
        self.assertTrue(read_4_1[1].startswith('TCATGACAGT'))
        self.assertTrue(read_4_1[1].endswith('CGAAGTAGGG'))
        self.assertTrue(read_4_2[1].startswith('ACGTGCTCGA'))
        self.assertTrue(read_4_2[1].endswith('CCTTCCTAGA'))
        self.assertTrue(read_4_3[1].startswith('GGCTCCCCCC'))
        self.assertTrue(read_4_3[1].endswith('AGTTCTCGTG'))
        if read_type == 'FASTQ':
            self.assertEqual(len(read_4_1[2]), 1990)
            self.assertEqual(len(read_4_2[2]), 2868)
            self.assertEqual(len(read_4_3[2]), 1878)
            self.assertTrue(read_4_1[2].startswith("45.2*25-',"))
            self.assertTrue(read_4_1[2].endswith("$**+,,&$%3"))
            self.assertTrue(read_4_2[2].startswith("]213438::5"))
            self.assertTrue(read_4_2[2].endswith(".3/4:23]/]"))
            self.assertTrue(read_4_3[2].startswith("0430/6202,"))
            self.assertTrue(read_4_3[2].endswith("+2+,853322"))

    def test_correct_read_correction(self):
        out, _ = self.run_command('porechop -i INPUT_FASTQ -o OUTPUT.fastq --custom_adapter GGTTGTTTCTGTTGGTGCTGATATTGCT GCAATATCAGCACCAACAGAAA "Custom Adapter 1" --custom_adapter AATGTACTTCGTTCAGTTACGTATTGCT GCAATACGTAACTGAACGAAGT "Custom Adapter 2 reverse" --correct_read_direction')

        self.assertTrue('Trimming adapters' in out)
        self.assertTrue('Custom Adapter 1_(start)' in out)
        self.assertTrue('Custom Adapter 1_(end)' in out)
        self.assertTrue('Custom Adapter 2 reverse_(start)' in out)
        self.assertTrue('Custom Adapter 2 reverse_(end)' in out)
        self.assertTrue('MAP006' not in out)
        self.assertTrue('NSK007' not in out)
        self.assertTrue('2 / 4 reads had adapters trimmed from their start (56 bp removed)' in out)
        self.assertTrue('1 / 4 reads had adapters trimmed from their end (22 bp removed)' in out)
        self.assertTrue('0 bp are removed due to cropping.' in out)
        self.assertTrue('1 sequences are reversed.' in out)
        self.assertTrue('0 sequences are discarded due not meeting filter requirements' in out)
        self.assertTrue('Splitting reads containing middle adapters' in out)
        self.assertTrue('1 / 4 reads were split based on middle adapters' in out)
        
        trimmed_reads, read_type = self.load_trimmed_reads()

        self.assertEqual(len(trimmed_reads), 6)
        
        read_1 = trimmed_reads[0]
        self.assertEqual(read_1[0], '1')
        self.assertEqual(len(read_1[1]), 10000)
        self.assertTrue(read_1[1].startswith('GCGGAAGTCA'))
        self.assertTrue(read_1[1].endswith('TCCAGCCCCA'))
        if read_type == 'FASTQ':
            self.assertEqual(len(read_1[2]), 10000)
            self.assertTrue(read_1[2].startswith("4,)4+311.'"))
            self.assertTrue(read_1[2].endswith("3-16('...,"))

        read_2 = trimmed_reads[1]
        self.assertEqual(read_2[0], '2')
        self.assertEqual(len(read_2[1]), 10972)
        self.assertTrue(read_2[1].startswith('GTTGCAGCAC'))
        self.assertTrue(read_2[1].endswith('AAAGAACTGG'))
        if read_type == 'FASTQ':
            self.assertEqual(len(read_2[2]), 10972)
            self.assertTrue(read_2[2].startswith("#$&)*-++$+"))
            self.assertTrue(read_2[2].endswith("(*'%)%.2+)"))

        read_3 = trimmed_reads[2]
        self.assertEqual(read_3[0], '3')
        self.assertEqual(len(read_3[1]), 7950)
        self.assertTrue(read_3[1].startswith('TCACACTCTC'))
        self.assertTrue(read_3[1].endswith('GCGAAATAGG'))
        if read_type == 'FASTQ':
            self.assertEqual(len(read_3[2]), 7950)
            self.assertTrue(read_3[2].startswith("++00+),*),"))
            self.assertTrue(read_3[2].endswith("%%%4.2:2.,"))

        read_4_1 = trimmed_reads[3]
        read_4_2 = trimmed_reads[4]
        read_4_3 = trimmed_reads[5]
        self.assertEqual(read_4_1[0], '4_1')
        self.assertEqual(read_4_2[0], '4_2')
        self.assertEqual(read_4_3[0], '4_3')
        self.assertEqual(len(read_4_1[1]), 1990)
        self.assertEqual(len(read_4_2[1]), 2868)
        self.assertEqual(len(read_4_3[1]), 1878)
        self.assertTrue(read_4_1[1].startswith('TCATGACAGT'))
        self.assertTrue(read_4_1[1].endswith('CGAAGTAGGG'))
        self.assertTrue(read_4_2[1].startswith('ACGTGCTCGA'))
        self.assertTrue(read_4_2[1].endswith('CCTTCCTAGA'))
        self.assertTrue(read_4_3[1].startswith('GGCTCCCCCC'))
        self.assertTrue(read_4_3[1].endswith('AGTTCTCGTG'))
        if read_type == 'FASTQ':
            self.assertEqual(len(read_4_1[2]), 1990)
            self.assertEqual(len(read_4_2[2]), 2868)
            self.assertEqual(len(read_4_3[2]), 1878)
            self.assertTrue(read_4_1[2].startswith("45.2*25-',"))
            self.assertTrue(read_4_1[2].endswith("$**+,,&$%3"))
            self.assertTrue(read_4_2[2].startswith("]213438::5"))
            self.assertTrue(read_4_2[2].endswith(".3/4:23]/]"))
            self.assertTrue(read_4_3[2].startswith("0430/6202,"))
            self.assertTrue(read_4_3[2].endswith("+2+,853322"))


    def test_everything_with_fasta(self):
        out, _ = self.run_command('porechop -i INPUT_FASTA -o OUTPUT.fasta --custom_adapter GGTTGTTTCTGTTGGTGCTGATATTGCT GCAATATCAGCACCAACAGAAA "Custom Adapter 1" --custom_adapter AATGTACTTCGTTCAGTTACGTATTGCT GCAATACGTAACTGAACGAAGT "Custom Adapter 2 reverse" --correct_read_direction --trimmed_only --min_length 2000 --max_length 10000 --head_crop 50 --tail_crop 100')

        self.assertTrue('Trimming adapters' in out)
        self.assertTrue('Custom Adapter 1_(start)' in out)
        self.assertTrue('Custom Adapter 1_(end)' in out)
        self.assertTrue('Custom Adapter 2 reverse_(start)' in out)
        self.assertTrue('Custom Adapter 2 reverse_(end)' in out)
        self.assertTrue('MAP006' not in out)
        self.assertTrue('NSK007' not in out)
        self.assertTrue('2 / 4 reads had adapters trimmed from their start (56 bp removed)' in out)
        self.assertTrue('1 / 4 reads had adapters trimmed from their end (22 bp removed)' in out)
        self.assertTrue('600 bp are removed due to cropping.' in out)
        self.assertTrue('1 sequences are reversed.' in out)
        self.assertTrue('3 sequences are discarded due not meeting filter requirements' in out)
        self.assertTrue('Splitting reads containing middle adapters' in out)
        self.assertTrue('1 / 4 reads were split based on middle adapters' in out)
        
        trimmed_reads, read_type = self.load_trimmed_reads()

        self.assertEqual(len(trimmed_reads), 2)

        read_3 = trimmed_reads[0]
        self.assertEqual(read_3[0], '3')
        self.assertEqual(len(read_3[1]), 7800)
        self.assertTrue(read_3[1].startswith('GTAATAACCC'))
        self.assertTrue(read_3[1].endswith('TCAATTGAGT'))

        read_4_2 = trimmed_reads[1]
        self.assertEqual(read_4_2[0], '4_2')
        self.assertEqual(len(read_4_2[1]), 2718)
        self.assertTrue(read_4_2[1].startswith('AGCACCCATT'))
        self.assertTrue(read_4_2[1].endswith('CACGATTACC'))
