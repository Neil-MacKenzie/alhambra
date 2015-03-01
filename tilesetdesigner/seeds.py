from tileutils import gettile
import copy

class seed_tileadapts:
    def create_adapter_sequences( self, tileset ):
        tset = copy.deepcopy( tileset )

        for adapter in tset['seed']['adapters']:
            tilelong = gettile( tset, adapter['tilestrand'][0] )['fullseqs'][adapter['tilestrand'][1]-1]
            tileshort = gettile( tset, adapter['tilestrand'][0] )['fullseqs'][adapter['tilestrand'][1]+1-1]

            adapter['seqs'] = [ tilelong[0:8]+self.cores[adapter['loc']-1]+tilelong[40:48], tileshort ]

        return tset

class longrect(seed_tileadapts):
    cores = ['CAGGAACGGTACGCCATTAAAGGGATTTTAGA',
             'CTACATTTTGACGCTCACGCTCATGGAAATAC',
             'CCAGCAGAAGATAAAAAATACCGAACGAACCA',
             'GCCGTCAATAGATAATCAACTAATAGATTAGA',
             'ACTTCTGAATAATGGATGATTGTTTGGATTAT',
             'GAAGATGATGAAACAAAATTACCTGAGCAAAA',
             'CATAGGTCTGAGAGACGTGAATTTATCAAAAT',
             'GAAAAAGCCTGTTTAGGGAATCATAATTACTA',
             'ACGCGCCTGTTTATCAGTTCAGCTAATGCAGA',
             'GCTTATCCGGTATTCTAAATCAGATATAGAAG',
             'AACGTCAAAAATGAAAAAACGATTTTTTGTTT',
             'GCAGATAGCCGAACAATTTTTAAGAAAAGTAA',
             'AGACAAAAGGGCGACAGGTTTACCAGCGCCAA',
             'GCGTCAGACTGTAGCGATCAAGTTTGCCTTTA',
             'GTCAGACGATTGGCCTCAGGAGGTTGAGGCAG',
             'TGAAAGTATTAAGAGGCTATTATTCTGAAACA']
    
class triangle_side2(seed_tileadapts):
    cores = ['AGAGAGTACCTTTAATCCAACAGGTCAGGATT',
             'TAAGAGGAAGCCCGAAATTGCATCAAAAAGAT',
             'CCCCCTCAAATGCTTTTAAATATTCATTGAAT',
             'TACCAGACGACGATAATATCATAACCCTCGTT',
             'CGTTAATAAAACGAACTGGGAAGAAAAATCTA',
             'TAACAAAGCTGCTCATATTACCCAAATCAACG',
             'ATTGTGTCGAAATCCGTATCATCGCCTGATAA',
             'CGAGGGTAGCAACGGCAAAGACAGCATCGGAA',
             'TCTCCAAAAAAAAGGCTTTTTCACGTTGAAAA',
             'TTTCGTCACCAGTACAGTACCGTAACACTGAG',
             'TCAGTACCAGGCGGATATTAGCGGGGTTTTGC',
             'GATACAGGAGTGTACTATACATGGCTTTTGAT',
             'ACCAGAGCCGCCGCCACGCCACCAGAACCACC',
             'GCGTTTGCCATCTTTTCATAGCCCCCTTATTA',
             'ATGAAACCATCGATAGGCCGGAAACGTCACCA',
             'ATCACCGTCACCGACTTCATTAAAGGTGAATT']

seedtypes = { 'longrect': longrect(), 'triangle_side2': triangle_side2() }
