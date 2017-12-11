from .tilestructures import gettile, check_edotparen_sequence, check_edotparen_consistency
import warnings
import copy

class seed_base:
    def create_adapter_sequence_diagram( self, adapterdef ):
        # At the moment, there are only three potential structures, so I'm just
        # going to put them here.
        ttype = 'tile_adapter_5up'
        if 'extra' in adapterdef.keys():
            ttype += '_'+adapterdef['extra']
        
        from lxml import etree
        import pkg_resources
        import os
        base_svg = etree.parse( pkg_resources.resource_stream(__name__, os.path.join('seqdiagrambases',ttype+'.svg') ) )

        try:
            s = adapterdef['seqs'][1]
            l = adapterdef['seqs'][0]
        except KeyError:
            s = adapterdef['fullseqs'][1]
            l = adapterdef['fullseqs'][0]
        name = adapterdef['name']
        if 'extra' not in adapterdef.keys():
            strings = [(s[0:5]+'--'+s[5:13])[::-1],s[13:21]+'--'+s[21:26], \
                      (l[0:8]+'--'+l[8:16])[::-1],l[16:24][::-1],l[24:32],l[32:40]+'--'+l[40:48], \
                      adapterdef['ends'][0],adapterdef['ends'][1],name]
        elif adapterdef['extra'] == '1h':
            hp = s[0:13]; s = s[13:]
            strings = [hp[0:5]+'-'+hp[5:7]+'-'+hp[7:9], \
                      (hp[9:11]+'-'+hp[11:13]+'-'+s[0:5]+'--'+s[5:13])[::-1],s[13:21]+'--'+s[21:26], \
                      (l[0:8]+'--'+l[8:16])[::-1],l[16:24][::-1],l[24:32],l[32:40]+'--'+l[40:48], \
                      adapterdef['ends'][1],name]
        elif adapterdef['extra'] == '2h':
            strings = [(s[0:5]+'--'+s[5:13])[::-1],s[13:21]+'--'+s[21:26]+'-'+s[26:28]+'-'+s[28:30], \
                      (s[30:32]+'-'+s[32:34]+'-'+s[34:39])[::-1], \
                      (l[0:8]+'--'+l[8:16])[::-1],l[16:24][::-1],l[24:32],l[32:40]+'--'+l[40:48], \
                      adapterdef['ends'][0],name]

        texts = base_svg.findall( "//{http://www.w3.org/2000/svg}text" )

        for text, string in zip(texts,strings):
            text.text = string

        return base_svg.xpath('/svg:svg/svg:g',namespaces={'svg':'http://www.w3.org/2000/svg'})[0]

class seed_tileadapts(seed_base):
    needspepper = False
    def create_adapter_sequences( self, tileset ):
        tset = copy.deepcopy( tileset )

        for adapter in tset['seed']['adapters']:
            mimic = self.getmimic_full( tset, adapter )
            if mimic['tile']['structure'] == 'tile_daoe_5up':
                # Tile is single. If we have ends, check that they match.
                if (('ends' in adapter.keys()) and (adapter['ends'] != mimic['tile']['ends'][1:3])):
                    raise ValueError("adapter {} and tile base {} ends don't match: adapter has {}, tile has {}".format(
                        adapter['name'], mimic['tile']['name'], adapter['ends'], mimic['tile']['ends'][1:3] ))
                adapter_strand_short = mimic['tile']['fullseqs'][-1]
                tile_long_strand = mimic['tile']['fullseqs'][-2]
            elif mimic['tile']['structure'] == 'tile_daoe_doublehoriz_35up' \
                    or mimic['tile']['structure'] == 'tile_daoe_doublevert_35up':
                if ( 'ends' in adapter.keys() ) and ( adapter['ends'] != mimic['tile']['ends'][2:4] ):
                    raise ValueError("adapter {} and tile base {} ends don't match: adapter has {}, tile has {}".format(
                        adapter['name'], mimic['tile']['name'], adapter['ends'], mimic['tile']['ends'][2:4] ))
                adapter_strand_short = mimic['tile']['fullseqs'][-1]
                tile_long_strand = mimic['tile']['fullseqs'][-2]

            adapter['seqs'] = [ tile_long_strand[0:8]+self.cores[adapter['loc']-1]+tile_long_strand[40:48], 
                                adapter_strand_short ]
            if 'extra' in mimic.keys():
                adapter['extra'] = mimic['extra']
        return tset

    def getmimic( self, tileset, adapter ):
        try:
            mimic = gettile( tileset, adapter['tilebase'] )
        except KeyError:
            warnings.warn("tilestrand format is deprecated; strand number will be ignored")
            mimic = gettile( tileset, adapter['tilestrand'][0] )
        return mimic

    def getmimic_full( self, tileset, adapter ):
        import re
        mimic = self.getmimic( tileset, adapter )
        out = {}
        e1 = None
        e2 = None
        if mimic['structure'] == 'tile_daoe_5up':
            out['tilestrand'] = mimic['fullseqs'][3]
            out['tileends'] = mimic['ends'][1:3]
            if 'extra' in mimic.keys():
                e1 = re.search(r'2([A-Za-z]+)',mimic['extra'])
                e2 = re.search(r'3([A-Za-z]+)',mimic['extra'])
        elif mimic['structure'] == 'tile_daoe_doublehoriz_35up' or \
           mimic['structure'] == 'tile_daoe_doublevert_35up':
            out['tilestrand'] = mimic['fullseqs'][5]
            out['tileends'] = mimic['ends'][2:4]
            if 'extra' in mimic.keys():
                e1 = re.search(r'3([A-Za-z]+)',mimic['extra'])
                e2 = re.search(r'4([A-Za-z]+)',mimic['extra'])
        else:
            raise ValueError("Tile type {} not known for adapter.".format(
                mimic['structure']), adapter, mimic)
        if e1:
            if e1[1] == 'h':
                out['extra']='1h'
            else:
                raise ValueError("Tile has unknown extra {}".format(e1),
                                 adapter, mimic)
        if e2:
            if e2[1] == 'h':
                out['extra']='2h'
            else:
                raise ValueError("Tile has unknown extra {}".format(e2),
                                 adapter, mimic)
        out['tile'] = mimic
        return out
            
                
    
    def check_consistent( self, tileset ):
        import re
        for adapt in tileset['seed']['adapters']:
            mimic = self.getmimic_full( tileset, adapt )
            assert mimic.get('extra') == adapt.get('extra')
            assert mimic['tileends'] == adapt['ends']
            

    edotparen_adapt = "8(32.8(+5.16)5."
    edotparen_adapt_2h = "8(32.8(+5.16)7(4.7)"
    edotparen_adapt_1h = "8(32.8(+7(4.7)16)5."
    
    def check_sequence( self, tileset ):
        import re
        check_edotparen_consistency( self.edotparen_adapt )
        check_edotparen_consistency( self.edotparen_adapt_1h )
        check_edotparen_consistency( self.edotparen_adapt_2h )
        self.check_consistent( tileset )
        for adapt in tileset['seed']['adapters']:
            if 'seqs' not in adapt.keys():
                continue
            mimic = self.getmimic_full(tileset, adapt)
            seqs = "+".join(adapt['seqs'])
            if 'extra' not in adapt.keys():
                check_edotparen_sequence(
                    self.edotparen_adapt, seqs )
            elif adapt['extra'] == '1h':
                check_edotparen_sequence(
                    self.edotparen_adapt_1h, seqs )
            elif adapt['extra'] == '2h':
                check_edotparen_sequence(
                    self.edotparen_adapt_2h, seqs )
            else:
                raise ValueError("Unknown extra", adapt['extra'])
            assert adapt['seqs'][0][8:(8+32)] == self.cores[adapt['loc']-1]
                
        
class tallrect_tileadapts(seed_tileadapts):
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
    # Note: this list is currently reversed because
    # the cadnano file for the triangle has 5' ends of
    # the last staple above 3' ends, while our adapter
    # scheme has 3' ends above 5' ends. 
    cores = list(reversed( \
            ['AGAGAGTACCTTTAATCCAACAGGTCAGGATT',
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
             'ATCACCGTCACCGACTTCATTAAAGGTGAATT']))

class tallrect_endadapts(seed_base):
    needspepper = True
    cores = tallrect_tileadapts.cores
    def create_pepper_input_files( self, seeddef, importlist, compstring ):
        for adapter in seeddef['adapters']:
            adaptertype = 'tile_adapter_5up'
            e = [[],[]]
            e[0].append('origamicore_{0}'.format(adapter['loc']))
            for end in adapter['ends']:
                if (end == 'hp'):
                    continue # skip hairpins, etc that aren't designed by stickydesign
                e[0].append( 'e_'+end.replace('/','*') )
                if end[-1]=='/':
                    a = 'c_'+end[:-1]+'*'
                else:
                    a = 'a_'+end
                e[1].append( a )
            s1 = " + ".join(e[0])
            s2 = " + ".join(e[1])
            if 'extra' in adapter.keys():
                adaptertype+='_'+adapter['extra']
            compstring += "component {} = {}: {} -> {}\n".format(adapter['name'],adaptertype,s1,s2)
            importlist.add( adaptertype )
        return (importlist, compstring)

seedtypes = { 
        'tallrect_tileadapts': tallrect_tileadapts(),
        'tallrect_endadapts': tallrect_endadapts(),
        'triangle_side2': triangle_side2()
        }
seedtypes['longrect'] = seedtypes['tallrect_tileadapts']
seedtypes['longrect_old'] = seedtypes['tallrect_endadapts']
