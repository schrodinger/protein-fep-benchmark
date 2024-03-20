'''
mutation.py - Handle protein residue mutation information.

'''
import re

from schrodinger import structure
from schrodinger.application.scisol.packages.fep import graph
from schrodinger.utils import log

# Local
from . import pka


# Configure logging
logger = log.get_output_logger(__name__)
# logger.format = '%(message)s'
logger.level = log.WARNING


WT_NODE_TITLE = "WT"
DEFAULT_MUTATIONS1_FILENAME = 'mutations1.txt'
DEFAULT_MUTATIONS_FILENAME = 'mutations.txt'


# {chain}:{start_aa1}{resnum}{optional_inscode}{end_aa1}
# DEFAULT_MUTATION_FORMAT = 'c:sri?e'  # e.g. A:S30AT  for A:SER30A->THR
MUTSTR1_REGEX = re.compile(r'''
    (?P<chain>[A-Z])        # 1-letter chain ID
    [:-]                    # either colon or hyphen separator
    (?P<start1>[A-Z])?      # optional starting 1-letter residue code
    (?P<resnum>\d+)         # numeric residue number
    (?P<inscode>[A-Z]?)     # optional insertion code
    (?:->)?                 # optional mutation arrow
    (?P<end1>[A-Z])         # ending 1-letter residue code
    ''', re.VERBOSE)


def parse_single_mutstr1(mutstr):
    logger.debug(f'single mut: `{mutstr}`')
    match = re.match(MUTSTR1_REGEX, mutstr)
    try:
        return match.groupdict()
    except AttributeError:
        logger.warn(f'unable to parse malformed mutation string: `{mutstr}`')
        return None


def sort_mut1_dict_list(mut1_dict_list):
    '''
    Sensibly sort a list of mutation dicts.
    '''
    return sorted(
        mut1_dict_list,
        key=lambda d:
            (d['chain'], d['resnum'], d['inscode'], d['start1'], d['end1'])
    )


def parse_mutations1_file(mutfile, sep=','):
    '''
    Parse a mutations1.txt file to a list of mutation dicts.
    '''
    with open(mutfile, 'r') as f:
        lines = [l.rstrip('\n') for l in f.readlines()]
    mutations = []
    for l in lines:
        muts = [parse_single_mutstr1(m) for m in l.split(sep)]
        sorted_muts = sort_mut1_dict_list(muts)
        mutations.append(sorted_muts)
    return mutations


_RS_MUTSTR_REGEX = re.compile(r'''
    (?P<chain>[A-Z])        # 1-letter chain ID
    :                       # colon separator
    (?P<resnum>\d+)         # numeric residue number
    (?P<inscode>[A-Z]?)     # optional insertion code
    \(                      # literal opening parenthesis
    (?P<start>[A-Z]{3})     # starting 3-letter residue code
    ->                      # mutation arrow
    (?P<end>[A-Z]{3})       # ending 3-letter residue code
    \)                      # literal closing parenthesis
    ''', re.VERBOSE)


def parse_residue_scanning_single_mutstr(mutstr):
    '''
    Parse a single mutation string from Residue Scanning.
    '''
    match = re.match(RS_MUTSTR_REGEX, mutstr)
    try:
        return match.groupdict()
    except AttributeError:
        logger.warn(f'unable to parse malformed mutation string: `{mutstr}`')
        return None


def parse_residue_scanning_mutstr(mutstr):
    '''
    Parse a mutations string (`Mutations` column value) from Residue Scanning.
    '''
    return [
        parse_residue_scanning_single_mutstr(m)
        for m in mutstr.split(",")
    ]


class ProteinMutation():
    '''
    Details about a protein mutation.
    '''
    def __init__(self, chain, resnum, inscode, start, end):
        self.chain = chain
        self.resnum = int(resnum)
        self.inscode = inscode
        self.start = start
        self.end = end

    @property
    def start1(self):
        try:
            return structure.RESIDUE_MAP_3_TO_1_LETTER[self.start]
        except KeyError:
            return self.start

    @property
    def end1(self):
        try:
            return structure.RESIDUE_MAP_3_TO_1_LETTER[self.end]
        except KeyError:
            return self.end

    @property
    def full_resnum(self):
        return f'{self.resnum}{self.inscode}'

    @property
    def biolum_mutstr(self):
        return f'{self.chain}:{self.full_resnum}({self.start}->{self.end})'

    @property
    def one_letter(self):
        if self.start1 == self.end1:
            return None
        return f'{self.chain}-{self.start1}{self.full_resnum}{self.end1}'

    @property
    def position(self):
        return f'{self.chain}:{self.full_resnum}'

    @property
    def site(self):
        return self.position

    # TODO add tests for this one
    @property
    def is_charged(self):
        # TODO: get charged list somewhere else
        charged = ['ASP', 'GLU', 'HIP', 'LYS', 'ARG']
        endpoints = [self.start, self.end]
        return len(set(charged) & set(endpoints)) > 0

    def __eq__(self, rhs):
        if isinstance(rhs, ProteinMutation):
            return all([
                self.__getattribute__(x) == rhs.__getattribute__(x)
                for x in ['chain', 'resnum', 'inscode', 'start', 'end']
            ])
        return False

    def __add__(self, rhs):
        if isinstance(rhs, ProteinMutation):
            if all([
                self.__getattribute__(x) == rhs.__getattribute__(x)
                for x in ['chain', 'resnum', 'inscode', 'start']
            ]):
                # both are mutations of the same position, i.e.
                # A->B + A->C = A->C
                return ProteinMutation(self.chain, self.resnum, self.inscode,
                                       self.start, rhs.end)
            elif all([
                self.__getattribute__(x) == rhs.__getattribute__(x)
                for x in ['chain', 'resnum', 'inscode']
            ]) and (self.end == rhs.start):
                # 2 mutations of the same position, applied sequentially, i.e.
                # A->B + B->C = A->C
                return ProteinMutation(self.chain, self.resnum, self.inscode,
                                       self.start, rhs.end)
            else:
                # mutations at different positions
                return [self, rhs]
        elif rhs is None:
            return self
        else:
            return NotImplemented

    def __radd__(self, lhs):
        if lhs is None:
            return self
        else:
            return NotImplemented

    def __sub__(self, rhs):
        if isinstance(rhs, ProteinMutation):
            if all([
                self.__getattribute__(x) == rhs.__getattribute__(x)
                for x in ['chain', 'resnum', 'inscode', 'start']
            ]):
                # 2 mutations of the same position, return the actual mutation
                return ProteinMutation(self.chain, self.resnum, self.inscode,
                                       rhs.end, self.end)
            else:
                # mutations at different positions
                return [reversed(rhs), self]
        elif rhs is None:
            return self
        else:
            return NotImplemented

    def __rsub__(self, lhs):
        if lhs is None:
            return reversed(self)
        else:
            return NotImplemented

    def __reversed__(self):
        return ProteinMutation(self.chain, self.resnum, self.inscode,
                               self.end, self.start)

    def __len__(self):
        return 1

    def __str__(self):
        return f'{self.chain}-{self.start}{self.full_resnum}{self.end}'

    def __repr__(self):
        return self.__str__()


_FEP_NODE_TITLE_REGEX = re.compile(r'''
    (?P<chain>[A-Z])-       # 1-letter chain ID + literal hyphen
    (?P<start>\w{3})        # starting residue 3-letter code
    (?P<resnum>\d+)         # numeric resnum
    (?P<inscode>[A-Za-z]?)  # optional 1-letter insertion code
    (?P<end>\w{3})          # end residue 3-letter code
    ''', re.VERBOSE)


class ProteinMutationNode():
    '''
    Protein-specific node information.

    '''
    def __init__(self, title, format='fep'):
        format_map = {
            'fep': _FEP_NODE_TITLE_REGEX,
            'rs': _RS_MUTSTR_REGEX,
        }
        self.title = title
        self._mutations = [
            ProteinMutation(**m.groupdict())
            for m in format_map[format].finditer(self.title)
        ]

    @property
    def mutations(self):
        return [
            m for m in self._mutations
            if m is not None
        ]

    @property
    def mutated_sites(self):
        return [
            m.position for m in self._mutations
            if m is not None
        ]

    @property
    def mutations1(self):
        return [
            m.one_letter for m in self._mutations
            if m.one_letter is not None
        ]

    @property
    def mutstr1(self):
        muts1 = [
            m.one_letter for m in self.mutations
            if m.one_letter is not None
        ]
        if len(muts1):
            return ','.join(muts1)
        else:
            return WT_NODE_TITLE


    @property
    def mutstr(self):
        return ','.join([str(m) for m in self.mutations])

    def __eq__(self, rhs):
        if isinstance(rhs, ProteinMutationNode):
            return all([
                self.__getattribute__(x) == rhs.__getattribute__(x)
                for x in ['title', 'mutations']
            ])
        return False

    def __str__(self):
        return self.mutstr

    def get_mutation_at_site(self, site):
        for m in self.mutations:
            if m.site == site:
                return m


class ProteinMutationEdge():
    '''
    Protein-specific edge information.

    init with FEP node titles, e.g. "WT" and "A-ASP100ASER", and optional
    associated graph.Edge instance.

    '''
    def __init__(self, n0_title=None, n1_title=None, edge=None):
        if n0_title and n1_title:
            self.nodes = [
                ProteinMutationNode(n0_title),
                ProteinMutationNode(n1_title),
            ]
        elif edge:
            self.nodes = [
                ProteinMutationNode(n.struc.title)
                for n in edge
            ]

        self.graph_edge = edge
        self._cached_mutation = None
        self._cached_bg_mutations = None

    @property
    def node_titles(self):
        return [n.title for n in self.nodes]

    @property
    def title(self):
        return '{0} -> {1}'.format(*self.node_titles)

    @property
    def mutations(self):
        return [n.mutations for n in self.nodes]

    @property
    def mutation(self):
        if self._cached_mutation is None:
            self._cached_mutation = self._calculate_mutation()
        return self._cached_mutation

    @property
    def bg_mutations(self):
        if self._cached_bg_mutations is None:
            self._cached_bg_mutations = self._calculate_bg_mutations()
        return self._cached_bg_mutations

    @property
    def bg_mutations1(self):
        return [
            m.one_letter for m in self.bg_mutations
            if not m.start1 == m.end1  # remove self/titration mutations
        ]

    @property
    def is_titration(self):
        if self.mutation is NotImplemented:
            return False
        elif len(self.mutation) == 1:
            return (self.mutation.start, self.mutation.end) in pka.PKA_LOOKUP.keys()
        else:
            return False

    @property
    def complex_pka(self):
        return self._calculate_leg_pka('complex')

    @property
    def solvent_pka(self):
        return self._calculate_leg_pka('solvent')

    def __str__(self):
        return self.title

    def __repr__(self):
        return str(self)

    def _calculate_leg_pka(self, leg):
        if self.mutation is NotImplemented or len(self.mutation) != 1:
            return None
        if not self.is_titration:
            return None
        e = self.graph_edge
        n0, n1 = e.nodes
        mdl_sol_shift = n1.pred_solvation_dg - n0.pred_solvation_dg
        sol_cpx_shift = n1.pred_dg - n0.pred_dg
        mdl_cpx_shift = mdl_sol_shift + sol_cpx_shift
        shifts = {
            'complex': mdl_cpx_shift,
            'solvent': mdl_sol_shift,
        }
        res_pair = (self.mutation.start, self.mutation.end)
        return pka.calculate_pka(shifts[leg], 0., res_pair)

    def _calculate_bg_mutations(self):
        '''
        Return a list of ProteinMutation instances common to both endpoints.
        '''
        muts0, muts1 = self.mutations
        return sorted(
            [m for m in muts0 if m in muts1],
            key=lambda m: (m.chain, m.resnum, m.inscode)
        )

    def _calculate_mutation(self):
        '''
        Return a ProteinMutation instance representing the edge perturbation.
        '''
        muts0, muts1 = self.mutations
        # Reverse direction for starting node muts to preserve edge orientation
        uniq_muts0 = [reversed(m) for m in muts0 if m not in muts1]
        uniq_muts1 = [m for m in muts1 if m not in muts0]

        if len(uniq_muts0) == len(uniq_muts1) == 1:
            # There is a different unique mutation at each end; return the
            # connecting mutation.
            edge_mut = uniq_muts0[0] + uniq_muts1[0]

        elif len(uniq_muts0) + len(uniq_muts1) == 1:
            # There is only one unique mutation in total between the two nodes.
            # This is the edge mutation.
            combined_uniq_muts = uniq_muts0 + uniq_muts1
            edge_mut = combined_uniq_muts[0]

        else:
            # Non-single mutations not currently handled
            # TODO: handle multiple mutations
            return NotImplemented
        return edge_mut

    def is_colocalized_with(self, other):
        '''
        Is the other ProteinMutationEdge a mutation of the same position in the
        same mutational background?
        '''
        if self.bg_mutations1 != other.bg_mutations1:
            # print(self.bg_mutations1, other.bg_mutations1)
            return False
        return self.mutation.position == other.mutation.position
