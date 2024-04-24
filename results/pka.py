'''
pka.py - Protein residue pKa calculations.

'''
import math

from schrodinger.utils import units


GAS_CONSTANT = units.constants.R / (1000 * units.JOULES_PER_CALORIE)

# Values from:
# [1] Thurlkill, R. L., Grimsley, G. R., Scholtz, J. M. & Pace, C. N.  "pK
#       values of the ionizable groups of proteins."  Protein Sci 15, 1214-1218
#       (2006).
# [2] Tanokura, M.  "1H-NMR study on the tautomerism of the imidazole ring of
#       histidine residues I. Microscopic pK values and molar ratios of
#       tautomers in histidine-containing peptides."  Biochimica Et Biophysica
#       Acta Bba - Protein Struct Mol Enzym 742, 576–585 (1983).
# Also in scisol-src/modules/lambda_dynamics/constant_ph.py:
PKA_ASH_ASP = 3.67   # [1]
PKA_GLH_GLU = 4.25   # [1]
PKA_HIP_HID = 6.92   # [2]
PKA_HIP_HIE = 6.53   # [2]
# Not yet in the constant_ph module:
PKA_CYS_CYM = 8.55   # [1]
PKA_LYS_LYN = 10.40  # [1]
PKA_TYR_TYN = 9.84   # [1]

# pKa values for titratable residues
PKA_LOOKUP = {
    ('CYM', 'CYS'): PKA_CYS_CYM,
    ('CYS', 'CYM'): PKA_CYS_CYM,
    ('ASH', 'ASP'): PKA_ASH_ASP,
    ('ASP', 'ASH'): PKA_ASH_ASP,
    ('GLH', 'GLU'): PKA_GLH_GLU,
    ('GLU', 'GLH'): PKA_GLH_GLU,
    ('HID', 'HIP'): PKA_HIP_HID,
    ('HIE', 'HIP'): PKA_HIP_HIE,
    ('HIP', 'HID'): PKA_HIP_HID,
    ('HIP', 'HIE'): PKA_HIP_HIE,
    ('LYN', 'LYS'): PKA_LYS_LYN,
    ('LYS', 'LYN'): PKA_LYS_LYN,
}

PROTONATED_RESIDUES = ['ASH', 'CYS', 'GLH', 'HIP', 'LYS']


def calculate_pka_shift(dG, dG_model, direction, temp=298.15):
    '''
    Calculate the pKa shift defined by dG - dG_model.

    '''
    ddG = dG - dG_model
    RT = GAS_CONSTANT * temp
    return (-1 * direction * ddG) / (math.log(10) * RT)


def calculate_pka(dG, dG_model, res_pair, temp=298.15):
    '''
    Calculate the (micro) pKa from the given dG values for this residue.

    '''
    try:
        model_pka = PKA_LOOKUP[res_pair]
    except KeyError:
        print(f'{res_pair} is not a supported titration residue pair.')
        raise
    # Deprotonation direction is -1, e.g. ASH->ASP or HIP->HID
    direction = -1 if res_pair[0] in PROTONATED_RESIDUES else 1
    pka_shift = calculate_pka_shift(dG, dG_model, direction, temp)
    pka = model_pka + pka_shift
    # print(res_pair)
    # print(f'  model pKa: {model_pka}')
    # print(f'  ddG model->prot: {dG - dG_model}')
    # print(f'  pka_shift: {pka_shift}')
    return pka
