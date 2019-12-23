### CB8 ternary snippets
import rdkit 
from rdkit import Chem
from rdkit.Chem import AllChem

def get_energy_contributions(rdkitmol, label='Guest:', fout='./sample.sdf'):
	"""
	PRE: Takes in a rdkit mol
	POST: Assigns to it energy properties of MMFF94
	"""
	mp = AllChem.MMFFGetMoleculeProperties(rdkitmol)
	for i in range(7):
		termList = [['BondStretch', False], ['AngleBend', False],
		['StretchBend', False], ['OopBend', False], ['Torsion', False],
		['VdW', False], ['Electrostatic', False]]
		termList[i][1] = True
		mp.SetMMFFBondTerm(termList[0][1])
		mp.SetMMFFAngleTerm(termList[1][1])
		mp.SetMMFFStretchBendTerm(termList[2][1])
		mp.SetMMFFOopTerm(termList[3][1])
		mp.SetMMFFTorsionTerm(termList[4][1])
		mp.SetMMFFVdWTerm(termList[5][1])
		mp.SetMMFFEleTerm(termList[6][1])
		ff = AllChem.MMFFGetMoleculeForceField(rdkitmol, mp)
		rdkitmol.SetProp(label+termList[i][0], '{0:12.4f}'.format(ff.CalcEnergy()))
		# print '{0:>16s} energy: {1:12.4f} kcal/mol'.format(termList[i][0],ff.CalcEnergy())
	ff = AllChem.MMFFGetMoleculeForceField(rdkitmol, pyMMFFMolProperties=AllChem.MMFFGetMoleculeProperties(rdkitmol, mmffVariant='MMFF94', mmffVerbosity = 1), ignoreInterfragInteractions=False, nonBondedThresh=100.0) 
	rdkitmol.SetProp(label+'Total', '{0:12.4f}'.format(ff.CalcEnergy()))
	w=Chem.SDWriter(fout)
	w.write(rdkitmol)
	w.close()
	return rdkitmol